#!/usr/bin/env python3
"""
Расчёт TMRCA (Walsh 2001) для таблицы попарных разностей гаплотипов.

Принимает на вход:
  --diffs   Файл с разностями (TSV или пробел-разделённый, с заголовком):
              sample1  sample2  LOCUS1  LOCUS2  LOCUS3 ...
              S0001    S0002    0       1       2      ...
              S0003    S0004    0       .       5      ...
            Значения — разности в числе STR-повторов (целые числа).
            Точка (.) означает, что локус не определён и пропускается.

  --bed     BED-файл с 7 колонками (без заголовка), разделитель — табуляция:
              chrom start end unit_size ref_repeats marker rate
            Используется последний столбец (rate). NA означает отсутствие
            мутационной частоты: такой локус участвует в расчёте постериора,
            но не влияет на среднюю μ.

  --output  Выходной TSV с результатами TMRCA для каждой пары.

Пример запуска:
  python tmrca_from_diffs.py \\
      --diffs pairwise_diffs.tsv \\
      --bed   panel/merged_hg19_full.bed \\
      --output TMRCA_results.tsv \\
      --lam 0.0 --min-loci 15 --threads 4
"""

import argparse
import multiprocessing
from collections import Counter

import numpy as np
import pandas as pd
from scipy.special import iv as bessel_iv
from scipy.signal import savgol_filter
from tqdm import tqdm


# ─────────────────────────────────────────────────────────────────────────────
# Параметры численной сетки для постериора
# ─────────────────────────────────────────────────────────────────────────────
_N_GRID = 200_000   # число точек сетки (разрешение ≈ 0.005 поколений)
_TT_MAX = 1_000.0   # максимальное время, поколений


def _make_grid():
    """Равномерная сетка по t ∈ (0, _TT_MAX]."""
    return np.linspace(1e-3, _TT_MAX, _N_GRID)


# ─────────────────────────────────────────────────────────────────────────────
# Walsh (2001): IAM-постериор
# Модель бесконечного числа аллелей: учитывается только факт наличия/отсутствия
# мутации, но не её размер.
# pi(t) ∝ (1 − e^{−2μt})^{n_diff} · e^{−t(2μk + λ)}
# ─────────────────────────────────────────────────────────────────────────────
def _log_weights_iam(tt, diffs, lam, mu):
    diffs  = np.asarray(diffs, dtype=float)
    n      = len(diffs)
    k      = int(np.sum(diffs == 0.0))    # локусов без разности
    n_diff = n - k                          # локусов с разностью
    if n_diff == 0:
        return -lam * tt
    z     = 2.0 * mu * tt
    # log1p(-exp(-z)) устойчиво при малых z
    log_p = np.log1p(-np.exp(-np.clip(z, 1e-12, None)))
    return n_diff * log_p - tt * (2.0 * mu * k + lam)


# ─────────────────────────────────────────────────────────────────────────────
# Walsh (2001): SMM-постериор
# Пошаговая модель мутаций: размер разности учитывается через функцию Бесселя.
# pi(t) ∝ e^{−(λ+2nμ)t} · ∏_i [I_{|d_i|}(2μt)]^{c_i}
# ─────────────────────────────────────────────────────────────────────────────
def _log_weights_smm(tt, diffs, lam, mu):
    diffs = np.asarray(diffs, dtype=float)
    n     = len(diffs)
    tab   = Counter(float(d) for d in diffs)   # частоты каждой уникальной разности
    arg   = 2.0 * mu * tt                        # аргумент для функции Бесселя
    log_w = -(lam + 2.0 * n * mu) * tt
    for val, count in tab.items():
        bess  = bessel_iv(abs(val), arg)          # I_{|d|}(2μt)
        log_w = log_w + count * np.log(np.where(bess > 1e-300, bess, 1e-300))
    return log_w


# ─────────────────────────────────────────────────────────────────────────────
# Статистики постериора из лог-весов на сетке
# ─────────────────────────────────────────────────────────────────────────────
def _stats_from_log_weights(tt, log_w):
    log_w  = log_w - log_w.max()    # нормировка для числовой устойчивости
    w      = np.exp(log_w)
    w     /= w.sum()

    mean   = float(np.dot(tt, w))
    cumsum = np.cumsum(w)
    median = float(tt[np.searchsorted(cumsum, 0.500)])
    ci_lo  = float(tt[np.searchsorted(cumsum, 0.025)])
    ci_hi  = float(tt[np.searchsorted(cumsum, 0.975)])

    # Мода после лёгкого сглаживания (убирает артефакты сетки у пика)
    win  = min(101, (_N_GRID // 1000) * 2 + 1)
    w_sm = savgol_filter(w, win, polyorder=3) if win >= 5 else w
    mode = float(tt[np.argmax(np.clip(w_sm, 0, None))])

    return {"mean": mean, "median": median, "mode": mode,
            "CI_2.5pct": ci_lo, "CI_97.5pct": ci_hi}


def walsh_posterior(diffs, lam, mu, model):
    """Постериор Walsh (2001) для одной пары, одной модели."""
    tt = _make_grid()
    if model == "IAM":
        log_w = _log_weights_iam(tt, diffs, lam, mu)
    elif model == "SMM":
        log_w = _log_weights_smm(tt, diffs, lam, mu)
    else:
        raise ValueError(f"Неизвестная модель: {model}")
    return _stats_from_log_weights(tt, log_w)


# ─────────────────────────────────────────────────────────────────────────────
# Глобальные переменные рабочих процессов (инициализируются один раз через fork)
# ─────────────────────────────────────────────────────────────────────────────
_PAIRS_DF    = None    # DataFrame: строки = пары, колонки = sample1/sample2/локусы
_RATED_LOCI  = None    # set: локусы с известной μ из BED-файла
_RATES       = None    # dict {locus: μ}
_LAM         = None    # параметр λ
_MIN_LOCI    = None    # минимальный порог локусов для расчёта

_NAN_KEYS = ("mean", "median", "mode", "CI_2.5pct", "CI_97.5pct")


def _init_pool(pairs_df, rated_loci, rates, lam, min_loci):
    """Инициализатор пула: записывает общие данные в глобальные переменные."""
    global _PAIRS_DF, _RATED_LOCI, _RATES, _LAM, _MIN_LOCI
    _PAIRS_DF   = pairs_df
    _RATED_LOCI = rated_loci
    _RATES      = rates
    _LAM        = lam
    _MIN_LOCI   = min_loci


def _process_row(idx):
    """
    Обрабатывает одну строку входного файла (одну пару сэмплов).

    Логика:
    1. Отбираем локусы без точки (валидные разности).
    2. Средняя μ считается только по локусам с известной мутационной частотой.
    3. Локусы без μ всё равно входят в модель, используя уже посчитанную среднюю μ.
    4. Запускаем Walsh (2001) постериор для IAM и SMM.
    """
    row = _PAIRS_DF.iloc[idx]
    s1  = row["sample1"]
    s2  = row["sample2"]

    locus_cols = [c for c in _PAIRS_DF.columns if c not in ("sample1", "sample2")]

    # Локусы с известной разностью (не точка)
    valid_loci = [c for c in locus_cols if not pd.isna(row[c])]
    n_loci     = len(valid_loci)

    result = {"sample1": s1, "sample2": s2, "n_loci": n_loci}

    # Слишком мало локусов — пропускаем
    if n_loci < _MIN_LOCI:
        for model in ("IAM", "SMM"):
            for k in _NAN_KEYS:
                result[f"{model}_{k}"] = float("nan")
        return result

    # Средняя μ — только по локусам, у которых есть мутационная частота в BED
    rated_here = [c for c in valid_loci if c in _RATED_LOCI]
    if not rated_here:
        # Ни одного локуса с известной μ — расчёт невозможен
        for model in ("IAM", "SMM"):
            for k in _NAN_KEYS:
                result[f"{model}_{k}"] = float("nan")
        return result

    avg_mu = float(np.mean([_RATES[c] for c in rated_here]))

    # Векторы разностей для всех валидных локусов (в т.ч. без μ в BED)
    iam_diffs = []
    smm_diffs = []
    for c in valid_loci:
        d = float(row[c])
        iam_diffs.append(0 if d == 0.0 else 1)   # IAM: 0/1
        smm_diffs.append(abs(d))                   # SMM: |Δ|

    # Постериор Walsh (2001) для IAM и SMM
    for model, diffs in [("IAM", iam_diffs), ("SMM", smm_diffs)]:
        try:
            stats = walsh_posterior(diffs, _LAM, avg_mu, model)
        except Exception:
            stats = {k: float("nan") for k in _NAN_KEYS}
        for k, v in stats.items():
            result[f"{model}_{k}"] = v

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции
# ─────────────────────────────────────────────────────────────────────────────
def load_bed(bed_path):
    """
    Загружает BED-файл (7 колонок, без заголовка).
    Возвращает dict {marker: rate} только для локусов с числовой μ (не NA).
    """
    cols = ["chrom", "start", "end", "unit_size", "ref_repeats", "marker", "rate"]
    bed  = pd.read_csv(bed_path, sep="\t", header=None, names=cols)
    bed["rate"] = pd.to_numeric(bed["rate"], errors="coerce")
    rates = {
        row["marker"]: float(row["rate"])
        for _, row in bed.iterrows()
        if not pd.isna(row["rate"])
    }
    return rates


def load_diffs(diffs_path):
    """
    Загружает таблицу разностей.
    Поддерживает TSV и файлы с разделением пробелами.
    Точку (.) заменяет на NaN; числа приводит к float.
    """
    # Пробуем автоматически определить разделитель
    df = pd.read_csv(diffs_path, sep=r"\s+", dtype=str, engine="python")

    # Проверяем наличие обязательных колонок
    if "sample1" not in df.columns or "sample2" not in df.columns:
        raise ValueError(
            "Первые два столбца должны называться 'sample1' и 'sample2'. "
            "Убедитесь, что в файле есть заголовок."
        )

    locus_cols = [c for c in df.columns if c not in ("sample1", "sample2")]
    df[locus_cols] = df[locus_cols].replace(".", np.nan).astype(float)
    return df


def parse_args():
    p = argparse.ArgumentParser(
        description="Расчёт TMRCA (Walsh 2001) для таблицы пар сэмплов.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--diffs",    required=True,
                   help="Файл с разностями: sample1 sample2 LOCUS1 LOCUS2 ...")
    p.add_argument("--bed",      required=True,
                   help="BED-файл с мутационными частотами (7 колонок, без заголовка)")
    p.add_argument("--output",   required=True,
                   help="Выходной TSV с результатами TMRCA")
    p.add_argument("--lam",      type=float, default=0.0,
                   help="Параметр λ (prior; по умолчанию 0.0)")
    p.add_argument("--min-loci", type=int,   default=15,
                   help="Минимальное число локусов для расчёта (по умолчанию 15)")
    p.add_argument("--threads",  type=int,   default=1,
                   help="Число параллельных процессов (по умолчанию 1)")
    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# Точка входа
# ─────────────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    # Загрузка мутационных частот из BED
    rates      = load_bed(args.bed)
    rated_loci = set(rates.keys())
    print(f"[BED]   Загружено {len(rates)} локусов с мутационными частотами.")

    # Загрузка таблицы разностей
    pairs_df   = load_diffs(args.diffs)
    n_pairs    = len(pairs_df)
    locus_cols = [c for c in pairs_df.columns if c not in ("sample1", "sample2")]
    n_rated_in_input = sum(1 for c in locus_cols if c in rated_loci)
    print(f"[INPUT] {n_pairs} пар, {len(locus_cols)} локусов "
          f"(из них с μ: {n_rated_in_input}).")
    print(f"[PARAMS] λ={args.lam}, min_loci={args.min_loci}, threads={args.threads}")

    # Параллельный расчёт TMRCA для каждой пары
    ctx     = multiprocessing.get_context("fork")
    indices = list(range(n_pairs))
    records = []

    with ctx.Pool(
        processes=args.threads,
        initializer=_init_pool,
        initargs=(pairs_df, rated_loci, rates, args.lam, args.min_loci),
    ) as pool:
        with tqdm(total=n_pairs, unit="pair", desc="TMRCA", ncols=80) as pbar:
            for res in pool.imap_unordered(_process_row, indices, chunksize=10):
                records.append(res)
                pbar.update(1)

    # Формируем и сохраняем результат
    col_order = [
        "sample1", "sample2", "n_loci",
        "IAM_mean", "IAM_median", "IAM_mode", "IAM_CI_2.5pct", "IAM_CI_97.5pct",
        "SMM_mean", "SMM_median", "SMM_mode", "SMM_CI_2.5pct", "SMM_CI_97.5pct",
    ]
    df_out    = pd.DataFrame(records)
    col_order = [c for c in col_order if c in df_out.columns]
    df_out[col_order].to_csv(args.output, sep="\t", index=False, float_format="%.2f")
    print(f"[OUT]   Результаты сохранены: {args.output}")


if __name__ == "__main__":
    main()
