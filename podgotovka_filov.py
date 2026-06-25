import pandas as pd
import numpy as np
import os
import csv
import shutil

# --- 1. ОЧИСТКА ПАПОК ---
# Удаляем старую папку, если она есть, чтобы не было "мусорных" файлов
if os.path.exists('mutation_matrices'):
    shutil.rmtree('mutation_matrices')
os.makedirs('mutation_matrices', exist_ok=True)

# --- 2. СЛОВАРЬ СКОРОСТЕЙ (BED) ---
# Здесь мы используем точно такие же имена, как в твоем CSV после замены '/' на '_'
bed_raw = """
DYS504.1_DYS660.1	0.009284
DYS504.2_DYS660.2	0.009284
DYS393	0.002598
DYS446	0.0006776
DYS490	0.0004406
DYS505	0.002988
DYS572	0.002419
DYS456	0.003266
DYS570	0.004203
DYS455	0.00214
DYS576	0.004184
DYS525	0.002356
DYS522	0.002766
DYS575	0.002162
DYS463	0.0006822
DYS520_DYS654	0.002318
DYS458	0.004777
DYS450.1	0.0004696
DYS450.2	0.0004696
DYS449.1	0.009642
DYS449.2	0.009642
DYS454_DYS639	0.002182
DYS532.1	0.004167
DYS532.2	0.004167
DYS481	0.006937
DYS531_DYS600	0.002301
DYS590	0.0004325
DYS568	0.0023
DYS487_DYS698	0.000457
DYS19_DYS394	0.002836
DYS716	0.00100
DYS726.1	0.00027
DYS726.2	0.00027
DYS552.1	0.002838
DYS552.2	0.002838
DYS391	0.002016
DYS635	0.002832
DYS434	0.002584
DYS437	0.00233
DYS435	0.002283
DYS439.2	0.002895
DYS439.1	0.002895
DYS389I	0.002196
DYS389II.1	0.002196
DYS389II.2	0.002196
DYS388	0.0004587
DYS442	0.001926
DYS438	0.0007527
DYS441	0.003709
DYS495	0.0005594
DYS436	0.0004414
DYS447	0.0007414
DYS712	0.0467
DYS632	0.00007
DYS413_2	0.00096
DYS641	0.002181
DYS413_1	0.00096
DYS472	0.0004094
DYS565	0.002415
DYS561	0.001827
DYS390_DYS708	0.004661
DYS510.1	0.002407
DYS510.2	0.002407
DYS510.3	0.002407
DYS511	0.002386
DYS717	0.00070
DYS492_DYS604	0.00044
DYS643	0.027
DYS513_DYS605	0.002275
DYS638	0.002245
DYS715	0.00416
DYS606_DYS640	0.002201
DYS587	0.0005036
DYS497	0.0004502
DYS534	0.002851
DYS533	0.002572
DYS607	0.003733
DYS540	0.002314
DYS593	0.0004387
Y-GATA-A10	0.002899
Y-GATA-H4	0.002194
DYS617	0.0005135
DYS426_DYS483	0.0004579
DYS444_DYS542	0.002003
DYS537	0.002278
DYS710	0.01828
YCAII_1	0.00123
DYS425_1	0.00024
DYS425_2	0.00024
YCAII_2	0.00123
DYS385_1	0.00286
DYS385_2	0.00286
DYS461	0.002972
DYS460	0.002488
DYS462	0.002774
DYS494	0.0004187
DYS549	0.002471
DYS452	0.0004213
DYS594	0.0005083
DYS445	0.002467
DYS485	0.0005591
DYS714	0.00773
DYS578	0.002548
DYS556	0.002505
DYS392	0.0004755
DYS636	0.002309
DYS557	0.003315
DYF406S1_DYS555	0.004728
DYS448.1	0.001653
DYS448.2	0.001653
DYS589	0.0006248
DYF386S1_4_DYS504_4	0.009284
DYS464_1.1	0.00566
DYS464_1.2	0.00566
DYS464_2.1	0.00566
DYS464_2.2	0.00566
DYS459_1	0.00132
DYS425_3	0.00024
DYS464_3.1	0.00566
DYS464_3.2	0.00566
DYS464_4.1	0.00566
DYS464_4.2	0.00566
DYS425_4	0.00024
DYS724	0.00245
DYS459_2	0.00132
"""
base_rates = {line.split('\t')[0].strip(): float(line.split('\t')[1]) for line in bed_raw.strip().split('\n')}

# --- 3. ПОДГОТОВКА HAPLOTYPES ---
df = pd.read_csv('Combined_Y-STR_Analysis.csv', sep=';')
df.columns = [c.replace('/', '_').strip() for c in df.columns]

# Уникализация имен колонок (те же DYS389I)
new_cols = []
seen = {}
for col in df.columns:
    if col in seen:
        seen[col] += 1
        new_cols.append(f"{col}_CPY{seen[col]}")
    else:
        seen[col] = 0
        new_cols.append(col)
df.columns = new_cols

# Чистка данных
for col in df.columns[1:]:
    df[col] = pd.to_numeric(df[col].astype(str).str.replace(',', '.'), errors='coerce')

df.to_csv('haplotypes.tsv', sep='\t', index=False)

# --- 4. ГЕНЕРАЦИЯ МАТРИЦ И TSV ---
rates_list = []


def generate_matrix(name, rate):
    alleles = np.arange(0, 60.25, 0.25)
    n = len(alleles)
    matrix = np.zeros((n + 1, n + 1))
    matrix[0, 1:] = alleles
    matrix[1:, 0] = alleles
    for i in range(1, n + 1):
        matrix[i, i] = 1 - rate
        if i > 1: matrix[i, i - 1] = rate / 2
        if i < n: matrix[i, i + 1] = rate / 2
    # Сохраняем без кавычек
    pd.DataFrame(matrix).to_csv(f'mutation_matrices/{name}.tsv', sep='\t', index=False, header=False)


for col in df.columns[1:]:
    # Ищем базовое имя (без _CPY)
    search_name = col.split('_CPY')[0] if '_CPY' in col else col

    if search_name in base_rates:
        rate = base_rates[search_name]
        rates_list.append([col, 0, rate])
        generate_matrix(col, rate)
    else:
        print(f"Warning: Нет скорости для {col}")

# Сохраняем итоговый TSV файл максимально чисто
with open('mutation_rates.tsv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    writer.writerows(rates_list)

print("Подготовка завершена! Теперь запускай ymrca.py")