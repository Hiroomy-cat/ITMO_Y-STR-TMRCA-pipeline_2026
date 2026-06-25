import pandas as pd
import numpy as np
import os
import re

# --- 1. ФУНКЦИЯ НОРМАЛИЗАЦИИ ---
def normalize(name):
    """Убирает всё лишнее для сопоставления имен"""
    name = str(name).split('.')[0] if '.' in str(name) and 'DYS' in str(name) else str(name)
    return re.sub(r'[^a-zA-Z0-9]', '', name).lower()

# --- 2. ЗАГРУЗКА ---
# Загружаем частоты мутаций (из твоего BED или файла фильтров)
filters = pd.read_csv('proposed_locus_filters_dp5_minimal.tsv', sep='\t')
rates_dict = {normalize(row['locus']): row['mutation_rate'] for _, row in filters.iterrows()}

# Загружаем гаплотипы Platinum
df = pd.read_csv('Combined_Y-STR_Analysis.csv', sep=';')
df['name'] = df['name'].astype(str).str.strip()

# --- 3. СОПОСТАВЛЕНИЕ И ГЕНЕРАЦИЯ ---
os.makedirs('matrices', exist_ok=True)
alleles = np.arange(0, 60.25, 0.25) # Сетка для Platinum с дробными шагами
n_alleles = len(alleles)

valid_loci = []
haplo_matrix = df[['name']].copy()

print("Начинаю сопоставление локусов...")

# Проходим по всем колонкам кроме 'name'
for col in df.columns:
    if col == 'name': continue
    
    norm_col = normalize(col)
    if norm_col in rates_dict:
        rate = rates_dict[norm_col]
        # Обрабатываем значения (замена , на . и перевод в числа)
        clean_vals = pd.to_numeric(df[col].astype(str).str.replace(',', '.'), errors='coerce')
        
        # Если в колонке есть хотя бы одно число
        if clean_vals.notna().any():
            valid_loci.append(col)
            haplo_matrix[col] = clean_vals
            
            # Генерируем матрицу мутаций для этого локуса
            matrix = np.zeros((n_alleles, n_alleles))
            for i in range(n_alleles):
                matrix[i, i] = 1 - rate
                if i > 0: matrix[i, i-1] = rate / 2
                if i < n_alleles - 1: matrix[i, i+1] = rate / 2
            
            # Сохраняем матрицу (заменяем / на _ для Windows-совместимости имен файлов)
            file_name = col.replace('/', '_')
            np.savetxt(f'matrices/{file_name}.csv', matrix, delimiter=',')
    else:
        # Для диагностики: если локус не найден в словаре скоростей
        pass

# Сохраняем финальные файлы
haplo_matrix.to_csv('haplo_matrix.csv', index=False)
np.savetxt('alleles_scale.txt', alleles)

print(f"ГОРЯЧО! Сопоставлено и сохранено локусов: {len(valid_loci)}")
print(f"Файлы haplo_matrix.csv и папка matrices/ готовы.")
