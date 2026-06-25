import pandas as pd
import numpy as np
import plotly.express as px
import os

# 1. ЗАГРУЗКА МЕТАДАННЫХ
df_meta = pd.read_csv('grishin_18.06/pairs_metadata.tsv', sep='\t')
# Порядок пар в метаданных соответствует номерам файлов couple1...15

results = []
for i in range(15):
    file_num = i + 1
    filename = f'grishin_18.06/couple{file_num}.tsv'
    
    if os.path.exists(filename):
        probs = pd.read_csv(filename, sep='\t', header=None).values.flatten()
        if np.sum(probs) > 0:
            mode_gen = np.argmax(probs) + 1
            # Считаем CI
            cum = np.cumsum(probs) / np.sum(probs)
            ci_low = np.where(cum >= 0.025)[0][0] + 1
            ci_high = np.where(cum >= 0.975)[0][0] + 1
        else:
            mode_gen, ci_low, ci_high = None, None, None
            
        results.append({
            'YMrCA_Mode': mode_gen,
            'CI': f"{ci_low}-{ci_high}" if mode_gen else "N/A"
        })

# Объединяем метаданные с результатами
df_final = pd.concat([df_meta, pd.DataFrame(results)], axis=1)

# --- ГРАФИК: IBD vs TMRCA (Главное доказательство) ---
fig = px.scatter(df_final.dropna(subset=['YMrCA_Mode']), 
                 x="ibd_sum", y="YMrCA_Mode", color="degree",
                 hover_data=["tubeid1", "tubeid2", "haplogroup"],
                 trendline="ols",
                 title="<b>Валидация: Связь IBD (весь геном) и TMRCA (Y-хромосома)</b>",
                 labels={"ibd_sum": "Сумма IBD сегментов (cM)", "YMrCA_Mode": "TMRCA (поколения)"},
                 template="plotly_white")

fig.write_html("COUSINS_IBD_VALIDATION.html")
print("Отчет по кузенам готов: COUSINS_IBD_VALIDATION.html")
fig.show()
