import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('FINAL_COUSINS_VALIDATION.tsv', sep='\t')
# Фильтруем те, что посчитались
df_plot = df[df['YMrCA_Mode'] != ">250"].copy()
df_plot['YMrCA_Mode'] = pd.to_numeric(df_plot['YMrCA_Mode'])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

# График А: IBD vs Y-TMRCA
sns.regplot(data=df_plot, x='ibd_sum', y='YMrCA_Mode', ax=ax1, 
            scatter_kws={'s':100, 'alpha':0.6}, line_kws={'color':'red'})
ax1.set_title('Валидация: Корреляция IBD (аутосомы) и TMRCA (Y-STR)')
ax1.set_xlabel('IBD Sum (cM)')
ax1.set_ylabel('Y-TMRCA (поколения)')

# График Б: Ожидание vs Реальность
sns.scatterplot(data=df_plot, x='Expected_Gen', y='YMrCA_Mode', hue='degree', s=150, ax=ax2)
ax2.plot([0, 10], [0, 10], 'k--', alpha=0.5, label='Ideal')
ax2.set_title('Точность: Предсказание YMrCA vs Степень родства')
ax2.set_xlabel('Ожидаемое поколение (по фамилии/degree)')
ax2.set_ylabel('Предсказанное поколение (YMrCA)')
ax2.legend()

plt.tight_layout()
plt.savefig('COUSINS_FINAL_ANALYSIS.png', dpi=300)
plt.show()
