import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Ler os arquivos CSV
df = pd.read_csv('model_results.csv')
data = {}
data['cd4'] = pd.read_csv("CD4_Critical.csv")
data['cd8'] = pd.read_csv("CD8_Critical.csv")
data['virus'] = pd.read_csv("Virus_Sev_Crit.csv")
data['igg'] = pd.read_csv("SNIgG_Critical.csv")

CD4 = data['cd4']
CD8 = data['cd8']
virus = data['virus']
IgG = data['igg']

# Plotar os dados
plt.figure(figsize=(10, 6))
plt.plot(df['Time'], df['The'], label="Resultado do modelo")
plt.title('Células T CD4+ maduras')
plt.xlabel('t (dias)')
plt.ylabel('$T_{he}$ (células/mL)')
plt.grid(True)

x = CD4[CD4.type == 'mean']['x']
y = CD4[CD4.type == 'mean']['y']

CD4_up = CD4[CD4.type == 'up']['y']
CD4_down = CD4[CD4.type == 'down']['y']
CD4_mean = CD4[CD4.type == 'mean']['y']

y_error = [CD4_mean.to_numpy() - CD4_down.to_numpy(), CD4_up.to_numpy() - CD4_mean.to_numpy()]

plt.errorbar(x,y, yerr = y_error, linestyle='None', label='Dados', fmt='o', color='orange', capsize=4, elinewidth=1)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.savefig('cd4_critical.png', bbox_inches='tight')
#plt.savefig('cd4_critical.pdf', bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(df['Time'], df['Tke'], label="Resultado do modelo")
plt.title('Células T CD8+ maduras')
plt.xlabel('t (dias)')
plt.ylabel('$T_{ke}$ (células/mL)')
plt.grid(True)

x = CD8[CD8.type == 'mean']['x']
y = CD8[CD8.type == 'mean']['y']

CD8_up = CD8[CD8.type == 'up']['y']
CD8_down = CD8[CD8.type == 'down']['y']
CD8_mean = CD8[CD8.type == 'mean']['y']

y_error = [CD8_mean.to_numpy() - CD8_down.to_numpy(), CD8_up.to_numpy() - CD8_mean.to_numpy()]

plt.errorbar(x,y, yerr = y_error, linestyle='None', label='Dados', fmt='o', color='orange', capsize=4, elinewidth=1)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.savefig('cd8_critical.png', bbox_inches='tight')
#plt.savefig('cd8_critical.pdf', bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(df['Time'], np.log10(df['V'] + 1), label="Resultado do modelo")
plt.title('Vírus')
plt.xlabel('t (dias)')
plt.ylabel('V ($\\log_{10}$ (cópias/mL + 1))')
plt.grid(True)

x = virus[virus.type == 'mean']['x']
y = np.log10(10**virus[virus.type == 'mean']['y'] + 1)

virus_up = np.log10(10**virus[virus.type == 'up']['y'] + 1)
virus_down = np.log10(10**virus[virus.type == 'down']['y'] + 1)
virus_mean = np.log10(10**virus[virus.type == 'mean']['y'] + 1)

y_error = [virus_mean.to_numpy() - virus_down.to_numpy(), virus_up.to_numpy() - virus_mean.to_numpy()]

plt.errorbar(x,y, yerr = y_error, linestyle='None', label='Dados', fmt='o', color='orange', capsize=4, elinewidth=1)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.savefig('virus_critical.png', bbox_inches='tight')
#plt.savefig('virus_critical.pdf', bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(df['Time'], np.log10(df['A'] + 1), label="Resultado do modelo")
plt.title('Anticorpos')
plt.xlabel('t (dias)')
plt.ylabel('A ($\\log_{10}$ (AU/mL + 1))')
plt.grid(True)

x = IgG[IgG.type == 'mean']['x']
y = np.log10(10**IgG[IgG.type == 'mean']['y'] + 1)

ant_up = np.log10(10**IgG[IgG.type == 'up']['y'] + 1)
ant_down = np.log10(10**IgG[IgG.type == 'down']['y'] + 1)
ant_mean = np.log10(10**IgG[IgG.type == 'mean']['y'] + 1)

y_error = [ant_mean.to_numpy() - ant_down.to_numpy(), ant_up.to_numpy() - ant_mean.to_numpy()]

plt.errorbar(x,y, yerr = y_error, linestyle='None', label='Dados', fmt='o', color='orange', capsize=4, elinewidth=1)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.savefig('ant_critical.png', bbox_inches='tight')
#plt.savefig('ant_critical.pdf', bbox_inches='tight')
plt.show()