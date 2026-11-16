import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Ler os arquivos CSV
df = pd.read_csv('model_results_uq.csv')
dh = pd.read_csv('param_values_uq.csv')
data = {}
data['cd4'] = pd.read_csv("CD4_Critical.csv")
data['cd8'] = pd.read_csv("CD8_Critical.csv")
data['virus'] = pd.read_csv("Virus_Sev_Crit.csv")
data['igg'] = pd.read_csv("SNIgG_Critical.csv")

CD4 = data['cd4']
CD8 = data['cd8']
virus = data['virus']
IgG = data['igg']

t = df['Time'].unique()

percl_cd4 = []
perch_cd4 = []
mean_cd4 = []

percl_cd8 = []
perch_cd8 = []
mean_cd8 = []

percl_virus = []
perch_virus = []
mean_virus = []

percl_ant = []
perch_ant = []
mean_ant = []

for i in range(len(t)):
    aux_cd4 = (df[df.Time == t[i]]['The']).to_numpy()
    mean_cd4.append(aux_cd4.mean())
    percl_cd4.append(np.percentile(aux_cd4, 2.5))
    perch_cd4.append(np.percentile(aux_cd4, 97.5))

    aux_cd8 = (df[df.Time == t[i]]['Tke']).to_numpy()
    mean_cd8.append(aux_cd8.mean())
    percl_cd8.append(np.percentile(aux_cd8, 2.5))
    perch_cd8.append(np.percentile(aux_cd8, 97.5))

    aux_virus = (df[df.Time == t[i]]['V']).to_numpy()
    mean_virus.append(aux_virus.mean())
    percl_virus.append(np.percentile(aux_virus, 2.5))
    perch_virus.append(np.percentile(aux_virus, 97.5))

    aux_ant = (df[df.Time == t[i]]['A']).to_numpy()
    mean_ant.append(aux_ant.mean())
    percl_ant.append(np.percentile(aux_ant, 2.5))
    perch_ant.append(np.percentile(aux_ant, 97.5))


# Plotar os dados
plt.figure(figsize=(10, 6))
plt.plot(t, mean_cd4, '-')
plt.fill_between(t, percl_cd4, perch_cd4, alpha=0.2)
plt.title('Células T CD4+ maduras - UQ Parâmetros ED (mais influentes para o caso)')
plt.xlabel('t (dias)')
plt.ylabel('$T_{he}$ (células/mL)')
plt.grid(True)

x = CD4[CD4.type == 'mean']['x']
y = CD4[CD4.type == 'mean']['y']

CD4_up = CD4[CD4.type == 'up']['y']
CD4_down = CD4[CD4.type == 'down']['y']
CD4_mean = CD4[CD4.type == 'mean']['y']

y_error = [CD4_mean.to_numpy() - CD4_down.to_numpy(), CD4_up.to_numpy() - CD4_mean.to_numpy()]

plt.errorbar(x,y, yerr = y_error, linestyle='None', label='Data', fmt='o', color='orange', capsize=4, elinewidth=1)

plt.savefig('cd4uq_influED.png', bbox_inches='tight')
plt.savefig('cd4uq_influED.pdf', bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, mean_cd8, '-')
plt.fill_between(t, percl_cd8, perch_cd8, alpha=0.2)
plt.title('Células T CD8+ maduras - UQ Parâmetros ED (mais influentes para o caso)')
plt.xlabel('t (dias)')
plt.ylabel('$T_{ke}$ (células/mL)')
plt.grid(True)

x = CD8[CD8.type == 'mean']['x']
y = CD8[CD8.type == 'mean']['y']

CD8_up = CD8[CD8.type == 'up']['y']
CD8_down = CD8[CD8.type == 'down']['y']
CD8_mean = CD8[CD8.type == 'mean']['y']

y_error = [CD8_mean.to_numpy() - CD8_down.to_numpy(), CD8_up.to_numpy() - CD8_mean.to_numpy()]

plt.errorbar(x,y, yerr = y_error, linestyle='None', label='Data', fmt='o', color='orange', capsize=4, elinewidth=1)

plt.savefig('cd8uq_influED.png', bbox_inches='tight')
plt.savefig('cd8uq_influED.pdf', bbox_inches='tight')
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(t, np.log10(np.asarray(mean_virus) + 1), '-')
plt.fill_between(t, np.log10(np.asarray(percl_virus) + 1), np.log10(np.asarray(perch_virus) + 1), alpha=0.2)
plt.title('Vírus - UQ Parâmetros ED (mais influentes para o caso)')
plt.xlabel('t (dias)')
plt.ylabel('V ($\\log_{10}$ (cópias/mL + 1))')
plt.grid(True)

x = virus[virus.type == 'mean']['x']
y = np.log10(10**virus[virus.type == 'mean']['y'] + 1)

virus_up = np.log10(10**virus[virus.type == 'up']['y'] + 1)
virus_down = np.log10(10**virus[virus.type == 'down']['y'] + 1)
virus_mean = np.log10(10**virus[virus.type == 'mean']['y'] + 1)

y_error = [virus_mean.to_numpy() - virus_down.to_numpy(), virus_up.to_numpy() - virus_mean.to_numpy()]

plt.errorbar(x,y, yerr = y_error, linestyle='None', label='Data', fmt='o', color='orange', capsize=4, elinewidth=1)

plt.savefig('virusuq_influED.png', bbox_inches='tight')
plt.savefig('virusuq_influED.pdf', bbox_inches='tight')
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(t, np.log10(np.asarray(mean_ant) + 1), '-')
plt.fill_between(t, np.log10(np.asarray(percl_ant) + 1), np.log10(np.asarray(perch_ant) + 1), alpha=0.2)
plt.title('Anticorpos - UQ Parâmetros ED (mais influentes para o caso)')
plt.xlabel('t (dias)')
plt.ylabel('A ($\\log_{10}$ (AU/mL + 1))')
plt.grid(True)

x = IgG[IgG.type == 'mean']['x']
y = np.log10(10**IgG[IgG.type == 'mean']['y'] + 1)

ant_up = np.log10(10**IgG[IgG.type == 'up']['y'] + 1)
ant_down = np.log10(10**IgG[IgG.type == 'down']['y'] + 1)
ant_mean = np.log10(10**IgG[IgG.type == 'mean']['y'] + 1)

y_error = [ant_mean.to_numpy() - ant_down.to_numpy(), ant_up.to_numpy() - ant_mean.to_numpy()]

plt.errorbar(x,y, yerr = y_error, linestyle='None', label='Data', fmt='o', color='orange', capsize=4, elinewidth=1)

plt.savefig('antuq_influED.png', bbox_inches='tight')
plt.savefig('antuq_influED.pdf', bbox_inches='tight')
plt.show()

# histograma
tot = dh['alpha_ap']
suc = dh[dh.Err == 0]['alpha_ap']
err = dh[dh.Err == 1]['alpha_ap']

comb = [tot, suc]


porcentagem = 0.2
auxMin = 1-porcentagem
auxMax = 1+porcentagem

alpha_ap = 0.389196
low = alpha_ap*auxMin
high = alpha_ap*auxMax

plt.figure(figsize=(10, 6))
height, bins, patches = plt.hist(comb, range=[low, high], edgecolor="black", bins = 10, label = ['Total', 'Sucesso'], histtype='bar')
plt.xticks(bins, rotation=45)
plt.title('alpha_ap')
plt.legend()
plt.savefig('histuq_tot_alpha_ap.png', bbox_inches='tight')
plt.savefig('histuq_tot_alpha_ap.pdf', bbox_inches='tight')
plt.show() 

plt.figure(figsize=(10, 6))
height, bins, patches = plt.hist(err, range=[low, high], edgecolor="black", bins = 10, histtype='bar', label = 'Erro', color = 'firebrick')
plt.xticks(bins, rotation=45)
plt.title('alpha_ap')
plt.legend()
plt.savefig('histuq_err_alpha_ap.png', bbox_inches='tight')
plt.savefig('histuq_err_alpha_ap.pdf', bbox_inches='tight')
plt.show() 

