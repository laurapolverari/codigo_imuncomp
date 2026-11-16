import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Ler os arquivos CSV
var = "pi_ps"
df = pd.read_csv(f'model_results_min1_{var}.csv')
df1 = pd.read_csv(f'model_results_min2_{var}.csv')
df2 = pd.read_csv(f'model_results_base_{var}.csv')
df3 = pd.read_csv(f'model_results_max2_{var}.csv')
df4 = pd.read_csv(f'model_results_max1_{var}.csv')

t = df['Time'].unique()

pi_ps = 4.00409160e4

porc1 = 0.2
auxMin1 = 1-porc1
auxMax1 = 1+porc1

porc2 = 0.1
auxMin2 = 1-porc2
auxMax2 = 1+porc2

# Plotar os dados
plt.figure(figsize=(10, 6))
plt.plot(df['Time'], df['The'], label="pi_ps - 20%")
plt.plot(df1['Time'], df1['The'], label="pi_ps - 10%")
plt.plot(df2['Time'], df2['The'], label=f"pi_ps = {pi_ps}")
plt.plot(df3['Time'], df3['The'], label="pi_ps + 10%")
plt.plot(df4['Time'], df4['The'], label="pi_ps + 20%")
plt.title('Células T CD4+ maduras')
plt.xlabel('t (dias)')
plt.ylabel('$T_{he}$ (células/mL)')
plt.grid(True)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.savefig('cd4_moderate_salocal_pi_ps.png', bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(df['Time'], df['Tke'], label="pi_ps - 20%")
plt.plot(df1['Time'], df1['Tke'], label="pi_ps - 10%")
plt.plot(df2['Time'], df2['Tke'], label=f"pi_ps = {pi_ps}")
plt.plot(df3['Time'], df3['Tke'], label="pi_ps + 10%")
plt.plot(df4['Time'], df4['Tke'], label="pi_ps + 20%")
plt.title('Células T CD8+ maduras')
plt.xlabel('t (dias)')
plt.ylabel('$T_{ke}$ (células/mL)')
plt.grid(True)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.savefig('cd8_moderate_salocal_pi_ps.png', bbox_inches='tight')
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(df['Time'], np.log10(df['V'] + 1), label="pi_ps - 20%")
plt.plot(df1['Time'], np.log10(df1['V'] + 1), label="pi_ps - 10%")
plt.plot(df2['Time'], np.log10(df2['V'] + 1), label=f"pi_ps = {pi_ps}")
plt.plot(df3['Time'], np.log10(df3['V'] + 1), label="pi_ps + 10%")
plt.plot(df4['Time'], np.log10(df4['V'] + 1), label="pi_ps + 20%")
plt.title('Vírus')
plt.xlabel('t (dias)')
plt.ylabel('V ($\\log_{10}$ (cópias/mL + 1))')
plt.grid(True)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.savefig('virus_moderate_salocal_pi_ps.png', bbox_inches='tight')
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(df['Time'], np.log10(df['A'] + 1), label="pi_ps - 20%")
plt.plot(df1['Time'], np.log10(df1['A'] + 1), label="pi_ps - 10%")
plt.plot(df2['Time'], np.log10(df2['A'] + 1), label=f"pi_ps = {pi_ps}")
plt.plot(df3['Time'], np.log10(df3['A'] + 1), label="pi_ps + 10%")
plt.plot(df4['Time'], np.log10(df4['A'] + 1), label="pi_ps + 20%")
plt.title('Anticorpos')
plt.xlabel('t (dias)')
plt.ylabel('A ($\\log_{10}$ (AU/mL + 1))')
plt.grid(True)

plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.savefig('ant_moderate_salocal_pi_ps.png', bbox_inches='tight')
plt.show()


dp = porc1 * pi_ps

dThedp = (df4['The'] - df['The']) / (2 * dp)
dTkedp = (df4['Tke'] - df['Tke']) / (2 * dp)
dVdp = (df4['V'] - df['V']) / (2 * dp)
dAdp = (df4['A'] - df['A']) / (2 * dp)

eps = 1e-12

# ========================================================
# SENSIBILIDADE RELATIVA
#
# S = (p/Y)*(dY/dp)
#
# equivalente a:
#
# S = (ΔY/Y)/(Δp/p)
# ========================================================

SThe = (
    (pi_ps / (df2['The'] + eps))
    * dThedp
)

STke = (
    (pi_ps / (df2['Tke'] + eps))
    * dTkedp
)

SV = (
    (pi_ps / (df2['V'] + eps))
    * dVdp
)

SA = (
    (pi_ps / (df2['A'] + eps))
    * dAdp
)

plt.figure(figsize=(10, 6))
plt.plot(t, SThe, label="pi_ps")
plt.xlabel('t (dias)')
plt.ylabel('Sensibilidade relativa')
plt.title('Sensibilidade local de The')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.grid()
plt.savefig('cd4_moderate_sens_pi_ps.png', bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, STke, label="pi_ps")
plt.xlabel('t (dias)')
plt.ylabel('Sensibilidade relativa')
plt.title('Sensibilidade local de Tke')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.grid()
plt.savefig('cd8_moderate_sens_pi_ps.png', bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, SV, label="pi_ps")
plt.xlabel('t (dias)')
plt.ylabel('Sensibilidade relativa')
plt.title('Sensibilidade local de V')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.grid()
plt.savefig('virus_moderate_sens_pi_ps.png', bbox_inches='tight')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, SA, label="pi_ps")
plt.xlabel('t (dias)')
plt.ylabel('Sensibilidade relativa')
plt.title('Sensibilidade local de A')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.grid()
plt.savefig('ant_moderate_sens_pi_ps.png', bbox_inches='tight')
plt.show()

print('Sensibilidade máxima absoluta:')

max_SThe = np.max(np.abs(SThe))
max_STke = np.max(np.abs(STke))
max_SV = np.max(np.abs(SV))
max_SA = np.max(np.abs(SA))

print(f'\nParametro: pi_ps')

print(f'Max |S_The| = {max_SThe:.4f}')
print(f'Max |S_Tke| = {max_STke:.4f}')
print(f'Max |S_V| = {max_SV:.4f}')
print(f'Max |S_A| = {max_SA:.4f}')

rows = [ ['S_The', 'S_Tke', 'S_V', 'S_A'],
         [max_SThe, max_STke, max_SV, max_SA] ]

np.savetxt('sens_max_abs_pi_ps.csv', rows, delimiter=",", fmt='%s')
