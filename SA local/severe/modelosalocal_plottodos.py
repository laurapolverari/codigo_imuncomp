import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# params = ['V0', 'pi_v', 'beta_tk', 'k_v1',
#           'k_v2', 'beta_ap', 'delta_apm', 'pi_tk',
#           'delta_tk', 'pi_pl', 'delta_a', 'alpha_tk',
#           'beta_th', 'pi_th', 'delta_th']

# vals = [18602.2,0.0030666,45.2902,9.99415e-05,
#         0.00484721,16.3662,2.83713e+08,1836.21,
#         0.0659293,20018.7,11738.9,0.0106936,
#         144.577,27119.2,0.0697364]

params = ['V0', 'pi_v', 'beta_tk', 'k_v1',
          'k_v2', 'beta_ap', 'delta_apm', 'pi_tk',
          'delta_tk', 'pi_pl', 'delta_a', 'alpha_tk',
          'beta_th', 'pi_th', 'delta_th', 'Ap0',
          'Apm0', 'Thn0', 'The0', 'Tkn0',
          'Tke0', 'B0', 'Ps0', 'Pl0',
          'Bm0', 'A0', 'c_v1', 'c_v2',
          'alpha_ap', 'c_ap1', 'c_ap2', 'alpha_th',
          'alpha_b', 'pi_b1', 'pi_b2', 'beta_ps',
          'beta_pl', 'beta_bm', 'delta_ps', 'delta_pl',
          'gamma_bm', 'pi_bm1', 'pi_bm2', 'pi_ps']

           

vals = [18602.2,0.0030666,45.2902,9.99415e-05,
        0.00484721,16.3662,2.83713e+08,1836.21,
        0.0659293,20018.7,11738.9,0.0106936,
        144.577,27119.2,0.0697364, 1e4,
        0.0, 1e3, 0.0, 5.0e3,
        0.0, 2.5e5, 0.0, 0.0,
        0.0, 0.0, 1.12583762e-9, 6.0e-1,
        2.59993413e-1, 4.08257620e0, 4.85488215e5, 2.17e-1,
        3.58e2, 8.98e-5, 1.27e-8, 6.0e-6,
        5.0e-6, 1.0e-6, 2.5e0, 3.5e-1,
        9.75e-4, 8.11167858e0, 3.79652355e3, 4.00409160e4]

porc1 = 0.2
auxMin1 = 1-porc1
auxMax1 = 1+porc1

porc2 = 0.1
auxMin2 = 1-porc2
auxMax2 = 1+porc2

eps = 1e-12

# para guardar as sensibilidades máximas absolutas de cada parâmetro
maxSThe = []
maxSTke = []
maxSV = []
maxSA = []

plt.figure(1, figsize=(10, 6))
plt.figure(2, figsize=(10, 6))
plt.figure(3, figsize=(10, 6))
plt.figure(4, figsize=(10, 6))

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',
          'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
          'tab:olive', 'tab:cyan', 'midnightblue', 'springgreen',
          'darkred', 'magenta', 'indigo']

linestyles = ['-', '--', '-.']

f = open("sensmaxabsparams.txt", "w")

for i in range(len(params)):

    df = pd.read_csv(f'Arq/model_results_min1_{params[i]}.csv')
    df1 = pd.read_csv(f'Arq/model_results_min2_{params[i]}.csv')
    df2 = pd.read_csv(f'Arq/model_results_base_{params[i]}.csv')
    df3 = pd.read_csv(f'Arq/model_results_max2_{params[i]}.csv')
    df4 = pd.read_csv(f'Arq/model_results_max1_{params[i]}.csv')

    t = df['Time'].unique()

    dp = porc1 * vals[i]

    dThedp = (df4['The'] - df['The']) / (2 * dp)
    dTkedp = (df4['Tke'] - df['Tke']) / (2 * dp)
    dVdp = (df4['V'] - df['V']) / (2 * dp)
    dAdp = (df4['A'] - df['A']) / (2 * dp)

    SThe = (
        (vals[i] / (df2['The'] + eps))
        * dThedp
    )

    STke = (
        (vals[i] / (df2['Tke'] + eps))
        * dTkedp
    )

    SV = (
        (vals[i] / (df2['V'] + eps))
        * dVdp
    )

    SA = (
        (vals[i] / (df2['A'] + eps))
        * dAdp
    )

    j =  int(i / len(colors))
    print("j: ", j)

    # Gráfico com a sensibilidade local das quatro populações para o parâmetro atual
    plt.figure(figsize=(10, 6))
    plt.plot(t, SThe, label=f"The")
    plt.plot(t, STke, label=f"Tke")
    plt.plot(t, SV, label=f"V")
    plt.plot(t, SA, label=f"A")
    plt.xlabel('t (dias)')
    plt.ylabel('Sensibilidade relativa')
    plt.title(f'Sensibilidade local de {params[i]}')
    plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.grid()
    plt.savefig(f'severe_sens_{params[i]}.png', bbox_inches='tight')
    #plt.show()
    plt.close()

    plt.figure(1)
    plt.plot(t, SThe, color = colors[i - len(colors)*j], linestyle = linestyles[j], label=f"{params[i]}")

    plt.figure(2)
    plt.plot(t, STke, color = colors[i - len(colors)*j], linestyle = linestyles[j], label=f"{params[i]}")

    plt.figure(3)
    plt.plot(t, SV, color = colors[i - len(colors)*j], linestyle = linestyles[j], label=f"{params[i]}")

    plt.figure(4)
    plt.plot(t, SA, color = colors[i - len(colors)*j], linestyle = linestyles[j], label=f"{params[i]}")

    print('Sensibilidade máxima absoluta:')
    print('Sensibilidade máxima absoluta:', file=f)

    max_SThe = np.max(np.abs(SThe))
    max_STke = np.max(np.abs(STke))
    max_SV = np.max(np.abs(SV))
    max_SA = np.max(np.abs(SA))

    print(f'\nParâmetro: {params[i]}')

    print(f'Max |S_The| = {max_SThe:.4f}')
    print(f'Max |S_Tke| = {max_STke:.4f}')
    print(f'Max |S_V| = {max_SV:.4f}')
    print(f'Max |S_A| = {max_SA:.4f}\n\n')

    print(f'\nParâmetro: {params[i]}', file=f)

    print(f'Max |S_The| = {max_SThe:.4f}', file=f)
    print(f'Max |S_Tke| = {max_STke:.4f}', file=f)
    print(f'Max |S_V| = {max_SV:.4f}', file=f)
    print(f'Max |S_A| = {max_SA:.4f}\n\n', file=f)

    if max_SThe > 1.0:
        maxSThe.append((i, max_SThe))

    if max_STke > 1.0:
        maxSTke.append((i, max_STke))

    if max_SV > 1.0:
        maxSV.append((i, max_SV))

    if max_SA > 1.0:
        maxSA.append((i, max_SA))

f.close()

plt.figure(1)
plt.xlabel('t (dias)')
plt.ylabel('Sensibilidade relativa')
plt.title('Sensibilidade local de The')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", ncol=2)
plt.grid()
plt.savefig('cd4_severe_sens_tot.png', bbox_inches='tight')
# plt.show()


plt.figure(2)
plt.xlabel('t (dias)')
plt.ylabel('Sensibilidade relativa')
plt.title('Sensibilidade local de Tke')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", ncol=2)
plt.grid()
plt.savefig('cd8_severe_sens_tot.png', bbox_inches='tight')
# plt.show()


plt.figure(3)
plt.xlabel('t (dias)')
plt.ylabel('Sensibilidade relativa')
plt.title('Sensibilidade local de V')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", ncol=2)
plt.grid()
plt.savefig('virus_severe_sens_tot.png', bbox_inches='tight')
# plt.show()


plt.figure(4)
plt.xlabel('t (dias)')
plt.ylabel('Sensibilidade relativa')
plt.title('Sensibilidade local de A')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", ncol=2)
plt.grid()
plt.savefig('ant_severe_sens_tot.png', bbox_inches='tight')
# plt.show()

plt.figure(5, figsize=(10, 6))
plt.figure(6, figsize=(10, 6))
plt.figure(7, figsize=(10, 6))
plt.figure(8, figsize=(10, 6))

if len(maxSThe) > 0:    
    maxSThe = sorted(maxSThe,key=lambda x: x[1], reverse=True)

    for j in range(len(maxSThe)):
        if j == 5:
            break
        else:
            i = maxSThe[j][0]
            df = pd.read_csv(f'Arq/model_results_min1_{params[i]}.csv')
            df2 = pd.read_csv(f'Arq/model_results_base_{params[i]}.csv')
            df4 = pd.read_csv(f'Arq/model_results_max1_{params[i]}.csv')

            t = df['Time'].unique()

            dp = porc1 * vals[i]

            dThedp = (df4['The'] - df['The']) / (2 * dp)

            SThe = (
                (vals[i] / (df2['The'] + eps))
                * dThedp
            )

            plt.figure(5)
            if(params[i] == 'beta_tk'):
                plt.plot(t, SThe, '--', label=f"{params[i]}")
            else:
                plt.plot(t, SThe, label=f"{params[i]}")

    plt.figure(5)
    plt.xlabel('t (dias)')
    plt.ylabel('Sensibilidade relativa')
    plt.title('Sensibilidade local de The')
    plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.grid()
    plt.savefig('cd4_severe_sens_selec.png', bbox_inches='tight')
    # plt.show()
  

if len(maxSTke) > 0:    
    maxSTke = sorted(maxSTke,key=lambda x: x[1], reverse=True)

    for j in range(len(maxSTke)):
        if j == 5:
            break
        else:
            i = maxSTke[j][0]
            df = pd.read_csv(f'Arq/model_results_min1_{params[i]}.csv')
            df2 = pd.read_csv(f'Arq/model_results_base_{params[i]}.csv')
            df4 = pd.read_csv(f'Arq/model_results_max1_{params[i]}.csv')

            t = df['Time'].unique()

            dp = porc1 * vals[i]

            dTkedp = (df4['Tke'] - df['Tke']) / (2 * dp)
            
            STke = (
                (vals[i] / (df2['Tke'] + eps))
                * dTkedp
            )

            plt.figure(6)
            plt.plot(t, STke, label=f"{params[i]}")

    plt.figure(6)
    plt.xlabel('t (dias)')
    plt.ylabel('Sensibilidade relativa')
    plt.title('Sensibilidade local de Tke')
    plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.grid()
    plt.savefig('cd8_severe_sens_selec.png', bbox_inches='tight')
    # plt.show()


if len(maxSV) > 0:    
    maxSV = sorted(maxSV,key=lambda x: x[1], reverse=True)
    
    for j in range(len(maxSV)):
        if j == 5:
            break
        else:
            i = maxSV[j][0]
            df = pd.read_csv(f'Arq/model_results_min1_{params[i]}.csv')
            df2 = pd.read_csv(f'Arq/model_results_base_{params[i]}.csv')
            df4 = pd.read_csv(f'Arq/model_results_max1_{params[i]}.csv')
  
            t = df['Time'].unique()

            dp = porc1 * vals[i]

            dVdp = (df4['V'] - df['V']) / (2 * dp)
            
            SV = (
                (vals[i] / (df2['V'] + eps))
                * dVdp
            )

            plt.figure(7)
            plt.plot(t, SV, label=f"{params[i]}")

    plt.figure(7)
    plt.xlabel('t (dias)')
    plt.ylabel('Sensibilidade relativa')
    plt.title('Sensibilidade local de V')
    plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.grid()
    plt.savefig('virus_severe_sens_selec.png', bbox_inches='tight')
    # plt.show()


if len(maxSA) > 0:    
    maxSA = sorted(maxSA,key=lambda x: x[1], reverse=True)
    
    for j in range(len(maxSA)):
        if j == 5:
            break
        else:
            i = maxSA[j][0]
            df = pd.read_csv(f'Arq/model_results_min1_{params[i]}.csv')
            df2 = pd.read_csv(f'Arq/model_results_base_{params[i]}.csv')
            df4 = pd.read_csv(f'Arq/model_results_max1_{params[i]}.csv')

            t = df['Time'].unique()

            dp = porc1 * vals[i]

            dAdp = (df4['A'] - df['A']) / (2 * dp)

            SA = (
                (vals[i] / (df2['A'] + eps))
                * dAdp
            )

            plt.figure(8)
            if(params[i] == 'beta_tk'):
                plt.plot(t, SA, '--', label=f"{params[i]}")
            else:
                plt.plot(t, SA, label=f"{params[i]}")
    
    plt.figure(8)
    plt.xlabel('t (dias)')
    plt.ylabel('Sensibilidade relativa')
    plt.title('Sensibilidade local de A')
    plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.grid()
    plt.savefig('ant_severe_sens_selec.png', bbox_inches='tight')
    # plt.show()

