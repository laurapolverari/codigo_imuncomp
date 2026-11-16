import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Ler os arquivos CSV
df = pd.read_csv('model_results_uq.csv')

t = df['Time'].unique()
#simnum = df['SimNum'].unique()
 
N_list = [625, 1250, 2500, 5000, 10000]
N_max = max(N_list)
N_total = int(N_max)

# Organizando os resultados de cada simulação para cada população
CD4_all = np.zeros((N_total, len(t)))
CD8_all = np.zeros((N_total, len(t)))
V_all = np.zeros((N_total, len(t)))
A_all = np.zeros((N_total, len(t)))

for i in range(N_total):
    CD4_all[i] = (df[df.SimNum == (i+1)]['The']).to_numpy()
    CD8_all[i] = (df[df.SimNum == (i+1)]['Tke']).to_numpy()
    V_all[i] = (df[df.SimNum == (i+1)]['V']).to_numpy()
    A_all[i] = (df[df.SimNum == (i+1)]['A']).to_numpy()


# ============================================================
# CALCULA AS ESTATÍSTICAS UTILIZANDO SUBCONJUNTOS DAS AMOSTRAS
# ============================================================

results = {}

for N in N_list:


    results[N] = {

        "CD4_mean": np.mean(CD4_all[:N], axis=0),
        "CD8_mean": np.mean(CD8_all[:N], axis=0),
        "V_mean": np.mean(V_all[:N], axis=0),
        "A_mean": np.mean(A_all[:N], axis=0),

        "CD4_p2": np.percentile(CD4_all[:N], 2.5, axis=0),
        "CD4_p97": np.percentile(CD4_all[:N], 97.5, axis=0),

        "CD8_p2": np.percentile(CD8_all[:N], 2.5, axis=0),
        "CD8_p97": np.percentile(CD8_all[:N], 97.5, axis=0),

        "V_p2": np.percentile(V_all[:N], 2.5, axis=0),
        "V_p97": np.percentile(V_all[:N], 97.5, axis=0),

        "A_p2": np.percentile(A_all[:N], 2.5, axis=0),
        "A_p97": np.percentile(A_all[:N], 97.5, axis=0)

    }


# ============================================================
# COMPARAÇÃO ENTRE VALORES CONSECUTIVOS
# ============================================================

erro_CD4_mean = []
erro_CD8_mean = []
erro_V_mean = []
erro_A_mean = []

erro_CD4_p2 = []
erro_CD4_p97 = []

erro_CD8_p2 = []
erro_CD8_p97 = []

erro_V_p2 = []
erro_V_p97 = []

erro_A_p2 = []
erro_A_p97 = []

print("\n")
print("=" * 95)
print("CONVERGÊNCIA ENTRE VALORES CONSECUTIVOS")
print("=" * 95)

print(
    f"{'Comparação':<15}"
    f"{'The média':>12}"
    f"{'The P2.5':>12}"
    f"{'The P97.5':>12}"
    f"{'Tke média':>12}"
    f"{'Tke P2.5':>12}"
    f"{'Tke P97.5':>12}"
    f"{'V média':>12}"
    f"{'V P2.5':>12}"
    f"{'V P97.5':>12}"
    f"{'A média':>12}"
    f"{'A P2.5':>12}"
    f"{'A P97.5':>12}"
)

for i in range(len(N_list) - 1):

    N1 = N_list[i]
    N2 = N_list[i + 1]

    r1 = results[N1]
    r2 = results[N2]

    err_CD4_mean = np.linalg.norm(
        r1["CD4_mean"] - r2["CD4_mean"]
    ) / np.linalg.norm(r2["CD4_mean"])

    err_CD8_mean = np.linalg.norm(
        r1["CD8_mean"] - r2["CD8_mean"]
    ) / np.linalg.norm(r2["CD8_mean"])

    err_V_mean = np.linalg.norm(
        r1["V_mean"] - r2["V_mean"]
    ) / np.linalg.norm(r2["V_mean"])

    err_A_mean = np.linalg.norm(
        r1["A_mean"] - r2["A_mean"]
    ) / np.linalg.norm(r2["A_mean"])

    err_CD4_p2 = np.linalg.norm(
        r1["CD4_p2"] - r2["CD4_p2"]
    ) / np.linalg.norm(r2["CD4_p2"])

    err_CD4_p97 = np.linalg.norm(
        r1["CD4_p97"] - r2["CD4_p97"]
    ) / np.linalg.norm(r2["CD4_p97"])

    err_CD8_p2 = np.linalg.norm(
        r1["CD8_p2"] - r2["CD8_p2"]
    ) / np.linalg.norm(r2["CD8_p2"])

    err_CD8_p97 = np.linalg.norm(
        r1["CD8_p97"] - r2["CD8_p97"]
    ) / np.linalg.norm(r2["CD8_p97"])

    err_V_p2 = np.linalg.norm(
        r1["V_p2"] - r2["V_p2"]
    ) / np.linalg.norm(r2["V_p2"])

    err_V_p97 = np.linalg.norm(
        r1["V_p97"] - r2["V_p97"]
    ) / np.linalg.norm(r2["V_p97"])

    err_A_p2 = np.linalg.norm(
        r1["A_p2"] - r2["A_p2"]
    ) / np.linalg.norm(r2["A_p2"])

    err_A_p97 = np.linalg.norm(
        r1["A_p97"] - r2["A_p97"]
    ) / np.linalg.norm(r2["A_p97"])

    erro_CD4_mean.append(err_CD4_mean)
    erro_CD8_mean.append(err_CD8_mean)
    erro_V_mean.append(err_V_mean)
    erro_A_mean.append(err_A_mean)

    erro_CD4_p2.append(err_CD4_p2)
    erro_CD4_p97.append(err_CD4_p97)

    erro_CD8_p2.append(err_CD8_p2)
    erro_CD8_p97.append(err_CD8_p97)

    erro_V_p2.append(err_V_p2)
    erro_V_p97.append(err_V_p97)

    erro_A_p2.append(err_A_p2)
    erro_A_p97.append(err_A_p97)

    print(
        f"{str(N1)+'->'+str(N2):<15}"
        f"{err_CD4_mean:12.3e}"
        f"{err_CD4_p2:12.3e}"
        f"{err_CD4_p97:12.3e}"
        f"{err_CD8_mean:12.3e}"
        f"{err_CD8_p2:12.3e}"
        f"{err_CD8_p97:12.3e}"
        f"{err_V_mean:12.3e}"
        f"{err_V_p2:12.3e}"
        f"{err_V_p97:12.3e}"
        f"{err_A_mean:12.3e}"
        f"{err_A_p2:12.3e}"
        f"{err_A_p97:12.3e}"
    )


# ============================================================
# GRÁFICO
# ============================================================

plt.figure(figsize=(8, 6))

plt.semilogy(N_list[1:], erro_CD4_mean, 'o-', label='The média')
plt.semilogy(N_list[1:], erro_CD4_p2, 'o--', label='The P2.5')
plt.semilogy(N_list[1:], erro_CD4_p97, 'o:', label='The P97.5')

plt.semilogy(N_list[1:], erro_CD8_mean, 'D-', label='Tke média')
plt.semilogy(N_list[1:], erro_CD8_p2, 'D--', label='Tke P2.5')
plt.semilogy(N_list[1:], erro_CD8_p97, 'D:', label='Tke P97.5')

plt.semilogy(N_list[1:], erro_V_mean, '^-', label='V média')
plt.semilogy(N_list[1:], erro_V_p2, '^--', label='V P2.5')
plt.semilogy(N_list[1:], erro_V_p97, '^:', label='V P97.5')

plt.semilogy(N_list[1:], erro_A_mean, 's-', label='A média')
plt.semilogy(N_list[1:], erro_A_p2, 's--', label='A P2.5')
plt.semilogy(N_list[1:], erro_A_p97, 's:', label='A P97.5')

# plt.loglog(N_list[1:], erro_A_p97, 's:', label='A P97.5')
# reta de referência N^{-1/2}

ref = erro_CD4_mean[0] * np.sqrt(
    N_list[1] / np.array(N_list[1:])
)

# plt.loglog(
#     N_list[1:],
#     ref,
#     'k--',
#     linewidth=2,
#     label=r'$N^{-1/2}$'
# )

plt.semilogy(
    N_list[1:],
    ref,
    'k--',
    linewidth=2,
    label=r'$N^{-1/2}$'
)
#plt.xscale("log", base=2)
plt.yscale("log")
plt.grid(True, which="both")
plt.title('UQ Parâmetros ED (mais influentes para o caso)')
plt.xlabel("Número de amostras")
plt.ylabel("Erro relativo")
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")

plt.tight_layout()
plt.savefig('uqconv_influED.png', bbox_inches='tight')
plt.savefig('uqconv_influED.pdf', bbox_inches='tight')
plt.show()