import pandas as pd
import matplotlib.pyplot as plt

# carregar os dados do CSV
df = pd.read_csv("resultados.csv")

# threads únicas
threads = sorted(df["threads"].unique())
sizes = sorted(df["matrix_size"].unique())

# --- 1. Tempo de execução ---
plt.figure(figsize=(8,6))
for n in sizes:
    subset = df[df["matrix_size"] == n]
    plt.plot(subset["threads"], subset["time"], marker='o', label=f"N={n}")
plt.xlabel("Número de threads")
plt.ylabel("Tempo de execução (s)")
plt.title("Tempo de execução vs Threads")
plt.legend()
plt.grid(True)
plt.savefig("tempo_execucao.png")

# --- 2. Speedup ---
plt.figure(figsize=(8,6))
for n in sizes:
    subset = df[df["matrix_size"] == n].sort_values("threads")
    t_seq = subset[subset["threads"]==1]["time"].values[0]
    plt.plot(subset["threads"], t_seq / subset["time"], marker='o', label=f"N={n}")
plt.xlabel("Número de threads")
plt.ylabel("Speedup")
plt.title("Speedup vs Threads")
plt.legend()
plt.grid(True)
plt.savefig("speedup.png")

# --- 3. Eficiência ---
plt.figure(figsize=(8,6))
for n in sizes:
    subset = df[df["matrix_size"] == n].sort_values("threads")
    t_seq = subset[subset["threads"]==1]["time"].values[0]
    eff = (t_seq / subset["time"]) / subset["threads"]
    plt.plot(subset["threads"], eff, marker='o', label=f"N={n}")
plt.xlabel("Número de threads")
plt.ylabel("Eficiência")
plt.title("Eficiência vs Threads")
plt.legend()
plt.grid(True)
plt.savefig("eficiencia.png")

plt.show()