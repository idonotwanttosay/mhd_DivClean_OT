import pandas as pd, numpy as np, matplotlib.pyplot as plt, glob, re, os
files = glob.glob("Result/out_rho_*.csv")
steps = sorted(int(re.findall(r"_rho_(\d+).csv",f)[0]) for f in files)
if not steps:
    raise SystemExit("no rho")
xs = np.unique(pd.read_csv(files[0],header=None)[0])
ys = np.unique(pd.read_csv(files[0],header=None)[1])
nx,ny = len(xs), len(ys)
for s in steps:
    df = pd.read_csv(f"Result/out_rho_{s}.csv", header=None)
    rho = df[2].values.reshape(nx, ny).T
    X, Y = np.meshgrid(xs, ys)
    plt.figure(figsize=(6,5))
    plt.contourf(X, Y, rho, levels=40, cmap='viridis')
    plt.colorbar()
    plt.title(f"rho step {s}")
    plt.gca().set_aspect('equal')
    out = f"Result/rho_{s}.png"
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved {out}")
