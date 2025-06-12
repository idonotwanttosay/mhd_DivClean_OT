"""Check mass and energy conservation from output CSV files."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt, glob, re

rho_files = glob.glob("Result/out_rho_*.csv")
if not rho_files:
    raise SystemExit("No output found")
steps = sorted(int(re.findall(r"_rho_(\d+).csv",f)[0]) for f in rho_files)

# grid info
sample = pd.read_csv(rho_files[0], header=None)
xs = np.unique(sample[0]); ys = np.unique(sample[1])
DX = xs[1]-xs[0]; DY = ys[1]-ys[0]
NX, NY = len(xs), len(ys)

def load(step, prefix):
    df = pd.read_csv(f"Result/out_{prefix}_{step}.csv", header=None)
    return df[2].values.reshape(NX, NY).T

mass=[]; energy=[]
for s in steps:
    rho = load(s, "rho")
    e   = load(s, "e")
    mass.append(rho.sum()*DX*DY)
    energy.append(e.sum()*DX*DY)

plt.plot(steps, mass, label='mass')
plt.plot(steps, energy, label='energy')
plt.xlabel('step'); plt.ylabel('total')
plt.title('Conservation check')
plt.legend()
plt.tight_layout()
plt.savefig('Result/conservation.png', dpi=200)
print('Saved Result/conservation.png')
