import numpy as np
import pandas as pd
import glob
import re
import matplotlib.pyplot as plt

# discover available time steps
rho_files = glob.glob("Result/out_rho_*.csv")
if not rho_files:
    raise SystemExit("No output found")
steps = sorted(int(re.findall(r"_rho_(\d+)", f)[0]) for f in rho_files)

# infer grid from first file
sample = pd.read_csv(rho_files[0], header=None)
xs = np.unique(sample[0]); ys = np.unique(sample[1])
DX = xs[1] - xs[0]
DY = ys[1] - ys[0]
NX, NY = len(xs), len(ys)


def load(step, prefix):
    df = pd.read_csv(f"Result/out_{prefix}_{step}.csv", header=None)
    return df[2].values.reshape(NX, NY).T

mass = []
energy = []
div_l2 = []
for s in steps:
    rho = load(s, "rho")
    e = load(s, "e")
    mass.append(rho.sum() * DX * DY)
    energy.append(e.sum() * DX * DY)

    bx = load(s, "bx")
    by = load(s, "by")
    div = (np.roll(bx, -1, 1) - np.roll(bx, 1, 1)) / (2 * DX) + \
          (np.roll(by, -1, 0) - np.roll(by, 1, 0)) / (2 * DY)
    div_l2.append(np.sqrt(np.mean(div ** 2)))

print(f"Mass change: {mass[0]:.6g} -> {mass[-1]:.6g}")
print(f"Energy change: {energy[0]:.6g} -> {energy[-1]:.6g}")
print(f"Final L2(divB): {div_l2[-1]:.3e}")

plt.figure(figsize=(6,4))
plt.plot(steps, mass, label='mass')
plt.plot(steps, energy, label='energy')
plt.xlabel('step'); plt.ylabel('total')
plt.legend()
plt.tight_layout()
plt.savefig('Result/summary_mass_energy.png', dpi=200)

plt.figure(figsize=(6,4))
plt.semilogy(steps, div_l2)
plt.xlabel('step'); plt.ylabel('L2(divB)')
plt.tight_layout()
plt.savefig('Result/summary_divB.png', dpi=200)
print('Saved summary_mass_energy.png and summary_divB.png')
