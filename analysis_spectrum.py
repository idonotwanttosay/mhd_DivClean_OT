"""Compute kinetic energy spectra eliminating mean flow and normalising.

Generates Result/energy_spectrum.png
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt, glob, re, os

u_files = glob.glob("Result/out_u_*.csv")
if not u_files:
    raise SystemExit("No velocity output found.")

steps = sorted(int(re.findall(r"_u_(\d+).csv",f)[0]) for f in u_files)

def grid():
    df = pd.read_csv(u_files[0],header=None)
    xs=np.unique(df[0]); ys=np.unique(df[1])
    return len(xs), len(ys), xs, ys
nx,ny,xs,ys = grid()

def load(step,comp):
    df = pd.read_csv(f"Result/out_{comp}_{step}.csv",header=None)
    return df.values[:,2].reshape(nx,ny).T

def spectrum(u,v):
    u = u - u.mean(); v = v - v.mean()
    uhat = np.fft.fft2(u); vhat = np.fft.fft2(v)
    E = 0.5*(np.abs(uhat)**2 + np.abs(vhat)**2)/(nx*ny)**2
    kx = np.fft.fftfreq(nx)*(nx)
    ky = np.fft.fftfreq(ny)*(ny)
    KX,KY = np.meshgrid(kx,ky)
    kmag = np.sqrt(KX**2+KY**2)
    k_bins = np.arange(0.5, min(nx,ny)//2 +1 ,1)
    Ek=np.zeros_like(k_bins)
    for i in range(len(k_bins)):
        mask = (kmag>=k_bins[i]-0.5)&(kmag<k_bins[i]+0.5)
        Ek[i]=E[mask].sum()
    return k_bins,Ek

plt.figure(figsize=(8,4))
for s in steps:
    u=load(s,"u"); v=load(s,"v")
    k,Ek = spectrum(u,v)
    plt.loglog(k,Ek,label=f"s{s}")
plt.loglog([1,10],[1e-2,1e-2*10**(-5/3)],'k--',label='k^-5/3')
plt.xlabel("k"); plt.ylabel("E(k)")
plt.title("Kinetic Energy Spectrum")
plt.legend(fontsize=7,ncol=2)
plt.tight_layout()
plt.savefig("Result/energy_spectrum.png",dpi=200)
print("Saved Result/energy_spectrum.png")
