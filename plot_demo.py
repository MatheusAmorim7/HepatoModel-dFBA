import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--output',type=Path,default=ROOT/'results/lobule.png')
p.add_argument('--result',type=Path)
args=p.parse_args()
rows=[line.split() for line in (ROOT/'data/hosts.txt').read_text().splitlines()]
assert len(rows)==3000
xyz=np.array([[float(v) for v in r[2:5]] for r in rows])
radii=np.cbrt(3*np.array([float(r[5]) for r in rows])/(4*np.pi))
zones=np.array([int(r[6].split('_')[-1]) for r in rows])
fig,ax=plt.subplots(1,2,figsize=(12,5),layout='constrained')
colors=['#258a75','#dbad38','#c45c85']
for zone,color in enumerate(colors,1):
    subset=np.flatnonzero(zones==zone)
    patches=[Circle(xyz[i,:2],radii[i]) for i in subset]
    ax[0].add_collection(PatchCollection(patches,facecolor=color,edgecolor='white',linewidth=.15,alpha=.8))
    ax[0].scatter([],[],color=color,label=f'Zone {zone} label')
ax[0].legend(loc='upper right',fontsize=8)
ax[0].set_title('3,000 agents: reference placement')
lines=(ROOT/'data/oxygen.txt').read_text().splitlines()
water=np.array([[float(x) for x in line.split()] for line in lines[49:2549]])
volume=water[:,0]
oxygen=water[:,2]
title='Initial oxygen (mM)'
if args.result:
    a=np.fromfile(args.result,dtype=np.float64)
    assert len(a)==589600
    oxygen=a[:50000].reshape(2500,20)[:,0]
    title='Oxygen after one 0.0075 s step (mM)'
im=ax[1].imshow((oxygen/volume).reshape(50,50),extent=[-500,500,-500,500],origin='lower',cmap='viridis')
ax[1].scatter(xyz[:,0],xyz[:,1],s=.6,c='black',alpha=.15)
fig.colorbar(im,ax=ax[1],label='Oxygen (mM)')
ax[1].set_title(title)
for a in ax:
    a.set(xlim=(-500,500),ylim=(-500,500),xlabel='x (um)',ylabel='y (um)',aspect='equal')
fig.suptitle('Exploratory lobule demo: labels do not establish biological zonation',fontsize=11)
args.output.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(args.output,dpi=160)
print(args.output)
