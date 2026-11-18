import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('result',type=Path)
args=parser.parse_args()
core='{http://www.sbml.org/sbml/level3/version1/core}'
model=ET.parse(ROOT/'data/host.xml').getroot().find(core+'model')
metabolites=list(model.find(core+'listOfSpecies'))
reactions=list(model.find(core+'listOfReactions'))
mi={m.attrib['id']:i for i,m in enumerate(metabolites)}
ri={r.attrib['id']:i for i,r in enumerate(reactions)}
S=np.zeros((len(mi),len(ri)))
for j,r in enumerate(reactions):
    for tag,sign in [('listOfReactants',-1),('listOfProducts',1)]:
        node=r.find(core+tag)
        if node is not None:
            for term in node:
                S[mi[term.attrib['species']],j]+=sign*float(term.attrib.get('stoichiometry',1))
assert S.shape==(175,179)
fields=list(ET.parse(ROOT/'data/fields.xml').getroot().find('finite_resource_control'))
cols=np.array([ri['R_'+f.attrib['reaction']] for f in fields])
atoms=np.array([[float(f.attrib['C']),float(f.attrib['N'])] for f in fields])
vox=np.array([int(r.split()[1]) for r in (ROOT/'data/hosts.txt').read_text().splitlines()])
with np.load(ROOT/'data/resources.npz') as z:
    q0=z['q'];bounds=z['bounds']
a=np.fromfile(args.result,dtype=np.float64)
assert len(a)==589600 and np.isfinite(a).all()
q=a[:50000].reshape(2500,20)
flux=a[50000:587000].reshape(3000,179)
sv=float(abs(S@flux.T).max())
bound=max(0.,float((bounds[:,:,0]-flux).max()),float((flux-bounds[:,:,1]).max()))
delta=flux[:,cols]*(.12*.0075/3600)
changes=np.zeros_like(q)
gross=np.zeros_like(q)
np.add.at(changes,vox,delta)
np.add.at(gross,vox,np.maximum(0.,-delta))
available=q0.copy()
available[:,0]=q[:,0]-changes[:,0]
excess=float((gross-available+np.minimum(1e-12,available*.5)).max())
update=float(abs(q[:,1:]-q0[:,1:]-changes[:,1:]).max())
cn=changes.sum(0)@atoms
assert sv<=2e-7 and bound<=2e-7 and excess<=1e-12 and update<=1e-12
assert abs(cn).max()<=1e-9 and (q>=0).all()
summary=dict(hosts=3000,reactions=179,metabolites=175,model_seconds=.0075,
             Sv_max=sv,bound_violation=bound,gross_stock_excess_pmol=excess,
             nonoxygen_update_error_pmol=update,C_error_pmol_atoms=float(cn[0]),N_error_pmol_atoms=float(cn[1]),
             biological_validation=False)
print(json.dumps(summary,indent=2))
