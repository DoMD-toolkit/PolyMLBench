import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdd

data = pd.read_csv('data.csv')

# rdkit descriptors
# chin

chi0n = []
chi1n = []
chi2n = []
chi2v = []
chi3n = []
chi3v = []
chi4n = []
chi4v = []
nrings = []
rtbonds = []
nnac = []
nhc = []
nab = []
nac = []
nahc = []
nhba = []
nhbd = []
nsc = []
nshc = []

bcut2d = []
for smi in data['smiles']:
    smi = smi.replace('[*]', '[H]')
    smi = smi.replace('*', '[H]')
    mol = Chem.MolFromSmiles(smi)
    chi0n_ = rdd.CalcChi0n(mol)
    chi0n.append(chi0n_)
    chi1n_ = rdd.CalcChi1n(mol)
    chi1n.append(chi1n_)
    chi2n_ = rdd.CalcChi2n(mol)
    chi2n.append(chi2n_)
    chi3n_ = rdd.CalcChi3n(mol)
    chi3n.append(chi3n_)
    chi4n_ = rdd.CalcChi4n(mol)
    chi4n.append(chi4n_)
    chi2v_ = rdd.CalcChi2v(mol)
    chi2v.append(chi2v_)
    chi3v_ = rdd.CalcChi3v(mol)
    chi3v.append(chi3v_)
    chi4v_ = rdd.CalcChi4v(mol)
    chi4v.append(chi4v_)
    nrings_ = rdd.CalcNumRings(mol)
    nrings.append(nrings_)
    rtbonds_ = rdd.CalcNumRotatableBonds(mol)
    rtbonds.append(rtbonds_)
    nnac_ = rdd.CalcNumAliphaticCarbocycles(mol)
    nnac.append(nnac_)
    nhc_ = rdd.CalcNumAliphaticHeterocycles(mol)
    nhc.append(nhc_)
    nab_ = rdd.CalcNumAmideBonds(mol)
    nab.append(nab_)
    nac_ = rdd.CalcNumAromaticCarbocycles(mol)
    nac.append(nac_)
    nahc_ = rdd.CalcNumAromaticHeterocycles(mol)
    nahc.append(nahc_)
    nhba_ = rdd.CalcNumHBA(mol)
    nhba.append(nhba_)
    nhbd_ = rdd.CalcNumHBD(mol)
    nhbd.append(nhbd_)
    nsc_ = rdd.CalcNumSaturatedCarbocycles(mol)
    nsc.append(nsc_)
    nshc_ = rdd.CalcNumSaturatedHeterocycles(mol)
    nshc.append(nshc_)
    try:
        bcut2d_ = rdd.BCUT2D(mol)
        bcut2d.append(bcut2d_)
    except ValueError:
        bcut2d.append([np.nan, ] * 8)
bcut2d = np.array(bcut2d)
for i in range(8):
    data[f'bcut{i + 1}'] = bcut2d.T[i]
data['chi0n'] = chi0n
data['chi1n'] = chi1n
data['chi2n'] = chi2n
data['chi2v'] = chi2v
data['chi3n'] = chi3n
data['chi3v'] = chi3v
data['chi4n'] = chi4n
data['chi4v'] = chi4v
data['NumRings'] = nrings
data['NumRotatableBonds'] = rtbonds
data['NumAliphaticCarbocycles'] = nnac
data['NumAliphaticHeterocycles'] = nhc
data['NumAromaticCarbocycles'] = nac
data['NumAromaticHeterocycles'] = nahc
data['NumHBondAcc'] = nhba
data['NumHBondDor'] = nhbd
data['NumAmideBonds'] = nab
data['NumSaturatedCarbocycles'] = nsc
data['NumSaturatedHeterocycles'] = nshc

data.to_csv('data_with_rdd.csv', index=False)
