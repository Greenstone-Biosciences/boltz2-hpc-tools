#!/usr/bin/env python3


import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import QED
from rdkit.Chem.QED import qed
import sys, os

smiles = str(sys.argv[1])
name = sys.argv[2]
outdir = sys.argv[3]

print(smiles)
print(name)
print(outdir + "/QED_adjust")
df = open(str(outdir) + "/QED_Adjust", "a")

mol = Chem.MolFromSmiles(str(smiles))
mol = Chem.AddHs(mol)
mol.SetProp("_Name", name)
cid = AllChem.EmbedMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol)
Chem.rdmolfiles.MolToMolFile(mol, os.path.join(outdir, '{n}.mol'.format(n=name)))
df.write(str(smiles) + "\t" + str(name) + "\t" + str(rdkit.Chem.QED.weights_mean(mol)) + "\n")
df.close()