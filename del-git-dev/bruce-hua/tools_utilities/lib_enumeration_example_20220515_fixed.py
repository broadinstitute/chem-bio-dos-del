#!/usr/bin/env python
# coding: utf-8

# In[80]:


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors3D
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen
import pandas as pd
import csv
import freesasa
import numpy as np
import itertools

dir_name = 'D:/Dropbox/Harvard/PhD/DELs/Screening/prototype/'

bb_meta = pd.read_csv(dir_name + 'meta/lib/bb_meta.csv')
agg_meta = pd.read_csv(dir_name + 'meta/an/agg_meta.csv')


# In[84]:


lib_id = 2

rxn1 = rdChemReactions.ReactionFromSmarts('[CX3:1](=[OX1])NC>>[CX3:1](=[OX1])NC[3H]')
rxn2 = rdChemReactions.ReactionFromSmarts('[NX3;H1;!$(NC=O):1].[C,c:2][SX4:3](=[OX1])(=[OX1])(Cl)>>[NX3;H0;!$(NC=O):1][SX4:3](=[OX1])(=[OX1])[C,c:2]')
rxn2sub1 = rdChemReactions.ReactionFromSmarts('[NX3;H1;!$(NC=O):1].[C,c:2][CX3H1:3](=O)>>[NX3;H0;!$(NC=O):1][C:3][C,c:2]')
rxn2sub2 = rdChemReactions.ReactionFromSmarts('[NX3;H1;!$(NC=O):1].[C,c:2][CX3:3](=O)[OX2H1]>>[NX3;H0;!$(NC=O):1][C:3](=O)[C,c:2]')
rxn3 = rdChemReactions.ReactionFromSmarts('[c:1]I.[C,c:2][#5X3](O)(O)>>[c:1][C,c:2]')
structures = []
emw = []
fsp3 = []
nrb = []
slogp = []

cy1 = bb_meta[(bb_meta['lib_id'] == lib_id) & (bb_meta['cyc'] == 1)]
cy2 = bb_meta[(bb_meta['lib_id'] == lib_id) & (bb_meta['cyc'] == 2)]
cy3 = bb_meta[(bb_meta['lib_id'] == lib_id) & (bb_meta['cyc'] == 3)]

for i in range(len(cy1)):
    r1 = cy1.iloc[i]['bb_smiles']
    products = rxn1.RunReactants((Chem.MolFromSmiles(r1),))
    p1 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
    for j in range(len(cy2)):
        r2 = cy2.iloc[j]['bb_smiles']
        if Chem.MolFromSmiles(r2) is None:
            p2 = p1
        else:
            products = rxn2.RunReactants((Chem.MolFromSmiles(p1),Chem.MolFromSmiles(r2)))
            if len(products) == 0:
                products = rxn2sub1.RunReactants((Chem.MolFromSmiles(p1),Chem.MolFromSmiles(r2)))
                if len(products) == 0:
                    products = rxn2sub2.RunReactants((Chem.MolFromSmiles(p1),Chem.MolFromSmiles(r2)))
            p2 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
        for k in range(len(cy3)):
            r3 = cy3.iloc[k]['bb_smiles']
            if Chem.MolFromSmiles(r3) is None:
                struct_out = p2
            else:
                products = rxn3.RunReactants((Chem.MolFromSmiles(p2),Chem.MolFromSmiles(r3)))
                struct_out = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
            structures.append(struct_out)
            cpd_out = Chem.MolFromSmiles(struct_out)
            try:
                emw.append(Descriptors.ExactMolWt(cpd_out))
            except:
                emw.append(np.nan)
            try:
                fsp3.append(Lipinski.FractionCSP3(cpd_out))
            except:
                fsp3.append(np.nan)
            try:
                nrb.append(Lipinski.NumRotatableBonds(cpd_out))
            except:
                nrb.append(np.nan)
            try:
                slogp.append(Crippen.MolLogP(cpd_out))
            except:
                slogp.append(np.nan)
            
            
pd.DataFrame({'structure': structures}).to_csv(dir_name + 'meta/lib/enum_struct/lib' + str(lib_id).zfill(3) + '.csv', index = False)
pd.DataFrame({'value': emw}).to_csv(dir_name + 'meta/lib/enum_prop/lib' + str(lib_id).zfill(3) + '_emw.csv', index = False)
pd.DataFrame({'value': fsp3}).to_csv(dir_name + 'meta/lib/enum_prop/lib' + str(lib_id).zfill(3) + '_fsp3.csv', index = False)
pd.DataFrame({'value': nrb}).to_csv(dir_name + 'meta/lib/enum_prop/lib' + str(lib_id).zfill(3) + '_nrb.csv', index = False)
pd.DataFrame({'value': slogp}).to_csv(dir_name + 'meta/lib/enum_prop/lib' + str(lib_id).zfill(3) + '_slogp.csv', index = False)

# Generate meta cols
tags = [cy1['tag_id'], cy2['tag_id'], cy3['tag_id']]
all_combinations = list(itertools.product(*tags))
enum_df = pd.concat([
    pd.DataFrame({'lib_id': [lib_id for x in all_combinations]}),
    pd.DataFrame(all_combinations, columns=['cycle1', 'cycle2', 'cycle3'])
], axis = 1)

for i in range(len(agg_meta)):
    agg_id = agg_meta.iloc[i]['agg_id']
    col_names = ['cycle' + x for x in agg_meta.iloc[i]['cyc'].split(', ')]
    col_names.insert(0, 'lib_id')
    
    enum_out = enum_df[col_names].drop_duplicates()
    enum_out.to_csv(dir_name + 'meta/lib/meta_cols/lib' + str(lib_id).zfill(3) + '_' +
                    'agg' + str(agg_id).zfill(2) + '_meta_cols.csv', index = False)


