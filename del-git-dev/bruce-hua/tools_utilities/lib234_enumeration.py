#!/usr/bin/env python
# coding: utf-8

# In[29]:


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors3D
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import pandas
import csv
import time

start = time.time()
df = pandas.read_csv('D:/Dropbox/Harvard/PhD/DELs/Screening/metadata/encoding_metadata_v3.csv')
dfe = pandas.read_csv('D:/Dropbox/Harvard/PhD/DELs/Screening/metadata/lib234_enumerated_20210209.csv')
print("load time:", round(time.time() - start,2), "seconds")


# In[30]:


start = time.time()
rxn1 = rdChemReactions.ReactionFromSmarts('[CX3:1](=[OX1])NC>>[CX3:1](=[OX1])NC[3H]')
rxn2 = rdChemReactions.ReactionFromSmarts('[NX3;H1;!$(NC=O):1].[C,c:2][SX4:3](=[OX1])(=[OX1])(Cl)>>[NX3;H0;!$(NC=O):1][SX4:3](=[OX1])(=[OX1])[C,c:2]')
rxn2sub1 = rdChemReactions.ReactionFromSmarts('[NX3;H1;!$(NC=O):1].[C,c:2][CX3H1:3](=O)>>[NX3;H0;!$(NC=O):1][C:3][C,c:2]')
rxn2sub2 = rdChemReactions.ReactionFromSmarts('[NX3;H1;!$(NC=O):1].[C,c:2][CX3:3](=O)[OX2H1]>>[NX3;H0;!$(NC=O):1][C:3](=O)[C,c:2]')
rxn3 = rdChemReactions.ReactionFromSmarts('[c:1]I.[C,c:2][#5X3](O)(O)>>[c:1][C,c:2]')
structures = []
libid = [2,3,4]
counter = 0
for lib in libid:
    cy1 = df[(df['library_id'] == lib) & (df['cycle'] == 1)]
    cy2 = df[(df['library_id'] == lib) & (df['cycle'] == 2)]
    cy3 = df[(df['library_id'] == lib) & (df['cycle'] == 3)]
    for i in range(len(cy1)):
        r1 = cy1.iloc[i]['structure_id']
        products = rxn1.RunReactants((Chem.MolFromSmiles(r1),))
        p1 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
        for j in range(len(cy2)):
            r2 = cy2.iloc[j]['structure_id']
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
                r3 = cy3.iloc[k]['structure_id']
                if Chem.MolFromSmiles(r3) is None:
                    structures.append(p2)
                else:
                    products = rxn3.RunReactants((Chem.MolFromSmiles(p2),Chem.MolFromSmiles(r3)))
                    structures.append(Chem.MolToSmiles(products[0][0], isomericSmiles = True))
                counter += 1
                if (counter > 0) & (counter % 10000 == 0):
                    print(counter, "completed, approx.", 
                          round((len(dfe)-counter)*(time.time() - start)/counter), "seconds left")

dfe['structure'] = structures
dfe.to_csv('D:/Dropbox/Harvard/PhD/DELs/Screening/metadata/lib234_enumerated_w_structures_20210210.csv', index = False)
print("time to enumerate structures:", round(time.time() - start,2), "seconds")


# In[18]:


Chem.MolFromSmiles('Brc1nc(Cl)nc(NC)n1')


# In[25]:


rxn0 = rdChemReactions.ReactionFromSmarts('[c:1]I.[C,c:2][#5X3](O)(O)>>[c:1][C,c:2]')
products = rxn0.RunReactants((Chem.MolFromSmiles('O=C([C@@H]1NCCC[C@H]1C2=CC=CC(I)=C2)NC'),Chem.MolFromSmiles('B(O)(O)c1ccc(cc1)C2(CC2)C#N')))
p0 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
Chem.MolFromSmiles(p0)


# In[23]:


len(products)


# In[ ]:




