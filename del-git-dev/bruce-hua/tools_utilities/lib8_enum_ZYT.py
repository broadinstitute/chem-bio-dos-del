from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors3D
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import SaltRemover
import pandas
import csv
import time


df = pandas.read_csv('Lib ID 8 (DOS-DEL-5) BB Structure.csv')
dfe = pandas.read_csv('lib8_enum.csv')

#reductive aminations with all sorts of requirements on amine
rxn_redam = rdChemReactions.ReactionFromSmarts('[O:4]=[CH:2][c:3].[#7;H1,H2!$(NS);!$(NC=O);!$(NC=C);!$(NC=N):1]>>[c:3][CH2:2][N:1]')
rxn_cuaac = rdChemReactions.ReactionFromSmarts('[#6:3][C:1]#[CH:2].[C:4][N:5]=[N:6]=[N:7]>>[C:4][N:5]1[C:2]=[C:1]([#6:3])[N+0:7]=[N+0:6]1')


remover = Chem.SaltRemover.SaltRemover(defnData="[Cl,Br]")

structures = []
libid = [8]
counter = 0

for lib in libid:
    cy1 = df[(df['master_lib_id'] == lib) & (df['cy_no'] == 1)]
    cy2 = df[(df['master_lib_id'] == lib) & (df['cy_no'] == 2)]
    cy3 = df[(df['master_lib_id'] == lib) & (df['cy_no'] == 3)]
    for i in range(len(cy1)):
        r1 = cy1.iloc[i]['bb_smiles']
        
        for j in range(len(cy2)):
            r2 = cy2.iloc[j]['bb_smiles']
            r2 = Chem.MolToSmiles(remover.StripMol(Chem.MolFromSmiles(r2)))
            products = rxn_cuaac.RunReactants((Chem.MolFromSmiles(r2),Chem.MolFromSmiles(r1)))
            p2 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
            
            for k in range(len(cy3)):
                r3 = cy3.iloc[k]['bb_smiles']
                r3 = Chem.MolToSmiles(remover.StripMol(Chem.MolFromSmiles(r3)))
                products = rxn_redam.RunReactants((Chem.MolFromSmiles(p2),Chem.MolFromSmiles(r3)))
                structures.append(Chem.MolToSmiles(products[0][0], isomericSmiles = True))
                
                counter += 1
                

dfe['structures'] = structures
dfe.to_csv('lib8_enum.csv', index = False)
len(structures)