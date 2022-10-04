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


#metadata here
df = pandas.read_csv('FKBP-CIP-DEL-BB structures.csv')
dfe = pandas.read_csv('lib300_enumerated.csv')

#define reactions
##adding FKBP part
rxn_amide = rdChemReactions.ReactionFromSmarts('[#7;A;X3;H1,H2!$(NS);!$(NC=O):1][2H]>>[3H]COC1=CC=C(CC[C@@H](NC(=O)C2CCCCN2C(=O)C(=O)C(C)(C)CC)C([N:1])=O)C=C1')
##adding FKBP part for shortest amine (Cy1 = 44)
rxn_amine2 = rdChemReactions.ReactionFromSmarts('[#7;!$(NC=O):1][2H]>>[3H]COC1=CC=C(CC[C@@H](NC(=O)C2CCCCN2C(=O)C(=O)C(C)(C)CC)C[N:1])C=C1')
##adding triazine part
rxn_triazine = rdChemReactions.ReactionFromSmarts('[N:1]S>>[N:1]C1=NC(I)=NC([At])=N1')
##adding triazine part for shortest amine (Cy1 = 44)
rxn_triazine2 = rdChemReactions.ReactionFromSmarts('[#7;A;X3;H1,H2!$(NS);!$(NC=O):1]>>[N:1]C1=NC(I)=NC([At])=N1')
##Cy2 sub
rxn_sub1 = rdChemReactions.ReactionFromSmarts('[#6:2]I.[#7;A;X3;H1,H2!$(NS);!$(NC=O):1]>>[#6:2]-[#7:1]')
##Cy3 sub
rxn_sub2 = rdChemReactions.ReactionFromSmarts('[#6:2][At].[#7;A;X3;H1,H2!$(NS);!$(NC=O):1]>>[#6:2]-[#7:1]')
##Cy3 Chlorine only
rxn_sub0 = rdChemReactions.ReactionFromSmarts('[#6:1][At]>>[#6:1][Cl]')
##desalt building blocks
remover = Chem.SaltRemover.SaltRemover(defnData="[Cl,Br]")


structures = []
libid = [300]
counter = 0



for lib in libid:
    cy1 = df[(df['master_lib_id'] == lib) & (df['cy_no'] == 1)]
    cy2 = df[(df['master_lib_id'] == lib) & (df['cy_no'] == 2)]
    cy3 = df[(df['master_lib_id'] == lib) & (df['cy_no'] == 3)]
    
    for i in range(len(cy1)):
        r1 = cy1.iloc[i]['bb_smiles']
        if rxn_amide.RunReactants((Chem.MolFromSmiles(r1),)) == ():
            #for shortest connector
            products = rxn_amine2.RunReactants((Chem.MolFromSmiles(r1),))
            p1 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)    
            products = rxn_triazine2.RunReactants((Chem.MolFromSmiles(p1),))
            p1 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
        else:
            products = rxn_amide.RunReactants((Chem.MolFromSmiles(r1),))
            p1 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
            products = rxn_triazine.RunReactants((Chem.MolFromSmiles(p1),))
            p1 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
        for j in range(50):
            r2 = cy2.iloc[j]['bb_smiles']
            r2 = Chem.MolToSmiles(remover.StripMol(Chem.MolFromSmiles(r2)))
            products = rxn_sub1.RunReactants((Chem.MolFromSmiles(p1),Chem.MolFromSmiles(r2)))
            p2 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)
            for k in range(50):
                r3 = cy3.iloc[k]['bb_smiles']
                r3 = Chem.MolToSmiles(remover.StripMol(Chem.MolFromSmiles(r3)))
                if rxn_sub2.RunReactants((Chem.MolFromSmiles(p2),Chem.MolFromSmiles(r3))) == ():
                    #for Cl
                    products = rxn_sub0.RunReactants((Chem.MolFromSmiles(p2),))
                    p3 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)    

                else:
                    products = rxn_sub2.RunReactants((Chem.MolFromSmiles(p2),Chem.MolFromSmiles(r3)))
                    p3 = Chem.MolToSmiles(products[0][0], isomericSmiles = True)

                structures.append(Chem.MolToSmiles(products[0][0], isomericSmiles = True))
                counter += 1 
                
dfe['structures'] = structures
dfe.to_csv('lib300_enumerated.csv', index = False)
