import rdkit.Chem as Chem
from rdkit.Chem import Draw
import rdkit.Chem.AllChem as AllChem

import numpy as np
import pandas as pd
import os, sys

file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 
    'encoding_metadata.csv')

df_tags = pd.read_csv(file)
df_tags = df_tags.loc[df_tags['library_id'].isin([2,4,5,6,7])]

# DD1S - need to replace structure_id names with SMILES
df_tags['structure_id'].replace(inplace=True, to_replace={
    '2R,3S-CisAzet':'IC1=CC=C([C@@H](C2)[C@H](C(N[*])=O)N2C(OC(C)(C)C)=O)C=C1',
    '2S,3R-CisAzet':'IC1=CC=C([C@H](C2)[C@@H](C(N[*])=O)N2C(OC(C)(C)C)=O)C=C1',
    '2R,3R-TransAzet':'IC1=CC=C([C@H](C2)[C@H](C(N[*])=O)N2C(OC(C)(C)C)=O)C=C1',
    '2S,3S-TransAzet':'IC1=CC=C([C@@H](C2)[C@@H](C(N[*])=O)N2C(OC(C)(C)C)=O)C=C1',
    '2R,3R-CisPyrr':'IC1=CC=C([C@@H](CC2)[C@H](C(N[*])=O)N2C(OC(C)(C)C)=O)C=C1',
    '2S,3S-CisPyrr':'IC1=CC=C([C@H](CC2)[C@@H](C(N[*])=O)N2C(OC(C)(C)C)=O)C=C1',
    '2S,3R-TransPyrr':'IC1=CC=C([C@@H](CC2)[C@@H](C(N[*])=O)N2C(OC(C)(C)C)=O)C=C1',
    '2R,3S-TransPyrr':'IC1=CC=C([C@H](CC2)[C@H](C(N[*])=O)N2C(OC(C)(C)C)=O)C=C1',
})

tag_to_struct = {
    (row['library_id'], row['cycle'], row['cycle_id']): row['structure_id'] \
    for i,row in df_tags.iterrows()
}

def canon_smi(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

def run_reaction(smi, smarts, print_ambiguous=False):
    ambig = False
    rxn = AllChem.ReactionFromSmarts(smarts)
    prods = rxn.RunReactants([Chem.MolFromSmiles(smi_frag) for smi_frag in smi])
    prod_smi = set([canon_smi(Chem.MolToSmiles(x[0])) for x in prods])
    if len(prod_smi) > 1 and print_ambiguous:
        ambig = True
        print('Ambiguous result!')
        print(smi)
        print(smarts)
        print(prod_smi)
    return prod_smi.pop(), ambig

SULFONYL_CHLORIDE = Chem.MolFromSmarts('S(=O)(=O)Cl')
ALDEHYDE = Chem.MolFromSmarts('[*][C;H1]=O')
PRIMARY_AMINE = Chem.MolFromSmarts('[NH2][*]')

def run_cycle1(library_id, cycle1_smi):
    a1 = False
    
    if library_id in [2, 3]:
        smi = cycle1_smi

    elif library_id == 4:
        smi, a1 = run_reaction((cycle1_smi,), '[C:1][CH0](=O)[NH1][CH3]>>[*:1][CH0](=O)[NH1][*]')

    elif library_id == 5:
        smi, a1 = run_reaction((cycle1_smi,), '[N:1]C(=O)OC(C)(C)C>>[*:1][*]')

    elif library_id == 6:
        smi, a1 = run_reaction((cycle1_smi,), '[NH2:1][C:2][C:3]>>[*][CH2][NH1:1][*:2][*:3]')

    elif library_id == 7: # triazine - use nonsense elements in place of chlorides
        if cycle1_smi != 'blank':
            smi, a1a = run_reaction((cycle1_smi,), '[N:1](C(=O)OCC1c2ccccc2-c2ccccc21)[C:2]>>[Pt]C1=NC([*:1][*:2])=NC([Pd])=N1')
            smi, a1b = run_reaction((smi,), '[#6:2]C(=O)[O]>>[*:2]C(=O)[NH][*]')
            a1 = a1a or a1b
        else:
            smi = '[Pt]C1=NC([NH][*])=NC([Pd])=N1'

    else:
        raise ValueError('Unknown library id {}'.format(library_id))
        
    return smi, a1
        

def run_cycle2(library_id, cycle2_smi, smi):
    a2 = False
    smi_bkp = smi
    
    if library_id in [2, 3]:
        if Chem.MolFromSmiles(cycle2_smi).HasSubstructMatch(SULFONYL_CHLORIDE): # sulfonylation
            smi, a2 = run_reaction((smi, cycle2_smi), '[N:1]C(=O)OC(C)(C)C.[*:2]S(=O)(=O)Cl>>[N:1]S(=O)(=O)[*:2]')
        else: # aldehyde for reductive amination
            smi, a2 = run_reaction((smi, cycle2_smi), '[N:1]C(=O)OC(C)(C)C.[*:2][C;H1:3]=O>>[N:1][CH2:3][*:2]')

    elif library_id == 4:
        if cycle2_smi != 'BLANK':
            if Chem.MolFromSmiles(cycle2_smi).HasSubstructMatch(SULFONYL_CHLORIDE): # sulfonylation
                smi, a2 = run_reaction((smi, cycle2_smi), '[C:1][NH:2][C:3].[*:4]S(=O)(=O)Cl>>[*:1][N:2]([*:3])S(=O)(=O)[*:4]')
            elif Chem.MolFromSmiles(cycle2_smi).HasSubstructMatch(ALDEHYDE): # aldehyde for reductive amination
                smi, a2 = run_reaction((smi, cycle2_smi), '[C:1][NH:2][C:3].[*:4][C;H1:5]=O>>[*:1][N:2]([*:3])[CH2:5][*:4]')
            else: # acid
                smi, a2 = run_reaction((smi, cycle2_smi), '[C:1][NH:2][C:3].[*:4][C:5](=O)[OH]>>[*:1][N:2]([*:3])[C:5](=O)[*:4]')

    elif library_id == 5:
        if cycle2_smi != 'blank':
            smi, a2 = run_reaction((smi, cycle2_smi), '[c:1]I.[*:2]B(O)O>>[*:1][*:2]')

    elif library_id == 6:
        smi, a2a = run_reaction((smi, cycle2_smi), '[C:1][NH:2][C:3].[*:4]S(=O)(=O)Cl>>[*:1][NH0:2]([*:3])S(=O)(=O)[*:4]')
        smi, a2b = run_reaction((smi,), '([c;$(ccS):1]F.[OH:2][C:3])>>[*:1][O:2][*:3]')
        a2 = a2a or a2b
        
    elif library_id == 7: # triazine - use nonsense elements in place of chlorides
        if Chem.MolFromSmiles(cycle2_smi).HasSubstructMatch(PRIMARY_AMINE):
            smi, a2 = run_reaction((smi, cycle2_smi), '[n:1][c:2][Pt].[NH2:4][*;!$(C=O);!$(C=N):6]>>[*:1][c:2][NH1:4][*:6]')
            if a2: # ambiguous reactivity, force aliphatic amine
                smi, a2 = run_reaction((smi_bkp, cycle2_smi), '[n:1][c:2][Pt].[NH2:4][C;!$(C=O);!$(C=N):6]>>[*:1][c:2][NH1:4][*:6]')
        else: # secondary
            smi, a2 = run_reaction((smi, cycle2_smi), '[n:1][c:2][Pt].[NH1:4]([*;!S;!$(C=O);!$(C=N):5])[*;!S;!$(C=O);!$(C=N):6]>>[*:1][c:2][NH0:4]([*:5])[*:6]')
            if a2: # favor unsubstituted neighbors
                smi, a2 = run_reaction((smi_bkp, cycle2_smi), '[n:1][c:2][Pt].[NH1:4]([CH2:5])[CH2:6]>>[*:1][c:2][NH0:4]([*:5])[*:6]')
                  
    else:
        raise ValueError('Unknown library id {}'.format(library_id))
        
    return smi, a2
        
        
def run_cycle3(library_id, cycle3_smi, smi):
    a3 = False
    smi_bkp = smi
    
    if library_id in [2, 3]:
        if cycle3_smi != 'blank':
            smi, a3 = run_reaction((smi, cycle3_smi), '[c:1]I.[*:2]B(O)O>>[*:1][*:2]')
    
    elif library_id == 4:
        if cycle3_smi != 'blank':
            smi, a3 = run_reaction((smi, cycle3_smi), '[c:1]I.[*:2]B(O)O>>[*:1][*:2]')
            
    elif library_id == 5:
        if cycle3_smi != 'blank':
            if Chem.MolFromSmiles(cycle3_smi).HasSubstructMatch(PRIMARY_AMINE):
                smi, a3 = run_reaction((smi, cycle3_smi), '[N:5][C:6][C:1](=O)[OH].[NH2:2][*;!$(C=O);!$(C=N):4]>>[*:5][*:6][*:1](=O)[NH1:2][*:4]')
                if a3: # ambiguous reactivity, force aliphatic amine
                    smi, a3 = run_reaction((smi_bkp, cycle3_smi), '[N:5][C:6][C:1](=O)[OH].[NH2:2][C;!$(C=O);!$(C=N):4]>>[*:5][*:6][*:1](=O)[NH1:2][*:4]')
            else: # secondary
                smi, a3 = run_reaction((smi, cycle3_smi), '[N:5][C:6][C:1](=O)[OH].[NH1:2]([*;!S;!$(C=O);!$(C=N):3])[*;!S;!$(C=O);!$(C=N):4]>>[*:5][*:6][*:1](=O)[NH0:2]([*:3])[*:4]')
                if a3: # favor unsubstituted neighbors
                    smi, a3 = run_reaction((smi_bkp, cycle3_smi), '[N:5][C:6][C:1](=O)[OH].[NH1:2]([CH2:3])[CH2:4]>>[*:5][*:6][*:1](=O)[NH0:2]([*:3])[*:4]')
                    
    elif library_id == 6:
        if cycle3_smi != 'blank':
            smi, a3 = run_reaction((smi, cycle3_smi), '[c:1]Br.[*:2]B(O)O>>[*:1][*:2]')
        
    elif library_id == 7: # triazine - use nonsense elements in place of chlorides
        if Chem.MolFromSmiles(cycle3_smi).HasSubstructMatch(PRIMARY_AMINE):
            smi, a3 = run_reaction((smi, cycle3_smi), '[n:1][c:2][Pd].[NH2:4][*:6]>>[*:1][c:2][NH1:4][*:6]')
            if a3: # ambiguous reactivity, force aliphatic amine
                smi, a3 = run_reaction((smi_bkp, cycle3_smi), '[n:1][c:2][Pd].[NH2:4][C;!$(C=O);!$(C=N):6]>>[*:1][c:2][NH1:4][*:6]')
        else: # secondary
            smi, a3 = run_reaction((smi, cycle3_smi), '[n:1][c:2][Pd].[NH1:4]([*;!S;!$(C=O);!$(C=N):5])[*;!S;!$(C=O):6]>>[*:1][c:2][NH0:4]([*:5])[*:6]')
            if a3: # favor unsubstituted neighbors
                smi, a3 = run_reaction((smi_bkp, cycle3_smi), '[n:1][c:2][Pd].[NH1:4]([CH2:5])[CH2:6]>>[*:1][c:2][NH0:4]([*:5])[*:6]')

    else:
        raise ValueError('Unknown library id {}'.format(library_id))
        
    return smi, a3


def cpd_def_to_smi(library_id, cycle1, cycle2, cycle3):
    cycle1_smi = tag_to_struct[(library_id, 1, cycle1)]
    cycle2_smi = tag_to_struct[(library_id, 2, cycle2)]
    cycle3_smi = tag_to_struct[(library_id, 3, cycle3)]
    a1, a2, a3 = False, False, False, # ambiguous cycle
    
    smi, a1 = run_cycle1(library_id, cycle1_smi)
    smi, a2 = run_cycle2(library_id, cycle2_smi, smi)
    smi, a3 = run_cycle3(library_id, cycle3_smi, smi)
    
    return smi, (a1, a2, a3)


def cpd_def_tup_to_smi(tup):
    return (tup[0],) + cpd_def_to_smi(*tup[1:])

if __name__ == '__main__': 
    cpd_def_to_smi(2, 1, 10, 15)
    cpd_def_to_smi(4, 1, 5, 8)
    cpd_def_to_smi(5, 1, 5, 8)
    cpd_def_to_smi(6, 1, 2, 8)
    cpd_def_to_smi(7, 1, 2, 8)