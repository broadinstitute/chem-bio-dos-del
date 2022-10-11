import pandas as pd
import csv
from csv import writer

def main():
    '''merge structures.csv file with enumerated.csv file according to an agg07_meta_cols.csv'''
if __name__ == '__main__':
    dir_name = '../Downloads/From_Zher_Yin/' 
    agg07_meta_cols_df = pd.read_csv(dir_name + 'meta/lib/meta_cols/lib008_agg07_meta_cols.csv')
    BB_structures_df   = pd.read_csv(dir_name + 'FKBP-CIP-DEL-BB structures.csv')
    lib008_enumerated_df  = pd.read_csv(dir_name + 'meta/lib/enum_struct/lib008.csv')
    new_enum_structure = dir_name + 'meta/new_enum_structure_008.csv'

    BB_structures_df[(BB_structures_df['cyc'] == 1) & (BB_structures_df['tag_id'] == 1)]
    # First, open the old CSV file in append mode, hence mentioned as 'a'
    # Then, for the CSV file, create a file object
    print('#####', new_enum_structure)
    with open(new_enum_structure, 'a', newline='') as f_object: 
        # Pass the CSV  file object to the writer() function
        writer_object = writer(f_object)
        # agg07_meta_cols_df
        # iterate through each row and select
        # 'cycle1', 'cycle2' and 'cycle3' column respectively.
        for i in range(len(agg07_meta_cols_df)):
            cycle1_tag = agg07_meta_cols_df.loc[i, "cycle1"]
            cycle2_tag = agg07_meta_cols_df.loc[i, "cycle2"]
            cycle3_tag = agg07_meta_cols_df.loc[i, "cycle3"]
            structure  = lib008_enumerated_df.loc[i, "structures"]
            df1 = BB_structures_df[(BB_structures_df['cyc'] == 1) & (BB_structures_df['tag_id'] == cycle1_tag)]
            df2 = BB_structures_df[(BB_structures_df['cyc'] == 2) & (BB_structures_df['tag_id'] == cycle2_tag)]
            df3 = BB_structures_df[(BB_structures_df['cyc'] == 3) & (BB_structures_df['tag_id'] == cycle3_tag)]
            #print(300,cycle1_tag,cycle2_tag,cycle3_tag, df1['cyc_seq'].item(), df2['cyc_seq'].item(), df3['cyc_seq'].item(), df1['bb_smiles'].item(), df2['bb_smiles'].item(), df3['bb_smiles'].item(), structure)

            # The data assigned to the list 
            list_data=['8',cycle1_tag,cycle2_tag,cycle3_tag, df1['cyc_seq'].item(), df2['cyc_seq'].item(), df3['cyc_seq'].item(), df1['bb_smiles'].item(), df2['bb_smiles'].item(), df3['bb_smiles'].item(), structure]
            
            # Result - a writer object
            # Pass the data in the list as an argument into the writerow() function
            writer_object.writerow(list_data) 
        
    
        # Close the file object
        f_object.close()
    main()