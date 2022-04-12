import pandas as pd
import os, sys
import argparse
import json

def main(config_tbl):
    comb_lib_meta = pd.read_csv('../data/comb_lib_meta.csv')

    config = {}
    config['project_dir'] = None
    config['enrichments'] = {}
    config['metadata'] = {}
    config['metadata']['library'] = '../lib_meta.csv'
    config['metadata']['encoding'] = '../bb_meta.csv'
    config['metadata']['encoding_enumerated'] = '../lib_enum_struct.csv'
    config['settings'] = {}
    config['settings']['fuzzy'] = True
    config['settings']['umi'] = "directional"
    config['settings']['minread'] = 1
    config['settings']['seq_slices'] = [[0,9],[9,18],[18,27],[27,37],[37,50]]
    config['settings']['read_subs'] = {}
    config['save_often'] = 2
    config['outfile'] = "baf_test"
    config['reprocess_reads'] = False
    config['reprocess_counts'] = False
    config['reprocess_enrichments'] = True
    config['njobs'] = 4
    config['metadata']['library_ids'] = comb_lib_meta[comb_lib_meta['comb_lib_id'].isin(config_tbl['comb_lib_id'])]['lib_id'].tolist()

    expts = {}
    for i,row in config_tbl.iterrows():
        expts[i] = {}
        expts[i]['id'] = "run{:03}_samp{:06}".format(row['run_id'], row['samp_id'])
        expts[i]['fname'] = row['fastq_name']
        expts[i]['library_ids'] = comb_lib_meta[comb_lib_meta['comb_lib_id'] == row['comb_lib_id']]['lib_id'].tolist()
    config['datasets'] = expts

    outname = args.config_file.name.replace('.csv', '.json')
    with open(outname, 'w') as outfile:
        json.dump(config, outfile)
    
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='create config json from csv')
    parser.add_argument('config_file', type=argparse.FileType('r'))
    args = parser.parse_args()
    
    config_tbl = pd.read_csv(args.config_file)
   
    main(config_tbl)


