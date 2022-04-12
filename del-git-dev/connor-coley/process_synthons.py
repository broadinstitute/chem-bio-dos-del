import os, sys
import argparse
import json
import logging 
import pandas as pd 
import numpy as np
from tqdm import tqdm 

from analysis import Analysis
from enrichments import R_range, R_ranges, R, R_lb, R_ub
from scipy.stats import rankdata
from itertools import combinations
from collections import OrderedDict


class SynthonAnalysis(object):
    def __init__(self, config):
        self.config = config
        self.outfile = config['outfile'] + '.csv'
        self.save_often = int(config['settings'].get('save_often', 12))

    def load_tags(self):
        self.df_tags = pd.read_csv(self.config['metadata']['encoding'])

    def save_results(self):
        self.results.to_csv(self.outfile)
        logging.info('Saved current results DF as csv at {}'.format(self.outfile))

class MonosynthonAnalysis(SynthonAnalysis):

    def try_restore_results(self):
        self.results = None
        if os.path.isfile(self.outfile):
            self.results = pd.read_csv(self.outfile)
            logging.info('Reloaded results file with columns ')
            logging.info('{}'.format(self.results.columns))
        else:
            logging.info('No existing results file found, creating from encoding metadata')
            results = pd.read_csv(self.config['metadata']['encoding'])
            self.results = results[results['library_id'].isin(self.config['metadata']['library_ids'])]

    def calculate_enrichments(self):
        '''Go through each data file and calcualte enrichments.

        Will add columns for R_{name}, Rlb_{name}, Rub_{name} as enrichment values.
        Will add columns for rank_{name}, where rank is the rank within this monosynthon series'''
        # TODO: using .groupby is almost certainly faster than this manual iteraction

        ctr = 0
        for dataset_key in self.config['datasets']:
            logging.info(f'Processing dataset {dataset_key}')
            counts = pd.read_csv(self.config['datasets'][dataset_key]['counts_file'])

            for enrichment_key in tqdm(self.config['datasets'][dataset_key]['enrichments']):
                logging.info(f'Processing enrichment {enrichment_key}')

                if f'R_{enrichment_key}' in self.results.columns and not self.config['settings']['reprocess']:
                    logging.info('Already found in results DF! Skipping')
                    continue

                # Get column names of counts to compare
                exp_tot, beads_tot = self.config['datasets'][dataset_key]['enrichments'][enrichment_key]

                all_Rs = []; all_Rs_lb = []; all_Rs_ub = [];
                all_ranks = []; all_ranks_lb = []; all_ranks_ub = [];

                for lib_id in self.results['library_id'].unique():
                    counts_lib = counts[counts['library_id'] == lib_id]
                    results_lib = self.results.loc[self.results['library_id']==lib_id]
                    
                    for cyclenum in results_lib['cycle'].unique():
                        cycle = f'cycle{int(cyclenum)}'

                        Rs = []; Rs_lb = []; Rs_ub = [];

                        exposure1 = counts[beads_tot].sum()
                        exposure2 = counts[exp_tot].sum()
                        for cycle_id in results_lib.loc[results_lib['cycle']==cyclenum]['cycle_id']:
                            count1 = counts_lib.loc[counts_lib[cycle]==cycle_id][beads_tot].sum()
                            count2 = counts_lib.loc[counts_lib[cycle]==cycle_id][exp_tot].sum()
                            R, R_lb, R_ub = R_range(count1, exposure1, count2, exposure2)

                            Rs.append(R)
                            Rs_lb.append(R_lb)
                            Rs_ub.append(R_ub)

                        # Calculate ranks within monosynthon series (fixed lib, fixed cycle #)
                        ranks = rankdata(-np.array(Rs), method='min')
                        ranks_lb = rankdata(-np.array(Rs_lb), method='min')
                        ranks_ub = rankdata(-np.array(Rs_ub), method='min')

                        all_Rs.extend(Rs)
                        all_Rs_lb.extend(Rs_lb)
                        all_Rs_ub.extend(Rs_ub)
                        all_ranks.extend(list(ranks))
                        all_ranks_lb.extend(list(ranks_lb))
                        all_ranks_ub.extend(list(ranks_ub))
                        
                # Save to DF
                self.results[f'R_{enrichment_key}'] = list(zip(all_Rs_lb, all_Rs, all_Rs_ub))
                self.results[f'rank_{enrichment_key}'] = list(zip(all_ranks_lb, all_ranks, all_ranks_ub))
                ctr += 1

                if (ctr % self.save_often) == 0:
                    self.save_results()
                    ctr = 1

        self.save_results()

class DisynthonAnalysis(SynthonAnalysis):

    def try_restore_results(self):
        self.results = None
        if os.path.isfile(self.outfile):
            self.results = pd.read_csv(self.outfile)
            logging.info('Reloaded results file with columns ')
            logging.info('{}'.format(self.results.columns))
        else:
            logging.info('No existing results file found, creating from encoding metadata')
            df_tags = pd.read_csv(self.config['metadata']['encoding'])
            df_tags = df_tags[df_tags['library_id'].isin(self.config['metadata']['library_ids'])]
            disynthons = {
                'library_id': [],
                'cycle_A': [],
                'cycle_B': [],
                'cycle_id_A': [],
                'cycle_id_B': [],
                'structure_A': [],
                'structure_B': [],
            }
            for lib_id, df_tags_lib in df_tags.groupby('library_id'):
                for cycle_A, cycle_B in combinations([1,2,3], 2):
                    df_tags_A = df_tags_lib[df_tags_lib['cycle']==cycle_A]
                    df_tags_B = df_tags_lib[df_tags_lib['cycle']==cycle_B]
                    disynthons['library_id'].extend([lib_id] * (len(df_tags_A)*len(df_tags_B)))
                    disynthons['cycle_A'].extend([cycle_A] * (len(df_tags_A)*len(df_tags_B)))
                    disynthons['cycle_B'].extend([cycle_B] * (len(df_tags_A)*len(df_tags_B)))
                                                 
                    for cycle_id_A, structure_A in zip(df_tags_A['cycle_id'], df_tags_A['structure_id']):
                        disynthons['cycle_id_A'].extend([cycle_id_A] * len(df_tags_B))
                        disynthons['structure_A'].extend([structure_A] * len(df_tags_B))
                                                        
                        for cycle_id_B, structure_B in zip(df_tags_B['cycle_id'], df_tags_B['structure_id']):
                            disynthons['cycle_id_B'].append(cycle_id_B)
                            disynthons['structure_B'].append(structure_B)

            self.results = pd.DataFrame.from_dict(disynthons)
            self.save_results()

    def calculate_enrichments(self):
        '''Go through each data file and calcualte enrichments.

        Will add columns for R_{name}, Rlb_{name}, Rub_{name} as enrichment values.
        Will add columns for rank_{name}, where rank is the rank within this disynthon series'''

        ctr = 0
        for dataset_key in self.config['datasets']:
            logging.info(f'Processing dataset {dataset_key}')
            counts = pd.read_csv(self.config['datasets'][dataset_key]['counts_file'])

            for enrichment_key in tqdm(self.config['datasets'][dataset_key]['enrichments']):
                logging.info(f'Processing enrichment {enrichment_key}')

                if f'R_{enrichment_key}' in self.results.columns and not self.config['settings']['reprocess']:
                    logging.info('Already found in results DF! Skipping')
                    continue

                # Get column names of counts to compare
                exp_tot, beads_tot = self.config['datasets'][dataset_key]['enrichments'][enrichment_key]
                exposure1 = counts[beads_tot].sum()
                exposure2 = counts[exp_tot].sum()

                all_Rs = []; all_Rs_lb = []; all_Rs_ub = [];
                all_ranks = []; all_ranks_lb = []; all_ranks_ub = [];

                for lib_id in self.results['library_id'].unique():
                    logging.info(f'....working on library {lib_id:g}')
                    counts_lib = counts[counts['library_id'] == lib_id]
                    
                    # for cyclenum_A in results_lib['cycle_A'].unique():
                    #     cycle_A = f'cycle{int(cyclenum_A)}'
                    #     results_lib_A = results_lib[results_lib['cycle_A']==cyclenum_A]

                    #     for cyclenum_B in results_lib_A['cycle_B'].unique():
                    #         cycle_B = f'cycle{int(cyclenum_B)}'
                    #         results_lib_A_B = results_lib_A[results_lib_A['cycle_B']==cyclenum_B]

                    #         # Rs = []; Rs_lb = []; Rs_ub = [];
                    #         count1s = np.zeros((len(results_lib_A_B),))
                    #         count2s = np.zeros((len(results_lib_A_B,)))

                    #         for i,(cycle_id_A, cycle_id_B) in tqdm(enumerate(
                    #                     zip(results_lib_A_B['cycle_id_A'], results_lib_A_B['cycle_id_B']))):
                    #             count1s[i] = counts_lib[
                    #                 (counts_lib[cycle_A]==cycle_id_A) & (counts_lib[cycle_B]==cycle_id_B)
                    #             ][beads_tot].sum()
                    #             count2s[i] = counts_lib[
                    #                 (counts_lib[cycle_A]==cycle_id_A) & (counts_lib[cycle_B]==cycle_id_B)
                    #             ][exp_tot].sum()

                    #         Rs, Rs_lb, Rs_ub = R_ranges(count1s, exposure1, count2s, exposure2)

                    #         # Calculate ranks within disynthon series (fixed lib, fixed cycle_A, fixed cycle_B)
                    #         ranks = rankdata(-np.array(Rs), method='min')
                    #         ranks_lb = rankdata(-np.array(Rs_lb), method='min')
                    #         ranks_ub = rankdata(-np.array(Rs_ub), method='min')

                    #         all_Rs.extend(list(Rs))
                    #         all_Rs_lb.extend(list(Rs_lb))
                    #         all_Rs_ub.extend(list(Rs_ub))
                    #         all_ranks.extend(list(ranks))
                    #         all_ranks_lb.extend(list(ranks_lb))
                    #         all_ranks_ub.extend(list(ranks_ub))

                    for cyclenum_A, cyclenum_B in combinations([1,2,3], 2):
                        cycle_A = f'cycle{int(cyclenum_A)}'
                        cycle_B = f'cycle{int(cyclenum_B)}'

                        count1s = []; count2s = []
                        for cyc_idxs,counts_cyc in tqdm(counts_lib.groupby([cycle_A, cycle_B])):
                            count1s.append(counts_cyc[beads_tot].sum())
                            count2s.append(counts_cyc[exp_tot].sum())

                        count1s = np.array(count1s)
                        count2s = np.array(count2s)

                        Rs, Rs_lb, Rs_ub = R_ranges(count1s, exposure1, count2s, exposure2)

                        # Calculate ranks within disynthon series (fixed lib, fixed cycle_A, fixed cycle_B)
                        ranks = rankdata(-np.array(Rs), method='min')
                        ranks_lb = rankdata(-np.array(Rs_lb), method='min')
                        ranks_ub = rankdata(-np.array(Rs_ub), method='min')

                        all_Rs.extend(list(Rs))
                        all_Rs_lb.extend(list(Rs_lb))
                        all_Rs_ub.extend(list(Rs_ub))
                        all_ranks.extend(list(ranks))
                        all_ranks_lb.extend(list(ranks_lb))
                        all_ranks_ub.extend(list(ranks_ub))
                        
                # Save to DF
                self.results[f'R_{enrichment_key}'] = list(zip(all_Rs_lb, all_Rs, all_Rs_ub))
                self.results[f'rank_{enrichment_key}'] = list(zip(all_ranks_lb, all_ranks, all_ranks_ub))
                ctr += 1

                if (ctr % self.save_often) == 0:
                    self.save_results()
                    ctr = 1

        self.save_results()

def main(config):
    if config['settings']['type'] == 'monosynthons':
        synthonanalysis = MonosynthonAnalysis(config)
    elif config['settings']['type'] == 'disynthons':
        synthonanalysis = DisynthonAnalysis(config)
    else:
        raise ValueError('Unknown analysis type in settings:type.')
    synthonanalysis.load_tags()
    synthonanalysis.try_restore_results()
    synthonanalysis.calculate_enrichments()
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='process count files from multiple DEL experiments')
    parser.add_argument('config_file', type=argparse.FileType('r'))
    args = parser.parse_args()

    config = json.load(args.config_file)
    os.chdir(config['project_dir'])

    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO,
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(config['outfile'] + '.log')
        ]
    )
    logging.info('Loaded config file')
    logging.info(config)

    main(config)