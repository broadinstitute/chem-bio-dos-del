import pandas as pd
import numpy as np
import logging
from collections import OrderedDict
from tqdm import tqdm
from umi_tools._dedup_umi import edit_distance
import time
import os
import json
from gzip import open as gzopen
from Bio import SeqIO, bgzf

from enrichments import R_ranges

class Library(object):
    '''Defines a DNA-encoded library.

    Useful for parsing a set of metadata files and preparing
    lists of expected tags, barcodes, etc.'''

    def __init__(self, config):
        # Initial loading
        self.df_cpds = pd.read_csv(config['metadata']['encoding_enumerated'])
        self.df_tags = pd.read_csv(config['metadata']['encoding'])
        self.df_libs = pd.read_csv(config['metadata']['library'])
        logging.info('Loaded library metadata')
        self.stats()

        # Truncation to libraries of interest
        library_ids = config['metadata']['library_ids']
        self.df_cpds = self.df_cpds[self.df_cpds['library_id'].isin(library_ids)]
        self.df_tags = self.df_tags[self.df_tags['library_id'].isin(library_ids)]
        self.df_libs = self.df_libs[self.df_libs['library_id'].isin(library_ids)]
        self.library_ids = config['metadata']['library_ids']
        logging.info('Cut library metadata to only include IDs {}'.format(config['metadata']['library_ids']))
        self.stats()

        # Create resolvers (map seq -> ID)
        self.tags = {lib_id: OrderedDict() for lib_id in self.library_ids}
        self.structures = {lib_id: OrderedDict() for lib_id in self.library_ids}
        for lib_id in self.tags:
            for cyc_idx in sorted(self.df_tags[self.df_tags['library_id']==lib_id]['cycle'].unique()):
                self.tags[lib_id][cyc_idx] = {}
                self.structures[lib_id][cyc_idx] = {}
        for i, row in self.df_tags.iterrows():
            self.tags[row['library_id']][row['cycle']][row['sequence']] = row['cycle_id']
            self.structures[row['library_id']][row['cycle']][row['cycle_id']] = row['structure_id']

        self.lib_tags = {}
        self.lib_structures = {}
        for i, row in self.df_libs.iterrows():
            self.lib_tags[row['library_sequence']] = row['library_id']
            self.lib_structures[row['library_id']] = row['library_name']

        # Check whether fuzzy matching is appropriate (requires highly unique tags!)
        self.use_fuzzy = {'lib': 0}
        for lib_id in self.tags:
            self.use_fuzzy[lib_id] = {cyc_idx: 0 for cyc_idx in self.tags[lib_id].keys()}
        if config['settings']['fuzzy']:
            for lib_id in self.tags:
                for cyc_idx in self.tags[lib_id].keys():
                    expected_tags = list(self.tags[lib_id][cyc_idx].keys())
                    min_dist = 5
                    for i,tag1 in enumerate(expected_tags):
                        for j,tag2 in enumerate(expected_tags):
                            if j <= i:
                                continue
                            min_dist = min(min_dist, edit_distance(tag1.encode('utf-8'), tag2.encode('utf-8')))
                    if min_dist >= 5: # only really possible for library tags
                        self.use_fuzzy[lib_id][cyc_idx] = 2 
                    elif min_dist >= 3:
                        self.use_fuzzy[lib_id][cyc_idx] = 1
            # Now check if library should be fuzzy
            expected_tags = list(self.lib_tags.keys())
            min_dist = 5
            for i,tag1 in enumerate(expected_tags):
                for j,tag2 in enumerate(expected_tags):
                    if j <= i:
                        continue
                    min_dist = min(min_dist, edit_distance(tag1.encode('utf-8'), tag2.encode('utf-8')))
            if min_dist >= 5: # only really possible for library tags
                self.use_fuzzy['lib'] = 2 
            elif min_dist >= 3:
                self.use_fuzzy['lib'] = 1

        # Report
        logging.info('Allowed fuzzy matching dists:?')
        for k,v in self.use_fuzzy.items():
            logging.info('    {}: {}'.format(k, v))
            

    def stats(self):
        logging.info('    {} unique libraries'.format(len(self.df_libs)))
        logging.info('    {} unique compounds'.format(len(self.df_cpds)))
        


    def resolve_seq(self, lib_id, cyc_idx, seq, cache=True):
        '''Resolve a sequence to a cycle_id. Must know the expected lib_id 
        If no match is found, returns None'''

        # Look for exact match
        if seq in self.tags[lib_id][cyc_idx]:
            return self.tags[lib_id][cyc_idx][seq]
        # If not using fuzzy, that was our only chance to match
        if not self.use_fuzzy[lib_id][cyc_idx]:
            if cache:
                self.tags[lib_id][cyc_idx][seq] = None
            return None
        # Look for match that is close enough
        for tag in self.tags[lib_id][cyc_idx]:
            if edit_distance(seq.encode('utf-8'), tag.encode('utf-8')) <= self.use_fuzzy[lib_id][cyc_idx]:
                if cache:
                    self.tags[lib_id][cyc_idx][seq] = self.tags[lib_id][cyc_idx][tag]
                return self.tags[lib_id][cyc_idx][tag]
        # Nothing found
        if cache:
            self.tags[lib_id][cyc_idx][seq] = None
        return None


    def resolve_seq_lib(self, seq, cache=True):
        '''Resolev a sequence to a library. If no match, return None'''

        # Exact
        if seq in self.lib_tags:
            return self.lib_tags[seq]
        # Was that our only shot?
        if not self.use_fuzzy['lib']:
            if cache:
                self.lib_tags[seq] = None 
            return None 
        # Find acceptable match
        for tag in self.lib_tags:
            if edit_distance(seq.encode('utf-8'), tag.encode('utf-8')) <= self.use_fuzzy['lib']:
                if cache:
                    self.lib_tags[seq] = self.lib_tags[tag]
                return self.lib_tags[tag]
        # Nothing found
        if cache:
            self.lib_tags[seq] = None
        return None

class Experiment(object):
    '''Defines the result of a single experiment/control.

    Kind of an unnecessary class for now, but may want to store counts more permanently'''

    def __init__(self, fpath):
        self.fpath = fpath # doubles as path
        self.fpath_reads = fpath + '.reads'
        self.fpath_counts = fpath + '.reads'

    def change_char(self, s, p, r):
        return s[:p]+r+s[p+1:]

    def get_reads(self, reprocess=False, read_subs={}):
        fpath = self.fpath
        fpath_reads = self.fpath_reads

        if os.path.isfile(fpath_reads):
            logging.info(f'Found already-processed read file {fpath_reads}')
            if not reprocess:
                with open(fpath_reads, 'r') as fid:
                    logging.info(f'Reloading...')
                    return json.load(fid)
            else:
                logging.info('Reprocessing according to "--reprocess" flag')

        logging.info(f'Preparing to read {fpath}')
        reads = {}
        if '.gz' in fpath:
            line_reader = SeqIO.parse(gzopen(fpath, "rt"), "fastq")
        else:
            line_reader = SeqIO.parse(open(fpath, "rt"), "fastq")
        for line in tqdm(line_reader):
            read = str(line.seq)
            # Fixed replacements (e.g., for conserved regions)
            for k,v in read_subs.items():
                if read[k] != v:
                    read = self.change_char(read, k, v)
            if read in reads:
                reads[read] += 1
            else:
                reads[read] = 1
        logging.info(f'Finished reading {fpath}')

        with open(fpath_reads, 'w') as fid:
            json.dump(reads, fid)
            logging.info(f'Saved reads to {fpath_reads}')

        return reads


class Analysis(object):
    '''Defines the total results from a DEL experiment'''
    def __init__(self):
        self.experiments = []
        self.enrichments = []
        self.beads = []
        self.library = None
        self.umi_mode = 'directional'
        self.clusterer = None
        self.seq_slices = []
        self.results = None
        self.outfile = ''
        self.fuzzy = True
        self.read_subs = {}


    def load(self, config):
        self.config = config
        self.umi_mode = config['settings']['umi']
        self.reprocess_reads = config['reprocess_reads']
        self.reprocess_counts = config['reprocess_counts']
        self.reprocess_enrichments = config['reprocess_enrichments']
        self.save_often = int(config.get('save_often', 1))
        self.seq_slices = config['settings']['seq_slices'] # cyc1, cyc2, ..., lib, umi
        self.outfile = config['outfile']+'.csv'
        self.fuzzy = config['settings']['fuzzy']
        self.minread = config['settings']['minread']
        self.enrichments = config['enrichments']
        self.read_subs = config['settings'].get('read_subs', {})
        self.read_subs = {int(k):v for k,v in self.read_subs.items()}
        

        if self.umi_mode == 'directional':
            from umi_tools.network import UMIClusterer
            self.clusterer = UMIClusterer()

        # Load library metadata 
        self.library = Library(config)

        # Initialize experiments
        self.experiments = OrderedDict()
        for i,exp_set in enumerate(config['datasets']):
            for j,fname in enumerate(config['datasets'][exp_set]):
                self.experiments['{}_r{}'.format(exp_set,j+1)] = Experiment(fname)
        
        # Initialize reporting dict
        self.record = OrderedDict()
        for exp in self.experiments:
            self.record[exp] = {}


    def run_counts(self):
        self.try_restore_results()
        if self.results is None:
            self.initialize_results()
            self.save_results()

        # Churn through experiments
        for i,exp in enumerate(self.experiments):
            if exp in self.results.columns and not self.reprocess_counts:
                logging.info(f'Already processed experiment {exp}, skipping')
                continue

            logging.info(f'Running analysis for experiment {exp}')

            # Raw reads
            reads = self.experiments[exp].get_reads(
                reprocess=self.reprocess_reads, read_subs=self.read_subs
            )
            self.read_stats(reads, record_exp=exp)

            # With sequence matching (no dealing with UMIs yet)
            valid_reads = self.validify_reads(reads)
            self.valid_read_stats(valid_reads, record_exp=exp)

            # Now convert full UMI read counts to barcode counts
            if isinstance(self.minread, int):
                minread = self.minread
            elif isinstance(self.minread, dict):
                if exp in self.minread:
                    minread = self.minread[exp]
                else:
                    minread = self.minread.get('default', 1)
            else:
                logging.error('Invalid minread encountered')

            logging.info(f'Using a minimum read cutoff of {minread} for {exp}')
            counts = self.convert_to_counts(valid_reads, umi_mode=self.umi_mode, 
                minread=minread, record_exp=exp)

            # Save
            self.update_results(exp, counts)
            if self.save_often and ((i+1) % self.save_often) == 0:
                self.save_results()

        # Assign totals
        logging.info('Defining total counts for different experiments')
        for i,exp_set in enumerate(self.config['datasets']):
            totals = np.zeros((len(self.results),), dtype=int)
            for j,fname in enumerate(self.config['datasets'][exp_set]):
                totals += np.array(self.results['{}_r{}'.format(exp_set,j+1)])
            self.results['{}_tot'.format(exp_set)] = totals

        # Save final DF
        self.save_results()

        # Report
        keys = ['total_reads', 'total_unique_reads', 'valid_reads',
            'valid_barcodes', 'total_counts']
        print('\t'.join(['experiment'] + keys))
        for exp in self.experiments:
            try:
                print('\t'.join([exp] + [str(self.record[exp][key]) for key in keys]))
            except KeyError:
                pass

    def run_enrichments(self):
        for enrichment in self.enrichments:
            if f'{enrichment}_R' in self.results.columns and not self.reprocess_enrichments:
                logging.info(f'Already calculated enrichments for {enrichment}, skipping')
                continue
            counts_exp = np.array(self.results[self.enrichments[enrichment][0] + '_tot'])
            counts_ref = np.array(self.results[self.enrichments[enrichment][1] + '_tot'])
            exposure_exp = counts_exp.sum() * np.ones(counts_exp.shape)
            exposure_ref = counts_ref.sum() * np.ones(counts_ref.shape)

            logging.info(f'Calculating R ranges for enrichment definition {enrichment}')
            Rs, Rs_lb, Rs_ub = R_ranges(counts_ref, exposure_ref, counts_exp, exposure_exp)

            self.results[f'{enrichment}_R'] = Rs
            self.results[f'{enrichment}_R_lb'] = Rs_lb
            self.results[f'{enrichment}_R_ub'] = Rs_ub
            if self.save_often and ((i+1) % self.save_often) == 0:
                self.save_results()
        # Final save
        self.save_results()

    def try_restore_results(self):
        self.results = None
        if os.path.isfile(self.outfile):
            self.results = pd.read_csv(self.outfile)
            logging.info('{}'.format(self.results.columns))

    def initialize_results(self):
        self.results = self.library.df_cpds.copy()
        logging.info('Creating copy of enumerated compound list to store results in')

        # # Do library first
        # lib_id_lookup = {v:k for k,v in self.library.lib_tags.items()}
        # self.results['library_sequence'] = [lib_id_lookup[tag] for tag in self.results['library_id']]
        # logging.info('Created new column in results DF: library_sequence')

        # # Then sequences
        # cyc_id_lookup = {lib_id: {} for lib_id in self.library.tags.keys()}
        # for lib_id in self.library.tags.keys():
        #     for cyc_idx in self.library.tags[lib_id].keys():
        #         cyc_id_lookup[lib_id][cyc_idx] = {v:k for k,v in self.library.tags[lib_id][cyc_idx].items()}

        # for cyc_idx in [1,2,3]: # don't want this to be hardcoded eventually!
        #     self.results[f'cycle{cyc_idx}_seq'] = [cyc_id_lookup[lib_id][cyc_idx][_id] for \
        #         lib_id, _id in zip(self.results['library_id'], self.results[f'cycle{cyc_idx}'])]
        #     logging.info(f'Created new column in results DF: cycle{cyc_idx}_seq')
            

    def update_results(self, exp, counts):
        '''Add the observed counts for this experiment to the results DF'''
        logging.info(f'Adding count data for {exp} to results DF')

        # There must be a better way to do this...
        count_list = np.zeros((len(self.results),), dtype=int)
        for i, barcode in tqdm(enumerate(zip(self.results['library_id'], 
                    self.results['cycle1'], self.results['cycle2'], self.results['cycle3']))):
            if barcode in counts:
                count_list[i] = counts[barcode]
        self.results[exp] = count_list
        logging.info('    found {}/{} tags ({:.2%})'.format((count_list>0).sum(), len(count_list), (count_list>0).mean()))


    def save_results(self):
        logging.info('Saving to {}...'.format(self.outfile))
        self.results.to_csv(self.outfile)
        logging.info('Saved current results DF as csv at {}'.format(self.outfile))


    def validify_reads(self, reads):
        '''Creates a new "reads" dictionary based on the expected tags.
        - Only keeps reads that match the expected library
        - Consolidates errors (modifies keys!) when using fuzzy matching
        First validates library, then cyc1, then cyc2, then cyc3'''
        logging.info('Validating reads...')
        
        lib_slice = self.seq_slices[-2]
        cyc_slices = self.seq_slices[:-2]
        umi_slice = self.seq_slices[-1]

        # If valid, convert raw sequence reads to a dictionary of dictionaries
        # containing a map of (lib, cyc1, cyc2, cyc3) -> (UMI) -> num_reads
        valid_reads = {}
        for seq,num_read in tqdm(reads.items()):
            invalid = False
            resolved_idxs = []

            # Validate library tag first
            lib_id = self.library.resolve_seq_lib(seq[slice(*lib_slice)])
            if lib_id is None:
                invalid = True 
                continue
            resolved_idxs.append(lib_id)

            # Then validate scaffold/BB tags
            for i,cyc_slice in enumerate(cyc_slices):
                cyc_idx = i + 1
                idx = self.library.resolve_seq(lib_id, cyc_idx, seq[slice(*cyc_slice)])
                if idx is None:
                    invalid = True
                    break
                resolved_idxs.append(idx)

            # Did we get a match?
            if not invalid:
                barcode = tuple(resolved_idxs)
                umi = seq[slice(*umi_slice)]

                if barcode in valid_reads:
                    if umi in valid_reads[barcode]:
                        valid_reads[barcode][umi] += num_read
                    else:
                        valid_reads[barcode][umi] = num_read
                else:
                    valid_reads[barcode] = {umi: num_read}

        return valid_reads

    def convert_to_counts(self, valid_reads, umi_mode='directional', minread=1,
            record_exp=None):
        '''Process UMIs to come up with a final count for each barcode combo'''
        logging.info('Converting valid reads to counts (consolidating UMIs)')
        counts = {}
        for barcode,umi_dict in tqdm(valid_reads.items()):
            if umi_mode == 'unique':
                counts[barcode] = len([v for v in umi_dict.values() if v >= minread])
            elif umi_mode == 'total':
                counts[barcode] = sum([v for v in umi_dict.values() if v >= minread])
            elif umi_mode == 'directional':
                filtered_umi_dict = {k.encode('utf-8'):v for k,v in umi_dict.items() if v >= minread}
                if filtered_umi_dict:
                    clustered_umis = self.clusterer(filtered_umi_dict, 1)
                    counts[barcode] = len(clustered_umis)
            else:
                raise ValueError('Invalid UMI mode specified')
        logging.info('Final total counts: {}'.format(sum(list(counts.values()))))

        if record_exp:
            self.record[record_exp]['total_counts'] = sum(list(counts.values()))

        return counts


    def read_stats(self, reads, record_exp=None):
        logging.info('    {} unique reads'.format(len(reads)))
        logging.info('    {} total reads'.format(sum(list(reads.values()))))
        if record_exp:
            self.record[record_exp]['total_reads'] = sum(list(reads.values()))
            self.record[record_exp]['total_unique_reads'] = len(reads)

    def valid_read_stats(self, valid_reads, record_exp=None):
        logging.info('    {} unique barcodes (out of {} possible)'.format(len(valid_reads), len(self.library.df_cpds)))
        tot = 0
        for umi_dict in valid_reads.values():
            tot += sum(list(umi_dict.values()))
        logging.info('    {} total reads'.format(tot))
        if record_exp:
            self.record[record_exp]['valid_reads'] = tot
            self.record[record_exp]['valid_barcodes'] = len(valid_reads)

