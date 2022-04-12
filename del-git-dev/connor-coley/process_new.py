import pandas as pd
import os, sys
import argparse
import json
import logging 

from analysis_new import Analysis

def main(config):
    '''Process a set of FastQ files according to a config json'''

    analysis = Analysis()
    analysis.load(config)
    analysis.run_counts()
    # analysis.run_enrichments()
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='process fastQ files for set of DEL experiments')
    parser.add_argument('config_file', type=argparse.FileType('r'))
    args = parser.parse_args()

    config = json.load(args.config_file)
    if config.get('project_dir', None) is None:
        config['project_dir'] = os.path.dirname(args.config_file.name)

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

