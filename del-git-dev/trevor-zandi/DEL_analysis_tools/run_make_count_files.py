import pandas as pd
import argparse

from DEL_analysis_functions import Enrich
def main(count_specs):
	enrich = Enrich()
	enrich.initializeMetadata(count_specs)
	enrich.loadCompoundInfo()
	enrich.createCountsFile(count_specs)
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='accept CSV file containing name of count files')
	parser.add_argument('count_specs_file', type=argparse.FileType('r'))
	args = parser.parse_args()
	count_specs = pd.read_csv(args.count_specs_file)
	main(count_specs)
