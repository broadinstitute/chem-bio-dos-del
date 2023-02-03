import pandas as pd
import argparse

from DEL_analysis_functions import Enrich
def main(count_specs):
	enrich = Enrich()
	enrich.loadCounts(count_specs)
	enrich.loadCompoundInfo()
	enrich.enrichment(2)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='accept CSV file containing name of count files')
	parser.add_argument('count_specs_file', type=argparse.FileType('r'))
	args = parser.parse_args()
	count_specs = pd.read_csv(args.count_specs_file)
	main(count_specs)
