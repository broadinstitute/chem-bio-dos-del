import pandas as pd
import numpy as np
import sys
import os
class Enrich(object):

	def __init__(self):
		self.cts = pd.DataFrame()
		self.lib = pd.DataFrame()
		self.hits = pd.DataFrame()
		self.lib_id = "0"
		self.exp_name = ""
		self.INPUT_PATH = '/home/zanditr1/DEL-analysis/inputs/'
		self.OUTPUT_PATH = '/home/zanditr1/DEL-analysis/outputs/'
		self.enrich_ub = []
		self.enrich_lb = []
		self.enrich_ml = []
		self.compound_codes = []
		self.targ_col = ""
		self.cont_col = ""
		self.count_file = pd.DataFrame()
		
	###read csv count files passed from main and extract the counts		
	def loadCounts(self, count_specs):
		### for now, assumes that we only have one pair of target and control counts
		self.exp_name = count_specs['experiment'][1]
		self.targ_col = self.exp_name+'_targ'
		self.cont_col = self.exp_name+'_cont'
		target_metadata = count_specs[count_specs['sample_type'] == "target"]
		control_metadata = count_specs[count_specs['sample_type'] == "control"]
		self.cts[self.targ_col] = self.consolidateCounts(target_metadata)
		self.cts[self.cont_col] = self.consolidateCounts(control_metadata)
		self.lib_id = str(count_specs['lib_id'][1])
		print(self.lib_id)	
	### lib_id must be assigned before calling
	def loadCompoundInfo(self):
		if(self.lib_id == "0"):
			print("No library chosen!")
			sys.exit()	
		#library_file_path = self.INPUT_PATH+ "lib_enum_struct"+self.lib_id+".csv"
		library_file_path = self.INPUT_PATH+"lib_enum_struct.csv"
		print(library_file_path)
		self.lib = pd.read_csv(library_file_path)
		self.lib.head(10)	
		self.lib = self.lib[self.lib['lib_id'] == int(self.lib_id)]
		self.lib['compound_code'] = self.lib['cycle1'].astype(str) +","+ self.lib['cycle2'].astype(str) + "," +self.lib['cycle3'].astype(str)
		self.lib.head(10)
		
	def enrichment(self,z):
		R0 = self.cts[self.targ_col]/self.cts[self.cont_col]
		Nratio = self.cts[self.cont_col].sum()/(self.cts[self.targ_col].sum())
		u1 = np.sqrt(self.cts[self.targ_col]+3/8)
		u2 = np.sqrt(self.cts[self.cont_col]+3/8)
		a = -z*z + 4*u2*u2
		aMaxL = 4*u2*u2
		b = -8*u1*u2
		c = -z*z + 4*u1*u1
		cMaxL = 4*u1*u1
		self.enrich_ub = Nratio*np.square((-b - np.sqrt(np.square(b)-4*a*c))/(2*a))
		self.enrich_lb = Nratio*np.square((-b + np.sqrt(np.square(b)-4*a*c))/(2*a))
		self.enrich_ml = Nratio*np.square((-b + np.sqrt(np.square(b)-4*aMaxL*cMaxL))/(2*aMaxL))
		topHits = pd.DataFrame(self.enrich_ub.sort_values(ascending=False), columns=[self.exp_name])
		topHits = pd.merge(topHits, self.cts, left_index=True, right_index=True)
		topHits.to_csv(os.path.join(self.OUTPUT_PATH,self.exp_name+".csv"))
		pd.merge(topHits, self.lib['compound_code'], left_index=True, right_index=True).to_csv(os.path.join(self.OUTPUT_PATH,self.exp_name+".csv"))

	def consolidateCounts(self, metadata):
		count_data = list()
		for f in metadata['file_name']:
			file_path = os.path.join(self.INPUT_PATH,f)
			if not os.path.exists(file_path):
				print("No such file! ",file_path)
				sys.exit()
			else:
				with open(os.path.join(self.INPUT_PATH,f)) as count_file:
					count_data.append([int(count) for count in count_file.read().splitlines()[1:]])
					
		counts = list(map(list, zip(*count_data)))
		count_tot = [sum(count) for count in counts]
		return count_tot
	
	def initializeMetadata(self, count_specs):
		self.lib_id = str(count_specs['lib_id'][1])
		self.exp_name = str(count_specs['experiment'][1])

	def createCountsFile(self, count_specs):
		self.count_file['compound_code'] = self.lib['compound_code']
		self.count_file['structure'] = self.lib['structure']
		self.count_file['bb1'] = self.lib['cycle1_bb']
		self.count_file['bb2'] = self.lib['cycle2_bb']
		self.count_file['bb3'] = self.lib['cycle3_bb']
		for fileNm in count_specs['file_name']:
			with open (os.path.join(self.INPUT_PATH, fileNm)) as cnt_file:
				print(count_specs[count_specs['file_name'] == fileNm]['sample_type'])
				self.count_file[count_specs[count_specs['file_name'] == fileNm]['sample_type'].values[0]] = cnt_file.read().splitlines()[1:]
		self.count_file.to_csv(os.path.join(self.OUTPUT_PATH,self.exp_name+"_counts.csv"))	
