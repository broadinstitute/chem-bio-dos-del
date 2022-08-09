import csv 
import pandas as pd

#########################################################
# 
# Convert list to string with separator
#
#########################################################
def listToString(stringaslist,delimiter): 
    # initialize an empty string
    str = "" 
    # traverse the list and concatenate to String variable
    for element in stringaslist: 
        str += (element + delimiter)   
    # return string  
    return str 

def read_Log_File():
    lines = ''
    types = ''
    sample_dict = {}
    column_names = ['sample','total_reads','total_unique_reads','valid_reads','valid_barcodes','total_counts']
    csv_filename = 'stats_table.csv'
    is_valid_total = False    
    is_data_reload = True
    total_reads	= 0
    total_unique_reads	= 0
    valid_reads	= 0
    valid_barcodes = 0
    total_counts = 0
    stats_df = pd.DataFrame(columns=column_names)

    with open(csv_filename, 'w') as outfile:
        outfile.write(listToString(column_names,',')[:-1])
        outfile.write('\n')
    with open('save_file.log', 'r') as infile:
        lines = [line.rstrip() for line in infile]

    for type in lines:
        if 'Saved reads to' in type:
            is_data_reload = False
        if is_data_reload == False:
            sample = Sample()
            if 'Validating reads...' in type:
                is_valid_total = True
            if type.endswith('unique reads'):
                total_unique_reads = (int(type[type.find('INFO:')+6:].strip().split(' ')[0]))
            if type.endswith('total reads') and is_valid_total == True :
                valid_reads = (int(type[type.find('INFO:')+6:].strip().split(' ')[0]))
                is_valid_total = False
            elif type.endswith('total reads'):
                total_reads = (int(type[type.find('INFO:')+6:].strip().split(' ')[0]))
            if 'unique barcodes' in type:
                valid_barcodes = (int(type[type.find('INFO:')+6:].strip().split(' ')[0]))
            if 'Final total counts:' in type:
                total_counts = (int(type[type.find('INFO:')+6:].strip().split(' ')[3]))
            if '_samp' in type and 'Processing count data' in type:
                sample_name = type[type.find('run'):].split('_')[1]
                is_data_reload = True
                if sample_name in sample_dict:
                    sample = sample_dict.get(sample_name)
                sample.add_total_reads(total_reads)
                sample.add_total_unique_reads(total_unique_reads)
                sample.add_valid_reads(valid_reads)
                sample.add_valid_barcodes(valid_barcodes)
                sample.add_total_counts(total_counts)
                sample_dict[sample_name] = sample
               
    for key in sample_dict.keys():
        stats_df.loc[0,['sample']] = key
        stats_df.loc[0,['total_reads']] = sample_dict[key].get_total_reads()
        stats_df.loc[0,['total_unique_reads']] = sample_dict[key].get_total_unique_reads()
        stats_df.loc[0,['valid_reads']] = sample_dict[key].get_valid_reads()
        stats_df.loc[0,['valid_barcodes']] = sample_dict[key].get_valid_barcodes()
        stats_df.loc[0,['total_counts']] = sample_dict[key].get_total_counts()
        stats_df.to_csv(csv_filename, sep=',', encoding='utf-8',mode='a', index=False, header=False)



class Sample(object):
    total_reads	= 0
    total_unique_reads	= 0
    valid_reads	= 0
    valid_barcodes = 0
    total_counts = 0
    def add_total_reads(self,value):
        self.total_reads +=  value
    def get_total_reads(self):
        return self.total_reads  
    def add_total_unique_reads(self,value):
        self.total_unique_reads += value
    def get_total_unique_reads(self):
        return self.total_unique_reads 
    def add_valid_reads(self,value):
        self.valid_reads += value
    def get_valid_reads(self):
        return self.valid_reads 
    def add_valid_barcodes(self,value):
        self.valid_barcodes += value
    def get_valid_barcodes(self):
        return self.valid_barcodes 
    def add_total_counts(self,value):
        self.total_counts += value
    def get_total_counts(self):
        return self.total_counts 


########################################################
#
# Main method to initiate chain of method calls 
# 
#########################################################
def main():
    read_Log_File()

if __name__ == '__main__':
    main()