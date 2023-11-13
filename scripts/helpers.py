import os
from itertools import zip_longest
import re
import pandas as pd



# check that the path ends with '/' if not then add it.
def is_correct_path(dir_path):
    if not dir_path.endswith('/'):
        return dir_path +'/' # for some reason os.path.join refused to add '/'
    else:
        return dir_path

# multiqc checks if jsons produced by fastp end with fast.json, 
# the function below renames fastp jsons incase they dont end with fastp.jsons
def fastp_jsons(path):
    for json in os.listdir(path):
        if not json.endswith('fastp.json'):
            base_name = json.split('.')[0]
            os.rename(path+json, path+base_name + '-fastp.json')
    return None



# select DeepVirFinder predictions with a given score and pval
def filterDeepVirFinder(report, score_cutoff=0.95, pval_cutoff=0.05):

    dvf_pd = pd.read_csv(report, sep="\t")
    
    query = f'score >= {score_cutoff} & pvalue < {pval_cutoff}'
    dvf_pd.query(query, inplace=True) 
    
    return dvf_pd


# # Analyse DeepvirFinder

# possibleViruses = filterDeepVirFinder("/lustre/shared/wfsr-mcfa/projects/internships/luka/viral_metagenomics_pipeline/scripts/test/deepvirfinder/before_rr.fasta_gt1bp_dvfpred.txt")

# print(possibleViruses.head())
######################################################################
####  NEVER MIND THE CODES BELOW; 
#### turns out glob_wildcards does what all the codes below perform
######################################################################

    

# reads [(forward, reverse), ...]
#https://docs.python.org/3/library/itertools.html

def grouper(iterable, n, *, incomplete='fill', fillvalue=None): 
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    if incomplete == 'fill':
        return zip_longest(*args, fillvalue=fillvalue)
    if incomplete == 'strict':
        return zip(*args, strict=True)
    if incomplete == 'ignore':
        return zip(*args)
    else:
        raise ValueError('Expected fill, strict, or ignore')


# create forward and reverse reads

def read_names(reads):
    SAMPLES = []
    FORWARDS = []
    REVERSES = []
    with open('fowards.txt', 'w') as fowards, \
         open('reverses.txt', 'w') as reverses, \
         open('samples.txt', 'w') as samples:
    
    
        for sample_no, read in enumerate(reads):

    # extract the base names for the reads
    
            base_name_r1 = re.split('_R1', read[0].split('.')[0], \
                            flags=re.IGNORECASE)[0]
    
            base_name_r2 = re.split('_R2', read[1].split('.')[0], \
                            flags=re.IGNORECASE)[0]
    
    #print(read[0].split('.')[0], read[1].split('.')[0])

    #ensure that r1 and r2 are correct.
            if base_name_r1 == base_name_r2:
               
          
               fowards.write(f'{os.path.join(os.getcwd(),read[0])}\n') # append any of r1/r2 to SAMPLES
               reverses.write(f'{read[1]}\n')
               samples.write(f'{base_name_r1}\n')  

               SAMPLES.append(base_name_r1)  
               FORWARDS.append(read[0])
               REVERSES.append(read[1])          
            
            else:
                raise NameError(f'forward and reverse reads do not have\
                the same base name for sample {sample_no + 1}')



    forward_path = os.path.abspath('forwards.txt')
    reverse_path = os.path.abspath('reverses.txt')
    samples_path= os.path.abspath('samples.txt')
    
    return  FORWARDS, REVERSES, SAMPLES#forward_path, reverse_path, samples_path


