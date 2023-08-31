import os
import sys
import argparse
import logging

import tqdm
import pandas as pd

FORMAT = '%(asctime)s:  %(message)s' 
logging.basicConfig(format=FORMAT,level=logging.INFO, handlers=[logging.StreamHandler()])

"""
A code to automatically download from AlphaFoldDB V4 models of protein 3D structure.
How to use : give to the script, as an arguments :  
	- a tsv file which has a column "Id_AF_Model", which contains the ID of the model to be downloaded
	- a directory where models will be stored (<outdir>/3D_structures/<models>)
"""


def download(url,outfile):
	#os.system('curl '+model_url+' -o '+AF_id_df['Id_AF_Model'][ind]+'.pdb')
	excode = os.system('curl -s {url} -o {out}'.format(url=url,out=outfile))
	if excode != 0:
		logging.error("\033[91m {}\033[00m".format('An error occured while downloading from {}'.format(url)))
		raise OSError('An error occured while downloading from {}'.format(url))

def download_models(tsv_file,outdir):

	AF_id_df = pd.read_csv(tsv_file, sep='\t')
	struc_dir = os.path.join(outdir, "3D_structures")
	os.makedirs(struc_dir,exist_ok=True)

	for ind in tqdm.tqdm(AF_id_df.index.tolist()) :
		if isinstance(AF_id_df['Id_AF_Model'][ind],str):
			model_outfile = os.path.join( struc_dir,  AF_id_df['Id_AF_Model'][ind]+'-model_v4.pdb')
			error_outfile = os.path.join( struc_dir,  AF_id_df['Id_AF_Model'][ind]+'-predicted_aligned_error_v4.json')
			model_url='https://alphafold.ebi.ac.uk/files/'+AF_id_df['Id_AF_Model'][ind]+'-model_v4.pdb'
			error_url='https://alphafold.ebi.ac.uk/files/'+AF_id_df['Id_AF_Model'][ind]+'-predicted_aligned_error_v4.json'
			
			download(model_url,model_outfile)
			download(error_url,error_outfile)
		else:
			logging.info("\033[93m {}\033[00m" .format("Missing Alphafold ID for "+AF_id_df['Id_Sequence'][ind]))

		

if __name__ == "__main__":
	tsv_file = sys.argv[1]
	outdir = sys.argv[2]
	logging.info('\033[92m {}\033[00m'.format("Start..."))
	download_models(tsv_file,outdir)
	logging.info('\033[92m {}\033[00m'.format("Done..."))