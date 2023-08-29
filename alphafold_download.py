import pandas as pd
import os
import sys

#A code to automatically download from AlphaFoldDB V4 models of protein 3D structure.
#How to use : give to the script, as an argument, a tsv file which has a column "Id_AF_Model", which contains the ID of the model to be downloaded

tsv_file = sys.argv[1]
AF_id_df = pd.read_csv(tsv_file, sep='\t')
base_path = os.getcwd()

for ind in AF_id_df.index :

	os.chdir(base_path)
	new_path = './3D_structures'
	isExist = os.path.exists(new_path)
	if not isExist :
		os.makedirs(new_path)
	os.chdir(new_path)

	model_url='https://alphafold.ebi.ac.uk/files/'+AF_id_df['Id_AF_Model'][ind]+'-model_v4.pdb'
	error_url='https://alphafold.ebi.ac.uk/files/'+AF_id_df['Id_AF_Model'][ind]+'-predicted_aligned_error_v4.json'

	os.system('curl '+model_url+' -o '+AF_id_df['Id_AF_Model'][ind]+'.pdb')
	os.system('curl '+error_url+' -o '+AF_id_df['Id_AF_Model'][ind]+'.json')
