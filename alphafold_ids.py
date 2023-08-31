import os
import csv
import sys
import logging

import tqdm

FORMAT = '%(asctime)s:  %(message)s' 
logging.basicConfig(format=FORMAT,level=logging.INFO, handlers=[logging.StreamHandler()])


"""
	A script to find wheter or not Uniprot ids have a model available in AlphaFoldDB.
	The file accession_ids.csv need to be downloaded from http://ftp.ebi.ac.uk/pub/databases/alphafold/

	How to use : give to the script, as an arguments :  
		- path to the file accession_ids.csv
		- path to a text file containing uniprot accession, one per line.
		- path to output file
		
"""
def parse_uniprot_id(id):
	id = id.strip().replace('>','').split()[0] # remove description and '>' if any
	if '|' in id: #UniProt case
		return id.split('|')[1]
	elif id.startswith('UPI'): #UniParc case - note that UniParc sequences are not avalaible in alphafold.
		return id
	else: #UniRef_case
		return id.split('_')[-1]

def uniprot_to_alphafold_id(accession_ids_file, identifiant_uniprot, outfile):
	
	with open (accession_ids_file, 'r') as accession_file :
		accession_reader =  csv.reader(accession_file)
		accession_dict = {}
		logging.info('\033[92m {}\033[00m'.format('Parsing alphafold accession file...'))
		for row in tqdm.tqdm(accession_reader):
			accession_dict[row[0]] = row[3]
	logging.info('\033[92m {}\033[00m'.format('Done.'))
	logging.info('\033[92m {}\033[00m'.format('Start converting uniprot ids to alphafold ids ...'))
	with open (outfile, 'w') as out_handler :
		with open (identifiant_uniprot,'r') as in_handler :
			identifiants = in_handler.readlines()
			out_handler.write('Id_Sequence\tId_AF_Model\n')	
			for ids in tqdm.tqdm(identifiants, total =len(identifiants)) :
				ids_redu = parse_uniprot_id(ids)									
				ids_alphafold = ''
				if ids_redu in accession_dict:
					ids_alphafold = accession_dict[ids_redu]	
				else:
					logging.warning("\033[93m {}\033[00m" .format("Model unavailable for "+ids_redu))
				out_handler.write("{}\t{}\n".format(ids_redu,ids_alphafold))
	logging.info('\033[92m {}\033[00m'.format('Done.'))


if __name__ == "__main__":
	accession_ids_file = sys.argv[1]  # Path to accession_ids.csv
	identifiant_uniprot = sys.argv[2] # Path to a txt file containing the Uniprot ID. One ID per line.
	outfile = sys.argv[3]
	if not os.path.exists(accession_ids_file):		
		raise FileExistsError('\033[91mAlphafold accession file not found : {}\033[00m'.format(accession_ids_file))		
	if not os.path.exists(identifiant_uniprot):		
		raise FileExistsError('\033[Uniprot accession file not found : {}\033[00m'.format(identifiant_uniprot))
	
	os.makedirs(os.path.abspath(os.path.dirname(outfile)),exist_ok=True)

	uniprot_to_alphafold_id(accession_ids_file,identifiant_uniprot,outfile)
