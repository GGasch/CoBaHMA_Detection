import csv
import tqdm
import sys

#A script to find wheter or not Uniprot ids have a model available in AlphaFoldDB. The file accession_ids.csv need to be downloaded from http://ftp.ebi.ac.uk/pub/databases/alphafold/

Accession_Ids = sys.argv[1] # Pathway to accession_ids.csv
Identifiant_Uniprot = sys.argv[2] # Pathway to a txt file containing the Uniprot ID. One ID per line.

with open (Accession_Ids", 'r') as accession_file :
	accession_reader =  csv.reader(accession_file)
	accession_dict = {}
	for row in accession_reader:
		accession_dict[row[0]] = row[3]

with open ("Correspondance_AF.tsv", 'w') as AF_file :
	with open (Identifiant_Uniprot,'r') as identifiant_file :
		identifiants = identifiant_file.readlines()
		AF_file.write('Id_Sequence\tId_AF_Model\n')	
		for ids in tqdm.tqdm(identifiants, total =len(identifiants)) :
			if ids.find('>') >=0 :
				ids_redu = ids[ids.find('>')+1:].strip() # We remove the >
			else :
				ids_redu = ids.strip()
			AF_file.write(ids_redu+'\t')
			if ids_redu.replace('UniRef100_','') in accession_dict:
				AF_file.write(accession_dict[ids_redu.replace('UniRef100_','')] + '\n')
			else :
				AF_file.write('None\n')
