import os
import sys
import difflib
import logging 

import tqdm
import plotly.figure_factory as ff
from Bio import SeqIO


FORMAT = '%(asctime)s:  %(message)s' 
logging.basicConfig(format=FORMAT,level=logging.INFO, handlers=[logging.StreamHandler()])

# A script to remove the redundancy of HHBLITS output, based on the ID of sequence as well as the sequence identity

def parse_uniprot_id(id):
	id = id.replace('>','').split()[0] # remove description and '>' if any 
	if '|' in id: #UniProt case
		return id.split('|')[1]
	elif id.startswith('UPI'): #UniParc case - note that UniParc sequences are not avalaible in alphafold.
		return id
	else: #UniRef_case
		return id.split('_')[-1]

def write_accession(dict_record,outfile):
	with open(outfile,'w') as fh:
		for k in dict_record:
			k = parse_uniprot_id(k)
			fh.write('{}\n'.format(k))


def record_fusion(a3m_file):
	dict_record = {}
	# imply to read fasta twice...
	record= SeqIO.parse(a3m_file, "fasta")
	lr = len(list(record))
	record= SeqIO.parse(a3m_file, "fasta")
	list_cov_ref = []
	list_cov_seq = []	
	# get generator length 
	for entry in tqdm.tqdm(record,total=lr) : # For each seq ID, the sequence that can be found in the file
		sequence_redu = str(entry.seq).replace('-','').upper()
		id_redu = entry.id.split('/')[0]
		
		if id_redu not in dict_record.keys() :
			dict_record[id_redu]=[sequence_redu]
		else :
			added = False #Was the sequence similar to one that was already described
			for i in range(len(dict_record[id_redu])):
				ref = dict_record[id_redu][i]
				sm = difflib.SequenceMatcher(None, sequence_redu, ref, autojunk=False)
				start_seq, start_ref, length = sm.find_longest_match(0,len(sequence_redu),0,len(ref))
				cov_seq = (start_seq+length)/len(sequence_redu)
				cov_ref = (start_ref+length)/len(ref)
				if cov_seq > 0.5 and cov_ref > 0.5: # If the sequences are identical at least at 50%
					if start_seq == 0 or start_ref == 0: # One sequence has to have its start in the other sequence
						seq_fusion =''
						added = True
						if start_ref > 0 :
							seq_fusion = ref
							if start_ref + length == len(ref) : # check if the other sequence has some bonus on its Nter
								seq_fusion = seq_fusion + sequence_redu[start_seq+length:]
						else :
							seq_fusion = sequence_redu
							if start_seq + length == len(sequence_redu) :
								seq_fusion = seq_fusion + ref[start_ref+length:]
						dict_record[id_redu][i] = seq_fusion
						list_cov_seq.append(cov_seq)
						list_cov_ref.append(cov_ref)	
			if not added :
				dict_record[id_redu]=dict_record[id_redu]+[sequence_redu]
	logging.info("\033[92m {}\033[00m".format('Number of sequences at the end : {}'.format(len(dict_record))))
	return dict_record, list_cov_seq, list_cov_ref		


def plot_coverage(cov_seq,cov_ref,outfile):
	if len(cov_seq) and len(cov_ref):
		fig = ff.create_distplot(
			[cov_seq,cov_ref], 
			["Seq_coverage","Ref_coverage"],
			colors=["green","lightgrey"],
			show_hist=False,
		)	
		fig.write_image(outfile)
	else:
		logging.info('Nothing to plot ....')
	

def write_outfile(dict_record , outfile):
	with open (outfile,'w') as file_reduite :	
		for keys in dict_record : 
			for value in dict_record[keys] :
				file_reduite.write('>'+ keys+'\n')
				file_reduite.write(value+'\n')

if __name__=="__main__":
	a3m_file = sys.argv[1] # Alignment file in a3m format
	outbasename = sys.argv[2]
	
	
	os.makedirs(os.path.abspath(os.path.dirname(outbasename)),exist_ok=True)
	if not os.path.exists(a3m_file):		
		raise FileExistsError('\033[91mAlphafold accession file not found : {}\033[00m'.format(a3m_file))	
	logging.info('\033[92m {}\033[00m'.format('Start multiple hits fusion ...'))
	dict_record,list_cov_seq, list_cov_ref = record_fusion(a3m_file)
	logging.info('\033[92m {}\033[00m'.format('done ...'))
	logging.info('\033[92m {}\033[00m'.format('Write reduced a3M ...'))
	write_outfile(dict_record, outbasename +".a3m")
	logging.info('\033[92m {}\033[00m'.format('done ... [{}]'.format(outbasename+".a3m")))
	logging.info('\033[92m {}\033[00m'.format('Plot coverage ...'))
	plot_coverage(list_cov_seq, list_cov_ref , outbasename +".pdf")
	logging.info('\033[92m {}\033[00m'.format('done ... [{}]'.format(outbasename+".pdf")))
	logging.info('\033[92m {}\033[00m'.format('Write accession file ...'))
	write_accession(dict_record,outbasename+".txt")
	logging.info('\033[92m {}\033[00m'.format('done ... [{}]'.format(outbasename+".txt")))