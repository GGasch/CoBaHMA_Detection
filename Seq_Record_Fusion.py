from Bio import SeqIO
import difflib
import plotly.express as px
import tqdm
import sys

# A script to remove the redundancy of HHBLITS output, based on the ID of sequence as well as the sequence identity

A3m_File = sys.argv[1] #The a3m file to be reduced

Dict_record = {}
record= SeqIO.parse(A3m_File, "fasta")
List_Cov = []

for entry in tqdm.tqdm(record, total = 18000000) : # For each seq ID, the sequence that can be found in the file

	sequence = str(entry.seq)
	if entry.id.find('/') >0 :
		id_redu = entry.id[:entry.id.find('/')]
	else :
		id_redu = entry.id
	sequence_redu=sequence.replace('-','').upper()

	if id_redu not in Dict_record.keys() :
		Dict_record[id_redu]=[sequence_redu]
	else :
		added = False #Was the sequence similar to one that was already described
		for i in range(len(Dict_record[id_redu])):
			ref = Dict_record[id_redu][i]
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
					Dict_record[id_redu][i] = seq_fusion
					List_Cov.append(cov_seq)
					List_Cov.append(cov_ref)	
		if not added :
			Dict_record[id_redu]=Dict_record[id_redu]+[sequence_redu]
			
with open ('Iteration1_redundancy_removed.fasta','w') as file_reduite :	
	for keys in Dict_record : 
		for value in Dict_record[keys] :
			file_reduite.write('>'+ keys+'\n')
			file_reduite.write(value+'\n')

fig = px.histogram(List_Cov)
fig.show()
