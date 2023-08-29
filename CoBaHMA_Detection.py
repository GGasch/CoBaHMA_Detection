from Bio import PDB
from Bio import AlignIO
from Bio import SeqIO

import difflib
import statistics as stat
import os
import sys
import pandas as pd
import numpy as np
import plotly.express as px
import re
import tqdm

#Convention : FULL_MAJ = Function, Majthenmin = Global_Variable, full_min = local variable


d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


####pLDDT Computation zone####
def EXTRACT_PLDDT (pdb_name, pdb_file) :
	'''Extract pLDDT value as a list from pdb_file. pdb_name is a variable that we don't really care about, but is needed for the parser'''
	parser = PDB.PDBParser(pdb_name)
	structure =parser.get_structure(pdb_name, pdb_file)
	plddt=[]
	seq=[]
	for model in structure :
		for chain in model :
			for residue in chain:
				seq.append(d3to1[residue.resname]) #We want the sequence in the single letter format
				for atom in residue :
					plddt.append(atom.get_bfactor()) #The bfactor of Alphafold model has the pLDDT value
					break #We only need the first b value since every atom of teh residue have the same b value
	return(''.join(seq),plddt)

def PLDDT_STAT (list_plddt) :
	'''Take a list of list of pLDDT. Each list of pLDDT is the pLDDT for the structure of one sequence. For every position (hence a list of list) and return mean and standard deviation'''
	list_plddt_mean = []
	list_plddt_sd = []
	list_plddt_t = np.array(list_plddt).T.tolist() #We transpose the pLDDT matrix so as each line = a position in the msa
	list_plddt_t_clnd = []
	for lists in list_plddt_t :
		lists_clnd = []
		for val in lists :
			if val != None :
				lists_clnd.append(val)
		list_plddt_t_clnd.append(lists_clnd)
	for pos in range(len(list_plddt_t_clnd)):
		if len(list_plddt_t_clnd[pos]) != 0:
			list_plddt_mean.append(stat.mean(list_plddt_t_clnd[pos]))
		else : 
			list_plddt_mean.append(None) #If there is no information available, put None in the case
		if len(list_plddt_t_clnd[pos]) > 2 :
			list_plddt_sd.append(stat.stdev(list_plddt_t_clnd[pos]))
		else :
			list_plddt_sd.append(None) #If there is less than 2 sequences, no SD can be computed, so none here
	return(list_plddt_mean, list_plddt_sd)

####DSSP Computation zone####
def DSSP (pdb_name, pdb_file) :
	'''get DSSP for a 3D model, same inputs has pLDDT'''
	parser = PDB.PDBParser(pdb_name)
	structure = parser.get_structure(pdb_name, pdb_file)
	model=structure[0]
	dssp = PDB.DSSP(model, pdb_file)
	dssp_l =list(dssp)
	dssp_ll=[]
	for i in range(len(dssp_l)) :
		dssp_ll.append(list(dssp_l[i])[0:3]) # DSSP return a lot of data, but only the first 3 (position, amino acids, 2ndary structure) interest us
	return(dssp_ll)


def DSSP_STAT (dssp_list):
	'''Compute the percentage of alpha helices beta sheet and coil per amino acid. Take a list of dssp. Return a list of 9 data for each position :  H,B,E,G,I,T,S,-,Not_Predicted'''
	dssp_stat=[] # Will create a list of size 8 lists that will contain the proportion of each 2dary structure prediction for each amino acids : [H,B,E,G,I,T,S,-,Not_Predicted]
	dssp_df =  pd.DataFrame(dssp_list)
	for _, row in dssp_df.iterrows():
		dssprediction = [0,0,0,0,0,0,0,0,0]
		if row[2] == 'H' :
			dssprediction[0]=dssprediction[0]+1
		if row[2] == 'B' :
			dssprediction[1]=dssprediction[1]+1
		if row[2] == 'E' :
			dssprediction[2]=dssprediction[2]+1
		if row[2] == 'G' :
			dssprediction[3]=dssprediction[3]+1		
		if row[2] == 'I' :
			dssprediction[4]=dssprediction[4]+1
		if row[2] == 'T' :
			dssprediction[5]=dssprediction[5]+1
		if row[2] == 'S' :
			dssprediction[6]=dssprediction[6]+1
		if row[2] == '-' :
				dssprediction[7]=dssprediction[7]+1
		dssp_stat.append(dssprediction)
	return(dssp_stat)


def DSSP_8to3 (dssp_stat) :
	'''Take a dssp output at 8 state, and return the equivalent at 3 state, plus a column for NoData'''
	dssp_stat_t = np.array(dssp_stat).T.tolist()
	dssp_h=[i*100 for i in dssp_stat_t[0]]
	dssp_b=[i*100 for i in dssp_stat_t[1]]
	dssp_e=[i*100 for i in dssp_stat_t[2]]
	dssp_g=[i*100 for i in dssp_stat_t[3]]
	dssp_i=[i*100 for i in dssp_stat_t[4]]
	dssp_t=[i*100 for i in dssp_stat_t[5]]
	dssp_s=[i*100 for i in dssp_stat_t[6]]
	dssp_c=[i*100 for i in dssp_stat_t[7]]
	dssp_NoData=[i*100 for i in dssp_stat_t[8]]
	
	dssp_hg =np.array(dssp_h)+np.array(dssp_g)+np.array(dssp_i)
	dssp_eg=np.array(dssp_e)+np.array(dssp_b)
	dssp_cg=np.array(dssp_c)+np.array(dssp_s)+np.array(dssp_t)#See Eudes et al. 2007 for the grouping of the secondary prediction

	return(dssp_hg, dssp_eg, dssp_cg, dssp_NoData)

def IDENTIFICATION_BETA0 (dssp_stat):
	'''Take a list of dssp, with a list for each position, whether it is alpha, beta, coil...'''

	dssp_hg, dssp_eg,dssp_cg, dssp_NoData = DSSP_8to3 (dssp_stat)
	mask_2ndaire = []
	counter = 0
	savd_pos = '0'
	structure=''
	for pos1, pos2, pos3 in zip(dssp_hg, dssp_eg, dssp_cg) : #We translate the dssp output into a mask
		if pos1 == 100: #Alpha helices case
			if savd_pos != 'A':
				if counter >= 3	and savd_pos != 'C':	#For coil I have a specific requirement (4a.a minimal)				
					structure=structure+savd_pos
				elif  counter >= 4 and savd_pos == 'C' :
					structure=structure+savd_pos
				counter = 0
			savd_pos = 'A'
			counter=counter+1
			continue
		if pos2 == 100: #Beta sheet case
			if savd_pos != 'B':
				if counter >= 3	and savd_pos != 'C':	#For coil I have a specific requirement (4a.a minimal)				
					structure=structure+savd_pos
				elif  counter >= 4 and savd_pos == 'C' :
					structure=structure+savd_pos
				counter = 0
			savd_pos = 'B'
			counter=counter+1
			continue
		if pos3 == 100: #Coil case
			if savd_pos != 'C':
				if counter >= 3	and savd_pos != 'C':	#For coil I have a specific requirement (4a.a minimal)				
					structure=structure+savd_pos
				counter = 0
			savd_pos = 'C'
			counter=counter+1
			continue
	if counter >= 3	and savd_pos != 'C':	#Add the last position saved in the loop
		structure=structure+savd_pos
	elif counter >= 4 and savd_pos == 'C' :
		structure=structure+savd_pos
	pattern='BBC*A+C*A*C*BA*C*BC*A+'
	patho_pattern_1='BC*A' #some specific HMA match in a bad way, by starting with that pattern
	search=re.search(pattern,structure)
	patho_match = re.match(patho_pattern_1,structure)
	has_beta0 = False
	if search and not patho_match and structure.count('B') < 6: # Some beta barrel have more than 6 beta sheets, it should remove them
		has_beta0 = True
	return(has_beta0, structure)

def IDENTIFICATION_CxxC (seq) :
	c_motif='C.{2}C|C.{3}C'
	if re.search(c_motif, str(seq)):
		return(True)
	else :
		return(False)

####Crop/Gap pLDDT/DSSSP####
def ADAPT2MSA (ref_seq, af_seq, list_ungapped) :
	'''Crop the list_ungapped to match the length of ref_seq and add gap in the list_ungapped to match the gap in the ref_seq'''
	debut_seq,debut_af,length_ali = BORNE_ALI(ref_seq.ungap("-").strip(), af_seq.strip()) # Find where the AF model match the ref_seq
	list_cropd = list_ungapped[debut_af:debut_af+length_ali] # Extract the area that matches the CoBaHMA-like sequence wise
	list_cropd_gapd = GAP(ref_seq, list_cropd) #Add gap in the pLLDT/dssp to match the gap in the sequence in the msa
	return(list_cropd_gapd, debut_af, length_ali)


def BORNE_ALI (seq_short, seq_long):
	'''Find where two sequences matches EXACTLY. This is not mafft, this is purely 100% identity finder. Return : start seq_short, start seq_long, length'''
	s=difflib.SequenceMatcher(None, seq_short.upper(), seq_long.upper(), autojunk=False)
	start_short,start_long,length = s.find_longest_match(0,len(seq_short),0,len(seq_long))
	if length != len(seq_short): # A little check wether or not difflib worked fine	
		return ('Error')

	return(start_short,start_long,length)

def GAP (gapped_seq, list_ref):
	'''Add gap in the list to match the gap in the gapped_seq'''
	gapped_list = []
	pointy = 0
	for pos in range (len(gapped_seq)):
		if gapped_seq[pos] == '-' :
			gapped_list.append(None)
		else :
			gapped_list.append(list_ref[pointy])
			pointy = pointy +1
	return(gapped_list)
############
def UNIREF2AF (name_seq, af_file):
	'''Return the name of the 3D model associated with each sequence, based on an already downloaded file that matches AF model with Uniref100'''
	AF_df = pd.read_csv(af_file, sep='\t')
	AF_name = AF_df[AF_df['Id_Sequence']==name_seq]['Id_AF_Model'].values[0]
	return(AF_name)

####FIGURE####
def FIGURE (dssp_stat, name, has_beta0, base_path, has_cxxc) :
	'''Draw the fig'''
	if has_beta0 and not has_cxxc: # Save the picture in the correct directory
		save_path = base_path+'/CoBaHMA/'
	elif has_beta0 and has_cxxc:
		save_path = base_path+'/CxxC/'			
	else :
		save_path = base_path+'/No_Beta0/'
	if not os.path.exists(save_path):
		os.mkdir(save_path)

	title = 'DSSP prediction by position for sequence'+name
	dssp_hg, dssp_eg, dssp_cg, dssp_NoData =  DSSP_8to3 (dssp_stat) 
	df_dssp2 =pd.DataFrame(list(zip(dssp_hg,dssp_eg,dssp_cg,dssp_NoData)), columns=['%Hg','%Eg','%Cg','%NoData'])

	fig = px.bar(df_dssp2,x=df_dssp2.index,y=df_dssp2['%Hg'])
	#fig.add_bar(name = 'Helix', x=df_dssp2.index,y =df_dssp2['%Hg'], marker_color='lightskyblue')
	fig.add_bar(name='Extended',x=df_dssp2.index,y =df_dssp2['%Eg'], marker_color='lightsalmon')
	fig.add_bar(name='Coil',x=df_dssp2.index,y =df_dssp2['%Cg'], marker_color='lightgoldenrodyellow')
	fig.add_bar(name='No Data',x=df_dssp2.index,y =df_dssp2['%NoData'], marker_color='lightgrey')
	fig.update_layout(barmode='stack')
	
	
	differenciator = 0
	while os.path.isfile(save_path+name+'_'+str(differenciator)+".svg"):
		differenciator = differenciator +1
	fig.write_image(save_path+name+'_'+str(differenciator)+".svg")
	return(name+'_'+str(differenciator))

####WRITING########

def SAVING (base_path, has_beta0, record_clnd, af_model, structure, debut_cob, length_cob, record_id, record_seq, tsv_file, has_cxxc, output_id):
	'''Writing in the excel file and saving the fasta for beta0 structure'''
	with open (tsv_file, 'a') as table :
		if has_beta0 :
			if has_cxxc :
				table.write(record_clnd+'\t'+af_model+'\t'+output_id+'\t'+structure+'\tTrue\tTrue\tFalse\t'+str(debut_cob)+'\t'+str(debut_cob+length_cob)+'\n')
			else :
				table.write(record_clnd+'\t'+af_model+'\t'+output_id+'\t'+structure+'\tTrue\tFalse\tTrue\t'+str(debut_cob)+'\t'+str(debut_cob+length_cob)+'\n')
		else :
			if has_cxxc :
				table.write(record_clnd+'\t'+af_model+'\t'+output_id+'\t'+structure+'\tFalse\tTrue\tFalse\t'+str(debut_cob)+'\t'+str(debut_cob+length_cob)+'\n')
			else :
				table.write(record_clnd+'\t'+af_model+'\t'+output_id+'\t'+structure+'\tFalse\tFalse\tFalse\t'+str(debut_cob)+'\t'+str(debut_cob+length_cob)+'\n')
	if has_beta0 and not has_cxxc : # Save the fasta of sequence that are CoBaHMA
		save_path = base_path+'/CoBaHMA_Fasta/'
		if not os.path.exists(save_path):
			os.mkdir(save_path)
		with open (save_path+output_id+'.fasta', 'w') as fasta_file :
			fasta_file.write ('>'+str(record_id)+'\n'+str(record_seq))

	if has_beta0 and has_cxxc : # Let's keep them on the side, just in case
		save_path = base_path+'/CxxC_Fasta/'
		if not os.path.exists(save_path):
			os.mkdir(save_path)
		with open (save_path+output_id+'.fasta', 'w') as fasta_file :
			fasta_file.write ('>'+str(record_id)+'\n'+str(record_seq))

############################################################Core Code################################################################
#Input : directory with every msa, each labelled Cxxx.fasta and tsv file that make the connection between a Uniref100 and its AF identifier

Fasta_File = sys.argv[1]
Connector_File = sys.argv[2]
Tsv_File = sys.argv[3]

Base_Path = os.getcwd()
with open (Tsv_File, 'w') as Table :
	Table.write('Uniref_ID\tModel_Name\tOutput_ID\t2ndary_Structure\tHas_Beta_0_?\tHas_Cxx(x)C_?\tCoBaHMA\tBeginning_Domain\tEnd_Domain\n')

Fasta = SeqIO.parse(Fasta_File,'fasta')

for Record in tqdm.tqdm(Fasta, total = 34688) :
	Af_Plddt_Cropd_Gapd =''
	Record_Clnd = Record.id
	Seq_Clnd = Record.seq.replace('-','').upper()
	print(Record_Clnd)
	Af_Model = UNIREF2AF(Record_Clnd, Connector_File) # Retrieve the pdb file for the sequence
	if Af_Model != 'None' :
		if os.path.isfile('/run/media/guest-biom/data4/AlphafoldDB/3D_Structures/'+Af_Model+'.pdb') :
			Af_Seq, Af_Plddt = EXTRACT_PLDDT(Af_Model,'/run/media/guest-biom/data4/AlphafoldDB/3D_Structures/'+Af_Model+'.pdb')
			Af_Dssp = DSSP(Af_Model, '/run/media/guest-biom/data4/AlphafoldDB/3D_Structures/'+Af_Model+'.pdb')

			Af_Plddt_Cropd_Gapd, Debut_Cob, Length_Cob = ADAPT2MSA (Seq_Clnd, Af_Seq, Af_Plddt) #Legacy line, useful if we want to re run the code and display plDDT on the picture		
			Af_Dssp_Cropd_Gapd, Debut_Cob, Length_Cob = ADAPT2MSA (Seq_Clnd, Af_Seq, Af_Dssp)
				
			if Af_Plddt_Cropd_Gapd:
				Dssp_Stat = DSSP_STAT(Af_Dssp_Cropd_Gapd)
				Has_Beta0, Structure = IDENTIFICATION_BETA0 (Dssp_Stat)
				Has_CxxC = IDENTIFICATION_CxxC(Seq_Clnd)
				Output_Id = FIGURE(Dssp_Stat, Record_Clnd, Has_Beta0, Base_Path, Has_CxxC)
				SAVING (Base_Path, Has_Beta0, Record_Clnd, Af_Model, Structure, Debut_Cob, Length_Cob, Record.id, Seq_Clnd, Tsv_File, Has_CxxC, Output_Id)
				
			else :
				with open (Tsv_File, 'a') as Table :
					Table.write(Record_Clnd+'\t\t\tNone\n')
		else :
			with open (Tsv_File, 'a') as Table :
				Table.write(Record_Clnd+'\t'+Af_Model+'\t\tModel_Not_Downloaded\n')
	else :
		with open (Tsv_File, 'a') as Table :
			Table.write(Record_Clnd+'\t'+Af_Model+'\n')
