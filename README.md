# CoBaHMA_Detection

---

Seq_Record_Fusion.py : a script to reduce an a3m file to remove redundant entry in it

Usage: `python3 Seq_Record_Fusion.py <a3m> <path/basename>`

where <out_basename> is a path without extension.

---

In order to find which UniProt sequence has a prediction available for its 3D structure in AlphaFoldDB, please download [*accession_ids.csv*](http://ftp.ebi.ac.uk/pub/databases/alphafold/).

alphafold_ids.py : a script to find whether or not a list of Uniprot ID has a model available in AlphaFoldDB, and download the name of the model. This script require *accession_ids.csv*.

Usage: `python3 alphafold_ids.py accession_ids.csv <acc_txt_file>`

where <acc_txt_file> is a text file containing uniprot accessions, one per line.

---

alphafold_download.py : a script to download model of protein 3D structure from [AlphaFoldDB](https://alphafold.ebi.ac.uk/). This code is made for AlphaFoldDB V4, in case a new version of AlphaFoldDB is released, the code needs to be updated.

Usage: `python3 alphafold_download.py <tabular_file> <output_dir>`

where <tabular_file> is a tabular file with two columns Id_Sequence and Id_AF_Model with Id_Sequence stand for UniProt accessions and Id_AF_Model are the corresponding Alphafold accessions.

---

CoBaHMA_Detection.py : Based on the 3D model and the sequence, detect wheter or not a sequence is a CoBaHMA
This scipt use [DSSP](https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html) as part of the Bio library
To use this module, you need a local instalation of DSSP, available at https://swift.cmbi.umcn.nl/gv/dssp/ or with conda : https://anaconda.org/salilab/dssp

Usage: `python3 CoBaHMA_Detection.py <fasta_file> <3D_models_dir> <output_dir>`

where :
- <fasta_file> A file with the sequences you want to annotate with the help of their 3D model.
- <3D_models_dir> A directory containing corresponding AF-\<accession\>-F1.pdb files.
- <output_dir> A directory where output files will be stored.

## Sequence classification

retrieve_sequence_from_uniprot.py: Download full length sequence from UniProt API based on accession.

Usage: `python3 retrieve_sequence_from_uniprot.py -i accessions.txt -o uniprot.fasta` 

---

clustering_ASL_constrained.py: Use mmseqs2 to cluster a fasta file for a range of threshold with a constraint based on sequence length amplitude (ASL).

Usage: `python3 clustering_ASL_constrained.py -f uniprot.fasta -o uniprot_clustering -c 90 -m 50`

where -c 90 is the starting point for the range of coverage [90:100] and -m 50 is the maxium sequence length amplitude for each cluster.

---

Intermediates steps using mmseqs2:

0) convert clustering results into tabular file:

```
mmseqs convertalis \
  clustering/clustering_cov_0.97/cluDB clustering/clustering_cov_0.97/cluDB \
  clustering/seqDB \
  clustering/clustering_cov_0.97/clualn.tsv \
  --format-output query,target,evalue,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qaln,taln,mismatch,qcov,tcov
```

1) extract representative sequences from clustering:

```
mmseqs createsubdb \
  clustering/clustering_cov_0.97/cluDB \
  clustering/seqDB \
  clustering/clustering_cov_0.97/repDB
```

2) perform a self-vs-self search with representative sequences:

```
mmseqs search \
  clustering/clustering_cov_0.97/repDB \
  clustering/clustering_cov_0.97/repDB  \
  clustering/clustering_cov_0.97/searchDB \
  tmp -a
```

3) convert search results in a tabular file:

```
mmseqs convertalis \
  clustering/clustering_cov_0.97/repDB \
  clustering/clustering_cov_0.97/repDB \
  clustering/clustering_cov_0.97/searchDB \
  clustering/clustering_cov_0.97/search.tsv \
  --format-output query,target,evalue,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qaln,taln,mismatch,qcov,tcov
```

4) Filter clustering and search alignments.

---

community_definition.py: Build a network from a tabular file with at least two columns : query and target. Louvain community detection algorithm will be used to detect group of nodes. Community will be refined by removing edges between unvalid and valid nodes (unvalid nodes are nodes with the sequence size farthest from the median sequence size within the community) until the sequence length amplitude of the community is lower than ASL (here ASL=50).

Usage: `python3 community_definition.py -t table.tsv -f corresponding.fasta -o output_dir -w weight_column -m 50 -s 5`

where :
- table.tsv - A tabular file containing node source, node target and optionaly a weight.
- corresponding.fasta - Fasta file containing nodes' sequences.
- output_dir - Where results will be stored
- weight_column - Name of a column which will be used to weight edges.
- -m 50 - maximum amplitude sequence length value (ASL)
- -s 5 - minimum number of nodes within a community to consider it as a 'big' community (worth seeing in details).

# Taxonomic, functionnal and structural annotations

retrieve_taxonomy.py: Retrieve sequence taxonomy from the UniProt API based on accession.

Usage: `python3 retrieve_taxonomy.py -i accessions.txt -o uniprot-taxonomy.tsv` 

Note that you can translate the taxID into a NCBI-like taxonomic lineage using [TaxonKit](https://github.com/shenwei356/taxonkit).


---

We use interproscan to annotate sequences with protein and domain families.

`interproscan.sh -t p -i uniprot.fasta -f TSV -T tmpdir -dp -b interproscan
`

---

We also use domainMapper and the ECOD  classification to detect non-contiguous, insertional and circularly permetud domains. See the [domainMAPPER](https://github.com/FriedLabJHU/DomainMapper) page for usage.

---

We use [deepTMHMM](https://dtu.biolib.com/DeepTMHMM) to annotate transmembrane regions.

---

Finally, we use [pCALF](https://github.com/K2SOHIGH/pcalf/tree/main) with default HMM profiles to detect and annotate calcyanin proteins. 



# Supplementary datas: 

Il y a 1) les fasta/msa/pdb du representant par communauté (sup data 3) et les orga modulaire par communauté (sup data 4). Du coup, pourrais tu inclure les fichiers suivants (il y a aussi tes deux fichiers tsv) dans un repertoire dédié (nommé par exemple SupData) sur GitHub, je donnerai alors comme référence le lien GitHub dans les suppléments et pourrai continuer la soumission du papier ainsi actualisé. Un grand merci pour ton aide! Et il y aurait aussi la table SI-1
      
- `SupplementaryData3_FASTA_MSA_AF2.zip`Archive including fasta, msa and pdb files for all representative sequences (community)
- `SupplementaryData4_ModularOrganization.zip` Archive containing the modular organization of each community.
- `SupplementaryData10_Community_desc.tsv` Table describing each community with:
  - Community (Community identifier)
  - Community Rep (The sequence representative of the community)
  - Rep degree (The number of connexion involving the representative sequence)
  - Community Size (The size of the community, i.e number of nodes)
  - subgraph (The subgraph which the community belong)
  - mean_seq_len (The mean sequence length in aa within the community)
  - longest_seqs (The longest sequence of the community)
  - longest_len (The maximum sequence length in aa within the community)
  - shortest_seqs (The shortest sequence of the community)
  - shortest_len (The minimum sequence length in aa within the community)
  - IPRS (List of IPRs ID found associated to the community)
  - IPRS_Functions (List of IPRs function associated to the community)
  - DomainMapper_Ecod (List of ECOD associated to the community)
  - DomainMapper_Functions (List of ECOD function associated to the community)
- `SupplementaryData11_cobahma_annotation_taxo.tsv` Table summarizing the taxonomic lineage and the functionnal/structural annotation of all CoBaHMA sequences. One functionnal/structural definition per row.
  - Accession (Uniprot accession) 
  - E-Value (E-value of the functionnal/structural annotation)
  - Residue Range (Feature start and stop-
  - Property (Domain Property (NC, IS, CP), see [DomainMapper](https://github.com/FriedLabJHU/DomainMapper) for details)
  - Architecture (ECOD Architecture, see [DomainMapper](https://github.com/FriedLabJHU/DomainMapper) for details)
  - X-group (ECOD X-Group, see [DomainMapper](https://github.com/FriedLabJHU/DomainMapper) for details)
  - T-group (ECOD T-Group, see [DomainMapper](https://github.com/FriedLabJHU/DomainMapper) for details)
  - Desc (Functional description)
  - Code (ECOD/IPR codes)
  - Tool (Tool used for detection)
  - Community (The community which the sequence belong)
  - Community Size (The community size)
  - Community Rep (The representative sequence of the community)
  - Seq length (aa) (The sequence lenfth in residus)
  - Taxonomy (The phylum which the sequence belong - less representated phyla are grouped under the label 'Other')
  - Taxonomy detail (The phylum which the sequence belong)
 








