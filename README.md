# CoBaHMA_Detection

**HERE ADD DETAILS FOR transitive search with HHBLIT...**

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
mmseqs convertalis clustering/clustering_cov_0.97/cluDB clustering/clustering_cov_0.97/cluDB clustering/seqDB clustering/clustering_cov_0.97/clualn.tsv --format-output query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid
```

1) extract representative sequences from clustering:

```
mmseqs createsubdb clustering/clustering_cov_0.97/cluDB clustering/seqDB clustering/clustering_cov_0.97/repDB
```

2) perform a self-vs-self search with representative sequences:

```
mmseqs search clustering/clustering_cov_0.97/repDB clustering/clustering_cov_0.97/repDB  clustering/clustering_cov_0.97/searchDB tmp -a
```

3) convert search results in a tabular file:

```
mmseqs convertalis clustering/clustering_cov_0.97/repDB clustering/clustering_cov_0.97/repDB clustering/clustering_cov_0.97/searchDB clustering/clustering_cov_0.97/search.tsv --format-output query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid
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