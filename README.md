# CoBaHMA_Detection
In order to find which UniProt sequence has a prediction available for its 3D structure in AlphaFoldDB, please download accession_ids.csv from http://ftp.ebi.ac.uk/pub/databases/alphafold/

_Seq_Record_Fusion.py : a script to reduce an a3m file to remove redundant entry in it

_alphafold_ids.py : a script to find whether or not a list of Uniprot ID has a model available in AlphaFoldDB, and download the name of the model. This script require accession_ids.csv.

_alphafold_download.py : a script to download model of protein 3D structure from AlphaFoldDB (https://alphafold.ebi.ac.uk/). This code is made for AlphaFoldDB V4, in case a new version of AlphaFoldDB is released, the code needs to be updated.
