#!/usr/bin/env python3


import os
import sys
import shutil
import argparse
import logging
import pandas as pd
from Bio import SeqIO, SeqRecord


FORMAT = '%(asctime)s:  %(message)s' 
logging.basicConfig(format=FORMAT,level=logging.INFO, handlers=[logging.StreamHandler()])




def fasta_to_length(FASTA):
    records = SeqIO.to_dict( 
        SeqIO.parse(open(FASTA),format="fasta")
    )
    length_dict = {}
    for i,j in records.items() : 
        if "|" not in i: 
            h = i
        else:
            h = i.split("|")[1]
        length_dict[h] = len(j.seq)
    return length_dict

def get_longest_seq_in_cluster(cluster_df, d):
    l = 0
    k = ""
    for i in cluster_df.seq:
        if i in d:
            if d[i] > l:
                l = d[i]
                k = i
        else:
            raise KeyError("fuck {}".format(i))
    return l

def is_ok(cludf, d):
    longest = get_longest_seq_in_cluster(cludf, d = d)
    for i in cludf.seq:
        if i in d:
            l = d[i]            
            if longest-l > MAXDIFF:
                print(longest,l)
                return False            
        else:
            raise KeyError  
    return True
            
def clutable(table):
    df = pd.read_csv(table,sep="\t",header=None)
    df.columns = ["rep","seq"]
    return df


def checktable(table, d , c ):    
    for _, cludf in table.groupby("rep"):
        v = True
        if not is_ok(cludf, d):
            print("threshold {} : fail".format(c) )
            v = False
            break
    if v:    
        print("threshold {} : success".format(c) )
    return v # True if clustering is OK


def get_args():
    parser = argparse.ArgumentParser(
            prog=os.path.basename(__file__),
            description='Get taxonomy from uniprot accession.'
        )

    parser.add_argument(
        '-f',
        '--fasta',
        type = str,
        required=True,
        help="fasta file."
        )
    
    parser.add_argument(
        '-o',
        '--output',        
        type = str,
        default=os.getcwd(),
        help="output directory."
    )

    parser.add_argument(
        '-c',
        '--coverage',        
        type = int,
        default=90,
        help="minimum coverage."
    )

    parser.add_argument(
        '-m',
        '--max-res',        
        type = int,
        default=50,
        help="amplitude maximum of sequence length in a cluster."
    )


    args = parser.parse_args()
    return args


# python3 full-seq-clustering ../../datas/

if __name__ == "__main__":
    args = get_args()
    FASTA = args.fasta
    OUTDIR = args.output
    COVERAGE = [c/100 for c in range(int(args.coverage),100)]
    MMSEQS2 = shutil.which("mmseqs")
    MAXDIFF = int(args.max_res)

    if not MMSEQS2:
        logging.error("\033[91m {}\033[00m" .format('Error - mmseqs not found in path'))
        sys.exit(-1)
    else:
        logging.info('MMSEQS FOUND : {}'.format(MMSEQS2))

    # CREATE SEQ DB
    os.makedirs(OUTDIR)
    SEQDB = os.path.join(OUTDIR,"seqDB")
    cmd = "{} createdb {} {}".format(
        MMSEQS2,
        FASTA,
        SEQDB,
    )
    os.system(cmd)
    length_dict = fasta_to_length(FASTA)
    for cov in COVERAGE:
        logging.info("\033[92m {}\033[00m".format("Clustering at {} ...".format(cov)))
        itedir = os.path.join(OUTDIR , "clustering_cov_{}".format(cov) )
        CLUDB = os.path.join(itedir,"cluDB")
        TMP =  os.path.join(itedir,"tmp")
        
        if not os.path.isdir(itedir):
            os.makedirs(TMP)
            cmd = "{exe} cluster {seqdb} {cludb} {tmp}  -v 0 -c {cov} -s 7.5 -a --cluster-reassign && {exe} createtsv -v 0 {seqdb} {seqdb} {cludb} {cludb}.tsv ".format(
                exe = MMSEQS2,
                cov = cov,
                seqdb = SEQDB,
                cludb = CLUDB,
                tmp = TMP
            )
            os.system(cmd)
        if os.path.exists(CLUDB+".tsv"):
            cludf = clutable(CLUDB+".tsv")  #pd.read_csv(CLUDB+".tsv",sep="\t",header=None,index_col=None)
            valid = checktable(cludf,length_dict , cov)
            if not valid:
                shutil.rmtree(itedir, ignore_errors=False, onerror=None)
                logging.info("\033[93m {}\033[00m".format("clustering with c={}  do not respect AMP=={}".format(cov,MAXDIFF)))
            else:                
                logging.info("\033[92m {}\033[00m".format("clustering with c={} is valid - AMP<{}".format(cov,MAXDIFF)))
        logging.info("\033[92m {}\033[00m".format("Clustering done."))

    # shutil.rmtree(OUTDIR, ignore_errors=False, onerror=None)
            