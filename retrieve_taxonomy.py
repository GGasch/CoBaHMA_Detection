import os
import argparse
import logging
import urllib
import time
import concurrent.futures

import pandas as pd
from rdflib import Graph
from SPARQLWrapper import SPARQLWrapper, JSON
import tqdm


FORMAT = '%(asctime)s:  %(message)s' 
logging.basicConfig(format=FORMAT,level=logging.INFO, handlers=[logging.StreamHandler()])


def taxid(protid,bdd):    
    try:
        g = Graph()
        g.parse("https://www.uniprot.org/{}/{}.ttl".format(bdd,protid))        
        o = g.serialize(format="turtle")#.decode("utf-8")
        for line in o.split("\n"):
            #up:organism taxon:
            if ":organism" in line and "taxon:" in line:
                for i in line.split():
                    if "taxon" in i:                        
                        return i.split(":")[-1]
        return None
    except urllib.error.HTTPError:
        return None

def lineage(taxid):
    link = "https://www.uniprot.org/taxonomy/{}.tsv".format(taxid)
    c=pd.read_csv(link,sep="\t")
    l = c.Lineage.values[0].split(", ")
    return ";".join(l[::-1])
    
def bdd(protid):
    if "UPI" in protid:
        if "UniRef" in protid:
            return protid,"uniref"
        else:
            return "UniRef100_{}".format(protid), "uniref"
    else:
        if "UniRef" in protid:
            return "_".join(protid.split("_")[1:]) ,"uniprot"
        else:
            return protid, "uniprot"

def retrieve_lineage(protid):
    uid,bd = bdd(protid)
    taxon = taxid(uid,bd)
    l = None
    if taxon:
        l = lineage(taxon)
    
    return uid,bd,taxon,l



def get_args():
    parser = argparse.ArgumentParser(
            prog=os.path.basename(__file__),
            description='Get taxonomy from uniprot accession.'
        )

    parser.add_argument(
        '-i',
        '--input',
        dest="INPUT",
        type = str,
        required=True,
        help="text file containing uniprot sequence accession. One per line."
        )
    
    parser.add_argument(
        '-o',
        '--output',
        dest="OUTPUT",
        type = str,
        default="taxo.tsv",
        help="output table file"
    )

    args = parser.parse_args()
    return args


def multithreads_request(ACCS):
    start = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm.tqdm(executor.map(retrieve_lineage, ACCS), total=len(ACCS)))
    finish = time.perf_counter()
    logging.info(f'Finished in {round(finish-start, 2)} second(s)')
    return results

def main():
    args = get_args()
    file = args.INPUT
    l_acc = []
    # get all accession from input file
    logging.info('\033[92m {}\033[00m'.format('Get accessions from file ...'))
    with open(file,'r') as stream:
        for line in stream.readlines():
            l_acc.append(line.strip())
    l_acc = list(set(l_acc))            
    logging.info('\033[92m {}\033[00m'.format('done ...'))

    logging.info('\033[92m {}\033[00m'.format('Fetch UNIPROT for taxIDs and lineages ...'))
    datas = multithreads_request(l_acc)  
    logging.info('\033[92m {}\033[00m'.format('done ...'))
    taxodf = pd.DataFrame(datas)
    taxodf.columns = ["seqid","db","taxid","lineage"]
    taxodf[taxodf.taxid.isna()]
    taxodf.set_index('seqid',inplace=True)
    outdir = os.path.dirname(os.path.abspath(args.OUTPUT))

    os.makedirs(outdir,exist_ok=True)
    logging.info('\033[92m {}\033[00m'.format('Save file to {} ...').format(args.OUTPUT))
    taxodf.taxid = taxodf.apply(lambda x : 0 if x.taxid == "-" else x.taxid , axis=1)
    taxodf.to_csv(args.OUTPUT,sep="\t",header=True,index=True)
    logging.info('\033[92m {}\033[00m'.format('done.'))
    

if __name__=="__main__":
    main()