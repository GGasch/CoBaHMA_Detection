#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import urllib
import logging
import requests
import time
import concurrent.futures

import tqdm

FORMAT = '%(asctime)s:  %(message)s' 
logging.basicConfig(format=FORMAT,level=logging.INFO, handlers=[logging.StreamHandler()])


def search_in_uniprot_save(accession):    
    acc = accession.split("_")[-1]
    url = "https://rest.uniprot.org/unisave/{}".format(acc)   
    r = requests.get(url)
    if r.status_code == 200: 
        # try:
            if "results" in r.json():            
                version = r.json()["results"][0]['entryVersion']
                accession = r.json()["results"][0]['accession']                               
                url = "https://rest.uniprot.org/unisave/{}?format=fasta&versions={}".format(accession,version) 
                c = requests.get(url)
                if c.status_code == 200:   
                    return c._content.decode()                                  
        # except:
            return accession
    return accession
    
def make_uniprot_url(accession):
    """
    https://rest.uniprot.org/uniprotkb/W2RNG4.fasta
    https://rest.uniprot.org/uniref/UniRef100_W7YBH4.fasta
    https://rest.uniprot.org/uniparc/UPI001F54928F.fasta
    """
    if accession.startswith("Uni"):
        if "UPI" in accession:
            if accession.split("_")[-1].startswith("UPI"):
                accession = accession.split("_")[-1]
                baseurl = "https://rest.uniprot.org/uniparc/"
            else:
                baseurl = "https://rest.uniprot.org/uniref/"
        else:
            baseurl = "https://rest.uniprot.org/uniref/"
    else:
        baseurl = "https://rest.uniprot.org/uniprotkb/"
    return baseurl + accession + ".fasta"

def get_args():
    parser = argparse.ArgumentParser(
            prog=os.path.basename(__file__),
            description='Get sequences in fasta format from uniprot.com.'
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
        default="uniprot.fasta",
        help="output fasta file"
    )

    args = parser.parse_args()
    return args


    
def multithreads_request(URLS):
    start = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm.tqdm(executor.map(get_fasta_content, URLS), total=len(URLS)))
    finish = time.perf_counter()
    logging.info(f'Finished in {round(finish-start, 2)} second(s)')
    return results

def get_fasta_content(acc):  
    url = make_uniprot_url(acc)                
    try:
        #filename = wget.download( url, out=outdir , bar=None)    
        c = requests.get(url)
        if c.status_code == 200:                                       
            return c._content.decode()
        else:
            return search_in_uniprot_save(acc)
    except urllib.error.HTTPError:                    
        return None

def download_uniprot_sequences(l_accession):    
    records = []
    errors = []
    results = multithreads_request(l_accession)
    records = [i for i in results if i.startswith('>')]
    errors = [i for i in results if not i.startswith('>')]
    del results
    return records,errors
    

def write_fasta(records, outfile):
    with open(outfile,'w+') as streamout:
        streamout.write("\n".join(records))
  
def write_errors(errors,outfile):
    with open(outfile,'w+') as streamout:
        for e in errors:
            streamout.write('{}\n'.format(e)) 


if __name__ == "__main__":
    args = get_args()
    l_acc = []
    logging.info('\033[92m {}\033[00m'.format('Get accession from file ...'))
    with open(args.INPUT , 'r') as stream:
        for i in stream.readlines():        
            l_acc.append(i.strip())
    l_acc = list(set(l_acc))
    logging.info('\033[92m {}\033[00m'.format('done ...'))
    logging.info('\033[92m {}\033[00m'.format('Download {} records ...'.format(len(l_acc))))
    records,errors  = download_uniprot_sequences(l_acc)
    logging.info('\033[92m {}\033[00m'.format('done ...'))
    logging.info('\033[92m {}\033[00m'.format('Number of records actually downloaded : {}'.format(len(records))))
    logging.warning('\033[93m {}\033[00m'.format('Number of missing records : {}'.format(len(errors))))
    logging.info('\033[92m {}\033[00m'.format('Write records to fasta and errors to log ...'))
    write_fasta(records,args.OUTPUT)
    write_errors(errors,args.OUTPUT+".log")
    logging.info('\033[92m {}\033[00m'.format("done."))
    
    