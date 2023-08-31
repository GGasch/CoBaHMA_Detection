import logging
import os
import argparse
import json

import pandas as pd
import numpy as np

import networkx as nwx
import networkx.algorithms.community as nxcommu

import tqdm

from Bio import SeqIO, SeqRecord

FORMAT = '%(asctime)s:  %(message)s' 
logging.basicConfig(format=FORMAT,level=logging.INFO, handlers=[logging.StreamHandler()])

def yellowMSG(msg): return "\033[93m {}\033[00m" .format(msg)
def greenMSG(msg): return "\033[92m {}\033[00m" .format(msg)
def redMSG(msg): return "\033[91m {}\033[00m" .format(msg)




def louvain_communities_to_subgraph(G,partition):
    communities = {}
    for p in range(1,len(partition)+1):    
        cid="C{}".format(p)
        communities[cid] = list(partition[p-1])
    return communities

def estimate_partition(graph,louvain,big_co_threshold = 5):
    coverage,perfo = nxcommu.partition_quality(graph,louvain)
    sumseq = 0
    totalseq = 0
    for c in louvain:
        totalseq+=len(c)
        if len(c)>=big_co_threshold:
            sumseq+=len(c)
    return round(coverage,2), round(perfo), len(louvain) , len([i for i in louvain if len(i)>=big_co_threshold]) , sumseq, round(sumseq/totalseq*100,2)

def get_seq_mean_length(lseq , length_dict ):
    lens = []
    for seq in lseq:
        if seq not in length_dict:
            raise KeyError("{} can't find in your sequence length directory".format(seq))
        else:
            lens.append(length_dict[seq])
    if len(lens) > 1:
        return np.mean(lens),np.sem(lens)                
    return lens[0],0

def remove_until_homogeneity(seq_dict,NRES=50):
    llen=list(seq_dict.values())
    median = np.median(llen)
    median_deviation = []
    # 3.2.2 we need to remove the furthest sequence of the median 
    # until the community is valid according to the max_diff threshold
    #score,pvalue = scipy.stats.shapiro(list_seq_len)             
    for seq,length in seq_dict.items(): # for each sequence in the community 
        median_deviation.append( abs(median-length))            
    sorted_sequences = [(_,x) for _,x in sorted(zip(median_deviation,list(seq_dict.keys())))]
    removed = []
    amplitude = max(llen)-min(llen)
    while amplitude >= NRES: 
        _ , node_removed = sorted_sequences.pop()
        removed.append(
            node_removed
        )
        seq_dict.pop(node_removed)        
        amplitude = max(seq_dict.values())-min(seq_dict.values())
    return removed,list(seq_dict.keys())



def make_raw_network_from_table(df,weight=None,seed=13, big_communnity_threshold = 5):
    PERFO = []

    edges = []
    graph = nwx.Graph() 
    
    for _ , row  in tqdm.tqdm(df.iterrows()):
        n1 = row["query"]
        n2 = row["target"]    
        for node in [n1,n2]:        
            if node and node not in graph.nodes():
                graph.add_node(
                    node
                )
                
        if n1 and n2 and n1!=n2:
            attr = {}
            if weight:
                attr["weight"] = row[weight]

            edges.append(
                (n1, n2, attr)
            )   
    graph.add_edges_from(edges)

    louv = nxcommu.louvain_communities(graph,seed=seed)
    sumseq = 0
    totalseq = 0
    big_co_threshold = 5
    for c in louv:
        totalseq+=len(c)
        if len(c)>=big_co_threshold:
            sumseq+=len(c)

    # Estimate Louvain partition quality.    
    PERFO.append( ('Raw',)+estimate_partition(graph,louv, big_communnity_threshold))
    return graph, PERFO



def refine_network(raw_graph,length_dict,seed=13,NRES=50,big_communnity_threshold=5):
    PERFO = []
    graph = raw_graph.copy()
    iteration = 0
    convergence_not_reached = True

    while convergence_not_reached: # and iteration < max_ite and no_res > 0: 
        #while edges have been removed or commmunity with a diff > 50 exists
        iteration += 1

        # 1) get partition
        louv = nxcommu.louvain_communities(graph,seed=seed)
        communities =  louvain_communities_to_subgraph(graph,louv)    

        # 2) estimate perfo and keep it
        PERFO.append(("Ite_{}".format(iteration),)+estimate_partition(graph,louv,big_communnity_threshold))
        
        # 3) for each community get mean sequence length and remove edge that 
        #do not respect a threshold of (mean_len +/- sem_len) +/- No of residues
        convergence_not_reached=False 
        
        for community, sequences in communities.items():
            #print('community : ', community)
            # 3.1) get length of sequences for this community
            seq_length_dict = {s:length_dict[s] for s in sequences}
            # 3.2) get amplitude
            amplitude = max(list(seq_length_dict.values())) - min(list(seq_length_dict.values()))
            if amplitude <= NRES:
                #if len(seq_length_dict)>1:
                #    valid_commu_at_first += 1
                #else:
                #    valid_singl_at_first +=1
                # 3.2.1 nothing to do here because the community is already good 
                continue            
            else:
                #mean_len, sem_len = get_seq_mean_length(sequences,length_dict)   
                sd_len = np.std(list(seq_length_dict.values()))
                
                unvalid_nodes,valid_nodes = remove_until_homogeneity(seq_length_dict)
                
                for uvnode in unvalid_nodes:
                    for vnode in valid_nodes:
                        if graph.has_edge(uvnode, vnode):                            
                            graph.remove_edge(uvnode,vnode)
                            convergence_not_reached=True   
                        if graph.has_edge(vnode, uvnode):
                            graph.remove_edge(vnode,uvnode)
                            convergence_not_reached=True   
                
                  
    nodes_and_community = {}
    for c,seqs in communities.items():                            
        for s in seqs:
            nodes_and_community[s]=c
    return communities,graph,nodes_and_community,PERFO


def dump_communities(graph,communities,fasta_records, outdir):    
    os.makedirs(os.path.join(outdir, "community-graphs") , exist_ok=True )
    os.makedirs(os.path.join(outdir, "community-fastas") , exist_ok=True )
    
    # Extract subgraph and fasta for each community."
    for cid , community_members in communities.items():            
        if len(community_members)>=0:
            subgraph = graph.subgraph(community_members)
            subfasta = {}
            for nodes in subgraph.nodes():
                if nodes in fasta_records:
                    subfasta[nodes] = fasta_records[nodes]     
            nwx.write_gexf(
                subgraph, 
                os.path.join(outdir, "community-graphs","{}.gexf".format(cid))
            )
            SeqIO.write(
                list(subfasta.values()), 
                open(os.path.join(outdir, "community-graphs", "{}.fasta".format(cid)),'w'), 
                format="fasta"
            )


def get_community_seq_length(community,communities,length_dict):
    """
    return a tuple compose of mean_len, longest_len , longest_seq, shortest_len, shortest_seq
    """
    llen = []
    lid  = [] 
    for seq in communities[community]:
        llen.append(length_dict[seq])
        lid.append(seq)
    meanlen = np.mean(llen)
    shortest = []
    shortest_len = 10000000000
    longest = []
    longest_len = 0
    for seq,lens in zip(lid,llen):
        if lens <= shortest_len:
            if lens < shortest_len:
                shortest = []
            shortest_len = lens
            shortest.append(seq)
        if lens >= longest_len:
            if lens > longest_len:
                longest = []
            longest_len = lens
            longest.append(seq)        

    return (
        community,
        meanlen,
        ";".join(longest),
        longest_len,
        ";".join(shortest),
        shortest_len)

def get_by_subgraph(graph):
    subgraphs = [graph.subgraph(c).copy() for c in nwx.connected_components(graph)]
    subgraphs = sorted(subgraphs, key = lambda x : len(x.nodes()) , reverse=True)
    si = 0
    subgraph_map_dict = {}
    for sg in subgraphs:
        si+=1
        for node in sg.nodes():
            subgraph_map_dict[node] = "Subgraph_{}".format(si)
    return subgraph_map_dict

def make_community_table(communities , graph , length_dict, outdir):
    d = []
    subgraph_map_dict = get_by_subgraph(graph)
    
    graph_seqs_degree = dict(graph.degree())
    community_rep_map = {}
    for cid, cseqs in communities.items() : 
        max_degree = -1
        rep = None
        for cseq in cseqs:
            if graph_seqs_degree[cseq] > max_degree:
                rep = cseq
                max_degree = graph_seqs_degree[cseq]
        community_rep_map[cid] = {"Community Rep" : rep, "Rep degree":max_degree}

    for community in communities:
        d.append(get_community_seq_length(community,communities,length_dict))

    community_seqlen_df = pd.DataFrame(d,columns=["Community","mean_seq_len","longest_seqs","longest_len","shortest_seqs","shortest_len"])
    community_seqlen_df.set_index("Community",inplace=True)


    community_table = pd.DataFrame(community_rep_map).T
    community_table["Community Size"] = community_table.apply(lambda x:
                                                            len(communities[x.name]),
                                                            axis=1)
    community_table["subgraph"] = community_table.apply(lambda x : subgraph_map_dict[x["Community Rep"]],axis=1)
    community_table["Rep Degree"] = community_table.apply(lambda x : graph_seqs_degree[x["Community Rep"]],axis=1)
    community_table = pd.concat([community_table,community_seqlen_df],axis=1)

    community_table.index.name="Community"    
    community_table = community_table[['Community Rep', 'Rep degree', 'Community Size', 'subgraph',
        'Rep Degree', 'mean_seq_len', 'longest_seqs', 'longest_len',
        'shortest_seqs', 'shortest_len', ]]
    community_table.to_csv(os.path.join(outdir,"community_table.tsv"),
                        sep="\t",header=True,index=True)

                  
    
    
def get_args():
    parser = argparse.ArgumentParser(
            prog=os.path.basename(__file__),
            description='Create network and identify communities with sequence length constrain.'
        )

    parser.add_argument(
        '-t',
        '--table',
        type = str,
        required=True,
        help="Tabular file with at least two columns named query and target." 
        )

    parser.add_argument(
        '-f',
        '--fasta',
        type = str,
        required=True,
        help="Fasta file used for generating the tabular file (query and target columns must have a sequence in it)."
        )   
    
    parser.add_argument(
        '-w',
        '--weight',
        type = str,
        default= None,
        help="Column name to use as weight for edges."
        )  

    parser.add_argument(
        '-m',
        '--max-res',
        type = int,
        default= 50,
        help="The maximum amplitude of sequence length within a community."
        )  

    parser.add_argument(
        '-s',
        '--size',
        type = int,
        default= 5,
        help="The minimum number of nodes in a community to consider it as a big community."
        )  

    parser.add_argument(
        '-o',
        '--output',
        required=True,
        type = str,        
        help="Output directory."
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    fasta = args.fasta
    max_amp = args.max_res
    outdir = args.output

    logging.info(greenMSG('Parse fasta file ... '))
    records =  {}
    for record in SeqIO.parse(open(fasta),format="fasta"):
        if record.id not in records:
            records[record.id] = record

    length_dict = {i: len(j.seq) for i,j in records.items()}        
    logging.info(greenMSG('Number of sequences : {} '.format(len(records))))

    logging.info(greenMSG('Parse tabular file ... '))
    table_df = pd.read_csv(args.table,sep='\t',header=0)
    cols = table_df.columns 
    if args.weight and args.weight not in cols:
        logging.error(redMSG('"{}" is not a column name ... [{}]'.format(args.weight,";".join(cols)) ))
        exit(-1)

    if 'query' not in cols or 'target' not in cols:
        logging.error(redMSG("query and/or target columns are missing ..."))
        exit(-1)

    unique_ids = list(set(table_df['query'].unique().tolist() + table_df['target'].unique().tolist()))
    if len(unique_ids) > len(records):
        logging.error('Number of unique identifiers ({}) in table is greater that the number of sequences {} in the fasta file ... '.format(
            len(unique_ids),len(records)
            )
        )
        exit(-1)

    logging.info(greenMSG('Make network ... '))
    raw_graph, raw_perfo = make_raw_network_from_table(
        table_df, 
        args.weight, 
        big_communnity_threshold=args.size
    )
    logging.info(greenMSG('Done ... ({}% nodes recruited)'.format(raw_perfo[0][-1])))
    
    logging.info(greenMSG('Refine network with AMP = {} ... '.format(max_amp)))
    refined_communities, refined_graph, refined_nodes_and_community, refined_perfo = refine_network(
        raw_graph, 
        length_dict, 
        NRES = max_amp, 
        big_communnity_threshold=args.size
    )
    
    logging.info(greenMSG('Done ... (last iteration : {}% nodes recruited)'.format(refined_perfo[-1][-1])))
    perfo_df = pd.DataFrame(raw_perfo + refined_perfo)
    perfo_df.columns = ["iteration","coverage","performance","#community","#big_community","#seqinclusters","%recruitment"]

    # write output
    logging.info(greenMSG('Create output files ... '))
    os.makedirs(outdir, exist_ok=True )
    perfo_df.to_csv(os.path.join(outdir,"partition_metrics.tsv"),sep='\t',header=True)
    json.dump(refined_communities,open(os.path.join(outdir,"communities.json"),'w'))
    json.dump(refined_nodes_and_community,open(os.path.join(outdir,"map_community_to_node.json"),'w'))
    nwx.write_gexf(raw_graph, os.path.join(outdir,"raw.network.gexf"))
    nwx.write_gexf(refined_graph, os.path.join(outdir,"raw.network.gexf"))
    make_community_table(refined_communities, refined_graph, length_dict, outdir)
    
    dump_communities(refined_graph, refined_communities , records, outdir)
    logging.info(greenMSG('Done ... '))