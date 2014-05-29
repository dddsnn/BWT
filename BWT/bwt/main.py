'''
Created on 10.05.2014

@author: dddsnn
'''

from collections import namedtuple
import bwt.coder as cd
import bwt.analyzer as an
from openopt import TSP, oosolver
import networkx as nx
import numpy as np
import pickle
from bwt.coder import bw_encode

def make_transitions(in_path, out_path=None):
    '''Create the transition analysis for a file.'''
    with open(in_path) as in_file:
        text = in_file.read()
    trs = an.analyze_transitions(text)
    if out_path:
        with open(out_path, 'wb') as out_file:
            pickle.dump(trs, out_file)
    return trs

def make_graph(transitions, out_path=None):
    # TODO
    trs = pickle.load(open('/home/dddsnn/tmp/book1/transitions', 'rb'))
    g = nx.DiGraph()
    edges = []
    # encode the real names as simpler strings, because numpy can't handle
    # '\x00' and openopt can't handle integers as node names
    # write 'nx' for the normal character x and 'sx' for special characters
    # '\x00' -> 's0'
    for k in trs.keys():
        if k[0] == '\x00':
            a = 's0'
        else:
            a = 'n' + k[0]
        if k[1] == '\x00':
            b = 's0'
        else:
            b = 'n' + k[1]
        edges.append((a, b, {'cost':trs[k].mean.diff}))
    g.add_edges_from(edges)
    print('graph created')
    pickle.dump(g, open('/home/dddsnn/tmp/book1/graph', 'wb'))
    return g

def simulate_compression(in_path, order=None):
    '''Simulate compression of a file and print achieved compression ratio.'''
    with open(in_path, 'rb') as in_file:
        bytes_ = in_file.read()
    bw_code = cd.bw_encode(bytes_, order)
    mtf_code = cd.mtf_enc(bw_code.encoded)
    huff_code = cd.huffman_enc(mtf_code)
    res_text = 'file: {0}\nin size: {1}\nout_size: {2}\nratio: {3}'
    res_text = res_text.format(in_path, len(bytes_), len(huff_code),
                               len(huff_code) / len(bytes_))
    print(res_text)

if __name__ == '__main__':
    simulate_compression('/home/dddsnn/Downloads/calgary/geo')

#     g = pickle.load(open('/home/dddsnn/tmp/book1/graph', 'rb'))
#     solver = oosolver('glpk', name='asf', fTol=1.0)
#     problem = TSP(g, objective='cost')
#     result = problem.solve('glpk', maxTime=1)
#     pickle.dump(result.nodes, open('/home/dddsnn/tmp/book1/result-glpk-nodes', 'wb'))
#     pickle.dump(result.edges, open('/home/dddsnn/tmp/book1/result-glpk-edges', 'wb'))
#     pickle.dump(result.Edges, open('/home/dddsnn/tmp/book1/result-glpk-Edges', 'wb'))

#     nodes = pickle.load(open('/home/dddsnn/tmp/book1/result-interalg-nodes', 'rb'))
#     edges = pickle.load(open('/home/dddsnn/tmp/book1/result-interalg-edges', 'rb'))
#     Edges = pickle.load(open('/home/dddsnn/tmp/book1/result-interalg-Edges', 'rb'))
