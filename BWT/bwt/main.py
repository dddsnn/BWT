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

def make_transitions(in_path, out_path=None):
    with open(in_path) as in_file:
        text = in_file.read()
    trs = an.analyze_transitions(text)
    if out_path:
        with open(out_path, 'wb') as out_file:
            pickle.dump(trs, out_file)
    return trs

def make_graph(transitions, out_path=None):
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

if __name__ == '__main__':
    text = 'abracadabra'
    bw1 = cd.bw_encode(text)
    print(bw1.firsts)
    print(bw1.encoded)
    print()
    bw2 = cd.bw_encode(text, {'r':0, 'b':12, 'c':2, 'a':9, 'd':3})
    print(bw2.firsts)
    print(bw2.encoded)
#     f = open('/home/dddsnn/Downloads/calgary/book1')
#     text = f.read()
#     trs = analyze_transitions(text)
#     print('transitions done')
#     pickle.dump(trs, open('/home/dddsnn/tmp/book1/transitions', 'wb'))

#     trs = pickle.load(open('/home/dddsnn/tmp/book1/transitions', 'rb'))
#     g = nx.DiGraph()
#     edges = []
#     # encode the real names as simpler strings, because numpy can't handle
#     # '\x00' and openopt can't handle integers as node names
#     # write 'nx' for the normal character x and 'sx' for special characters
#     # '\x00' -> 's0'
#     for k in trs.keys():
#         if k[0] == '\x00':
#             a = 's0'
#         else:
#             a = 'n' + k[0]
#         if k[1] == '\x00':
#             b = 's0'
#         else:
#             b = 'n' + k[1]
#         edges.append((a, b, {'cost':trs[k].mean.diff}))
#     g.add_edges_from(edges)
#     print('graph created')
#     pickle.dump(g, open('/home/dddsnn/tmp/book1/graph', 'wb'))

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
#     print()

#     bw = cd.bw_encode(text)
# #     print(bw2)
#     mtf = cd.mtf_enc(bw.encoded)
#     hf = cd.huffman_enc(mtf)
#     dec = cd.huffman_dec(hf)
#     print('{0} / {1} : {2}'
#         .format(len(text), len(hf), len(hf) / len(text)))
