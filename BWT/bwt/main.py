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
    '''Create the transition analysis for a file.'''
    with open(in_path) as in_file:
        text = in_file.read()
    trs = an.analyze_transitions(text)
    if out_path:
        with open(out_path, 'wb') as out_file:
            pickle.dump(trs, out_file)
    return trs


def make_graph(transitions, out_path=None):
    # TODO switch to bytes
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

def write_tsplib_files(graph, out_dir_path, file_name):
    '''Create and write the .tsp and .par files for the LKH program as well as
    the file mapping the LKH node ids back to the node names.'''
    num_nodes = len(graph.nodes())
    # name mapping part
    # make dicts nodename->number and number->nodename for the tsp file
    names_to_numbers = numbers_to_names = {}
    for i, n in enumerate(graph.nodes(), 1):
        names_to_numbers[n] = i
        numbers_to_names[i] = n

    # tsp file part
    # start with the header
    tsp_text = 'NAME: bwt\n'
    tsp_text += 'TYPE: ATSP\n'
    tsp_text += 'DIMENSION: {0}\n'.format(num_nodes)
    tsp_text += 'EDGE_WEIGHT_TYPE: EXPLICIT\n'
    tsp_text += 'EDGE_WEIGHT_FORMAT: FULL_MATRIX\n'
    tsp_text += 'EDGE_WEIGHT_SECTION\n'

    MAX_INT = 2 * 10 ** 7  # the maximum value to write to the tsp file
    INFINITY = 2147483648  # value to signify infinity in the tsp file
    # need to scale up all floats to integers
    # dictionary edge->cost
    costs = nx.get_edge_attributes(graph, 'cost')
    # min and max edge costs, and the ratio
    min_cost = min([abs(c) for c in costs.values() if c != 0])
    max_cost = max([abs(c) for c in costs.values()])
    ratio = max_cost / min_cost
    print('ratio between smallest and largest value: {0}'.format(ratio))
    # choose the scaling factor so that the max cost is close to MAX_INT
    factor = int((1 / max_cost) * MAX_INT)
    print('scaling with factor {0}'.format(factor))

    max_len = len(str(INFINITY))  # length of the longest number
    # append the scaled cost value matrix to the tsp text
    for irow in range(1, num_nodes + 1):
        for icol in range(1, num_nodes + 1):
            if irow == icol:
                # node to itself -> infinity
                write_num = INFINITY
            else:
                edge = (numbers_to_names[irow], numbers_to_names[icol])
                write_num = int(costs[edge] * factor)
            # leave 2 characters spacing
            space = ' ' * (max_len - len(str(write_num)) + 2)
            # append the string
            tsp_text += space + str(write_num)
        # newline at the end of a row
        tsp_text += '\n'
    # EOF at the end of the file
    tsp_text += 'EOF'

    # par file part
    par_text = 'PROBLEM_FILE = {0}\n'.format(file_name + '.atsp')
    par_text += 'TOUR_FILE = ' + file_name + '.tour'

    # write the files
    with open(out_dir_path + file_name + '.atsp', 'wt') as tsp_file:
        tsp_file.write(tsp_text)
    with open(out_dir_path + file_name + '.par', 'wt') as par_file:
        par_file.write(par_text)
    with open(out_dir_path + file_name + '.nodenames', 'wb') as names_file:
        pickle.dump(numbers_to_names, names_file)

def solve_tsp():
    # TODO
    g = pickle.load(open('/home/dddsnn/tmp/book1/graph', 'rb'))
    solver = oosolver('glpk', name='asf', fTol=1.0)
    problem = TSP(g, objective='cost')
    result = problem.solve('glpk', maxTime=1)
    pickle.dump(result.nodes, open('/home/dddsnn/tmp/book1/result-glpk-nodes', 'wb'))
    pickle.dump(result.edges, open('/home/dddsnn/tmp/book1/result-glpk-edges', 'wb'))
    pickle.dump(result.Edges, open('/home/dddsnn/tmp/book1/result-glpk-Edges', 'wb'))

    nodes = pickle.load(open('/home/dddsnn/tmp/book1/result-interalg-nodes', 'rb'))
    edges = pickle.load(open('/home/dddsnn/tmp/book1/result-interalg-edges', 'rb'))
    Edges = pickle.load(open('/home/dddsnn/tmp/book1/result-interalg-Edges', 'rb'))

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
    g = pickle.load(open('/home/dddsnn/tmp/book1/graph', 'rb'))
    e = nx.get_edge_attributes(g, 'cost')
    print(g.nodes())
    write_tsplib_files(g, '/home/dddsnn/tmp/', 'tsp_text')
#     order = b'aeioubcdgfhrlsmnpqjktwvxyzAEIOUBCDGFHRLSMNPQJKTWVXYZ'
#     simulate_compression('/home/dddsnn/Downloads/calgary/pic')
