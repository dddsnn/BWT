'''
Created on 10.05.2014

@author: dddsnn
'''

from bwt import *
import bwt.coder as cd
import bwt.analyzer as an
import networkx as nx
import pickle

def make_transitions(in_path, out_path=None):
    '''Create the transition analysis for a file.'''
    with open(in_path, 'rb') as in_file:
        text = in_file.read()
    trs = an.analyze_transitions(text)
    if out_path:
        with open(out_path, 'xb') as out_file:
            pickle.dump(trs, out_file)
    return trs


def make_graph(transitions, metric, out_path=None):
    g = nx.DiGraph()
    edges = []
    for a, b in transitions.keys():
        cost = getattr(transitions[(a, b)], metric)
        edges.append((a, b, {'cost':cost}))
    g.add_edges_from(edges)
    if out_path:
        with open(out_path, 'xb') as out_file:
            pickle.dump(g, out_file)
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
    par_text += 'RUNS = 100\n'
    par_text += 'TOUR_FILE = ' + file_name + '.tour'

    # write the files
    with open(out_dir_path + file_name + '.atsp', 'xt') as tsp_file:
        tsp_file.write(tsp_text)
    with open(out_dir_path + file_name + '.par', 'xt') as par_file:
        par_file.write(par_text)
    with open(out_dir_path + file_name + '.nodenames', 'xb') as names_file:
        pickle.dump(numbers_to_names, names_file)

def read_tsplib_files(in_path_tour, in_path_names):
    with open(in_path_tour, 'rt') as tour_file:
        line = tour_file.readline().rstrip()
        # skip until the TOUR_SECTION
        while line != 'TOUR_SECTION':
            line = tour_file.readline().rstrip()
        node_number_list = []
        line = tour_file.readline().rstrip()
        while line != '-1':
            node_number_list.append(int(line))
            line = tour_file.readline().rstrip()

    # now translate the numbers to names
    tour = []
    with open(in_path_names, 'rb') as names_file:
        names_dict = pickle.load(names_file)
        for number in node_number_list:
            tour.append(names_dict[number])
    return tour

def simulate_compression(in_path, title, order=None):
    '''Simulate compression of a file and print achieved compression ratio.'''
    with open(in_path, 'rb') as in_file:
        bytes_ = in_file.read()
    bw_code = cd.bw_encode(bytes_, order)
    mtf_code = cd.mtf_enc(bw_code.encoded)
    huff_code = cd.huffman_enc(mtf_code)
    res_text = '===================================================\n'
    res_text += title + '\n'
    res_text += 'file: {0}\n'.format(in_path)
    res_text += 'in size: {0}\n'.format(len(bytes_))
    res_text += 'out_size: {0}\n'.format(len(huff_code))
    res_text += 'ratio: {0}\n'.format(len(huff_code) / len(bytes_))
    res_text += '==================================================\n'
    print(res_text)

if __name__ == '__main__':
    wd = '/home/dddsnn/tmp/book1/'
    metrics = ['mean', 'median', 'num_chars', 'chapin1']

#     make_transitions('/home/dddsnn/Downloads/calgary/book1', wd + 'transitions')

#     with open(wd + 'transitions', 'rb') as trs_file:
#         trs = pickle.load(trs_file)
#     for metric in metrics:
#         g = make_graph(trs, metric)
#         write_tsplib_files(g, wd, metric)

    simulate_compression('/home/dddsnn/Downloads/calgary/book1', 'standard')
    simulate_compression('/home/dddsnn/Downloads/calgary/book1', 'aeiou...',
                         b'aeioubcdgfhrlsmnpqjktwvxyzAEIOUBCDGFHRLSMNPQJKTWVXYZ')
    for metric in metrics:
        tour = read_tsplib_files(wd + metric + '.tour',
                                 wd + metric + '.nodenames')
        simulate_compression('/home/dddsnn/Downloads/calgary/book1',
                             metric, tour)
