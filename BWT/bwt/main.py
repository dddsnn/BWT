from bwt import *
import bwt.coder as cd
import bwt.analyzer as an
import bwt.huffman as hf
import os.path
import networkx as nx
import pickle
import numpy as np
import time

def make_aux_data(in_path, out_path=None):
    """Create data commonly used in transition analysis."""
    with open(in_path, 'rb') as in_file:
        bs = in_file.read()
    raw = bs
    bw_code = cd.bw_encode(bs)
    mtf_code = cd.mtf_enc(bw_code.encoded)
    firsts = list(set([bytes([x[0]]) for x in bw_code.firsts]))
    firsts.extend(an.select_sequences(bs, 2))
    specializations = cd.specializations(firsts)
    bw_subcodes = {x:an.bw_block(bw_code, x, specializations) for x in firsts}
    # first add subcodes for individual bytes
    partial_mtf_subcodes = {x:cd.mtf_partial_enc(bw_subcodes[x])
                            for x in bw_subcodes}
    # now make and add the transitions
    partial_mtf_transitions = {(a, b):cd.mtf_partial_enc(bw_subcodes[a]
                                                         + bw_subcodes[b])
                               for a in bw_subcodes for b in bw_subcodes
                               if a != b}
    partial_mtf_subcodes.update(partial_mtf_transitions)
    partial_mtf_analyses = {x:an.analyze_partial_mtf(partial_mtf_subcodes[x])
                            for x in partial_mtf_subcodes}
    bw_subhistograms = {x:an.make_histogram(bw_subcodes[x])
                        for x in bw_subcodes}
    huffcode_len_complete = an.huffman_codeword_lengths(mtf_code, 'complete')
    huffcode_len_sparse = an.huffman_codeword_lengths(mtf_code, 'sparse')
    mtf_mean_steps = an.mtf_mean_steps(bw_code.encoded, mtf_code)
    freq_lists = {}
    for f in bw_subhistograms:
        # turn the histograms into lists and sort in decreasing order of frequency
        freq_list = sorted(bw_subhistograms[f].items(),
                           key=lambda x:x[1], reverse=True)
        # now just take the corresponding first symbols
        freq_list = [x[0] for x in freq_list]
        # add to the dict
        freq_lists[f] = freq_list

    result = AuxData(raw, firsts, bw_code, mtf_code, bw_subcodes,
                     partial_mtf_subcodes, partial_mtf_analyses,
                     bw_subhistograms, huffcode_len_complete,
                     huffcode_len_sparse, mtf_mean_steps, freq_lists)
    if out_path:
        with open(out_path, 'xb') as out_file:
            pickle.dump(result, out_file)
    return result

def make_transitions(in_path, aux_data, metric, out_path=None, **metric_opts):
    """Create the transition analysis for a file."""
    with open(in_path, 'rb') as in_file:
        bs = in_file.read()
    trs = an.analyze_transitions(bs, aux_data, metric, **metric_opts)
    if out_path:
        with open(out_path, 'xb') as out_file:
            pickle.dump(trs, out_file)
    return trs

def make_graph(transitions, out_path=None):
    g = nx.DiGraph()
    edges = []
    for a, b in transitions.keys():
        edges.append((a, b, {'cost':transitions[(a, b)]}))
    g.add_edges_from(edges)
    if out_path:
        with open(out_path, 'xb') as out_file:
            pickle.dump(g, out_file)
    return g

def write_tsplib_files(graph, out_dir_path, file_name):
    """Create and write the .tsp and .par files for the LKH program as well as
    the file mapping the LKH node ids back to the node names."""
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
    costs_nonzero = [abs(c) for c in costs.values() if c != 0]
    if costs_nonzero:
        min_cost = min(costs_nonzero)
    else:
        min_cost = 0
    max_cost = max([abs(c) for c in costs.values()])
    if min_cost:
        ratio = max_cost / min_cost
    else:
        ratio = float('inf')
    print('min value: {0}'.format(min_cost))
    print('max value: {0}'.format(max_cost))
    print('ratio: {0}'.format(ratio))
    # choose the scaling factor so that the max cost is close to MAX_INT
    if max_cost:
        factor = int((1 / max_cost) * MAX_INT)
    else:
        # in case everything is 0
        factor = 0
    print('scaling with factor {0}'.format(factor))
    print()

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

def very_greedy_tsp(transitions):
    # length of the tour is the number of different symbols
    length = len(set([x[0]for x in transitions.keys()]))
    # make a list of transitions from lowest cost to highest
    sorted_transitions = sorted(transitions.keys(),
                                key=lambda k: transitions[k])
    # list of subtours (lists of nodes)
    parts = []
    for t in sorted_transitions:
        # check that neither source nor destination of the transition are
        # already part of a partial tour
        inner_nodes = [n for p in parts for n in p[1:-1]]
        sources = [p[0] for p in parts]
        dests = [p[-1] for p in parts]
        # dicts mapping the source of a part to its destination
        source_dest = {p[0]:p[-1] for p in parts}
        if t[0] in inner_nodes or t[1] in inner_nodes or t[0] in sources \
                or t[1] in dests:
            # already part of a tour, continue
            continue
        # transitions are allowed to build onto existing parts, but only in the
        # right order
        # two parts being linked together
        if t[0] in dests and t[1] in sources:
            # check that it's not making a cycle
            if t[0] == source_dest[t[1]]:
                # cycle, continue
                continue
            left = next(p for p in parts if p[-1] == t[0])
            right = next(p for p in parts if p[0] == t[1])
            left.extend(right)
            parts.remove(right)
            continue

        # just one append/prepend to an existing part
        if t[0] in dests:
            # find the part and append
            part = next(p for p in parts if p[-1] == t[0])
            part.append(t[1])
            continue
        if t[1] in sources:
            # find the part and prepend
            part = next(p for p in parts if p[0] == t[1])
            part.insert(0, t[0])
            continue

        # transition not part of anything, make a new part
        parts.append([t[0], t[1]])

        # if all the symbols are in parts[0], we're done
        if len(parts[0]) == length:
            break
    return parts[0]

def simulate_compression(in_path, title, order=None):
    """Simulate compression of a file and print achieved compression ratio."""
    with open(in_path, 'rb') as in_file:
        bs = in_file.read()
    bw_code = cd.bw_encode(bs, order)
    mtf_code = cd.mtf_enc(bw_code.encoded)
    huff_code = hf.encode_to_bits_static(mtf_code)
    res_text = '======================================\n'
    res_text += title + '\n'
    res_text += 'file: {0}\n'.format(in_path)
    res_text += 'in size: {0}\n'.format(len(bs) * 8)
    res_text += 'out_size: {0}\n'.format(len(huff_code))
    res_text += 'ratio: {0}\n'.format(len(huff_code) / (len(bs) * 8))
    res_text += '======================================\n'
    print(res_text)

if __name__ == '__main__':
    def metric_file_name(metric):
        elems = [metric[0]]
        for opt in metric[1]:
            elems.append(opt)
            if opt != True:
                elems.append(metric[1][opt])
        return '_'.join(elems)

    start_time = time.time()
    in_dir = '/home/dddsnn/Dokumente/Studium/BA/calgary/'
    in_file_names = ['bib', 'book1', 'book2', 'geo', 'news', 'obj1', 'obj2',
                     'paper1', 'paper2', 'paper3', 'paper4', 'paper5',
                     'paper6', 'pic', 'progc', 'progl', 'progp', 'trans']
    in_file_names = ['book1']
    base_work_dir = '/home/dddsnn/tmp/'
    metrics = [('chapin_hst_diff', {}), ('chapin_inv', {}),
               ('chapin_inv', {'log':True}),
               ('badness', {'weighted':False, 'new_penalty':False,
                            'entropy_code_len':False}),
               ('badness', {'weighted':False, 'new_penalty':False,
                            'entropy_code_len':'complete'}),
               ('badness', {'weighted':False, 'new_penalty':False,
                            'entropy_code_len':'sparse'}),
               ('badness', {'weighted':False, 'new_penalty':'generic',
                            'entropy_code_len':False}),
               ('badness', {'weighted':False, 'new_penalty':'generic',
                            'entropy_code_len':'complete'}),
               ('badness', {'weighted':False, 'new_penalty':'generic',
                            'entropy_code_len':'sparse'}),
               ('badness', {'weighted':False, 'new_penalty':'specific',
                            'entropy_code_len':False}),
               ('badness', {'weighted':False, 'new_penalty':'specific',
                            'entropy_code_len':'complete'}),
               ('badness', {'weighted':False, 'new_penalty':'specific',
                            'entropy_code_len':'sparse'}),
               ('badness', {'weighted':True, 'new_penalty':False,
                            'entropy_code_len':False}),
               ('badness', {'weighted':True, 'new_penalty':False,
                            'entropy_code_len':'complete'}),
               ('badness', {'weighted':True, 'new_penalty':False,
                            'entropy_code_len':'sparse'}),
               ('badness', {'weighted':True, 'new_penalty':'generic',
                            'entropy_code_len':False}),
               ('badness', {'weighted':True, 'new_penalty':'generic',
                            'entropy_code_len':'complete'}),
               ('badness', {'weighted':True, 'new_penalty':'generic',
                            'entropy_code_len':'sparse'}),
               ('badness', {'weighted':True, 'new_penalty':'specific',
                            'entropy_code_len':False}),
               ('badness', {'weighted':True, 'new_penalty':'specific',
                            'entropy_code_len':'complete'}),
               ('badness', {'weighted':True, 'new_penalty':'specific',
                            'entropy_code_len':'sparse'}), ]

    # make directories
    for in_file_name in in_file_names:
        if not os.path.exists(base_work_dir + in_file_name):
            os.mkdir(base_work_dir + in_file_name)

    # make aux data
    for in_file_name in in_file_names:
        in_path = in_dir + in_file_name
        wd = base_work_dir + in_file_name + '/'
        make_aux_data(in_path, wd + 'aux')

    # make transitions
#     for in_file_name in in_file_names:
#         in_path = in_dir + in_file_name
#         wd = base_work_dir + in_file_name + '/'
#         with open(wd + 'aux', 'rb') as aux_file:
#             aux_data = pickle.load(aux_file)
#         for metric in metrics:
#             file_name = metric_file_name(metric)
#             make_transitions(in_path, aux_data, metric[0],
#                             wd + file_name + '.transitions', **metric[1])
#
# #     # write tsplib files
#     for in_file_name in in_file_names:
#         in_path = in_dir + in_file_name
#         wd = base_work_dir + in_file_name + '/'
#         for metric in metrics:
#             file_name = metric_file_name(metric)
#             with open(wd + file_name + '.transitions', 'rb') as trs_file:
#                 trs = pickle.load(trs_file)
#             g = make_graph(trs)
#             write_tsplib_files(g, wd, file_name)
#
#     # simulate compression
#     for in_file_name in in_file_names:
#         in_path = in_dir + in_file_name
#         wd = base_work_dir + in_file_name + '/'
#         handpicked_str = b'aeioubcdgfhrlsmnpqjktwvxyzAEIOUBCDGFHRLSMNPQJKTWVXYZ'
#         handpicked_order = [[bytes([c]) for c in handpicked_str]]
#         simulate_compression(in_path, 'aeiou...', handpicked_order)
#         simulate_compression(in_path, 'standard')
#
#         for metric in metrics:
#             file_name = metric_file_name(metric)
#             tsplib_tour = [read_tsplib_files(wd + file_name + '.tour',
#                                      wd + file_name + '.nodenames')]
#             with open(wd + file_name + '.transitions', 'rb') as trs_file:
#                 trs = pickle.load(trs_file)
#             very_greedy_tour = [very_greedy_tsp(trs)]
#             simulate_compression(in_path, file_name + ' tsplib', tsplib_tour)
# #             simulate_compression(in_path, file_name + ' very greedy',
# #                                 very_greedy_tour)

    print('time: {0}s'.format(time.time() - start_time))
