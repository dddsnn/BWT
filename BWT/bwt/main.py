from bwt import *
import bwt.coder as cd
import bwt.analyzer as an
import bwt.huffman as hf
import os.path
import networkx as nx
import pickle
import numpy as np
import time
import math

def make_aux_data(in_path, out_path=None):
    """Create data commonly used in transition analysis."""
    with open(in_path, 'rb') as in_file:
        bs = in_file.read()
    raw = bs
    num_symbols = len(set(bs))
    bw_code = cd.bw_encode(bs)
    mtf_code = cd.mtf_enc(bw_code.encoded)
    firsts = list(set([bytes([x[0]]) for x in bw_code.firsts]))
    firsts.extend(an.select_sequences(bs, 2))
    specializations = cd.specializations(firsts)
    bw_subcodes = {x:an.bw_block(bw_code.encoded, bw_code.firsts, x,
                                 specializations) for x in firsts}
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

    result = AuxData(raw, num_symbols, firsts, bw_code, mtf_code, bw_subcodes,
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
        with open(out_path + '.transitions', 'xb') as out_file:
            pickle.dump(trs, out_file)
    # write the new penalty log if it exists
    if out_path and 'new_penalty_log' in metric_opts:
        with open(out_path + '.new_penalty_log', 'xb') as out_file:
            pickle.dump(metric_opts['new_penalty_log'], out_file)
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

def simulate_compression(in_path, title, orders=None):
    """Simulate compression of a file and print achieved compression ratio."""
    with open(in_path, 'rb') as in_file:
        bs = in_file.read()
    bw_code = cd.bw_encode(bs, orders)
    mtf_code = cd.mtf_enc(bw_code.encoded)
    huff_code = hf.encode_to_bits_static(mtf_code)
    res_text = '======================================\n'
    res_text += title + '\n'
    res_text += 'file: {0}\n'.format(in_path)
    res_text += 'in size: {0}\n'.format(len(bs) * 8)
    res_text += 'out size: {0}\n'.format(len(huff_code))
    res_text += 'ratio: {0}\n'.format(len(huff_code) / (len(bs) * 8))
    res_text += '======================================\n'
    print(res_text)

def compare_new_penalty_predictions(aux_data, orders, new_penalty_log):
    """Compare the predictions for new codes with actual values.

    Args:
        aux_data: The bwt.AuxData for the input that the transitions were
            computed for.
        orders: A list of orders that is the result of the TSP on the computed
            transitions. Same format as in coder.bw_encode().
        new_penalty_log: The log dictionary generated by the badness metric when
            provided with the appropriate option.

    Returns:
        A list of tuples. For every transition between context blocks according
            to the given orders, for each new symbol on the right side of the
            transition a tuple is created containing the context of the block,
            the MTF code of the first appearance of the symbol in the actual
            code, and the logged prediction for that symbol.
            Symbols in the first block can't be compared, because the first
            block is not the right side of any transition.
    """
    bw_code = cd.bw_encode(aux_data.raw, orders)
    mtf = cd.mtf_enc(bw_code.encoded)
    bw_mtf = list(zip(bw_code.encoded, mtf))
    result = []
    # can only handle the first order for the moment
    order = orders[0]
    transitions = zip(order[:-1], order[1:])
    spec = cd.specializations(aux_data.firsts)
    for trs in transitions:
        if trs not in new_penalty_log:
            # there were not predictions necessary for this transition, skip
            continue
        # get the block of mtf-bw tuples corresponding to the right side of the
        # transition
        block = an.bw_block(bw_mtf, bw_code.firsts, trs[1], spec)
        predictions = new_penalty_log[trs]
        # for every prediction, record the current sequence, the actual code
        # and the prediction
        for s in sorted(predictions.keys()):
            # actual code is the mtf code of the first appearance of the symbol
            actual = next(b[1] for b in block if b[0] == s)
            result.append((trs[1], actual, predictions[s]))
    return result

def compare_entropy_len_predictions(aux_data, orders, predictions):
    """Compare the predictions for entropy code lengths with actual values.

    Args:
        aux_data: The bwt.AuxData for the input.
        orders: The computed sort orders yielding the actual entropy code
            lengths.
        predictions: A dictionary mapping MTF codes to the predicted entropy
            code length, as returned by bwt.analyzer.huffman_codeword_lengths().

    Returns:
        A list of tuples. Each tuple contains an MTF code, its actual entropy
        codeword length when the file is BW encoded with the given orders, and
        the predicted length for that MTF code.
    """
    mtf_code = cd.mtf_enc(cd.bw_encode(aux_data.raw, orders).encoded)
    freqs = hf.symbol_frequencies(mtf_code)
    actual = hf.codeword_lengths(freqs)
    result = []
    for s in actual:
        result.append((s, actual[s], predictions[s]))
    return result

if __name__ == '__main__':
    def metric_unique_name(metric):
        """Create a unique name for a metric with options."""
        elems = [metric[0]]
        opts = metric[1]
        for opt in sorted([o for o in opts if opts[o]]):
            if opt == 'new_penalty_log':
                # don't record the log as an option
                continue
            elems.append(opt)
            if opts[opt] != True:
                elems.append(opts[opt])
        return '_'.join(elems)

    start_time = time.time()
    in_dir = '/home/dddsnn/Dokumente/Studium/BA/calgary/'
    in_file_name = 'book1'
    base_work_dir = '/home/dddsnn/tmp/'
    metrics = [('chapin_hst_diff', {}), ('chapin_inv', {}),
               ('chapin_inv', {'log':True})]
    for w in [True, False]:
        for entr_len in [False, 'complete', 'sparse']:
            for new_pen in [False, 'generic', 'specific']:
                opts = {'weighted':w, 'entropy_code_len':entr_len,
                        'new_penalty':new_pen}
                if new_pen != False:
                    opts['new_penalty_log'] = {}
                metrics.append(('badness', opts))
    metrics = [('badness', {'weighted':False, 'new_penalty':False,
                            'entropy_code_len':False, 'new_penalty_log':{}}),
               ('badness', {'weighted':False, 'new_penalty':False,
                            'entropy_code_len':'complete', 'new_penalty_log':{}}),
               ('badness', {'weighted':False, 'new_penalty':False,
                            'entropy_code_len':'sparse', 'new_penalty_log':{}}),
               ('badness', {'weighted':False, 'new_penalty':'generic',
                            'entropy_code_len':False, 'new_penalty_log':{}}),
               ('badness', {'weighted':False, 'new_penalty':'generic',
                            'entropy_code_len':'complete',
                            'new_penalty_log':{}}),
               ('badness', {'weighted':False, 'new_penalty':'generic',
                            'entropy_code_len':'sparse', 'new_penalty_log':{}}),
               ('badness', {'weighted':False, 'new_penalty':'specific',
                            'entropy_code_len':False, 'new_penalty_log':{}}),
               ('badness', {'weighted':False, 'new_penalty':'specific',
                            'entropy_code_len':'complete',
                            'new_penalty_log':{}}),
               ('badness', {'weighted':False, 'new_penalty':'specific',
                            'entropy_code_len':'sparse', 'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':False,
                            'entropy_code_len':False, 'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':False,
                            'entropy_code_len':'complete', 'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':False,
                            'entropy_code_len':'sparse', 'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':'generic',
                            'entropy_code_len':False, 'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':'generic',
                            'entropy_code_len':'complete',
                            'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':'generic',
                            'entropy_code_len':'sparse', 'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':'specific',
                            'entropy_code_len':False, 'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':'specific',
                            'entropy_code_len':'complete',
                            'new_penalty_log':{}}),
               ('badness', {'weighted':True, 'new_penalty':'specific',
                            'entropy_code_len':'sparse', 'new_penalty_log':{}})]
#     metrics = [('badness', {'weighted':False, 'new_penalty':False,
#                             'entropy_code_len':False, 'new_penalty_log':{}}),
#                ('badness', {'weighted':True, 'new_penalty':'generic',
#                             'entropy_code_len':'sparse', 'new_penalty_log':{}})]

    # make directories
    if not os.path.exists(base_work_dir + in_file_name):
        os.mkdir(base_work_dir + in_file_name)

    # make aux data
    in_path = in_dir + in_file_name
    wd = base_work_dir + in_file_name + '/'
    make_aux_data(in_path, wd + 'aux')

    # make transitions
    in_path = in_dir + in_file_name
    wd = base_work_dir + in_file_name + '/'
    with open(wd + 'aux', 'rb') as aux_file:
        aux_data = pickle.load(aux_file)
    for metric in metrics:
        file_name = metric_unique_name(metric)
        make_transitions(in_path, aux_data, metric[0],
                        wd + file_name, **metric[1])

    # write tsplib files
    in_path = in_dir + in_file_name
    wd = base_work_dir + in_file_name + '/'
    for metric in metrics:
        file_name = metric_unique_name(metric)
        with open(wd + file_name + '.transitions', 'rb') as trs_file:
            trs = pickle.load(trs_file)
        g = make_graph(trs)
        write_tsplib_files(g, wd, file_name)

    # simulate compression
    in_path = in_dir + in_file_name
    wd = base_work_dir + in_file_name + '/'
    handpicked_str = b'aeioubcdgfhrlsmnpqjktwvxyzAEIOUBCDGFHRLSMNPQJKTWVXYZ'
    handpicked_order = [[bytes([c]) for c in handpicked_str]]
    simulate_compression(in_path, 'aeiou...', handpicked_order)
    simulate_compression(in_path, 'standard')

    for metric in metrics:
        file_name = metric_unique_name(metric)
        tsplib_tour = [read_tsplib_files(wd + file_name + '.tour',
                                 wd + file_name + '.nodenames')]
#         with open(wd + file_name + '.transitions', 'rb') as trs_file:
#             trs = pickle.load(trs_file)
#         very_greedy_tour = [very_greedy_tsp(trs)]
        simulate_compression(in_path, file_name + ' tsplib', tsplib_tour)
#         simulate_compression(in_path, file_name + ' very greedy',
#                             very_greedy_tour)

    # compare new penalty predictions with actual values
    print('evaluating new symbol penalties for file {0}'.format(in_file_name))
    print()
    in_path = in_dir + in_file_name
    wd = base_work_dir + in_file_name + '/'
    with open(wd + 'aux', 'rb') as aux_file:
        aux_data = pickle.load(aux_file)
    for metric in metrics:
        file_name = metric_unique_name(metric)
        with open(wd + file_name + '.new_penalty_log', 'rb') as in_file:
            new_penalty_log = pickle.load(in_file)
        print('metric: {0}'.format(file_name))
        tsplib_tour = read_tsplib_files(wd + file_name + '.tour',
                                         wd + file_name + '.nodenames')
        natural_order = [bytes([x]) for x in range(256)]
        orders = [tsplib_tour, natural_order]
        comp = compare_new_penalty_predictions(aux_data, orders,
                                               new_penalty_log)
        mtf_code = cd.mtf_enc(cd.bw_encode(aux_data.raw, orders).encoded)
        freqs = hf.symbol_frequencies(mtf_code)
        hf_len = hf.codeword_lengths(freqs)

        # TODO mean deviation etc for guesses that were larger than the actual
        # value/smaller than the actual value
        diffs = [c[1] - c[2] for c in comp]
        print('mean differences: {0}'.format(np.mean(diffs)))
        diffs_no_new = [c[1] - c[2] for c in comp
                        if c[1] < aux_data.num_symbols]
        print('mean differences w/o new symbols: {0}'
              .format(np.mean(diffs_no_new)))
        devs = map(abs, diffs)
        print('mean deviation: {0}'.format(np.mean(devs)))
        devs_no_new = map(abs, diffs_no_new)
        print('mean deviation w/o new symbols: {0}'
              .format(np.mean(devs_no_new)))
        variance = np.mean(map(lambda x:x ** 2, diffs))
        variance_no_new = np.mean(map(lambda x:x ** 2, diffs_no_new))
        std_deviation = math.sqrt(variance)
        std_deviation_no_new = math.sqrt(variance_no_new)
        print('standard deviation: {0}'.format(std_deviation))
        print('standard deviation w/o new symbols: {0}'
              .format(std_deviation_no_new))

        # take huffman codeword length of the closest mtf code to the
        # predicted one
        hf_diffs = [hf_len[c[1]] - hf_len[min(hf_len.keys(),
                                              key=lambda x:abs(x - c[2]))]
                    for c in comp]
        print('mean differences of entropy codes: {0}'.
              format(np.mean(hf_diffs)))
        hf_diffs_no_new = [hf_len[c[1]] -
                            hf_len[min(hf_len.keys(),
                                       key=lambda x:abs(x - c[2]))]
                                     for c in comp
                                     if c[1] < aux_data.num_symbols]
        print('mean differences of entropy codes w/o new symbols: {0}'
              .format(np.mean(hf_diffs_no_new)))
        hf_devs = map(abs, hf_diffs)
        print('mean deviation of entropy codes: {0}'.format(np.mean(hf_devs)))
        hf_devs_no_new = map(abs, hf_diffs_no_new)
        print('mean deviation of entropy codes w/o new symbols: {0}'
              .format(np.mean(hf_devs_no_new)))
        hf_variance = map(lambda x:x ** 2, hf_diffs)
        hf_variance_no_new = map(lambda x:x ** 2, hf_diffs_no_new)
        hf_std_deviation = math.sqrt(hf_variance)
        hf_std_deviation_no_new = math.sqrt(hf_variance_no_new)
        print('standard deviation of entropy codes: {0}'
              .format(hf_std_deviation))
        print('standard deviation of entropy codes w/o new symbols: {0}'
              .format(hf_std_deviation_no_new))
        print()

    # compare entropy length predictions with actual values
    print('evaluating entropy code code length predictions for file {0}'
          .format(in_file_name))
    print()
    in_path = in_dir + in_file_name
    wd = base_work_dir + in_file_name + '/'
    with open(wd + 'aux', 'rb') as aux_file:
        aux_data = pickle.load(aux_file)
    for metric in metrics:
        if 'entropy_code_len' in metric[1]:
            if metric[1]['entropy_code_len'] == 'complete':
                prediction = aux_data.huffman_codeword_lengths_complete
            elif metric[1]['entropy_code_len'] == 'sparse':
                prediction = aux_data.huffman_codeword_lengths_sparse
            else:
                continue
        file_name = metric_unique_name(metric)
        print('metric: {0}'.format(file_name))
        tsplib_tour = read_tsplib_files(wd + file_name + '.tour',
                                         wd + file_name + '.nodenames')
        comp = compare_entropy_len_predictions(aux_data, [tsplib_tour],
                                               prediction)
        diffs = [c[1] - c[2] for c in comp]
        print('mean differences: {0}'.format(np.mean(diffs)))
        devs = map(abs, diffs)
        print('mean deviation: {0}'.format(np.mean(devs)))
        variance = np.mean(map(lambda x:x ** 2, diffs))
        std_deviation = math.sqrt(variance)
        print('standard deviation: {0}'.format(std_deviation))
        print()

    print('time: {0:.0f}s'.format(time.time() - start_time))
