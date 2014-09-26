from bwt import *
import bwt.coder as cd
import bwt.analyzer as an
import bwt.huffman as hf
import pickle
import numpy as np
import time
import math
import os
from itertools import chain

def make_aux_data(work_dir, in_file_path, col_depth=1):
    """Create data commonly used in transition analysis."""
    with open(in_file_path, 'rb') as in_file:
        bs = in_file.read()
    raw = bs
    num_symbols = len(set(bs))
    bw_code = cd.bw_encode(bs)
    mtf_code = cd.mtf_encode(bw_code.encoded)
    firsts = list(set([bytes([x[0]]) for x in bw_code.firsts]))
    for i in range(1, col_depth):
        for prefix in [x for x in firsts if len(x) == i]:
            symbols = list(set([x[i] for x in bw_code.firsts
                                if x[:i] == prefix]))
            new = [prefix + bytes([sym]) for sym in symbols]
            firsts.extend(new)
    bw_subcodes = {x:an.context_block(bw_code.encoded, bw_code.firsts, x)
                   for x in firsts}
    # first add subcodes for individual bytes
    partial_mtf_subcodes = {x:cd.mtf_partial_encode(bw_subcodes[x])
                            for x in bw_subcodes}
    # now make and add the transitions
    partial_mtf_transitions = {(a, b):cd.mtf_partial_encode(bw_subcodes[a]
                                                         + bw_subcodes[b])
                               for a in bw_subcodes for b in bw_subcodes
                               if a != b and a[:-1] == b[:-1]}
    partial_mtf_subcodes.update(partial_mtf_transitions)
    partial_mtf_analyses = {x:an.analyze_partial_mtf(partial_mtf_subcodes[x])
                            for x in partial_mtf_subcodes}
    bw_subhistograms = {x:an.make_histogram(bw_subcodes[x])
                        for x in bw_subcodes}
    huffcode_len_complete = an.huffman_codeword_lengths(mtf_code, 'complete')
    huffcode_len_sparse = an.huffman_codeword_lengths(mtf_code, 'sparse')
    mtf_mean_steps = an.mtf_avg_steps(bw_code.encoded, mtf_code, np.mean)
    mtf_median_steps = an.mtf_avg_steps(bw_code.encoded, mtf_code, np.median)
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
                     huffcode_len_sparse, mtf_mean_steps, mtf_median_steps,
                     freq_lists)
    with open(work_dir + 'aux', 'xb') as out_file:
        pickle.dump(result, out_file)
    return result

def make_transitions(work_dir, metrics, col_depth=1):
    """Create the transition analysis for a file."""
    with open(work_dir + 'aux', 'rb') as aux_file:
        aux_data = pickle.load(aux_file)
    prefixes = set([x[:i] for x in aux_data.firsts for i in range(col_depth)])
    for metric in metrics:
        metric_opts = metric[1]
        for prefix in prefixes:
            file_name = metric_unique_name(metric, prefix)
            if os.path.exists(work_dir + file_name + '.transitions'):
                print('skipping transitions for {0}, because file exists'.
                      format(file_name))
                continue
            trs = an.analyze_transitions(aux_data, metric, prefix)
            with open(work_dir + file_name + '.transitions', 'xb') as out_file:
                pickle.dump(trs, out_file)
            # write the new penalty log if it exists
            if 'new_penalty_log' in metric_opts:
                with open(work_dir + file_name + '.new_penalty_log',
                          'xb') as out_file:
                    pickle.dump(metric_opts['new_penalty_log'], out_file)


def write_trivial_tour(path, transitions):
    if os.path.exists(path + '.tour'):
        print('skipping {0}, because file {0}.tour exists'.format(path))
        return
    tour_text = 'NAME : bwt.trivial.tour\n'
    tour_text += 'TYPE : TOUR\n'
    tour_text += 'DIMENSION : 256\n'
    tour_text += 'TOUR_SECTION\n'
    numbers_to_names = {}
    # either there's 2 nodes and transitions between them, or no nodes
    if len(transitions) == 2:
        sorted_nodes = sorted(transitions.keys(), key=lambda x:transitions[x])
        tour_text += '1\n2\n'
        numbers_to_names[1] = sorted_nodes[0]
        numbers_to_names[2] = sorted_nodes[1]
    else:
        # write natural tour
        for i in range(1, 257):
            tour_text += '{0}\n'.format(i)
            numbers_to_names = {i + 1:bytes([i]) for i in range(256)}
    tour_text += '-1\n'
    tour_text += 'EOF'
    with open(path + '.tour', 'xt') as tour_file:
        tour_file.write(tour_text)
    with open(path + '.nodenames', 'xb') as names_file:
        pickle.dump(numbers_to_names, names_file)

def write_tsplib_files(work_dir, metrics, print_rel_error=False):
    """Create and write the .tsp and .par files for the LKH program as well as
    the file mapping the LKH node ids back to the node names."""
    file_names_list = [find_metric_file_names(work_dir, metric)
                             for metric in metrics]
    file_names = list(chain.from_iterable(file_names_list))
    for file_name in file_names:
        if os.path.exists(work_dir + file_name + '.par'):
            print('skipping {0}, because file {0}.par exists'.format(file_name))
            continue
        with open(work_dir + file_name + '.transitions', 'rb') as trs_file:
            transitions = pickle.load(trs_file)
        # name mapping part
        # make dict number->nodename for the tsp file
        nodes = set([t[0] for t in transitions])
        num_nodes = len(nodes)
        if num_nodes < 3:
            # trivial case, no tsp necessary (LKH actually can't handle this)
            write_trivial_tour(work_dir + file_name, transitions)
            continue
        numbers_to_names = {}
        for i, n in enumerate(sorted(nodes), 1):
            numbers_to_names[i] = n

        # tsp file part
        # start with the header
        tsp_text = 'NAME: bwt\n'
        tsp_text += 'TYPE: ATSP\n'
        tsp_text += 'DIMENSION: {0}\n'.format(num_nodes)
        tsp_text += 'EDGE_WEIGHT_TYPE: EXPLICIT\n'
        tsp_text += 'EDGE_WEIGHT_FORMAT: FULL_MATRIX\n'
        tsp_text += 'EDGE_WEIGHT_SECTION\n'

        MAX_INT = 10 ** 7  # the maximum value to write to the tsp file
        INFINITY = 2147483648  # value to signify infinity in the tsp file
        # need to scale up all floats to integers
        # min and max trs costs, and the ratio
        costs_nonzero = [abs(c) for c in transitions.values() if c != 0]
        if costs_nonzero:
            min_cost = min(costs_nonzero)
        else:
            min_cost = 0
        absolute_values = [abs(c) for c in transitions.values()]
        if absolute_values:
            max_cost = max(absolute_values)
        else:
            max_cost = 0
        if min_cost:
            ratio = max_cost / min_cost
        else:
            ratio = float('inf')
        print('making tsplib for {0}'.format(file_name))
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

        new_values = []
        max_len = len(str(INFINITY))  # length of the longest number
        # append the scaled cost value matrix to the tsp text
        for irow in range(1, num_nodes + 1):
            for icol in range(1, num_nodes + 1):
                if irow == icol:
                    # node to itself -> infinity
                    write_num = INFINITY
                else:
                    trs = (numbers_to_names[irow], numbers_to_names[icol])
                    write_num = int(transitions[trs] * factor)
                    new_values.append((transitions[trs], write_num))
                # leave 2 characters spacing
                space = ' ' * (max_len - len(str(write_num)) + 2)
                # append the string
                tsp_text += space + str(write_num)
            # newline at the end of a row
            tsp_text += '\n'
        # EOF at the end of the file
        tsp_text += 'EOF'
        if print_rel_error:
            rel_errors = ((abs((t1[1] / t2[1]) -  # @UndefinedVariable weird PyDev bug
                               (t1[0] / t2[0])) / (t1[0] / t2[0]))  # @UndefinedVariable
                            for t1 in (t for t in new_values if t != (0, 0.0))
                            for t2 in (t for t in new_values if t != (0, 0.0)))
            try:
                max_error = max(rel_errors)
                print('max relative error: {0}'.format(max_error))
            except ValueError:
                # argument to max is empty
                print('all elements zero, no scaling, no errors')
            print()

        # par file part
        par_text = 'PROBLEM_FILE = {0}\n'.format(file_name + '.atsp')
        par_text += 'RUNS = 100\n'
        par_text += 'TOUR_FILE = ' + file_name + '.tour'

        # write the files
        with open(work_dir + file_name + '.atsp', 'xt') as tsp_file:
            tsp_file.write(tsp_text)
        with open(work_dir + file_name + '.par', 'xt') as par_file:
            par_file.write(par_text)
        with open(work_dir + file_name + '.nodenames', 'xb') as names_file:
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
            tour.append(names_dict[number][-1:])
    return tour

def print_simulated_compression_results(work_dir, metrics, in_file_path,
                                        mtf_exc_min_length=None,
                                        mtf_exc_threshold=None):
    """Simulate compression of a file and print achieved compression."""
    def final_bit_len(orders):
        bw = cd.bw_encode(bs, orders)
        if mtf_exc_min_length is not None and mtf_exc_threshold is not None:
            mtf_exceptions = an.select_mtf_exceptions(bw, mtf_exc_min_length,
                                                      mtf_exc_threshold)
            # exclude mtf exceptions from the bw code
            excepted_blocks = [an.context_block(bw.encoded, bw.firsts, e)
                               for e in mtf_exceptions]
            # make the bw code without the exceptions
            bw_ints = []
            for f, s in zip(bw.firsts, bw.encoded):
                if not any(f[:len(e)] == e for e in mtf_exceptions):
                    bw_ints.append(s)
            bw_code = bytes(bw_ints)
            excepted_huff_codes = [hf.encode_to_bits_static(b)
                                   for b in excepted_blocks]
        else:
            mtf_exceptions = None
            excepted_huff_codes = []
            bw_code = bw.encoded
        mtf_code = cd.mtf_encode(bw_code)
        huff_code_mtf = hf.encode_to_bits_static(mtf_code)
        size = len(huff_code_mtf) + sum(len(c) for c in excepted_huff_codes)
        result = str(size)
        if mtf_exceptions:
            result += ' (mtf exceptions: {0})'.format(mtf_exceptions)
        return result
    with open(in_file_path, 'rb') as in_file:
        bs = in_file.read()
    handpicked_str = b'aeioubcdgfhrlsmnpqjktwvxyzAEIOUBCDGFHRLSMNPQJKTWVXYZ'
    handpicked_orders = [[bytes([c]) for c in handpicked_str]]
    natural_order = [bytes([x]) for x in range(256)]
    print('simulating compression for file {0}'.format(in_file_path))
    if mtf_exc_min_length is not None and mtf_exc_threshold is not None:
        print('minimum length of an mtf block for an mtf exception: {0}'
              .format(mtf_exc_min_length))
        print('threshold for the mean mtf value in a block for an mtf '
              'exception: {0}'.format(mtf_exc_threshold))
    print('in size: {0}'.format(len(bs) * 8))
    print()
    print('natural order: {0}'.format(final_bit_len([natural_order])))
    print()
    print('handpicked order aeiou...:')
    print('all columns      : {0}'.format(final_bit_len(handpicked_orders)))
    handpicked_orders.append(natural_order)
    print('first column only: {0}'.format(final_bit_len(handpicked_orders)))
    print()
    for metric in metrics:
        file_name = metric_unique_name(metric, b'')
        print('{0}:'.format(file_name))
        orders = assemble_multicol_orders(work_dir, metric)
        print('using the last specified as default: {0}'
              .format(final_bit_len(orders)))
        orders.append(natural_order)
        print('using the natural order as default : {0}'
              .format(final_bit_len(orders)))
        print()

def print_mtf_prediction_evaluations(work_dir, metrics):
    with open(work_dir + 'aux', 'rb') as aux_file:
        aux_data = pickle.load(aux_file)
    for metric in metrics:
        if metric[0] != 'badness':
            continue
        file_name = metric_unique_name(metric, b'')
        with open(work_dir + file_name + '.new_penalty_log', 'rb') as in_file:
            new_penalty_log = pickle.load(in_file)
        print('metric: {0}'.format(file_name))
        tsplib_tour = read_tsplib_files(work_dir + file_name + '.tour',
                                         work_dir + file_name + '.nodenames')
        natural_order = [bytes([x]) for x in range(256)]

        for all_cols in [False, True]:
            if all_cols:
                print('  using the computed order for all columns:')
                orders = [tsplib_tour]
            else:
                print('  using the computed order for only the first column:')
                orders = [tsplib_tour, natural_order]
            comp = an.compare_new_penalty_predictions(aux_data, orders,
                                                      new_penalty_log)
            mtf_code = cd.mtf_encode(cd.bw_encode(aux_data.raw, orders).encoded)
            freqs = hf.symbol_frequencies(mtf_code)
            hf_len = hf.codeword_lengths(freqs)
            for subst_hf in [False, True]:
                if subst_hf:
                    print('    error of huffman codes:')
                else:
                    print('    error of mtf codes:')
                for filter_new in[False, True]:
                    if filter_new:
                        print('      not considering errors where the actual '
                              'value fetches a new symbol to the front:')
                    else:
                        print('      considering all errors:')
                    for filter_sides in ['all', 'ge', 'lt']:
                        if filter_sides == 'all':
                            print('        considering both sides of the actual '
                                  'value:')
                        elif filter_sides == 'ge':
                            print('        considering values where the '
                                  'prediction is >= the actual value:')
                        elif filter_sides == 'lt':
                            print('        considering values where the '
                                  'prediction is < the actual value:')
                        filtered_comp = comp[:]
                        if filter_sides == 'ge':
                            filtered_comp = [c for c in filtered_comp
                                             if c[2] >= c[1]]
                        elif filter_sides == 'lt':
                            filtered_comp = [c for c in filtered_comp
                                             if c[2] < c[1]]
                        if filter_new:
                            filtered_comp = [c for c in filtered_comp
                                    if c[1] < aux_data.num_symbols]
                        if subst_hf:
                            filtered_comp = [(c[0], hf_len[c[1]],
                                     hf_len[min(hf_len.keys(),
                                           key=lambda x:abs(x - c[2]))])
                                         for c in filtered_comp]
                        if not filtered_comp:
                            print('          no data matches these criteria')
                            continue
                        diffs = [c[1] - c[2] for c in filtered_comp]
                        devs = list(map(abs, diffs))
                        variance = np.mean(list(map(lambda x:x ** 2, diffs)))
                        std_deviation = math.sqrt(variance)
                        print('          mean differences: {0}'
                              .format(np.mean(diffs)))
                        print('          mean deviation: {0}'
                              .format(np.mean(devs)))
                        print('          standard deviation: {0}'
                              .format(std_deviation))
        print()

def print_entropy_length_prediction_evaluations(work_dir, metrics):
    with open(work_dir + 'aux', 'rb') as aux_file:
        aux_data = pickle.load(aux_file)
    for metric in metrics:
        if 'entropy_code_len' in metric[1]:
            if metric[1]['entropy_code_len'] == 'complete':
                prediction = aux_data.huffman_codeword_lengths_complete
            elif metric[1]['entropy_code_len'] == 'sparse':
                prediction = aux_data.huffman_codeword_lengths_sparse
            else:
                continue
        else:
            continue
        file_name = metric_unique_name(metric, b'')
        print('metric: {0}'.format(file_name))
        tsplib_tour = read_tsplib_files(work_dir + file_name + '.tour',
                                         work_dir + file_name + '.nodenames')
        natural_order = [bytes([x]) for x in range(256)]
        for all_cols in [False, True]:
            if all_cols:
                print('  using the computed order for all columns:')
                orders = [tsplib_tour]
            else:
                print('  using the computed order for only the first column:')
                orders = [tsplib_tour, natural_order]
            comp = an.compare_entropy_len_predictions(aux_data, orders,
                                               prediction)
            for filter_sides in ['all', 'ge', 'lt']:
                if filter_sides == 'all':
                    print('    considering both sides of the actual '
                          'value:')
                elif filter_sides == 'ge':
                    print('    considering values where the '
                          'prediction is >= the actual value:')
                elif filter_sides == 'lt':
                    print('    considering values where the '
                          'prediction is < the actual value:')
                filtered_comp = comp[:]
                if filter_sides == 'ge':
                    filtered_comp = [c for c in filtered_comp if c[2] >= c[1]]
                elif filter_sides == 'lt':
                    filtered_comp = [c for c in filtered_comp if c[2] < c[1]]
                if not filtered_comp:
                    print('      no data matches these criteria')
                    continue
                diffs = [c[1] - c[2] for c in filtered_comp]
                devs = list(map(abs, diffs))
                variance = np.mean(list(map(lambda x:x ** 2, diffs)))
                std_deviation = math.sqrt(variance)
                print('      mean differences: {0}'
                      .format(np.mean(diffs)))
                print('      mean deviation: {0}'
                      .format(np.mean(devs)))
                print('      standard deviation: {0}'
                      .format(std_deviation))
        print()

def metric_unique_name(metric, prefix):
    """Create a unique name for a metric with options."""
    elems = [metric[0]]
    opts = metric[1]
    for opt in sorted([o for o in opts if opts[o]]):
        if opt == 'new_penalty_log':
            # don't record the log as an option
            continue
        elems.append(opt)
        if opts[opt] != True:
            elems.append(str(opts[opt]))
    result = '_'.join(elems)
    if prefix:
        result += '.' + '.'.join([str(b) for b in prefix])
    return result

def find_metric_file_names(work_dir, metric):
    """Find base file names for all prefixes for metric in work_dir."""
    metric_name = metric_unique_name(metric, b'')
    files = os.listdir(path=work_dir)
    metric_files = [f for f in files if f.split('.')[0] == metric_name]
    return list(set([f.rsplit(sep='.', maxsplit=1)[0] for f in metric_files]))

def assemble_multicol_orders(work_dir, metric):
    base_file_name = metric_unique_name(metric, b'')
    file_name_splits = [f.split('.')
                        for f in find_metric_file_names(work_dir, metric)]
    prefixes = [bytes(map(lambda x:int(x), p[1:])) for p in file_name_splits]
    col_depth = max((len(s) for s in file_name_splits))
    orders = [read_tsplib_files(work_dir + base_file_name + '.tour',
                                work_dir + base_file_name + '.nodenames')]
    for i in range(1, col_depth):
        order = {}
        for prefix in [p for p in prefixes if len(p) == i]:
            file_name = base_file_name + '.' + '.'.join([str(b)for b in prefix])
            tour = read_tsplib_files(work_dir + file_name + '.tour',
                                     work_dir + file_name + '.nodenames')
            order[prefix] = tour
        orders.append(order)
    return orders

if __name__ == '__main__':
    start_time = time.time()
    work_dir = '/home/dddsnn/tmp/book1/'
    in_file_path = '/home/dddsnn/Dokumente/Studium/BA/calgary/book1'
    metrics = [('chapin_hst_diff', {}), ('chapin_inv', {}),
               ('chapin_inv', {'log':True})]
#     metrics = []
    for w in [True, False]:
        for entr_len in [False, 'complete', 'sparse']:
            for new_pen in [False, 'generic_mean', 'generic_median',
                            'specific_mean', 'specific_median']:
                opts = {'weighted':w, 'entropy_code_len':entr_len,
                        'new_penalty':new_pen, 'new_penalty_log':{}}
                metrics.append(('badness', opts))
#     metrics = [('badness', {'new_penalty': False, 'entropy_code_len': False,
#                             'weighted': False, 'new_penalty_log': {},
#                             'mtf_prediction_correction':14.76774193548387}),
#                ('badness', {'new_penalty': 'generic_mean',
#                             'entropy_code_len': False, 'weighted': False,
#                             'new_penalty_log': {},
#                             'mtf_prediction_correction':7.016073054355908}),
#                ('badness', {'new_penalty': 'generic_median',
#                             'entropy_code_len': False,
#                             'weighted': False, 'new_penalty_log': {},
#                             'mtf_prediction_correction':9.033012379642367}),
#                ('badness', {'new_penalty': 'specific_mean',
#                             'entropy_code_len': False,
#                             'weighted': False,
#                             'new_penalty_log': {},
#                             'mtf_prediction_correction':5.2337720496925035}),
#                ('badness', {'new_penalty': 'specific_median',
#                             'entropy_code_len': False,
#                             'weighted': False, 'new_penalty_log': {},
#                             'mtf_prediction_correction':6.9992559523809526})]
#     metrics = [('badness', {'new_penalty': False, 'entropy_code_len': False,
#                             'weighted': False, 'new_penalty_log': {}})]

#     metrics = []

#     make_aux_data(work_dir, in_file_path, col_depth=1)

#     make_transitions(work_dir, metrics, col_depth=1)

#     write_tsplib_files(work_dir, metrics, print_rel_error=True)

#     print_simulated_compression_results(work_dir, metrics, in_file_path)

    print_mtf_prediction_evaluations(work_dir, metrics)

#     print_entropy_length_prediction_evaluations(work_dir, metrics)

#     natural_order = [bytes([x]) for x in range(256)]
#     orders = [natural_order, list(reversed(natural_order))]
    # [15,24,25,14,26,13,17,21,4,11,9,12,16,20,6,0,27,7,1,28,8,2,30,3,18,19,22,29,5,10,23]
#     in_bs = b'missishjkdgfhjkdasdasdasdjklsdg'
#     in_bs = b'Lorem ipsum dolor sit amet, consectetur adipiscing elit.'
#     in_bs = b'mississippi'
#     with open(in_file_path, 'rb') as in_file:
#         in_bs = in_file.read()
#     in_bs = in_bs[:10000]
#     orders = assemble_multicol_orders(work_dir, metrics[10])
#     bw = cd.bw_encode(in_bs, orders)
#     dec = cd.bw_decode_2_orders(bw.encoded, bw.start_index, orders)
#     print(in_bs == dec)
#     print(dec)

    print('time: {0:.0f}s'.format(time.time() - start_time))
