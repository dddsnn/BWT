from bwt import *
import bwt.coder as cd
import bwt.analyzer as an
import bwt.huffman as hf
import pickle
import numpy as np
import time
import math
import os

def make_aux_data(work_dir, in_file_path):
    """Create data commonly used in transition analysis."""
    with open(in_file_path, 'rb') as in_file:
        bs = in_file.read()
    raw = bs
    num_symbols = len(set(bs))
    bw_code = cd.bw_encode(bs)
    mtf_code = cd.mtf_enc(bw_code.encoded)
    firsts = list(set([bytes([x[0]]) for x in bw_code.firsts]))
    firsts.extend(an.select_sequences(bs, 2))
    specializations = cd.specializations(firsts)
    bw_subcodes = {x:an.context_block(bw_code.encoded, bw_code.firsts, x,
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

def make_transitions(work_dir, metrics):
    """Create the transition analysis for a file."""
    with open(work_dir + 'aux', 'rb') as aux_file:
        aux_data = pickle.load(aux_file)
    for metric in metrics:
        metric_name = metric_unique_name(metric)
        if os.path.exists(work_dir + metric_name + '.transitions'):
            print('skipping transitions for {0}, because file exists'.
                  format(metric_name))
            continue
        metric_opts = metric[1]
        trs = an.analyze_transitions(aux_data, metric)
        with open(work_dir + metric_name + '.transitions', 'xb') as out_file:
            pickle.dump(trs, out_file)
        # write the new penalty log if it exists
        if 'new_penalty_log' in metric_opts:
            with open(work_dir + metric_name + '.new_penalty_log',
                      'xb') as out_file:
                pickle.dump(metric_opts['new_penalty_log'], out_file)

def write_tsplib_files(work_dir, metrics):
    """Create and write the .tsp and .par files for the LKH program as well as
    the file mapping the LKH node ids back to the node names."""
    file_names = [metric_unique_name(metric) for metric in metrics]
    for file_name in file_names:
        if os.path.exists(work_dir + file_name + '.par'):
            print('skipping {0}, because file {0}.par exists'.format(file_name))
            continue
        with open(work_dir + file_name + '.transitions', 'rb') as trs_file:
            transitions = pickle.load(trs_file)
        # name mapping part
        # make dicts nodename->number and number->nodename for the tsp file
        nodes = set([t[0] for t in transitions])
        num_nodes = len(nodes)
        names_to_numbers = numbers_to_names = {}
        for i, n in enumerate(sorted(nodes), 1):
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
        # min and max trs costs, and the ratio
        costs_nonzero = [abs(c) for c in transitions.values() if c != 0]
        if costs_nonzero:
            min_cost = min(costs_nonzero)
        else:
            min_cost = 0
        max_cost = max([abs(c) for c in transitions.values()])
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
                    trs = (numbers_to_names[irow], numbers_to_names[icol])
                    write_num = int(transitions[trs] * factor)
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
            tour.append(names_dict[number])
    return tour

def print_simulated_compression_results(work_dir, metrics, in_file_path):
    """Simulate compression of a file and print achieved compression."""
    def final_bit_len(bs, orders):
        bw_code = cd.bw_encode(bs, orders)
        mtf_code = cd.mtf_enc(bw_code.encoded)
        huff_code = hf.encode_to_bits_static(mtf_code)
        return len(huff_code)
    with open(in_file_path, 'rb') as in_file:
        bs = in_file.read()
    handpicked_str = b'aeioubcdgfhrlsmnpqjktwvxyzAEIOUBCDGFHRLSMNPQJKTWVXYZ'
    handpicked_orders = [[bytes([c]) for c in handpicked_str]]
    natural_order = [bytes([x]) for x in range(256)]
    print('simulating compression for file {0}'.format(in_file_path))
    print('in size: {0}'.format(len(bs) * 8))
    print()
    print('natural order: {0}'.format(final_bit_len(bs, [natural_order])))
    print()
    print('handpicked order aeiou...:')
    print('all columns      : {0}'.format(final_bit_len(bs, handpicked_orders)))
    handpicked_orders.append(natural_order)
    print('first column only: {0}'.format(final_bit_len(bs, handpicked_orders)))
    print()
    for metric in metrics:
        file_name = metric_unique_name(metric)
        print('{0}:'.format(file_name))
        orders = [read_tsplib_files(work_dir + file_name + '.tour',
                                   work_dir + file_name + '.nodenames')]
        print('all columns      : {0}'.format(final_bit_len(bs, orders)))
        orders.append(natural_order)
        print('first column only: {0}'.format(final_bit_len(bs, orders)))
        print()

def print_mtf_prediction_evaluations(work_dir, metrics):
    with open(work_dir + 'aux', 'rb') as aux_file:
        aux_data = pickle.load(aux_file)
    for metric in metrics:
        if metric[0] != 'badness':
            continue
        file_name = metric_unique_name(metric)
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
            mtf_code = cd.mtf_enc(cd.bw_encode(aux_data.raw, orders).encoded)
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
                        if filter_sides == 'ge':
                            comp = [c for c in comp if c[2] >= c[1]]
                        elif filter_sides == 'lt':
                            comp = [c for c in comp if c[2] < c[1]]
                        if filter_new:
                            comp = [c for c in comp
                                    if c[1] < aux_data.num_symbols]
                        if subst_hf:
                            comp = [(c[0], hf_len[c[1]],
                                     hf_len[min(hf_len.keys(),
                                           key=lambda x:abs(x - c[2]))])
                                         for c in comp]
                        if not comp:
                            print('          no data matches these criteria')
                            continue
                        diffs = [c[1] - c[2] for c in comp]
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
        file_name = metric_unique_name(metric)
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
                if filter_sides == 'ge':
                    comp = [c for c in comp if c[2] >= c[1]]
                elif filter_sides == 'lt':
                    comp = [c for c in comp if c[2] < c[1]]
                if not comp:
                    print('      no data matches these criteria')
                    continue
                diffs = [c[1] - c[2] for c in comp]
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
            elems.append(str(opts[opt]))
    return '_'.join(elems)

if __name__ == '__main__':
    start_time = time.time()
    work_dir = '/home/dddsnn/tmp/book1/'
    in_file_path = '/home/dddsnn/Dokumente/Studium/BA/calgary/book1'
    metrics = [('chapin_hst_diff', {}), ('chapin_inv', {}),
               ('chapin_inv', {'log':True})]
    for w in [False]:  # [True, False]:
        for entr_len in [False]:  # [False, 'complete', 'sparse']:
            for new_pen in [False, 'generic_mean', 'generic_median',
                            'specific_mean', 'specific_median']:
                opts = {'weighted':w, 'entropy_code_len':entr_len,
                        'new_penalty':new_pen, 'new_penalty_log':{}}
                metrics.append(('badness', opts))

    tmp_metrics = []
    for metric in metrics:
        if not 'new_penalty' in metric[1]:
            continue
        if metric[1]['new_penalty'] == False:
            tmp = (metric[0], metric[1].copy())
            tmp[1]['mtf_prediction_correction'] = 20.629411764705882
            tmp_metrics.append(tmp)
        elif metric[1]['new_penalty'] == 'generic_mean':
            tmp = (metric[0], metric[1].copy())
            tmp[1]['mtf_prediction_correction'] = 14.011834428781349
            tmp_metrics.append(tmp)
        elif metric[1]['new_penalty'] == 'generic_median':
            tmp = (metric[0], metric[1].copy())
            tmp[1]['mtf_prediction_correction'] = 14.692065491183879
            tmp_metrics.append(tmp)
        elif metric[1]['new_penalty'] == 'specific_mean':
            tmp = (metric[0], metric[1].copy())
            tmp[1]['mtf_prediction_correction'] = 12.722633592206217
            tmp_metrics.append(tmp)
        elif metric[1]['new_penalty'] == 'specific_median':
            tmp = (metric[0], metric[1].copy())
            tmp[1]['mtf_prediction_correction'] = 13.904411764705882
            tmp_metrics.append(tmp)
    metrics.extend(tmp_metrics)

#     make_aux_data(work_dir, in_file_path)

#     make_transitions(work_dir, metrics)

#     write_tsplib_files(work_dir, metrics)

#     print_simulated_compression_results(work_dir, metrics, in_file_path)

    print_mtf_prediction_evaluations(work_dir, metrics)

#     print_entropy_length_prediction_evaluations(work_dir, metrics)

    print('time: {0:.0f}s'.format(time.time() - start_time))
