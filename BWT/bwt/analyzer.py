'''
Created on 28.05.2014

@author: dddsnn
'''

from bwt import *
import bwt.coder as cd
import math

def analyze_partial_mtf(code):
    raw = code
    length = len(code)
    sorted_code = sorted(code)
    max_code = sorted_code[-1]
    # remove -1s
    sorted_code = [c for c in sorted_code if c != -1]
    # number of all recurring symbols (those that aren't -1)
    l = len(sorted_code)
    # number of -1s is number of all symbols minus number of all except -1s
    num_chars = length - l
    # number of recurring characters (not -1)
    length_rec = length - num_chars
    if l != 0:
        if l % 2 == 0:
            median = (sorted_code[(l // 2) - 1] +
                                 sorted_code[l // 2]) / 2
        else:
            median = sorted_code[l // 2]
        mean = sum(sorted_code) / l
    else:
        median = 0
        mean = 0
    result = PartialMTFAnalysisResult(raw, length, length_rec, num_chars,
                                      max_code, median, mean)
    return result

def analyze_transition(bw_code, first_symbol_a, first_symbol_b):
    '''Analyze pt_mtf_a single transition between two BW encoded strings.'''
    # get the two bw subcodes corresponding to the first symbols
    left_a = bw_code.firsts.index(first_symbol_a)
    right_a = bw_code.firsts.rindex(first_symbol_a)
    bw_a = bw_code.encoded[left_a:right_a + 1]
    left_b = bw_code.firsts.index(first_symbol_b)
    right_b = bw_code.firsts.rindex(first_symbol_b)
    bw_b = bw_code.encoded[left_b:right_b + 1]
    # make the partial mtf codes
    pt_mtf_a = analyze_partial_mtf(cd.mtf_partial_enc(bw_a))
    pt_mtf_b = analyze_partial_mtf(cd.mtf_partial_enc(bw_b))
    pt_mtf_ab = analyze_partial_mtf(cd.mtf_partial_enc(bw_a + bw_b))

    # LENGTH
    # metric: difference between length of left vs. length of right
    length = TransitionDataSet(pt_mtf_a.length, pt_mtf_b.length,
                               pt_mtf_ab.length)
    length_metric = abs(pt_mtf_a.length - pt_mtf_b.length)

    # NUMBER OF CHARACTERS
    # metric: number of -1s in the second half of the combined code (the one
    # that's being transitioned to)
    num_chars = TransitionDataSet(pt_mtf_a.num_chars, pt_mtf_b.num_chars,
                                  pt_mtf_ab.num_chars)
    num_chars_metric = pt_mtf_ab.num_chars - pt_mtf_a.num_chars

    # MAX CODE
    max_code = TransitionDataSet(pt_mtf_a.max_code, pt_mtf_b.max_code,
                                 pt_mtf_ab.max_code)
    max_code_metric = pt_mtf_ab.max_code - max(pt_mtf_a.max_code,
                                               pt_mtf_b.max_code)
    # to avoid division by zero
    if pt_mtf_ab.length_rec == 0:
        ab_length_rec = 1
    else:
        ab_length_rec = pt_mtf_ab.length_rec

    # MEDIAN
    median = TransitionDataSet(pt_mtf_a.median, pt_mtf_b.median,
                               pt_mtf_ab.median)
    # metric: difference of expected median and the achieved median
    median_metric = ((pt_mtf_a.median * pt_mtf_a.length_rec +
                      pt_mtf_b.median * pt_mtf_b.length_rec)
                       / ab_length_rec) - pt_mtf_ab.median
    # TODO by ignoring new characters (-1s) in the mean, good transitions get pt_mtf_a
    # penalty, because pt_mtf_a character that was already there before the transition
    # probably affects the mean in pt_mtf_a bad way, while an entirely new character
    # won't affect it at all
    # need to give new characters pt_mtf_a "penalty" value in the partial mtf encode
    # (but only for the target of the transition, lots of new characters in the
    # source don't mean anything bad)
    # or maybe an entirely new metric that takes this into account

    # MEAN
    mean = TransitionDataSet(pt_mtf_a.mean, pt_mtf_b.mean, pt_mtf_ab.mean)
    # metric: difference of expected mean and the achieved mean
    mean_metric = ((pt_mtf_a.mean * pt_mtf_a.length_rec +
                    pt_mtf_b.mean * pt_mtf_b.length_rec)
                     / ab_length_rec) - pt_mtf_ab.mean

    # CHAPIN: sum of squares of differences of logs
    hst_a = make_histogram(bw_a)
    hst_b = make_histogram(bw_b)
    log_diffs = []
    for k in hst_a:
        # get the logs
        # TODO not sure replacing -inf with 0 doesn't affect result
        if hst_a[k] == 0:
            log_a = 0
        else:
            log_a = math.log(hst_a[k])
        if hst_b[k] == 0:
            log_b = 0
        else:
            log_b = math.log(hst_b[k])
        log_diffs.append(log_a - log_b)
    # squares of differences of logarithms in the histograms
    chapin_hst_diff_metric = sum([d ** 2 for d in log_diffs])

    # CHAPIN: kullback-leibler
    logterms = []
    for k in hst_a:
        if hst_a[k] == 0 or hst_b[k] == 0:
            # in case of zeros, ignore this term
            # TODO http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
            continue
        logterms.append(math.log(hst_b[k] / hst_a[k]) * hst_b[k])
    chapin_kl_metric = sum(logterms)
    # CHAPIN: number of inversion between ordered histograms
    # turn the histograms into lists and sort in decreasing order of frequency
    freq_list_a = sorted(hst_a.items(), key=lambda x:x[1], reverse=True)
    freq_list_b = sorted(hst_b.items(), key=lambda x:x[1], reverse=True)
    # now just take the corresponding first symbols
    freq_list_a = [x[0] for x in freq_list_a]
    freq_list_b = [x[0] for x in freq_list_b]
    # metric is the number of inversions between the two lists
    chapin_inv_metric = num_inversions(freq_list_a, freq_list_b)

    # CHAPIN: log of previous
    if chapin_inv_metric == 0:
        chapin_inv_log_metric = 0
    else:
        chapin_inv_log_metric = math.log(chapin_inv_metric)

    data = TransitionData(length, num_chars, max_code, median, mean)
    result = TransitionAnalysisResult(data, length_metric, num_chars_metric,
                                      max_code_metric, median_metric,
                                      mean_metric, chapin_hst_diff_metric,
                                      chapin_kl_metric, chapin_inv_metric,
                                      chapin_inv_log_metric)
    return result

def num_inversions(list_a, list_b):
    '''Compute the number of inversions between two lists.'''
    inv = 0
    for i, x in enumerate(list_a):
        # TODO need to recheck that i'm not counting anything multiple times
        try:
            i_b = list_b.index(x)
        except ValueError:
            # no inversions if value isn't in hst_b, continue with next
            continue
        # all the items coming before x in hst_b
        before_list = list_b[:i_b]
        # all the items coming after x in hst_a, but not equal to x
        after_list = list_a[i + 1:]
        after_list = [a for a in after_list if a != x]
        while after_list:
            try:
                # get the index of the first element of after_list in
                # before_list
                idx_in_before = before_list.index(after_list[0])
                # inversion found, increment inv and delete the element
                # from both lists
                inv += 1
                del before_list[idx_in_before]
                del after_list[0]
            except ValueError:
                # element not found in before_list, delete from after_list
                # and continue
                del after_list[0]
                continue
    return inv

def make_histogram(bytes_):
    '''Make a histogram of symbol appearances for a bytes object.'''
    # initialize 0 for all possible bytes
    histogram = {i:0 for i in range(256)}
    for b in bytes_:
        histogram[b] += 1
    return histogram

def analyze_transitions(bytes_):
    '''Analyze all the transitions between bytes of a byte string.'''
    bw_code = cd.bw_encode(bytes_)
    firsts = set(bw_code.firsts)
    transitions = {(a, b): analyze_transition(bw_code, a, b)
                   for a in firsts
                   for b in firsts if a != b}
    return transitions

def make_table_string(table):
    result = ''
    col_width = [max(len(str(x)) for x in col) for col in zip(*table)]
    for line in table:
        result += "| " + " | ".join("{:{}}".format(x, col_width[i])
                                for i, x in enumerate(line)) + " |\n"
    return result

def print_transition_analyses(bytes_):
    trs = analyze_transitions(bytes_)
    for k in sorted(trs.keys()):
        print('<' + k[0] + '-' + k[1] + '>\n')
        header = ['', 'left', 'right', 'together', 'metric']
        len_line = ['lenght', trs[k].data.length.left,
                    trs[k].data.length.right,
                    trs[k].data.length.together, trs[k].length]
        num_chars_line = ['num_chars', trs[k].data.num_chars.left,
                          trs[k].data.num_chars.right,
                          trs[k].data.num_chars.together,
                          trs[k].num_chars]
        max_line = ['max_code', trs[k].data.max_code.left,
                    trs[k].data.max_code.right,
                    trs[k].data.max_code.together,
                    trs[k].max_code]
        median_line = ['median', trs[k].data.median.left,
                       trs[k].data.median.right,
                       trs[k].data.median.together, trs[k].median]
        mean_line = ['mean', trs[k].data.mean.left,
                     trs[k].data.mean.right,
                     trs[k].data.mean.together, trs[k].mean]
        table = [header, len_line, num_chars_line, max_line, median_line,
                 mean_line]
        print(make_table_string(table))
