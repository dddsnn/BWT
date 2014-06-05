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

def bw_block(bw_code, first_symbol):
    '''Get the BW code corresponding to a specific symbol in the first column.
    '''
    left = bw_code.firsts.index(first_symbol)
    right = bw_code.firsts.rindex(first_symbol)
    bw_block = bw_code.encoded[left:right + 1]
    return bw_block

def metric_num_chars(bw_code, first_symbol_a, first_symbol_b):
    # get the two bw subcodes corresponding to the first symbols
    bw_a = bw_block(bw_code, first_symbol_a)
    bw_b = bw_block(bw_code, first_symbol_b)
    # make the partial mtf codes
    pt_mtf_a = cd.mtf_partial_enc(bw_a)
    pt_mtf_ab = cd.mtf_partial_enc(bw_a + bw_b)
    # get the analyses for the parts
    an_a = analyze_partial_mtf(pt_mtf_a)
    an_ab = analyze_partial_mtf(pt_mtf_ab)

    # metric: number of -1s in the second half of the combined code (the one
    # that's being transitioned to)
    metric = an_ab.num_chars - an_a.num_chars
    return metric

def metric_max_code(bw_code, first_symbol_a, first_symbol_b):
    # get the two bw subcodes corresponding to the first symbols
    bw_a = bw_block(bw_code, first_symbol_a)
    bw_b = bw_block(bw_code, first_symbol_b)
    # make the partial mtf codes
    pt_mtf_a = cd.mtf_partial_enc(bw_a)
    pt_mtf_b = cd.mtf_partial_enc(bw_b)
    pt_mtf_ab = cd.mtf_partial_enc(bw_a + bw_b)
    # get the analyses for the parts
    an_a = analyze_partial_mtf(pt_mtf_a)
    an_b = analyze_partial_mtf(pt_mtf_b)
    an_ab = analyze_partial_mtf(pt_mtf_ab)

    metric = an_ab.max_code - max(an_a.max_code, an_b.max_code)
    return metric

def metric_median(bw_code, first_symbol_a, first_symbol_b):
    # get the two bw subcodes corresponding to the first symbols
    bw_a = bw_block(bw_code, first_symbol_a)
    bw_b = bw_block(bw_code, first_symbol_b)
    # make the partial mtf codes
    pt_mtf_a = cd.mtf_partial_enc(bw_a)
    pt_mtf_b = cd.mtf_partial_enc(bw_b)
    pt_mtf_ab = cd.mtf_partial_enc(bw_a + bw_b)
    # get the analyses for the parts
    an_a = analyze_partial_mtf(pt_mtf_a)
    an_b = analyze_partial_mtf(pt_mtf_b)
    an_ab = analyze_partial_mtf(pt_mtf_ab)

    # to avoid division by zero
    if an_ab.length_rec == 0:
        ab_length_rec = 1
    else:
        ab_length_rec = an_ab.length_rec

    # metric: difference of expected median and the achieved median
    metric = ((an_a.median * an_a.length_rec + an_b.median * an_b.length_rec)
                       / ab_length_rec) - an_ab.median
    return metric

def metric_mean(bw_code, first_symbol_a, first_symbol_b):
    # TODO by ignoring new characters (-1s) in the mean, good transitions get an_a
    # penalty, because an_a character that was already there before the transition
    # probably affects the mean in an_a bad way, while an entirely new character
    # won't affect it at all
    # need to give new characters an_a "penalty" value in the partial mtf encode
    # (but only for the target of the transition, lots of new characters in the
    # source don't mean anything bad)
    # or maybe an entirely new metric that takes this into account

    # get the two bw subcodes corresponding to the first symbols
    bw_a = bw_block(bw_code, first_symbol_a)
    bw_b = bw_block(bw_code, first_symbol_b)
    # make the partial mtf codes
    pt_mtf_a = cd.mtf_partial_enc(bw_a)
    pt_mtf_b = cd.mtf_partial_enc(bw_b)
    pt_mtf_ab = cd.mtf_partial_enc(bw_a + bw_b)
    # get the analyses for the parts
    an_a = analyze_partial_mtf(pt_mtf_a)
    an_b = analyze_partial_mtf(pt_mtf_b)
    an_ab = analyze_partial_mtf(pt_mtf_ab)

    # to avoid division by zero
    if an_ab.length_rec == 0:
        ab_length_rec = 1
    else:
        ab_length_rec = an_ab.length_rec

    # metric: difference of expected mean and the achieved mean
    metric = ((an_a.mean * an_a.length_rec + an_b.mean * an_b.length_rec)
                     / ab_length_rec) - an_ab.mean
    return metric

def metric_chapin_hst_diff(bw_code, first_symbol_a, first_symbol_b):
    # get the two bw subcodes corresponding to the first symbols
    bw_a = bw_block(bw_code, first_symbol_a)
    bw_b = bw_block(bw_code, first_symbol_b)

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
    metric = sum([d ** 2 for d in log_diffs])
    return metric

def metric_chapin_kl(bw_code, first_symbol_a, first_symbol_b):
    # get the two bw subcodes corresponding to the first symbols
    bw_a = bw_block(bw_code, first_symbol_a)
    bw_b = bw_block(bw_code, first_symbol_b)
    # make the histograms
    hst_a = make_histogram(bw_a)
    hst_b = make_histogram(bw_b)
    # CHAPIN: kullback-leibler
    logterms = []
    for k in hst_a:
        if hst_a[k] == 0 or hst_b[k] == 0:
            # in case of zeros, ignore this term
            # TODO http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
            continue
        logterms.append(math.log(hst_b[k] / hst_a[k]) * hst_b[k])
    metric = sum(logterms)
    return metric

def metric_chapin_inv(bw_code, first_symbol_a, first_symbol_b):
    # get the two bw subcodes corresponding to the first symbols
    bw_a = bw_block(bw_code, first_symbol_a)
    bw_b = bw_block(bw_code, first_symbol_b)
    # make the histograms
    hst_a = make_histogram(bw_a)
    hst_b = make_histogram(bw_b)

    # CHAPIN: number of inversion between ordered histograms
    # turn the histograms into lists and sort in decreasing order of frequency
    freq_list_a = sorted(hst_a.items(), key=lambda x:x[1], reverse=True)
    freq_list_b = sorted(hst_b.items(), key=lambda x:x[1], reverse=True)
    # now just take the corresponding first symbols
    freq_list_a = [x[0] for x in freq_list_a]
    freq_list_b = [x[0] for x in freq_list_b]
    # metric is the number of inversions between the two lists
    metric = num_inversions(freq_list_a, freq_list_b)
    return metric

def metric_chapin_inv_log(bw_code, first_symbol_a, first_symbol_b):
    inv_metric = metric_chapin_inv(bw_code, first_symbol_a, first_symbol_b)
    # CHAPIN: log of previous
    if inv_metric == 0:
        metric = 0
    else:
        metric = math.log(inv_metric)
    return metric

def metric_huffman(bw_code, first_symbol_a, first_symbol_b):
    # get the two bw subcodes corresponding to the first symbols
    bw_a = bw_block(bw_code, first_symbol_a)
    bw_b = bw_block(bw_code, first_symbol_b)
    pt_mtf_ab = cd.mtf_partial_enc(bw_a + bw_b)
    # HUFFMAN ENCODE METRIC
    # TODO -1s in the second part need to be weighted, see todo in mean
    # just encode the partial mtf (without -1s) with huffman and take the length
    # strip -1s from the partial mtf
    stripped_partial_mtf = list(filter(lambda x: x != -1, pt_mtf_ab))
    huffman_length = huffman_enc_length(bytes(stripped_partial_mtf))
    if len(stripped_partial_mtf) == 0:
        # if the stripped partial mtf code has length 0 it means there are no
        # recurring characters in it. since that means no compression benefits
        # for this transition, give it a bad value to penalize it
        # good transitions will have lengths far below 1, so 10 should be enough
        # to discourage tsp from including this transition
        metric = 10
    else:
        metric = huffman_length / len(stripped_partial_mtf)
    return metric

def huffman_enc_length(bytes_):
    '''Return the length of the huffman coded input.
    Use this wrapper if you just need the length because it catches some rare
    cases that will crash the huffman library, but still gives the correct (or
    at least a correct enough) result.
    '''
    # the library can't handle an input that contains only one distinct byte
    # e.g. [2, 2, 2, 2, 2, 2]
    l = len(bytes_)
    for i, b in enumerate(bytes_):
        if b != bytes_[0]:
            # at least two different bytes, break
            break
        if i == l - 1:
            # reached the end and all bytes were the same: library will fail, so
            # return a good estimate
            # huffman should code each byte with one bit, plus a little overhead
            # for the tree
            return int(math.ceil(l / 8)) + 5

    huffman_code = cd.huffman_enc(bytes_)
    return len(huffman_code)

def num_inversions(list_a, list_b):
    '''Compute the number of inversions between two lists.
    The lists need to be permutations of each other.
    '''
    inv = 0
    for i_a, x in enumerate(list_a):
        i_b = list_b.index(x)
        # all the items coming before x in list_b
        before_list = list_b[:i_b]
        # all the items coming after x in hst_a
        after_list = list_a[i_a + 1:]
        for a in after_list:
            if a in before_list:
                # inverions found, increment inv
                inv += 1
    return inv

def num_inversions_general(list_a, list_b):
    '''Compute the number of inversions between two lists.
    Can deal with lists of different lengths and multiple occurences.
    '''
    inv = 0
    for i, x in enumerate(list_a):
        # TODO need to recheck that i'm not counting anything multiple times
        try:
            i_b = list_b.index(x)
        except ValueError:
            # no inversions if value isn't in list_b, continue with next
            continue
        # all the items coming before x in list_b
        before_list = list_b[:i_b]
        # all the items coming after x in lst_a, but not equal to x
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

def analyze_transitions(bytes_, metric):
    '''Analyze all the transitions between bytes of a byte string.'''
    if metric == 'num_chars':
        an_func = metric_num_chars
    elif metric == 'max_code':
        an_func = metric_max_code
    elif metric == 'median':
        an_func = metric_median
    elif metric == 'mean':
        an_func = metric_mean
    elif metric == 'chapin_hst_diff':
        an_func = metric_chapin_hst_diff
    elif metric == 'chapin_kl':
        an_func = metric_chapin_kl
    elif metric == 'chapin_inv':
        an_func = metric_chapin_inv
    elif metric == 'chapin_inv_log':
        an_func = metric_chapin_inv_log
    elif metric == 'huffman_metric':
        an_func = metric_huffman
    else:
        raise ValueError('Unknown metric.')

    bw_code = cd.bw_encode(bytes_)
    firsts = set(bw_code.firsts)
    transitions = {(a, b): an_func(bw_code, a, b)
                   for a in firsts
                   for b in firsts if a != b}
    return transitions
