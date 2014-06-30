'''
Created on 28.05.2014

@author: dddsnn
'''

import sys
import queue
from bwt import *
import bwt.coder as cd
import math
import warnings

def select_sequences(bytes_, max_len):
    '''Select sequences from a file that could benefit from being ordered more
    specially than by the same order for all characters.
    '''
    # TODO this is a primitive placeholder and needs to be replaced by a
    # function that's aware of the quality of the transitions of each
    # possible sequence when it's not specially reordered

    # make histograms for sequences of length up to max_len
    hst = {}
    for l in range(2, max_len + 1):
        for i in range(len(bytes_) - l):
            if bytes_[i:i + l] in hst:
                hst[bytes_[i:i + l]] += 1
            else:
                hst[bytes_[i:i + l]] = 1
    # make a list of all sequences found and sort by number of occurrences
    seq = [s for s in hst]
    seq.sort(key=lambda x:hst[x], reverse=True)
    # select the first 255 and return
    return []  # TODO
    return seq[:255]

def bw_block(bw_code, seq, specializations):
    '''Get the block of BW code corresponding to a specific sequence at the
    beginning of the line. Exclude symbols in the BW code belonging to a more
    special sequence indicated in the specializations dictionary.
    '''
    # first index matching seq
    left = next(i for i, s in enumerate(bw_code.firsts) if s[:len(seq)] == seq)
    if not specializations[seq]:
        # no specializaions for this sequence, straightforward
        right = left + 1
        while right < len(bw_code.firsts) and bw_code.firsts[right][:len(seq)] == seq:
            right += 1
        bw_block = bw_code.encoded[left:right]
        return bw_block

    # reverse specializations so they're sorted from most general to most
    # specific
    spec = reversed(specializations[seq])
    # find the first occurrence of seq in firsts that isn't a specialization
    while left < len(bw_code.firsts) and bw_code.firsts[left][:len(seq)] == seq:
        for sp in spec:
            if bw_code.firsts[left][:len(sp)] == sp:
                # this line is a specialization, so increment left and break to
                # find another beginning
                left += 1
                break
        else:
            # no break, so this line is not a specialization, break out of the
            # while and keep the value of left
            break
    if not left < len(bw_code.firsts) or bw_code.firsts[left][:len(seq)] != seq:
        # this means that all elements in firsts starting with seq were
        # specializations, so return an empty bytes object as the bw block
        return b''

    # start the bw block
    bw_block = bytes([bw_code.encoded[left]])

    # now find the end of the bw block
    right = left + 1
    while right < len(bw_code.firsts) and bw_code.firsts[right][:len(seq)] == seq:
        for sp in spec:
            if bw_code.firsts[right][:len(sp)] == sp:
                # this line is a specialization, so increment right and break
                right += 1
                break
        else:
            # no break, so this line is not a specialization, append to the
            # result and increment right
            bw_block += bytes([bw_code.encoded[right]])
            right += 1
    return bw_block

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

def huffman_codeword_lengths(mtf_code):
    '''For each byte value, give an estimate of the length in bit of the huffman
    code word if the bytes_ were encoded with huffman while also including
    values that usually wouldn't be encoded because they don't occur. Use for
    MTF codes.
    '''
    all_counts = {n:mtf_code.count(n) for n in range(256)}
    # remove trailing zeroes (not occuring values) and set intermediate zeroes
    # to 0.2 (to make sure they get longer codes)
    last = max([n for n in range(256) if all_counts[n] != 0])
    counts = {k:v for k, v in all_counts.items() if k <= last}
    # this is ugly
    Node = namedtuple('Node', ['weight', 'inner', 'left', 'right', 'value'])
    q = queue.PriorityQueue()
    for n in range(256):
        # replace zeroes coming after the last non-zero with 1/1000000 so they
        # get the longest codes
        if not n in counts:
            counts[n] = 0.000001
        # replace other zeroes with 1/5 so they usually get longer codes than
        # non-zeroes but shorter ones than the trailing zeroes
        if counts[n] == 0:
            counts[n] = 0.2
        q.put(Node(counts[n], False, (), (), n))
    while q.qsize() > 1:
        # get the two lowest items from the queue, make them into a node and put
        # the node back into the queue with appropriate weight
        left = q.get()
        right = q.get()
        weight = left.weight + right.weight
        q.put(Node(weight, True, left, right, -1))
    # now that the queue has only one element, it must be the root of the tree
    NodeDepthTuple = namedtuple('NodeDepthTuple', ['node', 'depth'])
    root = q.get()
    root_tuple = NodeDepthTuple(root, 0)
    tuples = [root_tuple]
    lengths = {}
    # go through the tree from the root and write the lengths for each value
    # to the dict
    while tuples:
        t = tuples.pop()
        if not t.node.inner:
            # leaf node: write length into dict and continue
            lengths[t.node.value] = t.depth
            if t.depth == 0:
                # just to make sure there isn't just one value that gets length
                # zero (that's impossible)
                lengths[t.node.value] = 1
        else:
            # inner node: put children into tuples and increment depth
            tuples.append(NodeDepthTuple(t.node.left, t.depth + 1))
            tuples.append(NodeDepthTuple(t.node.right, t.depth + 1))
    return lengths

def analyze_transitions(bytes_, metric, aux_data):
    '''Analyze all the transitions between bytes of a byte string.'''
    if aux_data.raw != bytes_:
        # not the correct data
        raise ValueError('Not the correct aux data.')

    an_func = getattr(sys.modules[__name__], metric)
    bw_code = aux_data.bw_code
    firsts = aux_data.firsts
    transitions = {(a, b): an_func(bw_code, a, b, aux_data)
                   for a in [bytes([x]) for x in firsts]
                   for b in [bytes([x]) for x in firsts] if a != b}
    return transitions

def metric_num_chars(bw_code, first_seq_a, first_seq_b, aux_data):
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_ab = aux_data.partial_mtf_analyses[(first_seq_a, first_seq_b)]

    # metric: number of -1s in the second half of the combined code (the one
    # that's being transitioned to)
    metric = an_ab.num_chars - an_a.num_chars
    return metric

def metric_max_code(bw_code, first_seq_a, first_seq_b, aux_data):
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    an_ab = aux_data.partial_mtf_analyses[(first_seq_a, first_seq_b)]

    metric = an_ab.max_code - max(an_a.max_code, an_b.max_code)
    return metric

def metric_median(bw_code, first_seq_a, first_seq_b, aux_data):
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    an_ab = aux_data.partial_mtf_analyses[(first_seq_a, first_seq_b)]

    # to avoid division by zero
    if an_ab.length_rec == 0:
        ab_length_rec = 1
    else:
        ab_length_rec = an_ab.length_rec

    # metric: difference of expected median and the achieved median
    metric = ((an_a.median * an_a.length_rec + an_b.median * an_b.length_rec)
                       / ab_length_rec) - an_ab.median
    return metric

def metric_mean(bw_code, first_seq_a, first_seq_b, aux_data):
    # TODO by ignoring new characters (-1s) in the mean, good transitions get an_a
    # penalty, because a character that was already there before the transition
    # probably affects the mean in a bad way, while an entirely new character
    # won't affect it at all
    # need to give new characters an_a "penalty" value in the partial mtf encode
    # (but only for the target of the transition, lots of new characters in the
    # source don't mean anything bad)
    # or maybe an entirely new metric that takes this into account

    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    an_ab = aux_data.partial_mtf_analyses[(first_seq_a, first_seq_b)]

    # to avoid division by zero
    if an_ab.length_rec == 0:
        ab_length_rec = 1
    else:
        ab_length_rec = an_ab.length_rec

    # metric: difference of expected mean and the achieved mean
    metric = ((an_a.mean * an_a.length_rec + an_b.mean * an_b.length_rec)
                     / ab_length_rec) - an_ab.mean
    return metric

def metric_mean_new_penalty(bw_code, first_seq_a, first_seq_b, aux_data):
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    an_ab = aux_data.partial_mtf_analyses[(first_seq_a, first_seq_b)]

    # bytes that will be counted are all non-(-1)s and the -1s from the second
    # half
    length = an_a.length_rec + an_b.length

    new_in_right = an_ab.num_chars - an_a.num_chars
    # penalty for every new value: minimum possible value the worst of themwill
    # actually be encoded with
    penalty = an_a.max_code + an_a.num_chars

    # add a penalty for -1s in the right side
    metric = ((an_a.mean * an_a.length_rec + an_b.mean * an_b.length_rec
               + new_in_right * penalty) / length)
    return metric

def metric_mean_right(bw_code, first_seq_a, first_seq_b, aux_data):
    an_b = aux_data.partial_mtf_analyses[first_seq_b]

    metric = an_b.mean
    return metric

def metric_chapin_hst_diff(bw_code, first_seq_a, first_seq_b, aux_data):
    hst_a = aux_data.bw_subhistograms[first_seq_a]
    hst_b = aux_data.bw_subhistograms[first_seq_b]

    # CHAPIN: sum of squares of differences of logs
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

def metric_chapin_kl(bw_code, first_seq_a, first_seq_b, aux_data):
    hst_a = aux_data.bw_subhistograms[first_seq_a]
    hst_b = aux_data.bw_subhistograms[first_seq_b]

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

def metric_chapin_inv(bw_code, first_seq_a, first_seq_b, aux_data):
    freq_list_a = aux_data.freq_lists[first_seq_a]
    freq_list_b = aux_data.freq_lists[first_seq_b]

    # CHAPIN: number of inversion between ordered histograms
    # metric is the number of inversions between the two lists
    metric = num_inversions(freq_list_a, freq_list_b)
    return metric

def metric_chapin_inv_log(bw_code, first_seq_a, first_seq_b, aux_data):
    inv_metric = metric_chapin_inv(bw_code, first_seq_a, first_seq_b,
                                   aux_data)
    # CHAPIN: log of previous
    if inv_metric == 0:
        metric = 0
    else:
        metric = math.log(inv_metric)
    return metric

def metric_huffman(bw_code, first_seq_a, first_seq_b, aux_data):
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]

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

def metric_huffman_new_penalty(bw_code, first_seq_a, first_seq_b,
                               aux_data):
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    bw_a = aux_data.bw_subcodes[first_seq_a]

    # replace -1s in the right part of the partial mtf with values that haven't
    # occurred before
    length_a = len(bw_a)
    stripped_pt_mtf_a = list(filter(lambda x: x != -1, pt_mtf_ab[:length_a]))
    pt_mtf_b = pt_mtf_ab[length_a:]
    max_code = max(pt_mtf_ab) + 1
    stripped_pt_mtf_b = []
    for b in pt_mtf_b:
        if b == -1:
            stripped_pt_mtf_b.append(max_code)
            max_code += 1
        else:
            stripped_pt_mtf_b.append(b)
    stripped_pt_mtf = stripped_pt_mtf_a + stripped_pt_mtf_b
    huffman_length = huffman_enc_length(bytes(stripped_pt_mtf))
    if len(stripped_pt_mtf) == 0:
        # if the stripped partial mtf code has length 0 it means there are no
        # recurring characters in it. since that means no compression benefits
        # for this transition, give it a bad value to penalize it
        # good transitions will have lengths far below 1, so 10 should be enough
        # to discourage tsp from including this transition
        metric = 10
    else:
        metric = huffman_length / len(stripped_pt_mtf)
    return metric

def metric_badness(bw_code, first_seq_a, first_seq_b, aux_data):
    '''Compares the ordering of the MTF alphabet the left side of the transition
    "leaves" for the right side with the ideal ordering that would yield the
    greatest compression benefit and gives a badness value for non-ideal
    orderings.
    '''
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]

    badness = 0
    # make list of tuples (index in partial mtf of right side, optimal code) for
    # every new symbol in the partial mtf of the right side
    ideal_codes = []
    min_possible = 0
    for i, s in enumerate(pt_mtf_b):
        if s != -1:
            # recurring symbol, ignore
            continue
        ideal_codes.append((i, min_possible))
        min_possible += 1
    # length of the mtf code of the left side
    length_left = len(pt_mtf_ab) - len(pt_mtf_b)
    # the minimal possible code for any new symbols in the combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            badness += min_possible - ideal
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            badness += actual - ideal
    return badness

def metric_badness_weighted(bw_code, first_seq_a, first_seq_b, aux_data):
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]

    badness = 0
    # make list of tuples (index in partial mtf of right side, optimal code) for
    # every new symbol in the partial mtf of the right side
    ideal_codes = []
    min_possible = 0
    for i, s in enumerate(pt_mtf_b):
        if s != -1:
            # recurring symbol, ignore
            continue
        ideal_codes.append((i, min_possible))
        min_possible += 1
    # length of the mtf code of the left side
    length_left = len(pt_mtf_ab) - len(pt_mtf_b)
    # the minimal possible code for any new symbols in the combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            badness += min_possible - ideal
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            badness += actual - ideal
    return badness / an_b.num_chars

def metric_badness_mean_penalty(bw_code, first_seq_a, first_seq_b,
                                aux_data):
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    mtf_means = aux_data.mtf_mean_steps

    badness = 0
    # make list of tuples (index in partial mtf of right side, optimal code) for
    # every new symbol in the partial mtf of the right side
    ideal_codes = []
    min_possible = 0
    for i, s in enumerate(pt_mtf_b):
        if s != -1:
            # recurring symbol, ignore
            continue
        ideal_codes.append((i, min_possible))
        min_possible += 1
    # length of the mtf code of the left side
    length_left = len(pt_mtf_ab) - len(pt_mtf_b)
    # the minimal possible code for any new symbols in the combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            badness += mtf_means[min_possible] - ideal
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            badness += actual - ideal
    return badness

def metric_badness_weighted_mean_penalty(bw_code, first_seq_a,
                                         first_seq_b, aux_data):
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    mtf_means = aux_data.mtf_mean_steps

    badness = 0
    # make list of tuples (index in partial mtf of right side, optimal code) for
    # every new symbol in the partial mtf of the right side
    ideal_codes = []
    min_possible = 0
    for i, s in enumerate(pt_mtf_b):
        if s != -1:
            # recurring symbol, ignore
            continue
        ideal_codes.append((i, min_possible))
        min_possible += 1
    # length of the mtf code of the left side
    length_left = len(pt_mtf_ab) - len(pt_mtf_b)
    # the minimal possible code for any new symbols in the combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            badness += mtf_means[min_possible] - ideal
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            badness += actual - ideal
    return badness / an_b.num_chars

def metric_badness_huff_len(bw_code, first_seq_a, first_seq_b, aux_data):
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    hf_len = aux_data.huffman_codeword_lengths

    badness = 0
    # make list of tuples (index in partial mtf of right side, optimal code) for
    # every new symbol in the partial mtf of the right side
    ideal_codes = []
    min_possible = 0
    for i, s in enumerate(pt_mtf_b):
        if s != -1:
            # recurring symbol, ignore
            continue
        ideal_codes.append((i, min_possible))
        min_possible += 1
    # length of the mtf code of the left side
    length_left = len(pt_mtf_ab) - len(pt_mtf_b)
    # the minimal possible code for any new symbols in the combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            huff_dist = hf_len[min_possible] - hf_len[ideal]
            if huff_dist < 0:
                # for high code values, higher codes can accidentally be smaller
                # than lower ones. set to 0 if this happens
                huff_dist = 0
            badness += huff_dist
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            huff_dist = hf_len[actual] - hf_len[ideal]
            if huff_dist < 0:
                # for high code values, higher codes can accidentally be smaller
                # than lower ones. set to 0 if this happens
                huff_dist = 0
            badness += huff_dist
    return badness

def metric_badness_huff_len_weighted(bw_code, first_seq_a, first_seq_b,
                                     aux_data):
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    hf_len = aux_data.huffman_codeword_lengths

    badness = 0
    # make list of tuples (index in partial mtf of right side, optimal code) for
    # every new symbol in the partial mtf of the right side
    ideal_codes = []
    min_possible = 0
    for i, s in enumerate(pt_mtf_b):
        if s != -1:
            # recurring symbol, ignore
            continue
        ideal_codes.append((i, min_possible))
        min_possible += 1
    # length of the mtf code of the left side
    length_left = len(pt_mtf_ab) - len(pt_mtf_b)
    # the minimal possible code for any new symbols in the combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            huff_dist = hf_len[min_possible] - hf_len[ideal]
            if huff_dist < 0:
                # for high code values, higher codes can accidentally be smaller
                # than lower ones. set to 0 if this happens
                huff_dist = 0
            badness += huff_dist
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            huff_dist = hf_len[actual] - hf_len[ideal]
            if huff_dist < 0:
                # for high code values, higher codes can accidentally be smaller
                # than lower ones. set to 0 if this happens
                huff_dist = 0
            badness += huff_dist
    return badness / an_b.num_chars

def metric_badness_huff_len_mean_penalty(bw_code, first_seq_a,
                                         first_seq_b, aux_data):
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    hf_len = aux_data.huffman_codeword_lengths
    mtf_means = aux_data.mtf_mean_steps

    badness = 0
    # make list of tuples (index in partial mtf of right side, optimal code) for
    # every new symbol in the partial mtf of the right side
    ideal_codes = []
    min_possible = 0
    for i, s in enumerate(pt_mtf_b):
        if s != -1:
            # recurring symbol, ignore
            continue
        ideal_codes.append((i, min_possible))
        min_possible += 1
    # length of the mtf code of the left side
    length_left = len(pt_mtf_ab) - len(pt_mtf_b)
    # the minimal possible code for any new symbols in the combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            huff_dist = hf_len[int(round(mtf_means[min_possible]))] - hf_len[ideal]
            if huff_dist < 0:
                # for high code values, higher codes can accidentally be smaller
                # than lower ones. set to 0 if this happens
                huff_dist = 0
            badness += huff_dist
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            huff_dist = hf_len[actual] - hf_len[ideal]
            if huff_dist < 0:
                # for high code values, higher codes can accidentally be smaller
                # than lower ones. set to 0 if this happens
                huff_dist = 0
            badness += huff_dist
    return badness

def metric_badness_huff_len_weighted_mean_penalty(bw_code, first_seq_a,
                                                  first_seq_b, aux_data):
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    hf_len = aux_data.huffman_codeword_lengths
    mtf_means = aux_data.mtf_mean_steps

    badness = 0
    # make list of tuples (index in partial mtf of right side, optimal code) for
    # every new symbol in the partial mtf of the right side
    ideal_codes = []
    min_possible = 0
    for i, s in enumerate(pt_mtf_b):
        if s != -1:
            # recurring symbol, ignore
            continue
        ideal_codes.append((i, min_possible))
        min_possible += 1
    # length of the mtf code of the left side
    length_left = len(pt_mtf_ab) - len(pt_mtf_b)
    # the minimal possible code for any new symbols in the combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            huff_dist = hf_len[int(round(mtf_means[min_possible]))] - hf_len[ideal]
            if huff_dist < 0:
                # for high code values, higher codes can accidentally be smaller
                # than lower ones. set to 0 if this happens
                huff_dist = 0
            badness += huff_dist
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            huff_dist = hf_len[actual] - hf_len[ideal]
            if huff_dist < 0:
                # for high code values, higher codes can accidentally be smaller
                # than lower ones. set to 0 if this happens
                huff_dist = 0
            badness += huff_dist
    return badness / an_b.num_chars
