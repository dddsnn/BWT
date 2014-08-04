import sys
import queue
from bwt import *
import math

def select_sequences(bs, max_len):
    """Select sequences from a file that could benefit from being ordered more
    specially than by the same order for all characters.
    """
    # TODO this is a primitive placeholder and needs to be replaced by a
    # function that's aware of the quality of the transitions of each
    # possible sequence when it's not specially reordered

    # make histograms for sequences of length up to max_len
    hst = {}
    for l in range(2, max_len + 1):
        for i in range(len(bs) - l):
            if bs[i:i + l] in hst:
                hst[bs[i:i + l]] += 1
            else:
                hst[bs[i:i + l]] = 1
    # make a list of all sequences found and sort by number of occurrences
    seq = [s for s in hst]
    seq.sort(key=lambda x:hst[x], reverse=True)
    # select the first 255 and return
    return []
    return seq[:20]  # TODO

def bw_block(bw_code, seq, specializations):
    """Get a specific BW block.

    Gets the block of BW code corresponding to a specific sequence at the
    beginning of the line. Omit symbols in the BW code belonging to a more
    special sequence indicated in the specializations dictionary.

    Args:
        bw_code: The BWEncodeResult to use as a source.
        seq: The sequence at the beginning of the line whose block should be
            returned.
        specializations: A dictionary listing specializations of seq as returned
            by bwt.coder.specializations().

    Returns:
        The block of BW code whose first lines begin with seq, excluding those
        beginning with any of specializations[seq].
    """
    # first index matching seq
    left = next(i for i, s in enumerate(bw_code.firsts) if s[:len(seq)] == seq)
    if not specializations[seq]:
        # no specializations for this sequence, straightforward
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

def analyze_partial_mtf(mtf_code):
    """Get various stats about a partial MTF mtf_code.

    Args:
        mtf_code: The partial MTF to be analyzed.

    Returns:
        A bwt.PartialMTFAnalysisResult of the mtf_code.
    """
    raw = mtf_code
    length = len(mtf_code)
    sorted_code = sorted(mtf_code)
    max_code = sorted_code[-1]
    # remove -1s
    sorted_code = [c for c in sorted_code if c != -1]
    # number of all recurring symbols (those that aren't -1)
    l = len(sorted_code)
    # number of -1s is number of all symbols minus number of all except -1s
    num_chars = length - l
    # number of recurring characters (not -1)
    length_rec = length - num_chars
    result = PartialMTFAnalysisResult(raw, length, length_rec, num_chars,
                                      max_code)
    return result

def num_inversions(list_a, list_b):
    """Compute the number of inversions between two lists.

    The number of inversions is the number of times any pair of elements occur
    in a different order in the two lists.
    The lists need to be permutations of each other.

    Args:
        list_a: The first list.
        list_b: The second list.

    Returns:
        The number of inverions between list_a and list_b.
    """
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

def make_histogram(bs):
    """Make a histogram of symbol appearances for a bytes object.

    Args:
        bs: The input bytes object.

    Returns:
        A histogram giving count of appearance in bs for every byte value.
        E.g.: bytes([1, 1, 2, 3]) -> {0: 0, 1: 2, 3: 1, 4: 0, 5: 0 ...}
    """
    # initialize 0 for all possible bytes
    histogram = {i:0 for i in range(256)}
    for b in bs:
        histogram[b] += 1
    return histogram

def huffman_codeword_lengths(mtf_code):
    """Estimate huffman code word lengths for data similar to the input.

    For each byte value, give an estimate of the length in bit of the huffman
    code word if the bytes_ were encoded with huffman while also including
    values that usually wouldn't be encoded because they don't occur. Use for
    MTF codes.
    """
    all_counts = {n:mtf_code.count(n) for n in range(256)}
    # remove trailing zeroes (not occuring values) and set intermediate zeroes
    # to 0.2 (to make sure they get longer codes)
    last = max([n for n in range(256) if all_counts[n] != 0])
    counts = {k:v for k, v in all_counts.items() if k <= last}
    # TODO this is ugly
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

def analyze_transitions(bs, aux_data, metric, **metric_opts):
    """Analyze all possible transitions between values of a byte string.

    Args:
        bs: The input bytes object.
        aux_data: The bwt.AuxData object corresponding to bs.
        metric: The name of the metric to be used. The string 'metric_' will
            be prepended to this name to get the actual function.
        metric_opts: A dictionary of options to be passed to the metric function
            as keyword arguments.

    Returns:
        A dictionary mapping each pair of byte values from bs to that
        transition's weight according to the given metric.
    """
    if aux_data.raw != bs:
        # not the correct data
        raise ValueError('Not the correct aux data.')

    an_func = getattr(sys.modules[__name__], 'metric_' + metric)
    bw_code = aux_data.bw_code
    firsts = aux_data.firsts
    transitions = {(a, b): an_func(bw_code, a, b, aux_data, **metric_opts)
                   for a in firsts
                   for b in firsts if a != b}
    return transitions

def metric_chapin_hst_diff(bw_code, first_seq_a, first_seq_b, aux_data,
                           **kwargs):
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

def metric_chapin_kl(bw_code, first_seq_a, first_seq_b, aux_data, **kwargs):
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

def metric_chapin_inv(bw_code, first_seq_a, first_seq_b, aux_data, **kwargs):
    freq_list_a = aux_data.freq_lists[first_seq_a]
    freq_list_b = aux_data.freq_lists[first_seq_b]

    # CHAPIN: number of inversion between ordered histograms
    # metric is the number of inversions between the two lists
    metric = num_inversions(freq_list_a, freq_list_b)
    return metric

def metric_chapin_inv_log(bw_code, first_seq_a, first_seq_b, aux_data, **kwargs):
    inv_metric = metric_chapin_inv(bw_code, first_seq_a, first_seq_b,
                                   aux_data)
    # CHAPIN: log of previous
    if inv_metric == 0:
        metric = 0
    else:
        metric = math.log(inv_metric)
    return metric

def metric_badness(bw_code, first_seq_a, first_seq_b, aux_data, **kwargs):
    """Compares the ordering of the MTF alphabet the left side of the transition
    "leaves" for the right side with the ideal ordering that would yield the
    greatest compression benefit and gives a badness value for non-ideal
    orderings.

    entropy_code_len assumes lower mtfs get shorter entropy codes (at least for
    low values)
    """
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    mtf_means = aux_data.mtf_mean_steps
    hf_len = aux_data.huffman_codeword_lengths

    # get options from kwargs
    if 'weighted' in kwargs:
        weighted = kwargs['weighted']
    else:
        weighted = False
    if 'new_penalty' in kwargs:
        new_penalty = kwargs['new_penalty']
    else:
        new_penalty = False
    if 'entropy_code_len' in kwargs:
        entropy_code_len = kwargs['entropy_code_len']
    else:
        entropy_code_len = False

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
    # the minimal possible code for any new symbols in the right side of the
    # combined code
    min_possible = an_a.num_chars
    # now compare the ideal codes with the actual ones and add the difference
    # to the badness
    for i, ideal in ideal_codes:
        actual = pt_mtf_ab[i + length_left]
        if actual == -1:
            # new symbol in the combined mtf
            if new_penalty:
                # give a penalty for a new symbol, if this was requested in the
                # options
                assumed_actual = mtf_means[min_possible]
            else:
                # otherwise, assume the best possible
                assumed_actual = min_possible
            if entropy_code_len:
                # add the approximation of the actual number of bits this
                # transition will cost in the entropy coder, if requested in
                # the options
                huff_dist = hf_len[assumed_actual] - hf_len[ideal]
                if huff_dist < 0:
                    # for high code values, higher codes can accidentally be
                    # smaller than lower ones. set to 0 if this happens
                    huff_dist = 0
                badness += huff_dist
            else:
                # otherwise just add the distance
                badness += assumed_actual - ideal
            # minimal possible code for new symbols is now increased
            min_possible += 1
        else:
            if entropy_code_len:
                huff_dist = hf_len[actual] - hf_len[ideal]
                if huff_dist < 0:
                    huff_dist = 0
                badness += huff_dist
            else:
                badness += actual - ideal
    if weighted:
        # weight the result by the number of symbols in the right side, if this
        # was requested in the options
        badness /= an_b.num_chars
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
    # the minimal possible code for any new symbols in the right side of the
    # combined code
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
