import bwt.huffman as hf
import sys
from bwt import *
import math
import numpy as np

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

def bw_block(code, firsts, seq, specializations):
    """Get a specific BW block.

    Gets the block of BW code corresponding to a specific sequence at the
    beginning of the line. Omit symbols in the BW code belonging to a more
    special sequence indicated in the specializations dictionary.

    Args:
        code: An iterable from which the result will be selected.
        firsts: The sequences at the beginning of the BW table. Must have the
            same length as code and correspond to it. At least one element must
            contain seq as a substring.
        seq: The sequence at the beginning of the line whose block should be
            returned.
        specializations: A dictionary listing specializations of seq as returned
            by bwt.coder.specializations().

    Returns:
        The block of BW code whose first lines begin with seq, excluding those
        beginning with any of specializations[seq].
    """
    # first index matching seq
    left = next(i for i, s in enumerate(firsts) if s[:len(seq)] == seq)
    if not specializations[seq]:
        # no specializations for this sequence, straightforward
        right = left + 1
        while right < len(firsts) and firsts[right][:len(seq)] == seq:
            right += 1
        bw_block = code[left:right]
        return bw_block

    # reverse specializations so they're sorted from most general to most
    # specific
    spec = reversed(specializations[seq])
    # find the first occurrence of seq in firsts that isn't a specialization
    while left < len(firsts) and firsts[left][:len(seq)] == seq:
        for sp in spec:
            if firsts[left][:len(sp)] == sp:
                # this line is a specialization, so increment left and break to
                # find another beginning
                left += 1
                break
        else:
            # no break, so this line is not a specialization, break out of the
            # while and keep the value of left
            break
    if not left < len(firsts) or firsts[left][:len(seq)] != seq:
        # this means that all elements in firsts starting with seq were
        # specializations, so return an empty bytes object as the bw block
        return b''

    # start the bw block
    bw_block = bytes([code[left]])

    # now find the end of the bw block
    right = left + 1
    while right < len(firsts) and firsts[right][:len(seq)] == seq:
        for sp in spec:
            if firsts[right][:len(sp)] == sp:
                # this line is a specialization, so increment right and break
                right += 1
                break
        else:
            # no break, so this line is not a specialization, append to the
            # result and increment right
            bw_block += bytes([code[right]])
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

def huffman_codeword_lengths(mtf_code, zero_compensation):
    """Estimate huffman code word lengths for data similar to the input.

    For each byte value, give an estimate of the length in bit of the huffman
    code word if the bytes_ were encoded with huffman while also including
    values that usually wouldn't be encoded because they don't occur. Use for
    MTF codes.
    """
    freqs = hf.symbol_frequencies(mtf_code)
    if zero_compensation == 'complete':
        # find the highest value symbol that appears in the mtf code
        max_nonzero = max(freqs)
        # give all mtf codes that don't appear in mtf_code and are lower than
        # max_nonzero weight 1/256, so they get a longer code than any of the
        # ones actually appearing
        for i in range(max_nonzero):
            if i not in freqs:
                freqs[i] = 1 / 256
        # give all mtf codes that don't appear and are higher than max_nonzero
        # weight 1/256**2 so they get even longer codes
        for i in range(max_nonzero + 1, 256):
            if i not in freqs:
                freqs[i] = 1 / (256 ** 2)
        return hf.codeword_lengths(freqs)
    elif zero_compensation == 'sparse':
        hf_len = hf.codeword_lengths(freqs)
        # go through the dict and give each non-occurring symbol the same
        # codeword length as its predecessor
        last_len = hf_len[min(hf_len)]
        for i in range(256):
            if i not in hf_len:
                hf_len[i] = last_len
            else:
                last_len = hf_len[i]
        return hf_len
    else:
        raise ValueError('{0} is not a valid value for zero_compensation. Must'
                         ' be one of False, \'complete\' or \'sparse\'.'
                         .format(zero_compensation))

def mtf_mean_steps(bw_code, mtf_code):
    """Make a dictionary of average MTF codes of a minimum value.

    Args:
        bw_code: BW code as a bytes object.
        mtf_code: The MTF code of bw_code, a list of integers between 0 and 255.

    Returns:
        A dictionary containing, for every integer n between 0 and 255, the
        average of all values in mtf_code that are greater or equal than n.
        The dictionary also maps tuples (sym, n) to the average of all values
        in mtf_code that encode the symbol sym from bw_code and are greater or
        equal than n.
    """
    def mean_steps_generic(mtf_code):
        result = {}
        for n in range(256):
            l = [x for x in mtf_code if x >= n]
            if l:
                result[n] = np.mean(l)
            else:
                result[n] = n
        return result
    # first the generic means, not taking into account the underlying symbol
    result = mean_steps_generic(mtf_code)
    # now for every underlying symbol, for every mtf code
    symbol_mtf_dict = {s:[] for s in set(bw_code)}
    for sym, mtf in zip(bw_code, mtf_code):
        symbol_mtf_dict[sym].append(mtf)
    # make the generic means dict for every symbol's list of mtf codes
    for sym in symbol_mtf_dict:
        for n, mean_dev in mean_steps_generic(symbol_mtf_dict[sym]).items():
            result[(sym, n)] = mean_dev
    return result

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
    firsts = aux_data.firsts
    transitions = {(a, b): an_func(a, b, aux_data, **metric_opts)
                   for a in firsts
                   for b in firsts if a != b}
    return transitions

def metric_chapin_hst_diff(first_seq_a, first_seq_b, aux_data, **kwargs):
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

def metric_chapin_kl(first_seq_a, first_seq_b, aux_data, **kwargs):
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

def metric_chapin_inv(first_seq_a, first_seq_b, aux_data, **kwargs):
    freq_list_a = aux_data.freq_lists[first_seq_a]
    freq_list_b = aux_data.freq_lists[first_seq_b]

    # get options from kwargs
    if 'log' in kwargs:
        log = kwargs['log']
    else:
        log = False

    # CHAPIN: number of inversion between ordered histograms
    # metric is the number of inversions between the two lists
    metric = num_inversions(freq_list_a, freq_list_b)
    if log:
        # calculate log, if requested in the options
        if metric == 0:
            return metric
        else:
            return math.log(metric)
    else:
        return metric

def metric_badness(first_seq_a, first_seq_b, aux_data, **kwargs):
    """Assign a "badness" value to a given transition.

    Compares the ordering of the MTF alphabet the left side of the transition
    "leaves" for the right side with the ideal ordering that would yield the
    greatest compression benefit, and gives a value for non-ideal orderings
    denoting how "bad" the transition is.

    Args:
        first_seq_a: The sequence of bytes whose context block is the left side
            of the transition, as a bytes object.
        first_seq_b: Same as first_seq_a, for the right side of the transition.
        aux_data: The bwt.AuxData for the input file.
        kwargs: Keyword Arguments:
            weighted: Boolean, whether the badness should be weighted with the
                number of new symbols in the right side. Defaults to False.
            new_penalty: Boolean, whether a penalty should be given for new
                symbols in the right side that don't appear in the left side.
                Defaults to False.
            entropy_code_len: Boolean, whether the difference in the number of
                bits used by the entropy coder should be used to compute the
                badness, rather than the plain difference between the codes.
                Defaults to False.
            new_penalty_log: A dictionary to which the predictions of mtf codes
                that aren't know are written. If it doesn't exist, a new key
                (first_seq_a, first_seq_b) is created and a dictionary is
                stored, mapping the underlying symbol of the BW code for which
                the MTF code is predicted to the prediction.

    Returns:
        A number denoting how "bad" the given transition is.
    """
    pt_mtf_b = aux_data.partial_mtf_subcodes[first_seq_b]
    pt_mtf_ab = aux_data.partial_mtf_subcodes[(first_seq_a, first_seq_b)]
    an_a = aux_data.partial_mtf_analyses[first_seq_a]
    an_b = aux_data.partial_mtf_analyses[first_seq_b]
    mtf_means = aux_data.mtf_mean_steps

    # get options from kwargs
    if 'weighted' in kwargs:
        weighted = kwargs['weighted']
    else:
        weighted = False
    if 'new_penalty' in kwargs:
        new_penalty = kwargs['new_penalty']
        if new_penalty not in [False, 'generic', 'specific']:
            raise ValueError('{0} is not a valid value for new_penalty. Must be'
                             ' one of False, \'generic\' or \'specific\'.'
                             .format(new_penalty))
        if new_penalty == 'specific':
            bw_subcode_b = aux_data.bw_subcodes[first_seq_b]
    else:
        new_penalty = False
    if 'entropy_code_len' in kwargs:
        entropy_code_len = kwargs['entropy_code_len']
        if entropy_code_len not in [False, 'complete', 'sparse']:
            raise ValueError('{0} is not a valid value for entropy_code_len. '
                             'Must be one of False, \'complete\' or \'sparse\'.'
                             .format(kwargs['entropy_code_len']))
        if kwargs['entropy_code_len'] == 'complete':
            hf_len = aux_data.huffman_codeword_lengths_complete
        elif kwargs['entropy_code_len'] == 'sparse':
            hf_len = aux_data.huffman_codeword_lengths_sparse
    else:
        entropy_code_len = False
    if 'new_penalty_log' in kwargs:
        new_penalty_log = kwargs['new_penalty_log']
        bw_subcode_b = aux_data.bw_subcodes[first_seq_b]
    else:
        new_penalty_log = False

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
            if new_penalty == 'generic':
                # give a generic penalty for a new symbol, if this was requested
                # in the options
                assumed_actual = mtf_means[min_possible]
            elif new_penalty == 'specific':
                # give a penalty specific to the underlying symbol from the
                # bw code
                assumed_actual = mtf_means[(bw_subcode_b[i], min_possible)]
            else:
                # otherwise, assume the best possible
                assumed_actual = min_possible
            # write to the new penalty log if it exists
            if new_penalty_log != False:
                if (first_seq_a, first_seq_b) not in new_penalty_log:
                    new_penalty_log[(first_seq_a, first_seq_b)] = {}
                new_penalty_log[(first_seq_a,
                                 first_seq_b)][bw_subcode_b[i]] = assumed_actual
            if entropy_code_len:
                # add the approximation of the actual number of bits this
                # transition will cost in the entropy coder, if requested in
                # the options
                if new_penalty:
                    # if new_penalty is also selected, assumed actual isn't
                    # necessarily an integer. round to the next int in this case
                    # so we still find it in the dict
                    assumed_actual = int(round(assumed_actual))
                    if assumed_actual > 255:
                        # in case it gets rounded out of the range of valid
                        # dict keys
                        assumed_actual = 255
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
