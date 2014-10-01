"""This module contains functions to encode and decode BWT and MTF (full and
    partial."""

from bwt import BWEncodeResult
from collections.abc import Mapping
from collections import namedtuple
import itertools as it

NUM_CHARS = 25

def bw_table(text):
    """Create the Burrows-Wheeler table.

    This is only for demo purposes with small inputs, as it has terribly
    terrible memory requirements.
    """
    table = [text]
    for i in range(len(text) - 1):
        table.append(table[i][1:] + table[i][0])
    table.sort()
    return table

def first_column(table):
    """Get the first column of the BW table."""
    return ''.join([s[0] for s in table])

def last_column(table):
    """Get the last column of the BW table (the encoded string)."""
    return ''.join([s[-1] for s in table])

def print_demo(text):
    """Print the sorted BW table, and the string of first and last letters."""
    table = bw_table(text)
    firsts = first_column(table)
    encoded = last_column(table)

    for i, row in enumerate(table):
        print('{0:2} {1}'.format(i, row))
    print()
    print(firsts)
    print(encoded)

def order_list_to_dict(order_list, distinct_syms):
    result = {}
    for i, s in enumerate(order_list):
        # ignore multiple occurrences, count the first one
        if not s in result:
            result[s[0]] = i
            max_order = i
    # check that order is complete (all bytes from bs are in it)
    for b in distinct_syms:
        if not b in result:
            # append any missing bytes
            result[b] = max_order + 1
            max_order += 1
    return result

def make_order_lists(bs, orders, depth, start_depth=0):
    """Make a list of lists by which the columns of the BW table can be ordered.

    Args:
        bs: The input bytes object.
        orders: A list of orders, see bw_encode().
        depth: How many lists should be returned. This is important for orders
            that aren't generic, as order lists in this case are only valid for
            one column of the BW table. If there aren't enough orders to make
            lists, the last order is used as a default.
        start_depth: The number of lists that will be omitted from the
            beginning, so not all lists have to be created again when more are
            needed.

    Returns:
        A list of as many lists as there are orders in the orders argument. Each
        list is as long as bs. If the n-th list is zipped up with bs and the
        result sorted with the list values as the key, the bs values will be
        sorted according to the n-th order in orders.
    """
    distinct_syms = set(bs)
    # turn the order lists of bytes objects into dicts {byte: value} for
    # ordering
    order_dicts = []
    for order in orders:
        order_dict = {}
        if isinstance(order, Mapping):
            for prefix, order_list in order.items():
                order_dict[prefix] = order_list_to_dict(order_list,
                                                        distinct_syms)
        else:
            # None to indicate it's a general order
            order_dict[None] = order_list_to_dict(order, distinct_syms)
        # append the order dict to the list
        order_dicts.append(order_dict)
    # make lists by which the strings will be sorted
    order_lists = []
    for d in range(start_depth, depth):
        try:
            order_dict = order_dicts[d]
        except IndexError:
            order_dict = order_dicts[-1]
            if None in order_dict and order_lists:
                # if the default order is a general order, we can just copy the
                # last order list
                order_lists.append(order_lists[-1][:])
                continue
        order_list = []
        if None in order_dict:
            # only one general order
            for b in bs:
                order_list.append(order_dict[None][b])
        else:
            # lenght of prefix, assume all are equal
            l = len(next(iter(order_dict.keys())))
            for i, b in enumerate(bs):
                if i >= l:
                    prefix = bs[i - l:i]
                else:
                    prefix = bs[i - l:] + bs[:i]
                order_list.append(order_dict[prefix][b])
        order_lists.append(order_list)
    return order_lists

def bw_encode(bs, orders=None):
    """BW encode a string of bytes according to specified sort orders.

    Args:
        bs: The data to be encoded as a bytes object.
        orders: A list of orders by which each of the columns of the BW table
            should be sorted. An order is either a list of bytes objects of
            length one, in the order in which they should be sorted, to indicate
            a general order valid for any subcontext. Alternatively it can be a
            dictionary mapping the symbols in the columns before (as a bytes
            object) to a specific order for that context. The keys of this
            dictionary must have the same length as the order's position in the
            orders list. An order must contain all byte values in appearing in
            the relevant context. If less orders are given than there are input
            bytes, the last order is used for all following columns. If None is
            given, the natural order is assumed.

    Returns:
        A bwt.BWTEncodeResult containing the BW code and the list of unique
        first sequences the table was sorted by, in order.
    """
    # if no order was given, go with natural
    if not orders:
        orders = [[bytes([x]) for x in range(256)]]

    l = len(bs)
    num_chars = min(NUM_CHARS, l)  # number of characters to save for ordering
    # take the bs twice so we can get long substrings from the end
    bs = bs + bs
    order_lists = make_order_lists(bs, orders, num_chars)

    # make tuples (pos in the table, ordering first bytes, last byte,
    # actual first bytes)
    tuples = []
    for i in range(l):
        # make the list the tuple will be sorted by, according to the order
        # lists for the appropriate columns
        order_firsts = [order_lists[j][i + j] for j in range(num_chars)]
        tuples.append((i, order_firsts, bs[i - 1],
                       bs[i:i + num_chars]))
    tuples.sort(key=lambda x: x[1])
    # check for duplicate first bytes and compare longer strings
    for i in range(len(tuples) - 1):
        if tuples[i][3] == tuples[i + 1][3]:
            num_affected = 2
            # there might be more
            if i + 2 < len(tuples):
                while tuples[i + num_affected][3] == tuples[i][3]:
                    num_affected += 1
                    # break if we reach the end of the list
                    if not i + num_affected < len(tuples):
                        break
            # make tuples of the affected indices, but with enough bytes
            # to sort them
            long_tuples = []
            # take the affected tuples, take their positions and sort again
            # using more characters
            for t in tuples[i:i + num_affected]:
                others = [x for x in tuples[i:i + num_affected] if x != t]
                # find j such that the string starting at t[0] with length j is
                # unique, then add the string to the long tuples
                for j in range(num_chars, l):
                    # exclude tuples from others where the symbol at offset j is
                    # different; these are already distinct
                    others = [o for o in others if bs[t[0] + j] == bs[o[0] + j]]
                    # if the others is empty, j symbols are enough
                    if not others:
                        # make the list the tuple will be sorted by, see above
                        num_chars = j + 1
                        if num_chars > len(order_lists):
                            # we need new order lists with more elements for
                            # more columns
                            if len(orders) > len(order_lists):
                                # there are orders we've not made lists for
                                extension = make_order_lists(bs, orders,
                                                             num_chars,
                                                             len(order_lists))
                                order_lists.extend(extension)
                            else:
                                # we're using the default order, can just copy
                                times = num_chars - len(order_lists)
                                order_lists.extend(order_lists[-1:] * times)
                        order_firsts = [order_lists[k][t[0] + k]
                                        for k in range(num_chars)]
                        long_tuples.append((t[0], order_firsts,
                                            bs[t[0] - 1],
                                            bs[t[0]:t[0] + j + 1]))
                        break
            long_tuples.sort(key=lambda x: x[1])
            # replace the short tuples in the original list
            tuples[i:i + num_affected] = long_tuples
    firsts = [t[3] for t in tuples]
    encoded = bytes([t[2] for t in tuples])
    start_index = next((i for i, t in enumerate(tuples) if t[0] == 1))
    result = BWEncodeResult(firsts, encoded, start_index)
    return result

def bw_decode(bs, start_idx, orders=None):
    """Decode a string transformed with the BWT using an arbitrary number of
    sort orders.

    Does not produce correct results in its current form. E.g., the input
    b'missishjkdgfhjkdasdasdasdjklsdg' will lead to a livelock when the orders
    are a list containing the natural order, the reverse of the natural order,
    and then the natural order again.

    Args:
        bs: Input bytes object.
        start_idx: The index in bs of the first symbol of the output.
        orders: The list of orders that was used to encode bs.

    Returns:
        Anything between the correct result, garbage and maximum recursion depth
        reached. Is supposed to return the reverse of the BWT applied to bs.
    """
    SeqTuple = namedtuple('SeqTuple', ['idx_seq', 'sym_seq', 'possible',
                                       'need_more'])
    def next_indices(idx, history):
        level = 0
        history.append(idx)
        sym = first_col[idx]
        # how many times the same symbol occurs in the first column before
        num_in_first_col = sum(1 for i, b in enumerate(first_col)
                               if b == sym and i < idx)
        possible_seq = [SeqTuple([i], [first_col[i]], True, True)
                        for i, b in enumerate(bs) if b == sym]
        while True:
            possible_seq.sort(key=lambda x:
                              [order_dicts[min(len(order_dicts) - 1,
                                               i)][x.sym_seq[i - 1]]
                               for i in range(1, min(level + 2,
                                                     len(x.sym_seq) + 1))])
            if not possible_seq[num_in_first_col].possible:
                # if the sequence at position num_in_first_col is not possible
                # after sorting, there is no correct continuation
                # TODO is this correct? necessary?
                return None
            if len(order_dicts) <= level + 2:
                # no need to get any more indices and sort any further, since
                # all following symbols are already sorted with the default
                # order and num_in_first_col is the correct continuation
                return possible_seq[num_in_first_col]
            next_sym = possible_seq[num_in_first_col].sym_seq[level]
            for i, t in enumerate(possible_seq):
                if not t.possible and not t.need_more:
                    # it's already impossible and doesn't need more symbols, no
                    # need to worry about it
                    continue
                l = t.idx_seq
                # keep indices that haven't already appeared and where the next
                # symbol matches
                if (t.possible and l[level] not in history
                    and first_col[l[level]] == next_sym):
                    continue
                # in case the history plus the sequence is as long as the input
                # and now wraps around to the correct symbol, that's also ok
                elif (t.possible and len(history) + len(l) >= len(bs) + 1
                      and history[0] == l[len(bs) - len(history)]
                      and first_col[l[level]] == next_sym):
                    continue
                elif (l[level] in history
                    and first_col[l[level]] == next_sym):
                    # not possible because it's in history, but we need more
                    # syms to sort the still possible ones correctly
                    # TODO can this lead to infinite loops? where, to get the
                    # idx for sym a you need more syms for the (impossible) idx
                    # b, and for b you need syms for a
                    possible_seq[i] = SeqTuple(t[0], t[1], False, True)
                else:
                    # not possible anymore, set possible to False
                    possible_seq[i] = SeqTuple(t[0], t[1], False, False)
            num_possible = sum(1 for t in possible_seq if t.possible)
            if num_possible == 0:
                # return None to indicate no successor index can be found with
                # the given history
                return None
            elif num_possible == 1:
                # find and return the possible one
                return next(t for t in possible_seq if t.possible)
            elif all(len(history) + len(l.idx_seq) >= len(bs) + 1
                     for l in (t for t in possible_seq if t.possible)):
                length = len(bs) - len(history)
                first_res = [bs[i] for i in possible_seq[0].idx_seq[:length]]
                if all([bs[i] == first_res for i in l[0][:length]]
                       for l in (t for t in possible_seq[1:] if t.possible)):
                    # all possible sequences have maximum length and yield the
                    # same result, it's safe to return
                    return SeqTuple(possible_seq[0][0][:length],
                            possible_seq[0][1][:length], True, True)
                else:
                    raise Exception  # TODO i don't think this can happen
            for i, t in enumerate(possible_seq):
                if t.need_more and len(t.idx_seq) <= level + 1:
                    if t.possible:
                        next_seq = next_indices(t.idx_seq[-1],
                                            history + t.idx_seq[:-1])
                    else:
                        # if the sequence is not possible but we need more syms,
                        # give an empty history
                        # TODO this can cause liveloops
                        next_seq = next_indices(t.idx_seq[-1], [])
                    if next_seq is None:
                        # indicates this sequence is impossible
                        # TODO what if we needed more syms, but can't get any?
                        #     can this happen?
                        if not t.possible and t.need_more:
                            raise Exception  # TODO
                        possible_seq[i] = SeqTuple(t[0], t[1], False, False)
                        continue
                    t.idx_seq.extend(next_seq.idx_seq)
                    t.sym_seq.extend(next_seq.sym_seq)
            num_possible = sum(1 for t in possible_seq if t.possible)
            if num_possible == 0:
                # no sequence has a correct continuation
                return None
            level += 1
    # if no order was given, assume natural
    if not orders:
        orders = [[bytes([x]) for x in range(256)]]
    distinct_syms = set(bs)
    order_dicts = [order_list_to_dict(o, distinct_syms) for o in orders]
    first_order_list = make_order_lists(bs, orders, 1)
    first_col = bytes((x[1] for x in sorted(zip(first_order_list[0], bs),
                                            key=lambda x:x[0])))
    result_idx = [start_idx]
    current_idx_seq = []
    last_idx = start_idx
    while len(result_idx) < len(bs):
        if not current_idx_seq:
            # find more indices
            idx_sym_tuple = next_indices(last_idx, result_idx[:-1])
            current_idx_seq = idx_sym_tuple.idx_seq
            last_idx = current_idx_seq[-1]
        result_idx.append(current_idx_seq[0])
        del current_idx_seq[0]
    return bytes([bs[i] for i in result_idx])

def bw_decode_2_orders(bs, start_idx, orders=None):
    """Decode a string encoded with the BWT using up to two orders.

    Args:
        bs: The encoded input, as a bytes object.
        start_idx: The index in bs containing the first symbol of the output.
        order: A list of at most two orders, as passed to bw_encode(). If more
            than two are provided, a ValueError is raised.

    Returns:
        The reverse BWT of bs using the given orders, as a bytes object.
    """
    # if no order was given, assume natural
    if not orders:
        orders = [[bytes([x]) for x in range(256)]]
    if len(orders) > 2:
        raise ValueError('Only up to 2 different orders are allowed.')
    distinct_syms = set(bs)
    if len(orders) == 2:
        second_order_dict = {}
        if isinstance(orders[1], Mapping):
            for prefix, order_list in orders[1].items():
                second_order_dict[prefix[0]] = order_list_to_dict(order_list,
                                                        distinct_syms)
        else:
            # copy the general order for all prefixes
            order_dict = order_list_to_dict(orders[1], distinct_syms)
            for s in distinct_syms:
                second_order_dict[s] = order_dict
    first_order_list = make_order_lists(bs, orders, 1)
    first_col = bytes((x[1] for x in sorted(zip(first_order_list[0], bs),
                                            key=lambda x:x[0])))
    result_idx = [start_idx]
    while len(result_idx) < len(bs):
        idx = result_idx[-1]
        next_sym = first_col[idx]
        num_in_first_col = sum(1 for i, b in enumerate(first_col)
                               if b == next_sym and i < idx)
        possible_idx = [i for i, b in enumerate(bs) if b == next_sym]
        if len(orders) == 2:
            # reordering is only necessary if there are 2 orders
            possible_idx.sort(key=lambda x:second_order_dict[next_sym]
                              [first_col[x]])
        result_idx.append(possible_idx[num_in_first_col])
    return bytes([bs[i] for i in result_idx])

def mtf_partial_encode(bs):
    """Do a partial MTF encode of a byte sequence.

    Partial MTF means that values that haven't occurred before are marked as -1.

    Args:
        bs: The bytes object to be encoded.

    Returns:
        The MTF code of the input, with new values marked as -1.
    """
    alphabet = []
    result = []
    for byte in bs:
        if byte in alphabet:
            # append index of the byte
            index = alphabet.index(byte)
            result.append(index)
            # shift alphabet
            alphabet.pop(index)
            alphabet.insert(0, byte)
        else:
            # append -1 to signal new character
            result.append(-1)
            alphabet.insert(0, byte)
    return result

def mtf_encode(bs):
    """Encode a byte sequence with Move-to-front.

    Args:
        bs: The bytes object to be encoded.

    Returns:
        The MTF code of the input, with an initial alphabet [0..255].
    """
    # initialize the alphabet in natural order
    alphabet = list(range(256))
    result = []
    for byte in bs:
        if byte in alphabet:
            # append index of the character
            index = alphabet.index(byte)
            result.append(index)
            # shift alphabet
            alphabet.pop(index)
            alphabet.insert(0, byte)
        else:
            raise ValueError(str(byte) + ' is not a byte.')
    return result
