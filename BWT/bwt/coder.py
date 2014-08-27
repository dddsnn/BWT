import warnings
from bwt import BWEncodeResult
from collections.abc import Mapping

def bw_table(text):
    """Create the Burrows-Wheeler table."""
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

    for row in table:
        print(row)
    print()
    print(firsts)
    print(encoded)

def make_order_lists(bs, orders):
    """Make a list of lists by which the columns of the BW table can be ordered.

    Args:
        bs: The input bytes object.
        orders: A list of orders, see bw_encode().

    Returns:
        A list of as many lists as there are orders in the orders argument. Each
        list is as long as bs. If the n-th list is zipped up with bs and the
        result sorted with the list values as the key, the bs values will be
        sorted according to the n-th order in orders.
    """
    def list_to_dict(order_list, distinct_syms):
        result = {}
        for i, s in enumerate(order):
            # ignore multiple occurrences, count the first one
            if not s in result:
                result[s] = i
                max_order = i
            else:
                warnings.warn('multiple occurence of symbol {0} in '
                              'order. ignoring.'.format(s))
        # check that order is complete (all bytes from bs are in it)
        for b in distinct_syms:
            if not bytes([b]) in result:
                warnings.warn('symbol {0} not in the custom order but '
                              'in the bs.  appending'.format(b))
                # append any missing bytes
                result[bytes([b])] = max_order + 1
                max_order += 1
        return result
    distinct_syms = set(bs)
    # turn the order lists of bytes objects into dicts {byte: value} for
    # ordering
    order_dicts = []
    for order in orders:
        order_dict = {}
        if isinstance(order, Mapping):
            for prefix, order_list in order.items():
                order_dict[prefix] = list_to_dict(order_list, distinct_syms)
        else:
            # None to indicate it's a general order
            order_dict[None] = list_to_dict(order, distinct_syms)
        # append the order dict to the list
        order_dicts.append(order_dict)
    # make lists by which the strings will be sorted
    order_lists = []
    for order_dict in order_dicts:
        order_list = []
        if None in order_dict:
            # only one general order
            order_dict = order_dict[None]
            for b in bs:
                order_list.append(order_dict[bytes([b])])
        else:
            # lenght of prefix, assume all are equal
            l = len(next(iter(order_dict.keys())))
            for i, b in enumerate(bs):
                if i >= l:
                    prefix = bs[i - l:i]
                else:
                    prefix = bs[i - l:] + bs[:i]
                order_list.append(order_dict[prefix][bytes([b])])
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
    # take the bs twice so we can get long substrings from the end
    bs = bs + bs
    order_lists = make_order_lists(bs, orders)

    NUM_CHARS = 25  # number of characters to save for ordering
    # make tuples (pos in the table, ordering first bytes, last byte,
    # actual first bytes)
    tuples = []
    for i in range(l):
        # make the list the tuple will be sorted by, according to the order
        # lists for the appropriate columns
        order_firsts = []
        # first the columns for which there is a special ordering
        for j in range(min(NUM_CHARS, len(order_lists))):
            order_firsts.append(order_lists[j][i + j])
        # now use the last ordering for the rest
        order_firsts.extend(order_lists[-1][i + min(NUM_CHARS,
                                                    len(order_lists)):i + NUM_CHARS])
        tuples.append((i, order_firsts, bs[i - 1],
                       bs[i:i + NUM_CHARS]))
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
                # if the j-th character in the string is different from all the
                # other affected strings, add the string of length j to the
                # long tuples
                for j in range(NUM_CHARS, l):
                    # make a list of chars the other strings have a position j
                    other_chars = [bs[o[0] + j] for o in others]
                    # if the current string's char doesn't appear in it,
                    # j characters are enough
                    if not bs[t[0] + j] in other_chars:
                        # make the list the tuple will be sorted by, see above
                        order_firsts = []
                        num_chars = j + 1
                        for k in range(min(num_chars, len(order_lists))):
                            order_firsts.append(order_lists[k][t[0] + k])
                        order_firsts.extend(order_lists[-1][t[0] + min(num_chars, len(order_lists)):
                                             t[0] + num_chars])
                        long_tuples.append((t[0], order_firsts,
                                            bs[t[0] - 1],
                                            bs[t[0]:t[0] + j + 1]))
                        break
            long_tuples.sort(key=lambda x: x[1])
            # replace the short tuples in the original list
            tuples[i:i + num_affected] = long_tuples
    firsts = [t[3] for t in tuples]
    encoded = bytes([t[2] for t in tuples])
    result = BWEncodeResult(firsts, encoded)
    return result

def mtf_partial_enc(bs):
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

def mtf_enc(bs):
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
