import ctypes as ct
import warnings
from bwt import BWEncodeResult

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

def specializations(string_list):
    """Find all specialized strings in a list.

    For every element in the list, find all other elements in the list that are
    specializations of that element, i.e. that start with the same symbols, but
    are longer.

    Args:
        string_list: A list of strings.

    Returns:
        A dictionary mapping every element of the input to a list of more
        specific strings in the input, sorted from most specific to most
        general, i.e. from longest to shortest.
        E.g.: {'a': ['asd', 'as'], 'as': ['asd'], 'asd': []}
    """
    spec_dict = {s:[] for s in string_list}
    for s1 in spec_dict:
        for s2 in spec_dict:
            if s1 != s2 and len(s2) > len(s1) and s2[:len(s1)] == s1:
                # more specific sequence found, add to dict
                spec_dict[s1].append(s2)
    # sort the lists in spec_dict from most specific to most general
    for s in spec_dict:
        spec_dict[s].sort(key=lambda x:len(x), reverse=True)
    return spec_dict

def bw_encode(bs, orders=None):
    """BW encode a string of bytes according to specified sort orders.

    Args:
        bs: The data to be encoded as a bytes object.
        orders: A list of bytes objects representing sort orders for each of the
            columns in the BW table. An order must contain all byte values in
            bs, in the order in which they should be sorted. If less orders are
            given than there are input bytes, the last order is used for all
            following columns. If None is given, the natural order is assumed.

    Returns:
        A bwt.BWTEncodeResult containing the BW code and the list of unique
        first sequences the table was sorted by, in order.
    """
    # if no order was given, go with natural
    if not orders:
        orders = [[bytes([x]) for x in range(256)]]
    # turn the order lists of bytes objects into dicts {byte: value} for
    # ordering
    order_dicts = []
    for order in orders:
        order_dict = {}
        for i, s in enumerate(order):
            # ignore multiple occurrences, count the first one
            if not s in order_dict:
                order_dict[s] = i
                max_order = i
            else:
                warnings.warn('multiple occurence of symbol or sequence {0} in '
                              'order. ignoring.'.format(s))
        # check that order is complete (all bytes from bs are in it)
        for b in bs:
            if not bytes([b]) in order_dict:
                warnings.warn('symbol or sequence {0} not in the custom order but '
                              'in the bs.  appending'.format(b))
                # append any missing bytes
                order_dict[bytes([b])] = max_order + 1
                max_order += 1
        # append the order dict to the list
        order_dicts.append(order_dict)
    # make lists by which the strings will be sorted
    order_lists = []
    for order_dict in order_dicts:
        order_list = []
        # dict with all keys in order_dict that lists more specific keys
        # (e.g. b'a' -> [b'adf', b'as']
        spec_dict = specializations(order_dict)
        for i, b in enumerate(bs + bs):
            # in case there is a more specific sequence
            for s in spec_dict[bytes([b])]:
                if bs[i:i + len(s)] == s:
                    # more specific sequence has been found, append the order for
                    # that sequence and break (so we don't get into the else)
                    order_list.append(order_dict[s])
                    break
            else:
                order_list.append(order_dict[bytes([b])])
        # append this list to the list of order lists
        order_lists.append(order_list)

    l = len(bs)
    # take the bs  twice so we can get long substrings from the end
    # (order list is already twice as long)
    bs = bs + bs

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

def huffman_enc(bs):
    """Huffman encode a list of bytes."""
    lib = ct.cdll.LoadLibrary('../libhuffman.so')
    enc_in_len = len(bs)
    if enc_in_len == 0:
        # the library segfaults for length 0. just return empty list
        return []
    enc_in = ct.create_string_buffer(enc_in_len)
    enc_in = ct.cast(enc_in, ct.POINTER(ct.c_ubyte))
    for i in range(enc_in_len):
        if bs[i] > 256 or bs[i] < 0:
            raise ValueError('Byte value not between 0 and 256')
        enc_in[i] = bs[i]
    enc_out = ct.pointer(ct.c_ubyte())
    enc_out_len = ct.c_uint()
    lib.huffman_encode_memory(enc_in, ct.c_uint(enc_in_len), ct.byref(enc_out),
                              ct.byref(enc_out_len))
    result = []
    for i in range(enc_out_len.value):
        result.append(enc_out[i])
    return bytes(result)

def huffman_dec(bs):
    lib = ct.cdll.LoadLibrary('../libhuffman.so')
    dec_in_len = len(bs)
    dec_in = ct.create_string_buffer(dec_in_len)
    dec_in = ct.cast(dec_in, ct.POINTER(ct.c_ubyte))
    for i in range(dec_in_len):
        if bs[i] > 256 or bs[i] < 0:
            raise ValueError('Byte value not between 0 and 256')
        dec_in[i] = bs[i]
    dec_out = ct.pointer(ct.c_ubyte())
    dec_out_len = ct.c_uint()
    lib.huffman_decode_memory(dec_in, ct.c_uint(dec_in_len), ct.byref(dec_out),
                              ct.byref(dec_out_len))
    result = []
    for i in range(dec_out_len.value):
        result.append(dec_out[i])
    return bytes(result)
