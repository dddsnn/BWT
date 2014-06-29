'''
Created on 23.05.2014

@author: dddsnn
'''

import ctypes as ct
import warnings
from bwt import BWEncodeResult

def bw_table(text):
    '''Create the Burrows-Wheeler table.'''
    table = [text]
    for i in range(len(text) - 1):
        table.append(table[i][1:] + table[i][0])
    table.sort()
    return table

def first_column(table):
    '''Get the first column of the BW table.'''
    return ''.join([s[0] for s in table])

def last_column(table):
    '''Get the last column of the BW table (the encoded string).'''
    return ''.join([s[-1] for s in table])

def print_demo(text):
    '''Print the sorted BW table, and the string of first and last letters.'''
    table = bw_table(text)
    firsts = first_column(table)
    encoded = last_column(table)

    for row in table:
        print(row)
        print()
    print()
    print(firsts)
    print(encoded)

def bw_encode(bytes_, order=None):
    '''BW encode a string of bytes.'''
    # if no order was given, go with natural
    if not order:
        order = [bytes([x]) for x in range(256)]
    # turn the order list of characters into a dict {byte: value} for ordering
    order_dict = {}
    for i, s in enumerate(order):
        # ignore multiple occurrences, count the first one
        if not s in order_dict:
            order_dict[s] = i
            max_order = i
        else:
            warnings.warn('multiple occurence of symbol or sequence {0} in '
                          'order. ignoring.'.format(s))
    # check that order is complete (all bytes from bytes_ are in it)
    for b in bytes_:
        if not bytes([b]) in order_dict:
            warnings.warn('symbol or sequence {0} not in the custom order but '
                          'in the bytes_.  appending'.format(b))
            # append any missing bytes
            order_dict[bytes([b])] = max_order + 1
            max_order += 1
    # dict with all keys in order_dict that lists more specific keys
    # (e.g. b'a' -> [b'adf', b'as']
    spec_dict = {s:[] for s in order_dict}
    for s1 in spec_dict:
        for s2 in spec_dict:
            if s1 != s2 and len(s2) > len(s1) and s2[:len(s1)] == s1:
                # more specific sequence found, add to dict
                spec_dict[s1].append(s2)
    # sort the lists in spec_dict from most specific to most general
    for s in spec_dict:
        spec_dict[s].sort(key=lambda x:len(x), reverse=True)
    # make a list by which the strings will be sorted
    order_list = []
    for i, b in enumerate(bytes_ + bytes_):
        # in case there is a more specific sequence
        for s in spec_dict[bytes([b])]:
            if bytes_[i:i + len(s)] == s:
                # more specific sequence has been found, append the order for
                # that sequence and break (so we don't get into the else)
                order_list.append(order_dict[s])
                break
        else:
            order_list.append(order_dict[bytes([b])])
    # take it twice so we can get long substrings from the end
    order_list = order_list + order_list

    NUM_CHARS = 25  # number of characters to save for ordering
    l = len(bytes_)
    # make tuples (pos in the table, first byte, last byte)
    tuples = []
    for i in range(l):
        tuples.append((i, order_list[i:i + NUM_CHARS], bytes_[i - 1]))
    tuples.sort(key=lambda x: x[1])
    # check for duplicate first bytes and compare longer strings
    for i in range(len(tuples) - 1):
        if tuples[i][1] == tuples[i + 1][1]:
            num_affected = 2
            # there might be more
            if i + 2 < len(tuples):
                while tuples[i + num_affected][1] == tuples[i][1]:
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
                    other_chars = [order_list[o[0] + j] for o in others]
                    # if the current string's char doesn't appear in it,
                    # j characters are enough
                    if not order_list[t[0] + j] in other_chars:
                        long_tuples.append((t[0], order_list[t[0]:t[0] + j + 1],
                                            bytes_[t[0] - 1]))
                        break
            long_tuples.sort(key=lambda x: x[1])
            # replace the short tuples in the original list
            tuples[i:i + num_affected] = long_tuples
    firsts = bytes([t[1][0] for t in tuples])
    encoded = bytes([t[2] for t in tuples])
    result = BWEncodeResult(firsts, encoded)
    return result

def mtf_partial_enc(bytes_):
    alphabet = []
    result = []
    for byte in bytes_:
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

def mtf_enc(bytes_):
    '''Encode bytes_ with mtf.'''
    # initialize the alphabet in natural order
    alphabet = list(range(256))
    result = []
    for byte in bytes_:
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

def huffman_enc(bytes_):
    '''Huffman encode a list of bytes.'''
    lib = ct.cdll.LoadLibrary('../libhuffman.so')
    enc_in_len = len(bytes_)
    if enc_in_len == 0:
        # the library segfaults for length 0. just return empty list
        return []
    enc_in = ct.create_string_buffer(enc_in_len)
    enc_in = ct.cast(enc_in, ct.POINTER(ct.c_ubyte))
    for i in range(enc_in_len):
        if bytes_[i] > 256 or bytes_[i] < 0:
            raise ValueError('Byte value not between 0 and 256')
        enc_in[i] = bytes_[i]
    enc_out = ct.pointer(ct.c_ubyte())
    enc_out_len = ct.c_uint()
    lib.huffman_encode_memory(enc_in, ct.c_uint(enc_in_len), ct.byref(enc_out),
                              ct.byref(enc_out_len))
    result = []
    for i in range(enc_out_len.value):
        result.append(enc_out[i])
    return bytes(result)

def huffman_dec(bytes_):
    lib = ct.cdll.LoadLibrary('../libhuffman.so')
    dec_in_len = len(bytes_)
    dec_in = ct.create_string_buffer(dec_in_len)
    dec_in = ct.cast(dec_in, ct.POINTER(ct.c_ubyte))
    for i in range(dec_in_len):
        if bytes_[i] > 256 or bytes_[i] < 0:
            raise ValueError('Byte value not between 0 and 256')
        dec_in[i] = bytes_[i]
    dec_out = ct.pointer(ct.c_ubyte())
    dec_out_len = ct.c_uint()
    lib.huffman_decode_memory(dec_in, ct.c_uint(dec_in_len), ct.byref(dec_out),
                              ct.byref(dec_out_len))
    result = []
    for i in range(dec_out_len.value):
        result.append(dec_out[i])
    return bytes(result)
