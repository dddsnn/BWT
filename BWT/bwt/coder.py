'''
Created on 23.05.2014

@author: dddsnn
'''

import ctypes as ct
from collections import namedtuple

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

def bw_encode(text):
    num_chars = 25  # number of characters to save for ordering
    l = len(text)
    # loop the text around so we can get long substrings from the end
    looptext = text + text
    # make tuples (pos in the table, first chars, last char)
    tuples = []
    for i in range(l):
        tuples.append((i, looptext[i:i + num_chars], text[i - 1]))
    tuples.sort(key=lambda x: x[1])
    # check for duplicates and compare more characters
    for i in range(len(tuples) - 1):
        if tuples[i][1] == tuples[i + 1][1]:
            affected_idx = 2
            # there might be more
            if i + 2 < len(tuples):
                while tuples[i + affected_idx][1] == tuples[i][1]:
                    affected_idx += 1
                    # break if we reach the end of the list
                    if not i + affected_idx < len(tuples):
                        break
            # make tuples of the affected indices, but with enough characters
            # to sort them
            long_tuples = []
            # take the affected tuples, take their positions and sort again
            # using more characters
            for t in tuples[i:i + affected_idx]:
                others = [x for x in tuples[i:i + affected_idx] if x != t]
                # if the j-th character in the string is different from all the
                # other affected strings, add the string of length j to the
                # long tuples
                for j in range(25, l):
                    # make a list of chars the other strings have a position j
                    other_chars = [looptext[o[0] + j] for o in others]
                    # if the current string's char doesn't appear in it,
                    # j characters are enough
                    if not looptext[t[0] + j] in other_chars:
                        long_tuples.append((t[0], looptext[t[0]:t[0] + j + 1],
                                            text[t[0] - 1]))
                        break
            long_tuples.sort(key=lambda x: x[1])
            # replace the short tuples in the original list
            tuples[i:i + affected_idx] = long_tuples
    firsts = ''.join([t[1][0] for t in tuples])
    encoded = ''.join([t[2] for t in tuples])
    BWEncodeResult = namedtuple('BWEncodeResult', ['firsts', 'encoded'])
    result = BWEncodeResult(firsts, encoded)
    return result

def print_demo(text):
    """Print the sorted BW table, and the string of first and last letters."""
    table = bw_table(text)
    firsts = first_column(table)
    encoded = last_column(table)

    for row in table:
        print(row)
        print()
    print()
    print(firsts)
    print(encoded)

def mtf_partial_enc(text):
    alphabet = []
    result = []
    for char in text:
        if char in alphabet:
            # append index of the character
            index = alphabet.index(char)
            result.append(index)
            # shift alphabet
            alphabet.pop(index)
            alphabet.insert(0, char)
        else:
            # append -1 to signal new character
            result.append(-1)
            alphabet.insert(0, char)
    return result

def mtf_enc(text):
    """Encode text with mtf. Only ASCII (for now)."""
    # initialize list of ascii characters
    alphabet = [chr(i) for i in range(128)]
    result = []
    for char in text:
        if char in alphabet:
            # append index of the character
            index = alphabet.index(char)
            result.append(index)
            # shift alphabet
            alphabet.pop(index)
            alphabet.insert(0, char)
        else:
            # not an ascii character
            raise ValueError(char + ' is not an ASCII character')
    return result

def huffman_enc(byte_list):
    """Huffman encode a list of bytes."""
    lib = ct.cdll.LoadLibrary('../libhuffman.so')
    enc_in_len = len(byte_list)
    enc_in = ct.create_string_buffer(enc_in_len)
    enc_in = ct.cast(enc_in, ct.POINTER(ct.c_ubyte))
    for i in range(enc_in_len):
        if byte_list[i] > 256 or byte_list[i] < 0:
            raise ValueError('Byte value not between 0 and 256')
        enc_in[i] = byte_list[i]
    enc_out = ct.pointer(ct.c_ubyte())
    enc_out_len = ct.c_uint()
    lib.huffman_encode_memory(enc_in, ct.c_uint(enc_in_len), ct.byref(enc_out),
                              ct.byref(enc_out_len))
    result = []
    for i in range(enc_out_len.value):
        result.append(enc_out[i])
    return result

def huffman_dec(byte_list):
    lib = ct.cdll.LoadLibrary('../libhuffman.so')
    dec_in_len = len(byte_list)
    dec_in = ct.create_string_buffer(dec_in_len)
    dec_in = ct.cast(dec_in, ct.POINTER(ct.c_ubyte))
    for i in range(dec_in_len):
        if byte_list[i] > 256 or byte_list[i] < 0:
            raise ValueError('Byte value not between 0 and 256')
        dec_in[i] = byte_list[i]
    dec_out = ct.pointer(ct.c_ubyte())
    dec_out_len = ct.c_uint()
    lib.huffman_decode_memory(dec_in, ct.c_uint(dec_in_len), ct.byref(dec_out),
                              ct.byref(dec_out_len))
    result = []
    for i in range(dec_out_len.value):
        result.append(dec_out[i])
    return result

if __name__ == '__main__':
    pass
#     for i in range(2, 100):
#         orig = list(range(i))
#         enc = huffman_enc(orig)
#         dec = huffman_dec(enc)
#         print(str(orig == dec) + ' : ' + str(len(enc) / len(orig)))
