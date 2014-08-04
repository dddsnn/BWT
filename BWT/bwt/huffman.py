import bitstring
from collections import namedtuple
import heapq

# TODO raise errors where needed

def symbol_frequencies(iterable):
    """Determine the frequencies of symbols in an iterable.

    Args:
        iterable: An iterable.

    Returns:
        A dictionary mapping element that appears in iterable to the number of
        its occurrences in iterable.
    """
    # initialize
    freq = {i:0 for i in set(iterable)}
    for s in iterable:
        freq[s] += 1
    return freq

def codeword_lengths(weights):
    """Determine ideal codeword lengths for symbols with given weigths.

    Args:
        weights: A dictionary mapping symbols to weights.

    Returns:
        A dictionary mapping each of the symbols in weights to the length of its
        huffman codeword.
    """
    # if no symbols, return empty codebook
    if not weights:
        return {}
    # nodes to build the huffman tree with. data can be a child node tuple
    # (inner node) or a symbol (leaf node)
    # the id field is there as a tie breaker, to avoid unorderable types
    Node = namedtuple('Node', ['weight', 'id', 'data'])
    ChildNodes = namedtuple('ChildNodes', ['left', 'right'])
    Symbol = namedtuple('Symbol', ['symbol'])
    # make the huffman tree
    q = []
    for symbol in weights:
        # put all the symbols into the queue
        symbol_tuple = Symbol(symbol)
        heapq.heappush(q, Node(weights[symbol], id(symbol_tuple), symbol_tuple))
    # now merge the two lowest weight nodes until there is only one
    while len(q) > 1:
        left = heapq.heappop(q)
        right = heapq.heappop(q)
        weight = left.weight + right.weight
        child_nodes = ChildNodes(left, right)
        heapq.heappush(q, Node(weight, id(child_nodes), child_nodes))
    # the last node is the root
    root = heapq.heappop(q)

    # if root is the only node, give its symbol code length 1 and return
    if hasattr(root.data, 'symbol'):
        return {root.data.symbol:1}
    # go through the tree and assign code lengths
    def traverse_child_nodes(node, current_length):
        # if it's a leaf node, return the current code length for the symbol
        if hasattr(node.data, 'symbol'):
            return {node.data.symbol:current_length}
        # go into child nodes and merge their codebooks
        lengths = traverse_child_nodes(node.data.left, current_length + 1)
        lengths.update(traverse_child_nodes(node.data.right,
                                            current_length + 1))
        return lengths
    # start recursion
    return traverse_child_nodes(root, 0)

def canonical_codebook(codeword_lengths):
    """Make a canonical codebook based on given codeword lengths.

    https://en.wikipedia.org/w/index.php?title=Canonical_Huffman_code&oldid=600211231

    Args:
        codeword_lengths: A dictionary mapping symbols to the lengths of their
            huffman code.

    Returns:
        A dictionary mapping each of the symbols to its canonical huffman code
        represented as a Bits object.
    """
    # if input is empty, return empty dict
    if not codeword_lengths:
        return {}
    max_length = max(codeword_lengths.values())
    # make a list containing, at the i-th position, all symbols whose codeword
    # length is i
    length_buckets = [[s for s in codeword_lengths if codeword_lengths[s] == i]
                      for i in range(max_length + 1)]
    # sort all the buckets
    for bucket in length_buckets:
        bucket.sort()
    codebook = {}
    # set the code for the first element of the first non-empty bucket to all
    # zeros and remove it from its bucket
    for i, bucket in enumerate(length_buckets):
        if not bucket:
            continue
        last_code = bitstring.Bits(i)
        codebook[bucket[0]] = last_code
        del bucket[0]
        break
    # go through buckets and assign codes
    for i, bucket in enumerate(length_buckets):
        for symbol in bucket:
            # increment the last_code to get the current one
            last_code = bitstring.Bits(uint=last_code.uint + 1,
                                       length=len(last_code))
            if i != len(last_code):
                # if we are in a new bucket for longer codes, append zeros
                last_code += bitstring.Bits(i - len(last_code))
            # now add to the codebook
            codebook[symbol] = last_code
    return codebook

def serialize_codebook(codebook):
    """Serialize a canonical codebook.

    The symbols used in the codebook have to be integers between 0 and 255
    (i.e. byte values).
    The format is as follows:
    The first 8 bits are an unsigned integer giving the number of different
    symbols in the codebook minus one (meaning 00000000 stands for one symbol,
    00000001 for two etc.).
    The next 4 bits are an unsigned integer giving the number of bits with which
    each of the following numbers are encoded (so 0000 means 1 bit etc.).
    After this comes a sequence of unsigned integers, the i-th of which gives
    the number of symbols whose codeword is i bits long, starting at i=1. The
    sequence is over when the sum of their values is equal to the number of
    symbols in the codebook.
    After it comes the sequence of symbols used in the codebook, each
    represented as an 8 bit unsigned integer. They are sorted by the length of
    their codewords, beginning with the lowest. Tie breaker is their natural
    order.

    Args:
        codebook: A dictionary mapping byte values given as integers to
            canonical huffman codewords.

    Returns:
        A BitStream representing a serialized version of the codebook.
    """
    # if input is empty, return empty BitStream
    if not codebook:
        return bitstring.BitStream()
    num_symbols = len(codebook)
    max_length = max((len(c) for c in codebook.values()))
    # list of how many symbols have length i, starting at 1
    lengths = [len([v for v in codebook.values() if len(v) == i])
               for i in range(1, max_length + 1)]
    # determine how many bits are needed to store the max
    lengths_bit_len = (max(lengths)).bit_length()
    # buckets for symbols of length i, starting at 1
    length_buckets = [[s for s in codebook if len(codebook[s]) == i]
                      for i in range(1, max_length + 1)]
    # sort the buckets
    for bucket in length_buckets:
        bucket.sort()
    # piece together the result
    # first the number of symbols - 1, 8 bits
    res = bitstring.BitStream(uint=num_symbols - 1, length=8)
    # now number of bits for lengths sequence, 4 bits
    res += bitstring.Bits(uint=lengths_bit_len, length=4)
    # now the lengths, each encoded with lengths_bit_len bits
    for l in lengths:
        res += bitstring.Bits(uint=l, length=lengths_bit_len)
    # now the symbols from the sorted buckets, 8 bits each
    for bucket in length_buckets:
        for symbol in bucket:
            res += bitstring.Bits(uint=symbol, length=8)
    return res

def deserialize_codebook(bits):
    """Deserialize a canonical codebook.

    The input must be in the format described in serialize_codebook().

    Args:
        bits: A BitStream representing a serialized codebook.

    Returns:
        A codebook, i.e. a dictionary mapping symbols (integers) to code words
            (Bits objects).
    """
    # if input is empty, return empty dict
    if len(bits) == 0:
        return {}
    # first 8 bits: number of symbols - 1
    num_symbols = bits.read('uint:8') + 1
    # next 4 bits: lengths of the numbers in the lengths sequence
    lengths_bit_len = bits.read('uint:4')
    # read numbers from the sequence until they add up to num_symbols
    count = 0
    lengths = []
    while count < num_symbols:
        num = bits.read('uint:{0}'.format(lengths_bit_len))
        lengths.append(num)
        count += num
    # now read num_symbols times 8 bits to get the symbols
    symbols_list = []
    for _ in range(num_symbols):
        symbols_list.append(bits.read('uint:8'))
    # make codeword length dict
    cw_lengths = {}
    sym_iter = iter(symbols_list)
    for length, num in enumerate(lengths, start=1):
        for _ in range(num):
            cw_lengths[next(sym_iter)] = length
    return canonical_codebook(cw_lengths)

def encode_data(bs, codebook):
    """Encode data with a given codebook.

    Args:
        bitstring: Data to be encoded as a bytes object.
        codebook: A dictionary mapping symbols (byte values as integers) to code
            words (Bits objects).

    Returns:
        bitstring encoded according to codebook, as a BitStream.
    """
    res = bitstring.BitStream()
    for b in bs:
        res += codebook[b]
    return res

def decode_data(bits, codebook):
    """Decode data with a given codebook.

    Args:
        bits: Data to be decoded as a BitStream.
        codebook: A dictionary mapping symbols (byte values as integers) to code
            words (Bits objects).

    Returns:
        The decoded bits as a bytes object.
    """
    # inverse the codebook
    cb_inv = {v:k for (k, v) in codebook.items()}
    res = b''
    # since it's a prefix free code, read bits until a sequence is in the
    # codebook
    current_word = bitstring.Bits()
    while True:
        # TODO this is inefficient and slow
        while current_word not in cb_inv:
            try:
                new_bit = bits.read(1)
                if len(current_word) == 0:
                    # can't append to an empty Bits, or it will be turned into
                    # a BitStream
                    current_word = bitstring.Bits(bin=new_bit.bin)
                else:
                    current_word += new_bit
            except bitstring.ReadError:
                break
        else:
            # no errors reading, and we have a code word
            res += bytes([cb_inv[current_word]])
            current_word = bitstring.Bits()
            # continue with the outer loop
            continue
        # read() in the inner loop raised an error
        if len(current_word) != 0:
            # some bits couldn't be matched to a symbol
            # TODO
            raise Exception('there\'s leftovers')
        break
    return res

def encode_to_bits_static(bs):
    """Encode data to a bit string with static Huffman coding.

    Writes a header containing the codebook and the data so that the result is
    all you need to recreate the input.

    Args:
        bs: A bytes object.

    Returns:
        The Huffman encoded input as a BitStream.
    """
    weights = symbol_frequencies(bs)
    lengths = codeword_lengths(weights)
    cb = canonical_codebook(lengths)
    ser_cb = serialize_codebook(cb)
    enc = encode_data(bs, cb)
    return ser_cb + enc

def decode_from_bits_static(bits):
    """Decode a bit string containing a codebook and Huffman coded data.

    Args:
        A BitStream.

    Returns:
        A bytes object containing the data encoded by the input.
    """
    cb = deserialize_codebook(bits)
    return decode_data(bits, cb)

if __name__ == '__main__':
    with open('/home/dddsnn/Downloads/calgary/book1', 'rb') as file:
        bs = file.read(10000)
    s = b'aaaaaaaaaaaaaaaaab'
    enc = encode_to_bits_static(bs)
    dec = decode_from_bits_static(enc)
    print(bs == dec)

#     weights = symbol_frequencies(s)
#     lengths = codeword_lengths(weights)
#     codebook = canonical_codebook(lengths)
#     ser = serialize_codebook(codebook)
#     des_cb = deserialize_codebook(ser)
#     print(des_cb == codebook)
#     enc = encode_data(s, codebook)
#     dec = decode_data(enc, codebook)
#     print(s == dec)
