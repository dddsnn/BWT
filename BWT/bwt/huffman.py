import bitstring as bs
from collections import namedtuple
import heapq

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
        last_code = bs.Bits(i)
        codebook[bucket[0]] = last_code
        del bucket[0]
        break
    # go through buckets and assign codes
    for i, bucket in enumerate(length_buckets):
        for symbol in bucket:
            # increment the last_code to get the current one
            last_code = bs.Bits(uint=last_code.uint + 1, length=len(last_code))
            if i != len(last_code):
                # if we are in a new bucket for longer codes, append zeros
                last_code += bs.Bits(i - len(last_code))
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
        codebook: A dictionary mapping byte values to canonical huffman
            codewords.

    Returns:
        A BitStream representing a serialized version of the codebook.
    """
    num_symbols = len(codebook)
    max_length = max((len(c) for c in codebook.values()))
    # list of how many symbols have length i, starting at 1
    lengths = [list(codebook.values()).count(i)
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
    res = bs.BitStream(uint=num_symbols - 1, length=8)
    # now number of bits - 1 for lengths sequence, 4 bits
    res += bs.Bits(uint=lengths_bit_len, length=4)
    # now the lengths, each encoded with lengths_bit_len bits
    for l in lengths:
        res += bs.Bits(uint=l, length=lengths_bit_len)
    # now the symbols from the sorted buckets, 8 bits each
    for bucket in length_buckets:
        for symbol in bucket:
            res += bs.Bits(uint=symbol, length=8)
    return res

def deserialize_codebook(bits):
    """Deserialize a canonical codebook.

    The input must be in the format described in serialize_codebook().
    This function returns the inverse of the input to serialize_codebook(), so
    the result is not a mapping symbol->code but code->symbol.

    Args:
        bits: A BitStream representing a serialized codebook.

    Returns:
        A codebook, i.e. a dictionary mapping code words to symbols.
    """
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
    for i in range(num_symbols):
        symbols_list.append(bits.read('uint:8'))
    # TODO

if __name__ == '__main__':
    weights = symbol_frequencies(b'aabc')
    lengths = codeword_lengths(weights)
    cb = canonical_codebook(lengths)
    ser = serialize_codebook(cb)
    print(len(ser))
    deserialize_codebook(ser)
