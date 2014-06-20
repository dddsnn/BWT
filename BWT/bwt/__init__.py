from collections import namedtuple


BWEncodeResult = namedtuple('BWEncodeResult', ['firsts', 'encoded'])

PartialMTFAnalysisResult = namedtuple('PartialMTFAnalysisResult',
                                          ['raw', 'length', 'length_rec',
                                           'num_chars', 'max_code', 'median',
                                           'mean'])
AuxData = namedtuple('AuxData',
                     ['raw', 'bw_code', 'mtf_code', 'bw_subcodes',
                       'partial_mtf_subcodes', 'partial_mtf_analyses',
                       'bw_subhistograms', 'huffman_codeword_lengths',
                       'mtf_mean_steps', 'freq_lists'])
