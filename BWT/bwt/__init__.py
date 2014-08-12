from collections import namedtuple


BWEncodeResult = namedtuple('BWEncodeResult', ['firsts', 'encoded'])

PartialMTFAnalysisResult = namedtuple('PartialMTFAnalysisResult',
                                          ['raw', 'length', 'length_rec',
                                           'num_chars', 'max_code'])
AuxData = namedtuple('AuxData',
                     ['raw', 'num_symbols', 'firsts', 'bw_code', 'mtf_code',
                      'bw_subcodes',
                      'partial_mtf_subcodes', 'partial_mtf_analyses',
                      'bw_subhistograms', 'huffman_codeword_lengths_complete',
                      'huffman_codeword_lengths_sparse', 'mtf_mean_steps',
                      'mtf_median_steps', 'freq_lists'])
