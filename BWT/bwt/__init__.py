from collections import namedtuple


BWEncodeResult = namedtuple('BWEncodeResult', ['firsts', 'encoded'])

PartialMTFAnalysisResult = namedtuple('PartialMTFAnalysisResult',
                                          ['raw', 'length', 'length_rec',
                                           'num_chars', 'max_code', 'median',
                                           'mean'])
