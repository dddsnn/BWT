from collections import namedtuple


BWEncodeResult = namedtuple('BWEncodeResult', ['firsts', 'encoded'])

PartialMTFAnalysisResult = namedtuple('PartialMTFAnalysisResult',
                                          ['raw', 'length', 'length_rec',
                                           'num_chars', 'max_code', 'median',
                                           'mean'])

TransitionAnalysisResult = namedtuple('TransitionAnalysisResult',
                                          ['data', 'length', 'num_chars',
                                           'max_code', 'median', 'mean',
                                           'chapin1'])

TransitionData = namedtuple('TransitionData',
                                    ['length', 'num_chars', 'max_code',
                                     'median', 'mean'])

TransitionDataSet = namedtuple('TransitionDataSet',
                              ['left', 'right', 'together'])
