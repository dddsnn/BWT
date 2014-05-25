'''
Created on 10.05.2014

@author: dddsnn
'''

from collections import namedtuple
import bwt.coder as cd
from openopt import TSP
import networkx as nx
import numpy as np

def analyze_partial_mtf(code):
    raw = code
    length = len(code)
    sorted_code = sorted(code)
    max_code = sorted_code[-1]
    # remove -1s
    sorted_code = [c for c in sorted_code if c != -1]
    # number of all recurring symbols (those that aren't -1)
    l = len(sorted_code)
    # number of -1s is number of all symbols minus number of all except -1s
    num_chars = length - l
    # number of recurring characters (not -1)
    length_rec = length - num_chars
    if l != 0:
        if l % 2 == 0:
            median = (sorted_code[(l // 2) - 1] +
                                 sorted_code[l // 2]) / 2
        else:
            median = sorted_code[l // 2]
        mean = sum(sorted_code) / l
    else:
        median = 0
        mean = 0
    PartialMTFAnalysisResult = namedtuple('PartialMTFAnalysisResult',
                                          ['raw', 'length', 'length_rec',
                                           'num_chars', 'max_code', 'median',
                                           'mean'])
    result = PartialMTFAnalysisResult(raw, length, length_rec, num_chars,
                                      max_code, median, mean)
    return result

def analyze_transition(bw1, bw2):
    a = analyze_partial_mtf(cd.mtf_partial_enc(bw1))
    b = analyze_partial_mtf(cd.mtf_partial_enc(bw2))
    ab = analyze_partial_mtf(cd.mtf_partial_enc(bw1 + bw2))
    TransitionResult = namedtuple('TransitionResult', ['left', 'right',
                                                         'together', 'diff'])
    length = TransitionResult(a.length, b.length, ab.length,
                               abs(a.length - b.length))
    num_chars = TransitionResult(a.num_chars, b.num_chars, ab.num_chars,
                                  a.num_chars + b.num_chars - ab.num_chars)
    max_code = TransitionResult(a.max_code, b.max_code, ab.max_code,
                             ab.max_code - max(a.max_code, b.max_code))
    # to avoid division by zero
    if ab.length_rec == 0:
        ab_length_rec = 1
    else:
        ab_length_rec = ab.length_rec
    weigthed_median = ((a.median * a.length_rec + b.median * b.length_rec)
                       / ab_length_rec) - ab.median
    median = TransitionResult(a.median, b.median, ab.median, weigthed_median)
    # TODO by ignoring new characters (-1s) in the mean, good transitions get a
    # penalty, because a character that was already there before the transition
    # probably affects the mean in a bad way, while an entirely new character
    # won't affect it at all
    # need to give new characters a "penalty" value in the partial mtf encode
    # (but only for the target of the transition, lots of new characters in the
    # source don't mean anything bad)
    # or maybe an entirely new metric that takes this into account
    weighted_mean = ((a.mean * a.length_rec + b.mean * b.length_rec)
                     / ab_length_rec) - ab.mean
    mean = TransitionResult(a.mean, b.mean, ab.mean, weighted_mean)
    TransitionAnalysisResult = namedtuple('TransitionAnalysisResult',
                                          ['length', 'num_chars', 'max_code',
                                           'median', 'mean'])
    result = TransitionAnalysisResult(length, num_chars, max_code, median, mean)
    return result

def make_table_string(table):
    result = ''
    col_width = [max(len(str(x)) for x in col) for col in zip(*table)]
    for line in table:
        result += "| " + " | ".join("{:{}}".format(x, col_width[i])
                                for i, x in enumerate(line)) + " |\n"
    return result

def print_transition_analyses(text):
    code = cd.bw_encode(text)
    first = code.firsts
    bw_code = code.encoded
    subcodes = {c: '' for c in first}
    for i in range(len(first)):
        subcodes[first[i]] += bw_code[i]
    transitions = {a + b: analyze_transition(subcodes[a], subcodes[b])
                   for a in subcodes.keys()
                   for b in subcodes.keys() if a != b}

    for k in sorted(transitions.keys()):
        print('<' + k + '>\n')
        header = ['', 'left', 'right', 'together', 'diff']
        len_line = ['lenght', transitions[k].length.left, transitions[k].length.right, transitions[k].length.together, transitions[k].length.diff]
        num_chars_line = ['num_chars', transitions[k].num_chars.left, transitions[k].num_chars.right, transitions[k].num_chars.together, transitions[k].num_chars.diff]
        max_line = ['max_code', transitions[k].max_code.left, transitions[k].max_code.right, transitions[k].max_code.together, transitions[k].max_code.diff]
        median_line = ['median', transitions[k].median.left, transitions[k].median.right, transitions[k].median.together, transitions[k].median.diff]
        mean_line = ['mean', transitions[k].mean.left, transitions[k].mean.right, transitions[k].mean.together, transitions[k].mean.diff]
        table = [header, len_line, num_chars_line, max_line, median_line, mean_line]
        print(make_table_string(table))

if __name__ == '__main__':
#     g = nx.Graph()
#     g.add_edges_from([(i, j, {'cost': i + j}) for i in range(50) for j in range(50)])
#     problem = TSP(g, objective='cost')
#     result = problem.solve('glpk')
    f = open('/home/dddsnn/Downloads/calgary/book2')
    text = f.read()
    text = text[0:50000]
    print_transition_analyses(text)
#     bw = cd.bw_encode(text)
# #     print(bw2)
#     mtf = cd.mtf_enc(bw[1])
#     hf = cd.huffman_enc(mtf)
#     dec = cd.huffman_dec(hf)
#     print('{0} / {1} : {2}'
#         .format(len(text), len(hf), len(hf) / len(text)))
