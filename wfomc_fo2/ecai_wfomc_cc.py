from __future__ import annotations
from collections import Counter

import sys

from enum import Enum
from functools import reduce
import math
import os
import argparse
import logging
import logzero
import pynauty

from sympy import factorint, factor_list

from logzero import logger
from typing import Callable
from contexttimer import Timer

from wfomc_fo2.problems import WFOMCSProblem
from wfomc_fo2.utils import RingElement, Rational
from wfomc_fo2.cell_graph import CellGraph
from wfomc_fo2.context import WFOMCContext_cc
from wfomc_fo2.parser import parse_input
from wfomc_fo2.fol.syntax import Const, Pred, QFFormula
from wfomc_fo2.utils.polynomial import coeff_dict, expand


class Algo(Enum):
    INCREMENTAL = 'inc'
    RECURSIVE = 'rec'
    RECURSIVE_real = 'rec_real'
    
    def __str__(self):
        return self.value

class IGCache(object):
    def init(self, domain_size: int):
        self.cache = []
        self.cache_hit_count = []
        for _ in range(domain_size):
            self.cache.append({})
            self.cache_hit_count.append(0)
    
    def get(self, level: int, color_kind: tuple[int], color_count: tuple[int], can_label):
        if color_kind not in self.cache[level]:
            self.cache[level][color_kind] = {}
            return None
        if color_count not in self.cache[level][color_kind]:
            self.cache[level][color_kind][color_count] = {}
            return None
        if can_label not in self.cache[level][color_kind][color_count]:
            return None
        self.cache_hit_count[level] += 1
        return self.cache[level][color_kind][color_count][can_label]
    
    def set(self, level: int, color_kind: tuple[int], color_count: tuple[int], can_label, value):
        if color_kind not in self.cache[level]:
            self.cache[level][color_kind] = {}
        if color_count not in self.cache[level][color_kind]:
            self.cache[level][color_kind][color_count] = {}
        self.cache[level][color_kind][color_count][can_label] = value

# cache for isormophic graph
IG_CACHE = IGCache()
# edge weight matrix of original cell graph
ORIGINAL_WEIGHT_EDGE = []
# how many layers is needed when convert edge-colored graph
COLOR_EDGE_d = 0
# the number of colors of edge
NUM_EDGE_COLOR = 0
# key: edge weight, value: edge color
EDGE_COLOR_MAP = {1:0}
# edge color matrix of original cell graph (ORIGINAL_WEIGHT_EDGE + EDGE_COLOR_MAP = COLOR_EDGE)
COLOR_EDGE = []
# the global no. of color of vertex
NO_VERTEX_COLOR = 0
# the number of vertries of the original cell graph
NUM_CELLS = 0
# the number of vertries of the extend cell graph
NUM_EXT_CELLS = 0
# key: vertex weight, value: vertex color
VERTEX_COLOR_MAP: dict[any, int] = {}
# adjacency_dict
ADJACENCY_DICT = {}
# reduce the call of pynauty.certificate
CACHE_FOR_NAUTY = {}
# key is factor and value is the index of factor
FACTOR_DICT = {}
# the index of factor 0
ZERO_FACTOR_INDEX = -1
# the factor set of edge weights
EDGE_FACTOR_SET = []


def update_factor_dict(factor):
    global FACTOR_DICT, ZERO_FACTOR_INDEX
    if FACTOR_DICT.get(factor) is None:
        FACTOR_DICT[factor] = len(FACTOR_DICT)
        if factor == 0:
            ZERO_FACTOR_INDEX = FACTOR_DICT[factor]
def prime_init_factors(cell_weights, edge_weights):
    '''
    prime init factors for the cell weights and edge weights (including expression with symbols)
    all factors are stored in FACTOR_DICT
    '''
    for w in cell_weights:
        factored_list = factor_list(w)
        coef = factored_list[0]
        syms = factored_list[1]
        for k,_ in factorint(coef).items():
            update_factor_dict(k)
        for sym in syms:
            update_factor_dict(sym[0])
    for rs in edge_weights:
        for r in rs:
            factored_list = factor_list(r)
            coef = factored_list[0]
            syms = factored_list[1]
            for k,_ in factorint(coef).items():
                update_factor_dict(k)
            for sym in syms:
                update_factor_dict(sym[0])

def get_init_factor_set(cell_weights, edge_weights):
    cell_factor_set = []
    for w in cell_weights:
        vector = [0] * len(FACTOR_DICT)
        factored_list = factor_list(w)
        coef = factored_list[0]
        syms = factored_list[1]
        for k,v in factorint(coef).items():
            vector[FACTOR_DICT[k]] = v
        for sym in syms:
            vector[FACTOR_DICT[sym[0]]] = int(sym[1])
        cell_factor_set.append(tuple(vector))
    global EDGE_FACTOR_SET
    for i in range(len(edge_weights)):
        rs = edge_weights[i]
        vecs = []
        for j in range(len(rs)):
            r = rs[j]
            vector = [0] * len(FACTOR_DICT)
            factored_list = factor_list(r)
            coef = factored_list[0]
            syms = factored_list[1]
            for k,v in factorint(coef).items():
                vector[FACTOR_DICT[k]] = v
            for sym in syms:
                vector[FACTOR_DICT[sym[0]]] = int(sym[1])
            vecs.append(tuple(vector))
        EDGE_FACTOR_SET.append(vecs)
    return cell_factor_set

def get_vertex_color(weight):
    global NO_VERTEX_COLOR
    if weight not in VERTEX_COLOR_MAP:
        VERTEX_COLOR_MAP[weight] = NO_VERTEX_COLOR
        NO_VERTEX_COLOR += 1
    return VERTEX_COLOR_MAP[weight]

def edgeWeight_To_edgeColor():
    global NUM_EDGE_COLOR, EDGE_COLOR_MAP, COLOR_EDGE
    NUM_EDGE_COLOR = 0
    EDGE_COLOR_MAP = {1:0}
    COLOR_EDGE = []
    
    for lst in ORIGINAL_WEIGHT_EDGE:
        tmp_list = []
        for w in lst:
            if w not in EDGE_COLOR_MAP:
                NUM_EDGE_COLOR += 1
                EDGE_COLOR_MAP[w] = NUM_EDGE_COLOR
            tmp_list.append(EDGE_COLOR_MAP[w])
        COLOR_EDGE.append(tmp_list)
    
    global COLOR_EDGE_d, NUM_EXT_CELLS
    COLOR_EDGE_d = math.ceil(math.log2(NUM_EDGE_COLOR+1))
    NUM_EXT_CELLS = NUM_CELLS * COLOR_EDGE_d
    
def cellWeight_To_vertexColor(cell_weights):
    vertex_colors = []
    for w in cell_weights:
        vertex_colors.append(get_vertex_color(w))
    color_dict = Counter(vertex_colors)
    color_kind = tuple(sorted(color_dict))
    color_count = tuple(color_dict[num] for num in color_kind)
    return vertex_colors, color_kind, color_count

def calculate_adjacency_dict():
    # Generate new edges
    adjacency_dict = {}
    for i in range(NUM_EXT_CELLS):
        adjacency_dict[i] = []
    
    c2layers = {}
    for k in range(NUM_EDGE_COLOR+1):
        bi = bin(k)
        bi = bi[2:]
        bi = bi[::-1]
        layers = [i for i in range(len(bi)) if bi[i] == '1']
        c2layers[k] = layers
    
    for i in range(NUM_CELLS):
        for j in range(NUM_CELLS):
            layers = c2layers[COLOR_EDGE[i][j]]
            for l in layers:
                adjacency_dict[l*NUM_CELLS+i].append(l*NUM_CELLS+j)
    
    # The vertical threads (each corresponding to one vertex of the original graph) 
    # can be connected using either paths or cliques.
    for i in range(NUM_CELLS):
        clique = [i + j*NUM_CELLS for j in range(COLOR_EDGE_d)]
        for ii in clique:
            for jj in clique:
                if ii == jj:
                    continue
                adjacency_dict[ii].append(jj)
    global ADJACENCY_DICT
    ADJACENCY_DICT = adjacency_dict

def create_graph():
    global GRAPH
    GRAPH = pynauty.Graph(NUM_EXT_CELLS, 
                          directed=False, 
                          adjacency_dict=ADJACENCY_DICT)

def adjust_vertex_coloring(colored_vertices):
    '''
    Adjust the color no. of vertices to make the color no. start from 0 and be continuous.
    
    Args:
        colored_vertices: list[int]
            The color no. of vertices.
            
        Returns:
            new_colored_vertices: list[int]
                The adjusted color no. of vertices.
            num_color: int
                The number of colors. 
    
    Example:
        colored_vertices = [7, 5, 7, 3, 5, 7]
        new_colored_vertices, num_color = adjust_vertex_coloring(colored_vertices)
        print(new_colored_vertices)  # [0, 1, 0, 2, 1, 0]
        print(num_color)  # 3
    '''
    color_map = {}
    new_colored_vertices = []
    num_color = 0
    for c in colored_vertices:
        if c not in color_map:
            color_map[c] = num_color
            num_color += 1
        new_colored_vertices.append(color_map[c])
    return new_colored_vertices, num_color

def extend_vertex_coloring(colored_vertices, no_color):
    '''
    Extend the vertex set to convert colored edge
    
    Args:
        colored_vertices: list[int]
            The color no. of vertices.
        no_color: int
            The number of colors.
    
    Returns:
        vertex_coloring: list[set[int]]
            The color set of vertices.
    
    Example:
        colored_vertices = [0, 1, 0, 2, 1, 0]
        no_color = 3
        vertex_coloring = extend_vertex_coloring(colored_vertices, no_color)
        print(vertex_coloring)  # [{0, 2, 5}, {1, 4}, {3}]
    '''
    # Extend the vertex set to convert colored edge
    ext_colored_vertices = []
    for i in range(COLOR_EDGE_d):
        ext_colored_vertices += [x + no_color * i for x in colored_vertices]
    
    # Get color set of vertices
    no_color *= COLOR_EDGE_d
    vertex_coloring = [ set() for _ in range(no_color)]
    for i in range(len(ext_colored_vertices)):
        c = ext_colored_vertices[i]
        vertex_coloring[c].add(i)
    
    return vertex_coloring

def update_graph(colored_vertices):
    # for speed up, we have modified the function 'set_vertex_coloring' in graph.py of pynauty
    GRAPH.set_vertex_coloring(colored_vertices)
    return GRAPH

def value_filter(value):
    expanded_expr = expand(value)
    res = Rational(0,1)
    for degrees, coef in coeff_dict(expanded_expr, [VAR_SYMBOL]):
        if degrees[0] <= SINGLE_VAR_LIMIT:
            res += coef*VAR_SYMBOL**degrees[0]
    return res

def dfs_wfomc(cell_weights, domain_size, cell_factor_tuple_list):  
    res = 0
    for l in range(NUM_CELLS):
        w_l = cell_weights[l]
        new_cell_weights = [cell_weights[i] * ORIGINAL_WEIGHT_EDGE[l][i] for i in range(NUM_CELLS)]
        if domain_size - 1 == 1:
            value = sum(new_cell_weights)
        else:
            new_cell_factor_tuple_list = []
            for i in range(NUM_CELLS):
                new_cell_factor_tuple_list.append([x+y for x, y in zip(cell_factor_tuple_list[i], EDGE_FACTOR_SET[l][i])])
                if ZERO_FACTOR_INDEX >= 0 and new_cell_factor_tuple_list[-1][ZERO_FACTOR_INDEX] > 0:
                    new_cell_factor_tuple_list[-1][ZERO_FACTOR_INDEX] = 1
                new_cell_factor_tuple_list[-1] = tuple(new_cell_factor_tuple_list[-1])
            original_vertex_colors, vertex_color_kind, vertex_color_count = cellWeight_To_vertexColor(new_cell_factor_tuple_list) # convert cell weights to vertex colors
            adjust_vertex_colors, no_color = adjust_vertex_coloring(original_vertex_colors) # adjust the color no. of vertices to make them start from 0 and be continuous
            if tuple(adjust_vertex_colors) not in CACHE_FOR_NAUTY:
                can_label = pynauty.certificate(update_graph(extend_vertex_coloring(adjust_vertex_colors, no_color)))
                CACHE_FOR_NAUTY[tuple(adjust_vertex_colors)] = can_label
            else:
                can_label = CACHE_FOR_NAUTY[tuple(adjust_vertex_colors)]

            value = IG_CACHE.get(domain_size-1, vertex_color_kind, vertex_color_count, can_label)
            if value is None:
                value = dfs_wfomc(new_cell_weights, domain_size - 1, new_cell_factor_tuple_list)
                value = value_filter(value)
                IG_CACHE.set(domain_size-1, vertex_color_kind, vertex_color_count, can_label, value)
                
        res += w_l * value
    return res

def dfs_wfomc_real(cell_weights, domain_size):
    raise NotImplementedError('Not implemented')

def dp_wfomc_init(formula: QFFormula,
                  domain: set[Const],
                  get_weight: Callable[[Pred], tuple[RingElement, RingElement]],
                  leq_pred: Pred = None,
                  algo = Algo.RECURSIVE_real) -> RingElement:
    domain_size = len(domain)
    cell_graph = CellGraph(formula, get_weight, leq_pred)
    cell_weights = cell_graph.get_all_weights()[0]
    edge_weights = cell_graph.get_all_weights()[1]
    
    IG_CACHE.init(domain_size)
    global ORIGINAL_WEIGHT_EDGE, NUM_CELLS
    # not change in the same problem
    NUM_CELLS = len(cell_weights)
    ORIGINAL_WEIGHT_EDGE = edge_weights
    edgeWeight_To_edgeColor()
    calculate_adjacency_dict()
    create_graph()
        
    if algo == Algo.RECURSIVE:
        prime_init_factors(cell_weights, edge_weights)
        cell_factor_tuple_list = get_init_factor_set(cell_weights, edge_weights)
        return dfs_wfomc(cell_weights, domain_size, cell_factor_tuple_list)
    elif algo == Algo.RECURSIVE_real:
        return dfs_wfomc_real(cell_weights, domain_size)
    else:
        raise ValueError('Invalid algorithm')

def incremental_wfomc(formula: QFFormula,
                      domain: set[Const],
                      get_weight: Callable[[Pred], tuple[RingElement, RingElement]],
                      leq_pred: Pred = None) -> RingElement:
    cell_graph = CellGraph(formula, get_weight, leq_pred)
    cells = cell_graph.get_cells()
    n_cells = len(cells)
    domain_size = len(domain)
    
    cell_weights = cell_graph.get_all_weights()[0]
    edge_weights = cell_graph.get_all_weights()[1]

    table = dict(
        (tuple(int(k == i) for k in range(n_cells)),
         cell_weights[i])
         for i in range(n_cells)
    )
    for _ in range(domain_size - 1):
        old_table = table
        table = dict()
        for j in range(n_cells):
            w = cell_weights[j]
            for ivec, w_old in old_table.items():
                w_new = w_old * w * reduce(
                    lambda x, y: x * y,
                    (edge_weights[j][l] ** int(ivec[l]) for l in range(n_cells)),
                    Rational(1, 1)
                )
                w_new = value_filter(w_new)
                
                ivec = list(ivec)
                ivec[j] += 1
                ivec = tuple(ivec)

                w_new = w_new + table.get(ivec, Rational(0, 1))
                table[tuple(ivec)] = w_new
    ret = sum(table.values())
    return ret

SINGLE_VAR_LIMIT = -1
VAR_SYMBOL = None
def wfomc(problem: WFOMCSProblem, algo: Algo = Algo.INCREMENTAL,
          use_partition_constraint: bool = False) -> Rational:

    context = WFOMCContext_cc(problem, use_partition_constraint)
    
    global SINGLE_VAR_LIMIT, VAR_SYMBOL
    VAR_SYMBOL = context.var_symbol
    if context.sentence_name == 'ba':
        SINGLE_VAR_LIMIT = len(context.domain) * 3
    
    leq_pred = Pred('LEQ', 2)
    with Timer() as t:
        if algo == Algo.INCREMENTAL:
            res = incremental_wfomc(
                context.formula, context.domain,
                context.get_weight, leq_pred
            )
        elif algo == Algo.RECURSIVE_real:
            res = dp_wfomc_init(
                context.formula, context.domain, context.get_weight, leq_pred, algo
            )
        elif algo == Algo.RECURSIVE:
            res = dp_wfomc_init(
                context.formula, context.domain, context.get_weight, leq_pred, algo
            )
    res = context.decode_result(res)
    logger.info('WFOMC time: %s', t.elapsed)
    return res

def parse_args():
    parser = argparse.ArgumentParser(
        description='WFOMC for MLN',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--output_dir', '-o', type=str,
                        default='./check-points')
    parser.add_argument('--algo', '-a', type=Algo,
                        choices=list(Algo), default=Algo.INCREMENTAL)
    parser.add_argument('--domain_recursive',
                        action='store_true', default=False,
                        help='use domain recursive algorithm '
                             '(only for existential quantified MLN)')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    sys.setrecursionlimit(int(1e6))
    args = parse_args()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if args.debug:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.CRITICAL)
        # logzero.loglevel(logging.INFO)
    logzero.logfile('{}/log.txt'.format(args.output_dir), mode='w')

    with Timer() as t:
        problem = parse_input(args.input)
    logger.info('Parse input: %ss', t)

    res = wfomc(problem, algo=args.algo)
    logger.critical('WFOMC (arbitrary precision): %s', res)
