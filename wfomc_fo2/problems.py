from __future__ import annotations

from wfomc_fo2.fol.sc2 import SC2, to_sc2
from wfomc_fo2.fol.syntax import AtomicFormula, Const, Pred, top, AUXILIARY_PRED_NAME, \
    Formula, QuantifiedFormula, Universal, Equivalence
from wfomc_fo2.fol.utils import new_predicate
from wfomc_fo2.network.constraint import CardinalityConstraint
from wfomc_fo2.utils.polynomial import Rational
from fractions import Fraction
import math


class WFOMCSProblem(object):
    """
    A weighted first-order model counting/sampling problem.
    """

    def __init__(self, sentence: Formula,
                 domain: set[Const],
                 weights: dict[Pred, tuple[Rational, Rational]],
                 cardinality_constraint: CardinalityConstraint = None,
                 unary_evidence: set[AtomicFormula] = None):
        self.domain: set[Const] = domain
        try:
            self.sentence: SC2 = to_sc2(sentence)
        except:
            raise ValueError('Sentence must be a valid SC2 formula.')
        self.weights: dict[Pred, tuple[Rational, Rational]] = weights
        self.cardinality_constraint: CardinalityConstraint = cardinality_constraint
        if self.cardinality_constraint is None:
            self.cardinality_constraint = CardinalityConstraint()
        self.unary_evidence = unary_evidence
        if self.unary_evidence is not None:
            # check if the evidence is unary and consistent with the domain
            for atom in self.unary_evidence:
                if len(atom.args) != 1:
                    raise ValueError('Evidence must be unary.')
                if atom.args[0] not in self.domain:
                    raise ValueError(f'Evidence must be consistent with the domain: {atom.args[0]} not in {self.domain}.')


    def __str__(self) -> str:
        s = ''
        s += 'Domain: \n'
        s += '\t' + str(self.domain) + '\n'
        s += 'Sentence: \n'
        s += '\t' + str(self.sentence) + '\n'
        s += 'Weights: \n'
        s += '\t' + str(self.weights) + '\n'
        if self.cardinality_constraint is not None:
            s += 'Cardinality Constraint: \n'
            s += '\t' + str(self.cardinality_constraint) + '\n'
        if self.unary_evidence is not None:
            s += 'Unary Evidence: \n'
            s += '\t' + str(self.unary_evidence) + '\n'
        return s

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        s = ''
        s += 'Domain: \n'
        s += '\t' + str(self.domain) + '\n'
        s += 'Sentence: \n'
        s += '\t' + str(self.sentence) + '\n'
        s += 'Weights: \n'
        s += '\t' + str(self.weights) + '\n'
        if self.cardinality_constraint is not None:
            s += 'Cardinality Constraint: \n'
            s += '\t' + str(self.cardinality_constraint) + '\n'
        return s

    def __repr__(self) -> str:
        return str(self)


class MLNProblem(object):
    """
    A Markov Logic Network problem.
    """

    def __init__(self, rules: tuple[list[tuple[Rational, Rational]], list[Formula]],
                 domain: set[Const],
                 cardinality_constraint: CardinalityConstraint):
        self.rules = rules
        # self.formulas: rules[1]
        # self.formula_weights: = dict(zip(rules[1], rules[0]))
        self.domain: set[Const] = domain
        self.cardinality_constraint: CardinalityConstraint = cardinality_constraint


def MLN_to_WFOMC(mln: MLNProblem):
    sentence = top
    weightings: dict[Pred, tuple[Rational, Rational]] = dict()
    for weighting, formula in zip(*mln.rules):
        free_vars = formula.free_vars()
        if weighting != float('inf'):
            aux_pred = new_predicate(len(free_vars), AUXILIARY_PRED_NAME)
            formula = Equivalence(formula, aux_pred(*free_vars))
            weightings[aux_pred] = (Rational(Fraction(math.exp(weighting)).numerator,
                                             Fraction(math.exp(weighting)).denominator), Rational(1, 1))
        for free_var in free_vars:
            formula = QuantifiedFormula(Universal(free_var), formula)
        sentence = sentence & formula

    # sc2 = to_sc2(sentence)
    return WFOMCSProblem(sentence, mln.domain, weightings, mln.cardinality_constraint)
