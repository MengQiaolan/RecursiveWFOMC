from __future__ import annotations
from logzero import logger
from wfomc_fo2.fol.sc2 import SC2
from wfomc_fo2.fol.utils import new_predicate, convert_counting_formula, unary_evidence_to_ccs, unary_evidence_to_pc

from wfomc_fo2.network.constraint import CardinalityConstraint, PartitionConstraint
from wfomc_fo2.fol.syntax import *
from wfomc_fo2.problems import WFOMCSProblem
from wfomc_fo2.fol.syntax import AUXILIARY_PRED_NAME, SKOLEM_PRED_NAME
from wfomc_fo2.utils.multinomial import MultinomialCoefficients
from wfomc_fo2.utils.third_typing import RingElement, Rational

from symengine import Symbol

class WFOMCContext_cc(object):
    """
    Context for WFOMC algorithm
    """

    def __init__(self, problem: WFOMCSProblem, use_partition_constraint: bool = False):
        self.domain: set[Const] = problem.domain
        self.sentence: SC2 = problem.sentence
        self.weights: dict[Pred, tuple[Rational, Rational]] = problem.weights
        self.cardinality_constraint: CardinalityConstraint = problem.cardinality_constraint
        if self.cardinality_constraint is None:
            self.cardinality_constraint = CardinalityConstraint()

        self.repeat_factor = 1
        self.unary_evidence = problem.unary_evidence
        self.use_partition_constraint = use_partition_constraint

        logger.info('sentence: \n%s', self.sentence)
        logger.info('domain: \n%s', self.domain)
        
        self.sentence_name = ''
        self.var_symbol = Symbol('a')
        for pred in self.sentence.preds():
            if pred.name == 'E1' or pred.name == 'E2' or pred.name == 'E3' or pred.name == 'E4' or pred.name == 'E5' or pred.name == 'E6' :
                self.sentence_name = 'match'
                self.weights[pred] = (self.var_symbol, Rational(1, 1))
            if pred.name == 'R1' or pred.name == 'R2' or pred.name == 'R3' or pred.name == 'R4' or pred.name == 'R5' or pred.name == 'R6' :
                self.sentence_name = 'ba'
                self.weights[pred] = (self.var_symbol, Rational(1, 1))
            if pred.name == 'F' :
                self.sentence_name = 'regular'
                self.weights[pred] = (self.var_symbol, Rational(1, 1))
            if pred.name == 'D' :
                self.sentence_name = 'derangement'
                self.weights[pred] = (self.var_symbol, Rational(1, 1))
        
        logger.info('weights:')
        for pred, w in self.weights.items():
            logger.info('%s: %s', pred, w)
        logger.info('cardinality constraint: %s', self.cardinality_constraint)
        logger.info('unary evidence: %s', self.unary_evidence)

        self.formula: QFFormula
        # for handling unary evidence
        self.partition_constraint: PartitionConstraint = None
        self._build()
        logger.info('Skolemized formula for WFOMC: \n%s', self.formula)
        logger.info('weights for WFOMC: \n%s', self.weights)
        logger.info('repeat factor: %d', self.repeat_factor)
        logger.info('partition constraint: %s', self.partition_constraint)

    def contain_cardinality_constraint(self) -> bool:
        return not self.cardinality_constraint.empty()

    def contain_partition_constraint(self) -> bool:
        return self.partition_constraint is not None

    def contain_existential_quantifier(self) -> bool:
        return self.sentence.contain_existential_quantifier()

    def get_weight(self, pred: Pred) -> tuple[RingElement, RingElement]:
        default = Rational(1, 1)
        if pred in self.weights:
            return self.weights[pred]
        return (default, default)

    def decode_result(self, res: RingElement):
        if not self.contain_cardinality_constraint():
            return res / self.repeat_factor
        res = self.cardinality_constraint.decode_poly(res)
        return res / self.repeat_factor

    def _skolemize_one_formula(self, formula: QuantifiedFormula) -> QFFormula:
        """
        Only need to deal with \forall X \exists Y: f(X,Y) or \exists X: f(X,Y)
        """
        quantified_formula = formula.quantified_formula
        quantifier_num = 1
        while(not isinstance(quantified_formula, QFFormula)):
            quantified_formula = quantified_formula.quantified_formula
            quantifier_num += 1

        formula: QFFormula = top
        ext_formula = quantified_formula
        if not isinstance(ext_formula, AtomicFormula):
            aux_pred = new_predicate(quantifier_num, AUXILIARY_PRED_NAME)
            aux_atom = aux_pred(X, Y) if quantifier_num == 2 else aux_pred(X)
            formula = formula & (ext_formula.equivalent(aux_atom))
            ext_formula = aux_atom

        if quantifier_num == 2:
            skolem_pred = new_predicate(1, SKOLEM_PRED_NAME)
            skolem_atom = skolem_pred(X)
        elif quantifier_num == 1:
            skolem_pred = new_predicate(0, SKOLEM_PRED_NAME)
            skolem_atom = skolem_pred()
        formula = formula & (skolem_atom | ~ext_formula)
        self.weights[skolem_pred] = (Rational(1, 1), Rational(-1, 1))
        return formula

    def _skolemize(self) -> QFFormula:
        """
        Skolemize the sentence
        """
        formula = self.sentence.uni_formula
        while(not isinstance(formula, QFFormula)):
            formula = formula.quantified_formula
        for ext_formula in self.sentence.ext_formulas:
            formula = formula & self._skolemize_one_formula(ext_formula)
        return formula

    def _build(self):
        self.formula = self.sentence.uni_formula
        while(not isinstance(self.formula, QFFormula)):
            self.formula = self.formula.quantified_formula

        if self.unary_evidence:
            if self.use_partition_constraint:
                logger.info('Use partition constraint to encode unary evidence')
                evi_formula, partition = unary_evidence_to_pc(self.unary_evidence)
                logger.info('formula to encode unary evidence: %s', evi_formula)
                logger.info('partition constraint: %s', partition)
                self.formula = self.formula & evi_formula
                self.partition_constraint = partition
            else:
                logger.info('Use cardinality constraint to encode unary evidence')
                evi_formula, ccs, repeat_factor = unary_evidence_to_ccs(
                    self.unary_evidence, self.domain
                )
                logger.info('formula to encode unary evidence: %s', evi_formula)
                logger.info('cardinality constraints: %s', ccs)
                self.formula = self.formula & evi_formula
                self.cardinality_constraint.extend_simple_constraints(ccs)
                self.repeat_factor *= repeat_factor

        self.ext_formulas = self.sentence.ext_formulas
        if self.sentence.contain_counting_quantifier():
            logger.info('translate SC2 to SNF')
            for cnt_formula in self.sentence.cnt_formulas:
                uni_formula, ext_formulas, cardinality_constraint, repeat_factor = \
                    convert_counting_formula(cnt_formula, self.domain)
                self.formula = self.formula & uni_formula
                self.ext_formulas = self.ext_formulas + ext_formulas
                self.cardinality_constraint.add_simple_constraint(*cardinality_constraint)
                self.repeat_factor *= repeat_factor

        if self.contain_cardinality_constraint():
            self.cardinality_constraint.build()

        for ext_formula in self.ext_formulas:
            self.formula = self.formula & self._skolemize_one_formula(ext_formula)

        # self.formula = self.formula.simplify()

        if self.contain_cardinality_constraint():
            self.weights.update(
                self.cardinality_constraint.transform_weighting(
                    self.get_weight,
                )
            )
