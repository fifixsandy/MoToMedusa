#ifndef INTERFACE_H
#define INTERFACE_H


#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>    
#include <stddef.h> 
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>
typedef struct LEAF_TYPE_IMPL LEAF_TYPE_IMPL;
#ifndef QBDD_TYPE_DEFINED
typedef void* qBDD;
#endif


typedef struct LEAF_TYPE{
    LEAF_TYPE_IMPL* pImpl;
} LEAF_TYPE;


/**
 * @brief Creates new qBDD
 * @param target Target variable number
 * @param low qBDD representing low branch of the new qBDD
 * @param high qBDD representing high of the new qBDD
 */
qBDD newqBDD(unsigned int target, qBDD low, qBDD high);
qBDD qBDD_protect(qBDD);
qBDD qBDD_unprotect(qBDD);

qBDD qBDD_false();
qBDD qBDD_true();

int qBDD_isFalse(qBDD);
int qBDD_isInternal(qBDD);
int qBDD_isTerminal(qBDD);
qBDD qBDD_getLow(qBDD);
qBDD qBDD_getHigh(qBDD);
size_t qBDD_getVar(qBDD);
LEAF_TYPE qBDD_getTerminalValue(qBDD);
qBDD qBDD_maketerminal(size_t, void*);

int qBDD_leafcount(qBDD);

double qBDD_calculateProb(qBDD);
double interface_prob_sum(qBDD t, uint32_t xt, char* curr_state, int n);

size_t qBDD_symbolicMapLType();
size_t qBDD_symbolicValLType();
size_t qBDD_classicLType();


/**
 * @brief Performs a operation on a target var node of a given qBDD
 * with optional control nodes
 */
qBDD bdd_operation(qBDD operand, size_t* targets, size_t controlNum, qBDD(*op)(size_t tgt, qBDD, qBDD));
qBDD bdd_operation_guarded(qBDD operand, size_t* targets, size_t controlNum, qBDD(*op)(size_t tgt, qBDD));

qBDD binary_apply(qBDD l, qBDD r, LEAF_TYPE(*op)(LEAF_TYPE, LEAF_TYPE));
qBDD binary_apply_guarded(qBDD l, qBDD r, qBDD(*op)(qBDD, qBDD));
qBDD binary_apply_guarded_param(qBDD l, qBDD r, qBDD(*op)(qBDD, qBDD, size_t), size_t param);
qBDD unary_apply(qBDD l, LEAF_TYPE(*op)(LEAF_TYPE));
qBDD unary_apply_guarded(qBDD l, qBDD(*op)(qBDD, size_t), size_t arg);
qBDD unary_apply_param(qBDD l, qBDD(*op)(qBDD, size_t), size_t arg);
qBDD cube(int value, int width, qBDD *variables, qBDD leaf1, qBDD leaf0);

void incInvSqrtCoeff();
void incInvSqrtCoeffSymb();
void initInvSqrtCoeffSymb();
void clearInvSqrtCoeffNormal();
void clearInvSqrtCoeffSymb();
void mulInvSqrtSymbCoeff(unsigned long mul);
void addInvSqrtCoeffs();
void resetInvSqrtCoeffSymb();
mpz_srcptr getInvSqrtCoeff(void);
mpz_srcptr getInvSqrtCoeffSymb(void);


void validateOperationResult();
void invalidateOperationResult();
void validateApplyResult();
void invalidateApplyResult();


LEAF_TYPE invertLeaf(LEAF_TYPE a);
LEAF_TYPE addLeaf(LEAF_TYPE a, LEAF_TYPE b);
LEAF_TYPE subLeaf(LEAF_TYPE a, LEAF_TYPE b);
LEAF_TYPE mulLeaf(LEAF_TYPE a, LEAF_TYPE b);
LEAF_TYPE divLeaf(LEAF_TYPE a, LEAF_TYPE b);
LEAF_TYPE sqrtLeaf(LEAF_TYPE a);
LEAF_TYPE invsqrtLeaf(LEAF_TYPE a);
LEAF_TYPE rotateCoef1(LEAF_TYPE a);
LEAF_TYPE rotateCoef2(LEAF_TYPE a);
LEAF_TYPE times2Leaf(LEAF_TYPE l);
void q_printdot(qBDD a);
void q_fprintdot(FILE* out, qBDD a);

qBDD initCircuit(unsigned varNum); 
void circuit_init_interface(qBDD *c, const uint32_t n);
void initPackage(unsigned cacheSize, unsigned nodeSize, unsigned varNum);
void deleteCircuit(qBDD *circ);
void freePackage();

void init_terminal_symb_val_i();
void init_terminal_symb_map_i();

void forceGC();

size_t qBDD_getTerminalType(qBDD terminal);

LEAF_TYPE mtbdd_symb_neg_i(LEAF_TYPE t);
LEAF_TYPE mtbdd_symb_coef_rot1_i(LEAF_TYPE t);
LEAF_TYPE mtbdd_symb_coef_rot2_i(LEAF_TYPE t);
LEAF_TYPE mtbdd_symb_plus_i(LEAF_TYPE a, LEAF_TYPE b);
LEAF_TYPE mtbdd_symb_minus_i(LEAF_TYPE a, LEAF_TYPE b);
#ifdef __cplusplus
}
#endif

#endif
