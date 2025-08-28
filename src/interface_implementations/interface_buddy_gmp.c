#include "bdd.h"
#include "mtbdd.h"
#define QBDD_TYPE_DEFINED
typedef BDD qBDD;
#include "../interface.h"
#include <math.h>
#include <gmp.h>
#include <stdint.h>
#include "../hash.h"
#include <string.h>
#include <assert.h>
#include "../symexp.h"
#include "../symexp_list.h"
#include "../mtbdd_symb_map.h"
#include "../mtbdd_symb_val.h"
#include "../mtbdd_out.h"
#include "../error.h"

#define MAX_LEAF_STR_LEN 250
/// Max. number of digits written in the .dot output file of a single number
#define MAX_NUM_LEN 50

/// Return min from two numbers
#define GET_MIN(a, b) ((a) < (b))? (a) : (b)



mtbdd_terminal_type lt_classic;
mtbdd_terminal_type lt_symb_map;
mtbdd_terminal_type lt_symb_val;

struct LEAF_TYPE_IMPL {
        mpz_t a;
        mpz_t b;
        mpz_t c;
        mpz_t d;
};

uint64_t combine_mpz_t_interface(mpz_t x) {
    uint64_t hash = 0;

    const mp_limb_t *limbs = mpz_limbs_read(x);
    size_t size = mpz_size(x);

    for (size_t i = 0; i < size; i++) {
        hash ^= limbs[i];
    }

    return hash;
}

mpz_t globalSquareRootCoeff;
mpz_t globalSquareRootCoeffSymb;

static inline mpz_t* leafA(const LEAF_TYPE* leaf) { return &leaf->pImpl->a; }
static inline mpz_t* leafB(const LEAF_TYPE* leaf) { return &leaf->pImpl->b; }
static inline mpz_t* leafC(const LEAF_TYPE* leaf) { return &leaf->pImpl->c; }
static inline mpz_t* leafD(const LEAF_TYPE* leaf) { return &leaf->pImpl->d; }

static inline void allocPimpl(LEAF_TYPE* result) {
    result->pImpl = malloc(sizeof(LEAF_TYPE_IMPL));
    if (!result->pImpl) exit(1); // TODO EXIT
    mpz_init(result->pImpl->a);
    mpz_init(result->pImpl->b);
    mpz_init(result->pImpl->c);
    mpz_init(result->pImpl->d);
}

void freePimpl(void* leafraw) {
    return;
    if (!leafraw) return;

    LEAF_TYPE *leaf = (LEAF_TYPE*) leafraw;

    if (!leaf->pImpl) return;

    mpz_clear(leaf->pImpl->a);
    mpz_clear(leaf->pImpl->b);
    mpz_clear(leaf->pImpl->c);
    mpz_clear(leaf->pImpl->d);

    free(leaf->pImpl);
    leaf->pImpl = NULL;
}

void freePackage(){
    bdd_done();
}

LEAF_TYPE (*applyOperationToConvert)(LEAF_TYPE, LEAF_TYPE);
LEAF_TYPE (*applyOperationToConvertUnary)(LEAF_TYPE);

void *convertedApplyOperation(void *a, void *b) {
    LEAF_TYPE *leaf_a = (LEAF_TYPE*)a;
    LEAF_TYPE *leaf_b = (LEAF_TYPE*)b;
    LEAF_TYPE leaf_a_val;
    LEAF_TYPE leaf_b_val;

    leaf_a_val.pImpl = NULL;
    leaf_b_val.pImpl = NULL;

    if(leaf_a == NULL) {
        leaf_a_val.pImpl = NULL;
    } else {
        leaf_a_val = *leaf_a;
    }
    if(leaf_b == NULL) {
        leaf_b_val.pImpl = NULL;
    } else {
        leaf_b_val = *leaf_b;
    }
    
    // apply operation on LEAF_TYPE values to get a new LEAF_TYPE (by value)
    LEAF_TYPE result = applyOperationToConvert(leaf_a_val, leaf_b_val);

    if (result.pImpl == NULL) {
        // no internal data, just copy pointer as NULL
        return NULL;
    }

        // allocate memory for the new LEAF_TYPE pointer to return
    LEAF_TYPE *result_ptr = (LEAF_TYPE*)malloc(sizeof(LEAF_TYPE));
    if (result_ptr == NULL) {
        fprintf(stderr, "ERROR: malloc failed in convertedApplyOperation\n");
        return NULL;
    }


    result_ptr->pImpl = result.pImpl;

    // allocate new LEAF_TYPE_IMPL for pImpl
    // LEAF_TYPE_IMPL *new_impl = (LEAF_TYPE_IMPL*)malloc(sizeof(LEAF_TYPE_IMPL));
    // if (new_impl == NULL) {
    //     fprintf(stderr, "ERROR: malloc failed for LEAF_TYPE_IMPL\n");
    //     free(result_ptr);
    //     return NULL;
    // }

    // // deep copy each mpz_t field
    // mpz_init_set(new_impl->a, ((LEAF_TYPE_IMPL*)result.pImpl)->a);
    // mpz_init_set(new_impl->b, ((LEAF_TYPE_IMPL*)result.pImpl)->b);
    // mpz_init_set(new_impl->c, ((LEAF_TYPE_IMPL*)result.pImpl)->c);
    // mpz_init_set(new_impl->d, ((LEAF_TYPE_IMPL*)result.pImpl)->d);

    // // assign the new_impl to result_ptr
    // result_ptr->pImpl = new_impl;

    return result_ptr;
}

void *convertedApplyOperationAdd(void *a, void *b) {
    return convertedApplyOperation(a, b);
}

void *convertedApplyOperationMask(void *a, void *b) {
    return convertedApplyOperation(a, b);
}

void *convertedApplyOperationSub(void *a, void *b) {
    return convertedApplyOperation(a, b);
}

void *convertedApplyOperationPlusI(void *a, void *b) {
    return convertedApplyOperation(a, b);
}

void *convertedApplyOperationMinusI(void *a, void *b) {
    return convertedApplyOperation(a, b);
}

void *convertedApplyOperationUnary(void *a){
    LEAF_TYPE first;
    if (a == NULL) {first.pImpl=NULL;}
    else {first = *(LEAF_TYPE*)a;}

    LEAF_TYPE result = applyOperationToConvertUnary(first);
    if (result.pImpl == NULL) {
        return NULL;
    }
    LEAF_TYPE *resultPointer = (LEAF_TYPE*)malloc(sizeof(LEAF_TYPE));
    if(resultPointer == NULL){
        printf("ERROR\n");
        return NULL;
        // some err
    }
 // Allocate new LEAF_TYPE_IMPL for pImpl
    // LEAF_TYPE_IMPL *new_impl = (LEAF_TYPE_IMPL*)malloc(sizeof(LEAF_TYPE_IMPL));
    // if (new_impl == NULL) {
    //     fprintf(stderr, "ERROR: malloc failed for LEAF_TYPE_IMPL\n");
    //     free(resultPointer);
    //     return NULL;
    // }

    // // Deep copy each mpz_t field
    // mpz_init_set(new_impl->a, ((LEAF_TYPE_IMPL*)result.pImpl)->a);
    // mpz_init_set(new_impl->b, ((LEAF_TYPE_IMPL*)result.pImpl)->b);
    // mpz_init_set(new_impl->c, ((LEAF_TYPE_IMPL*)result.pImpl)->c);
    // mpz_init_set(new_impl->d, ((LEAF_TYPE_IMPL*)result.pImpl)->d);

    // Assign the new_impl to result_ptr
    resultPointer->pImpl = result.pImpl;

    return resultPointer;
}

void *convertedApplyOperationInv(void *a) {
    return convertedApplyOperationUnary(a);
}

void *convertedApplyOperationRot1(void *a) {
    return convertedApplyOperationUnary(a);
}

void *convertedApplyOperationRot2(void *a) {
    return convertedApplyOperationUnary(a);
}

void *convertedApplyOperationTimes2(void *a) {
    return convertedApplyOperationUnary(a);
}

void *convertedApplyOperationNegI(void *a) {
    return convertedApplyOperationUnary(a);
}

void *convertedApplyOperationRot1I(void *a) {
    return convertedApplyOperationUnary(a);
}

void *convertedApplyOperationRot2I(void *a) {
    return convertedApplyOperationUnary(a);
}

qBDD qBDD_protect(qBDD toProtect) {
    return bdd_addref(toProtect);    
}

qBDD qBDD_unprotect(qBDD toUnprotect) {
    return bdd_delref(toUnprotect);
}

int qBDD_isFalse(qBDD toCheck) {
    return toCheck == bdd_false();
}

qBDD qBDD_false() {
    return bdd_false();
}

qBDD qBDD_true() {
    return bdd_true();
}

int qBDD_isTerminal(qBDD toCheck) {
    //printf("a%d\n", toCheck);
    fflush(stdout);
    return ISCONST(toCheck) || ISTERMINAL(toCheck);
}

int qBDD_isInternal(qBDD toCheck) {
    return !ISTERMINAL(toCheck) && !ISCONST(toCheck);
}

qBDD qBDD_getHigh(qBDD a) {
    return HIGH(a);
}

qBDD qBDD_getLow(qBDD a) {
    return LOW(a);
}

size_t qBDD_getVar(qBDD a) {
    return bdd_var(a);
}

qBDD qBDD_maketerminal(size_t type, void* valuep) {
    return mtbdd_maketerminal(valuep, type);
}

size_t qBDD_symbolicMapLType() {
    return lt_symb_map;
}

size_t qBDD_symbolicValLType() {
    return lt_symb_val;
}

size_t qBDD_classicLType() {
    return lt_classic;
}

size_t qBDD_getTerminalType(qBDD terminal) {
    return mtbdd_get_terminal_type(terminal);
}

double qBDD_calculateProb(qBDD terminal) {
    mpf_t a, b, c, d;
    LEAF_TYPE_IMPL* prob = ((LEAF_TYPE*)mtbdd_getTerminalValue(terminal))->pImpl;
    mpf_inits(a, b, c, d, NULL);
    mpf_set_z(a, prob->a);
    mpf_set_z(b, prob->b);
    mpf_set_z(c, prob->c);
    mpf_set_z(d, prob->d);

    double prob_re, prob_im;
    double c_a, c_b, c_c, c_d;

    assert(mpz_fits_uint_p(globalSquareRootCoeff));
    mp_bitcnt_t shift_cnt = mpz_get_ui(globalSquareRootCoeff);

    shift_cnt = shift_cnt >> 1; // for k odd actually: shift_cnt = (shift_cnt - 1) >> 1, but the result is same

    // k even, k+1 odd
    if (mpz_even_p(globalSquareRootCoeff) != 0) {
        mpf_div_2exp(a, a, shift_cnt); // k/2 right shifts
        mpf_div_2exp(b, b, shift_cnt);
        mpf_div_2exp(c, c, shift_cnt);
        mpf_div_2exp(d, d, shift_cnt);

        c_a = mpf_get_d(a);
        c_b = mpf_get_d(b);
        c_c = mpf_get_d(c);
        c_d = mpf_get_d(d);

        prob_re = pow(c_a + c_b * M_SQRT1_2 - c_d * M_SQRT1_2, 2);
        prob_im = pow(c_c + c_b * M_SQRT1_2 + c_d * M_SQRT1_2, 2);

    }
    // k odd, k+1 even
    else {
        mpf_div_2exp(a, a, shift_cnt); // k-1/2 right shifts
        mpf_div_2exp(b, b, shift_cnt + 1); // k+1/2 right shifts
        mpf_div_2exp(c, c, shift_cnt);
        mpf_div_2exp(d, d, shift_cnt + 1);

        c_a = mpf_get_d(a);
        c_b = mpf_get_d(b);
        c_c = mpf_get_d(c);
        c_d = mpf_get_d(d);

        prob_re = pow(c_a * M_SQRT1_2 + c_b - c_d, 2);
        prob_im = pow(c_c * M_SQRT1_2 + c_b + c_d, 2);
    }

    mpf_clears(a, b, c, d, NULL);
    return prob_re+prob_im;
}

qBDD bdd_operation(qBDD operand, size_t* targets, size_t controlNum, qBDD(*op)(size_t tgt, qBDD, qBDD)){
    BDD returned = (mtbdd_operation(operand, targets, controlNum, op));
    return returned;
}

qBDD bdd_operation_guarded(qBDD operand, size_t* targets, size_t controlNum, qBDD(*op)(size_t tgt, qBDD)) {
    BDD returned = (mtbdd_operation_guarded(operand, targets, controlNum, op));
    return returned;
}

qBDD binary_apply(qBDD l, qBDD r, LEAF_TYPE(*op)(LEAF_TYPE, LEAF_TYPE)){
    applyOperationToConvert = op;
    if (op == addLeaf) {
        return (mtbdd_apply(l, r, convertedApplyOperationAdd));
    }
    else if (op == subLeaf){
        return (mtbdd_apply(l, r, convertedApplyOperationSub));
    }
    else if (op == mtbdd_symb_minus_i) {
        return (mtbdd_apply(l, r, convertedApplyOperationMinusI));
    }
    else if (op == mtbdd_symb_plus_i) {
        return (mtbdd_apply(l, r, convertedApplyOperationPlusI));
    }
    else {
        return (mtbdd_apply(l, r, convertedApplyOperation));
    }
}

qBDD binary_apply_guarded(qBDD l, qBDD r, qBDD(*op)(qBDD, qBDD)) {
    return (mtbdd_apply_guarded(l, r, op));
}

qBDD binary_apply_guarded_param(qBDD l, qBDD r, qBDD(*op)(qBDD, qBDD, size_t), size_t param) {
    return (mtbdd_apply_guarded_param(l,r,op,param));
}

qBDD unary_apply(qBDD l, LEAF_TYPE(*op)(LEAF_TYPE)){
    applyOperationToConvertUnary = op;
    if (op == rotateCoef1) {
        return (mtbdd_apply_unary(l,convertedApplyOperationRot1));
    } else if (op == rotateCoef2) {
        return (mtbdd_apply_unary(l,convertedApplyOperationRot2));
    } else if (op == times2Leaf) {
        return (mtbdd_apply_unary(l,convertedApplyOperationTimes2));
    } else if (op == invertLeaf){
        return (mtbdd_apply_unary(l, convertedApplyOperationInv));
    } else if (op == mtbdd_symb_neg_i) {
        return (mtbdd_apply_unary(l, convertedApplyOperationNegI));
    } else if (op == mtbdd_symb_coef_rot1_i) {
        return (mtbdd_apply_unary(l, convertedApplyOperationRot1I));
    } else if (op == mtbdd_symb_coef_rot2_i) {
        return (mtbdd_apply_unary(l, convertedApplyOperationRot2I));
    }
    return (mtbdd_apply_unary(l, convertedApplyOperationUnary));
}

int qBDD_leafcount(qBDD a) {
    return mtbdd_leaf_count(a);
}

qBDD unary_apply_guarded(qBDD l, qBDD(*op)(qBDD, size_t), size_t arg) {
    return (mtbdd_apply_unary_guarded(l, op, arg));
}

qBDD unary_apply_param(qBDD l, qBDD(*op)(qBDD, size_t), size_t arg) {
    return (mtbdd_apply_unary_param(l, op, arg));
}

qBDD newqBDD(unsigned int target, qBDD lhs, qBDD rhs){
    return bdd_makenode(target, lhs, rhs); 
}

LEAF_TYPE invertLeaf(LEAF_TYPE a){
    if (a.pImpl == NULL) return a;
    LEAF_TYPE result;
    allocPimpl(&result);

    mpz_t *aa = leafA(&a);
    mpz_t *ab = leafB(&a);
    mpz_t *ac = leafC(&a);
    mpz_t *ad = leafD(&a);
    mpz_neg(result.pImpl->a, *aa);
    mpz_neg(result.pImpl->b, *ab);
    mpz_neg(result.pImpl->c, *ac);
    mpz_neg(result.pImpl->d, *ad);
    return result;
}


LEAF_TYPE addLeaf(LEAF_TYPE a, LEAF_TYPE b){
    if (a.pImpl == NULL) return b; // zero
    if (b.pImpl == NULL) return a; // zero
    LEAF_TYPE result;
    allocPimpl(&result);

    mpz_t *aa = leafA(&a);
    mpz_t *ab = leafB(&a);
    mpz_t *ac = leafC(&a);
    mpz_t *ad = leafD(&a);

    mpz_t *ba = leafA(&b);
    mpz_t *bb = leafB(&b);
    mpz_t *bc = leafC(&b);
    mpz_t *bd = leafD(&b);

    mpz_add(result.pImpl->a, *aa, *ba);
    mpz_add(result.pImpl->b, *ab, *bb);
    mpz_add(result.pImpl->c, *ac, *bc);
    mpz_add(result.pImpl->d, *ad, *bd);

        // Return zero leaf if all values are zero
    if (!mpz_sgn(result.pImpl->a) && !mpz_sgn(result.pImpl->b) &&
        !mpz_sgn(result.pImpl->c) && !mpz_sgn(result.pImpl->d)) {
        mpz_clears(result.pImpl->a, result.pImpl->b, result.pImpl->c, result.pImpl->d, NULL);
        free(result.pImpl);
        LEAF_TYPE zero_leaf = { .pImpl = NULL };
        return zero_leaf;
    }

    return result;
}

LEAF_TYPE subLeaf(LEAF_TYPE a, LEAF_TYPE b){
    if (a.pImpl == NULL) return invertLeaf(b);
    if (b.pImpl == NULL) return a;
    LEAF_TYPE result;
    allocPimpl(&result);

    mpz_t *aa = leafA(&a);
    mpz_t *ab = leafB(&a);
    mpz_t *ac = leafC(&a);
    mpz_t *ad = leafD(&a);

    mpz_t *ba = leafA(&b);
    mpz_t *bb = leafB(&b);
    mpz_t *bc = leafC(&b);
    mpz_t *bd = leafD(&b);

    mpz_sub(result.pImpl->a, *aa, *ba);
    mpz_sub(result.pImpl->b, *ab, *bb);
    mpz_sub(result.pImpl->c, *ac, *bc);
    mpz_sub(result.pImpl->d, *ad, *bd);

    // Return zero leaf if all values are zero
    if (!mpz_sgn(result.pImpl->a) && !mpz_sgn(result.pImpl->b) &&
        !mpz_sgn(result.pImpl->c) && !mpz_sgn(result.pImpl->d)) {
        mpz_clears(result.pImpl->a, result.pImpl->b, result.pImpl->c, result.pImpl->d, NULL);
        free(result.pImpl);
        LEAF_TYPE zero_leaf = { .pImpl = NULL };
        return zero_leaf;
    }

    return result;
}

LEAF_TYPE mulLeaf(LEAF_TYPE a, LEAF_TYPE b){
    LEAF_TYPE result;
    allocPimpl(&result);

    mpz_t *aa = leafA(&a);
    mpz_t *ab = leafB(&a);
    mpz_t *ac = leafC(&a);
    mpz_t *ad = leafD(&a);

    mpz_t *ba = leafA(&b);
    mpz_t *bb = leafB(&b);
    mpz_t *bc = leafC(&b);
    mpz_t *bd = leafD(&b);

    mpz_mul(result.pImpl->a, *aa, *ba);
    mpz_mul(result.pImpl->b, *ab, *bb);
    mpz_mul(result.pImpl->c, *ac, *bc);
    mpz_mul(result.pImpl->d, *ad, *bd);

    return result;
}
LEAF_TYPE divLeaf(LEAF_TYPE a, LEAF_TYPE b){
    LEAF_TYPE result;
    allocPimpl(&result);

    mpz_t *aa = leafA(&a);
    mpz_t *ab = leafB(&a);
    mpz_t *ac = leafC(&a);
    mpz_t *ad = leafD(&a);

    mpz_t *ba = leafA(&b);
    mpz_t *bb = leafB(&b);
    mpz_t *bc = leafC(&b);
    mpz_t *bd = leafD(&b);

    mpz_div(result.pImpl->a, *aa, *ba);
    mpz_div(result.pImpl->b, *ab, *bb);
    mpz_div(result.pImpl->c, *ac, *bc);
    mpz_div(result.pImpl->d, *ad, *bd);

    return result;
}

LEAF_TYPE sqrtLeaf(LEAF_TYPE a){
    // TODO
    return a;
}

LEAF_TYPE rotateCoef1(LEAF_TYPE l) {
    if (l.pImpl == NULL) return l;
    LEAF_TYPE result;
    allocPimpl(&result);

    mpz_t *a = leafA(&l);
    mpz_t *b = leafB(&l);
    mpz_t *c = leafC(&l);
    mpz_t *d = leafD(&l);

    mpz_init(result.pImpl->a);
    mpz_neg(result.pImpl->a, *d);               // a' = -d
    mpz_init_set(result.pImpl->b, *a);          // b' = a
    mpz_init_set(result.pImpl->c, *b);          // c' = b
    mpz_init_set(result.pImpl->d, *c);          // d' = c

    return result;
}

LEAF_TYPE rotateCoef2(LEAF_TYPE l) {
    if (l.pImpl == NULL) return l;
    LEAF_TYPE result;
    allocPimpl(&result);

    mpz_t *a = leafA(&l);
    mpz_t *b = leafB(&l);
    mpz_t *c = leafC(&l);
    mpz_t *d = leafD(&l);

    mpz_init(result.pImpl->a);
    mpz_neg(result.pImpl->a, *c);
    mpz_init(result.pImpl->b);
    mpz_neg(result.pImpl->b, *d);
    mpz_init_set(result.pImpl->c, *a);
    mpz_init_set(result.pImpl->d, *b);

    return result;
}
LEAF_TYPE times2Leaf(LEAF_TYPE l) {
    if (l.pImpl == NULL) {return l;}
    unsigned long cc = (unsigned long) 2;
    LEAF_TYPE result;
    allocPimpl(&result);

    mpz_t *a = leafA(&l);
    mpz_t *b = leafB(&l);
    mpz_t *c = leafC(&l);
    mpz_t *d = leafD(&l);

    mpz_mul_ui(result.pImpl->a, *a, cc);
    mpz_mul_ui(result.pImpl->b, *b, cc);
    mpz_mul_ui(result.pImpl->c, *c, cc);
    mpz_mul_ui(result.pImpl->d, *d, cc);

    return result;
}

void incInvSqrtCoeff() {
    mpz_add_ui(globalSquareRootCoeff, globalSquareRootCoeff, 1);
}

void incInvSqrtCoeffSymb() {
    mpz_add_ui(globalSquareRootCoeffSymb, globalSquareRootCoeffSymb, 1);
}

void initInvSqrtCoeffSymb() {
    mpz_init(globalSquareRootCoeffSymb);
}

void addInvSqrtCoeffs() {
    mpz_add(globalSquareRootCoeff, globalSquareRootCoeff, globalSquareRootCoeffSymb);
}

void mulInvSqrtCoeff(mpz_t res, mpz_t a, unsigned long mul) {
    mpz_mul_ui(res, a, mul);
}

void clearInvSqrtCoeffSymb() {
    mpz_clear(globalSquareRootCoeffSymb);
}

void clearInvSqrtCoeffNormal() {
    mpz_clear(globalSquareRootCoeff);
}

void resetInvSqrtCoeffSymb() {
    mpz_set_ui(globalSquareRootCoeffSymb, 0);
}

mpz_srcptr getInvSqrtCoeff(void) {
    return globalSquareRootCoeff;
}
mpz_srcptr getInvSqrtCoeffSymb(void) {
    return globalSquareRootCoeffSymb;
}

void mulInvSqrtSymbCoeff(unsigned long a) {
    mpz_mul_ui(globalSquareRootCoeffSymb, globalSquareRootCoeffSymb, a);
}

void validateOperationResult() {
    FLAG_VALID_OPERATION();
}

void invalidateOperationResult() {
    FLAG_INVALID_OPERATION();
}

void validateApplyResult() {
    FLAG_VALID_APPLY();
}

void invalidateApplyResult() {
    FLAG_INVALID_APPLY();
}

CUSTOM_COMPARE_DECLARE(myCompare);
CUSTOM_COMPARE_DEFINE_START(myCompare, a, b)

    LEAF_TYPE *ldata_a = (LEAF_TYPE *) a;
    LEAF_TYPE *ldata_b = (LEAF_TYPE *) b;

    if (ldata_a == NULL || ldata_a->pImpl == NULL) {
        return 0;
    }
    if (ldata_b == NULL || ldata_b->pImpl == NULL) {
        return 0;
    }

    return !mpz_cmp(*leafA(ldata_a), *leafA(ldata_b)) &&
        !mpz_cmp(*leafB(ldata_a), *leafB(ldata_b)) &&
        !mpz_cmp(*leafC(ldata_a), *leafC(ldata_b)) &&
        !mpz_cmp(*leafD(ldata_a), *leafD(ldata_b));
CUSTOM_COMPARE_DEFINE_END

CUSTOM_HASH_DECLARE(myHash);
CUSTOM_HASH_DEFINE_START(myHash, q)
    LEAF_TYPE *ldata = (LEAF_TYPE*) q;
//    if(ldata == NULL) printf("daldkalfowefiuber *******\n");
    uint64_t val = 1;
    val = MY_HASH_COMB_GMP(val, *leafD(ldata));
    val = MY_HASH_COMB_GMP(val, *leafC(ldata));
    val = MY_HASH_COMB_GMP(val, *leafB(ldata));
    val = MY_HASH_COMB_GMP(val, *leafA(ldata));

    return (unsigned)(val ^ (val >> 32));
CUSTOM_COMPARE_DEFINE_END

static int _leaf_to_str_output_i(char *buf, LEAF_TYPE_IMPL *ldata, mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t k, mp_bitcnt_t shift_cnt)
{
    char buf_a[MAX_NUM_LEN + 2] = {0}; // +2 for sign and '\0'
    char buf_b[MAX_NUM_LEN + 2] = {0};
    char buf_c[MAX_NUM_LEN + 2] = {0};
    char buf_d[MAX_NUM_LEN + 2] = {0};
    int chars_written;
    
    if (mpz_sizeinbase(a, 10) > MAX_NUM_LEN) {
        chars_written = snprintf(buf_a, MAX_NUM_LEN + 2, VAR_NAME_FMT, lnum_map_add(&(ldata->a), shift_cnt));
        // variable length shouldn't exceed the max length but check to make sure
        assert(chars_written < MAX_NUM_LEN + 2 && chars_written >= 0);
    }
    else {
        // will always fit
        gmp_snprintf(buf_a, MAX_NUM_LEN + 2, "%Zd", a);
    }

    if (mpz_sizeinbase(b, 10) > MAX_NUM_LEN) {
        chars_written = snprintf(buf_b, MAX_NUM_LEN + 2, "+"VAR_NAME_FMT, lnum_map_add(&(ldata->b), shift_cnt));
        assert(chars_written < MAX_NUM_LEN + 2 && chars_written >= 0);
    }
    else {
        gmp_snprintf(buf_b, MAX_NUM_LEN + 2, "%+Zd", b);
    }

    if (mpz_sizeinbase(c, 10) > MAX_NUM_LEN) {
        chars_written = snprintf(buf_c, MAX_NUM_LEN + 2, "+"VAR_NAME_FMT, lnum_map_add(&(ldata->c), shift_cnt));
        assert(chars_written < MAX_NUM_LEN + 2 && chars_written >= 0);
    }
    else {
        gmp_snprintf(buf_c, MAX_NUM_LEN + 2, "%+Zd", c);
    }

    if (mpz_sizeinbase(d, 10) > MAX_NUM_LEN) {
        chars_written = snprintf(buf_d, MAX_NUM_LEN + 2, "+"VAR_NAME_FMT, lnum_map_add(&(ldata->d), shift_cnt));
        assert(chars_written < MAX_NUM_LEN + 2 && chars_written >= 0);
    }
    else {
        gmp_snprintf(buf_d, MAX_NUM_LEN + 2, "%+Zd", d);
    }

    chars_written = gmp_snprintf(buf, MAX_LEAF_STR_LEN, "(1/√2)^(%Zd) * (%s%sω%sω²%sω³)", k, buf_a, buf_b, buf_c, buf_d);
    return chars_written;
}


char* my_leaf_to_str_i(void* ldata_raw, char *buddy_buf, size_t buddy_bufsize)
{
    LEAF_TYPE *leafValP = (LEAF_TYPE*) ldata_raw;
    LEAF_TYPE_IMPL *ldata = leafValP->pImpl;
    char ldata_string[MAX_LEAF_STR_LEN] = {0};
    int chars_written;
    // Result in leaf will be divided by GCD for better clarity
    if (mpz_even_p(ldata->a) && mpz_even_p(ldata->b) && mpz_even_p(ldata->c) && mpz_even_p(ldata->d)) {
        mp_bitcnt_t greatest_pow2_a = mpz_scan1(ldata->a, 0);
        mp_bitcnt_t greatest_pow2_b = mpz_scan1(ldata->b, 0);
        mp_bitcnt_t greatest_pow2_c = mpz_scan1(ldata->c, 0);
        mp_bitcnt_t greatest_pow2_d = mpz_scan1(ldata->d, 0);

        mp_bitcnt_t greatest_pow2 = GET_MIN(GET_MIN(GET_MIN(greatest_pow2_a, greatest_pow2_b), greatest_pow2_c), greatest_pow2_d);

        // Get the power of 2, by which the result will be divided
        // Will be k/2 for even k and (k-1)/2 for odd k
        assert(mpz_fits_uint_p(globalSquareRootCoeff));
        mp_bitcnt_t pow2 = mpz_get_ui(globalSquareRootCoeff) >> 1;

        mp_bitcnt_t shift_cnt = GET_MIN(greatest_pow2, pow2);

        mpz_t a,b,c,d,k;
        mpz_inits(a, b, c, d, k, NULL);
        mpz_fdiv_q_2exp(a, ldata->a, shift_cnt);
        mpz_fdiv_q_2exp(b, ldata->b, shift_cnt);
        mpz_fdiv_q_2exp(c, ldata->c, shift_cnt);
        mpz_fdiv_q_2exp(d, ldata->d, shift_cnt);
        unsigned long k_decrement = shift_cnt << 1; // need to subtract 2*shift_cnt from c_k
        mpz_sub_ui(k, globalSquareRootCoeff, k_decrement);

        chars_written = _leaf_to_str_output_i(ldata_string, ldata, a, b, c, d,k, shift_cnt);
        mpz_clears(a, b, c, d, k, NULL);
    }
    else {
        chars_written = _leaf_to_str_output_i(ldata_string, ldata, ldata->a, ldata->b, ldata->c, ldata->d, globalSquareRootCoeff, 0);
    }

    // Was string truncated?
    if (chars_written >= MAX_LEAF_STR_LEN) {
        error_exit("Allocated string length for leaf value output has not been sufficient.\n");
    }
    else if (chars_written < 0) {
        error_exit("An encoding error has occured when producing leaf value output.\n");
    }

    // Is buffer large enough?
    if (chars_written < buddy_bufsize) {
        memcpy(buddy_buf, ldata_string, chars_written * sizeof(char));
        buddy_buf[chars_written] = '\0';
        return buddy_buf;
    }
    
    // Else return newly allocated string
    char *new_buf = (char*)my_malloc((chars_written + 1) * sizeof(char));
    memcpy(new_buf, ldata_string, chars_written * sizeof(char));
    new_buf[chars_written] = '\0';
    return new_buf;
}


void q_fprintdot(FILE* out, qBDD a) {
   buddy_mtbdd_fprintdot(out, a);
}


void initPackage(unsigned cacheSize, unsigned nodeSize, unsigned varNum){
    mtbdd = 1;
    bdd_init(10000, 10000);
    bdd_setvarnum(1);
    lt_classic = mtbdd_new_terminal_type();
    mtbdd_register_compare_function(lt_classic, myCompare);
    mtbdd_register_free_function(lt_classic, freePimpl);
    mtbdd_register_hash_function(lt_classic, myHash);
    mtbdd_register_to_str_function(lt_classic, my_leaf_to_str_i);
    SETDOMAIN(CUSTOM);
}


qBDD cube(int value, int width, qBDD *variables, qBDD leaf1, qBDD leaf0) {
    return mtbdd_cube2(0x0, width, variables, leaf1, leaf0);
}

void circuit_init_interface(qBDD *c, const uint32_t n)
{
    unsigned varNum = n;
    if(bdd_varnum() < varNum) {
        bdd_setvarnum(varNum);
    }
    
    mpz_init(globalSquareRootCoeff);  
    LEAF_TYPE* onePtr = malloc(sizeof(LEAF_TYPE));
    if(onePtr == NULL){bdd_error(BDD_MEMORY);}
    onePtr->pImpl = NULL;
    allocPimpl(onePtr);
    mpz_init_set_ui(*leafA(onePtr), 1);
    mpz_inits(*leafB(onePtr), *leafC(onePtr), *leafD(onePtr), globalSquareRootCoeff, NULL);
    
    BDD* variables = malloc(varNum * sizeof(BDD));
    if (variables == NULL) { bdd_error(BDD_MEMORY); }

    for (unsigned i = 0; i < varNum; i++) {
        variables[i] = bdd_ithvar(i); 
    }
    
    BDD leaf1 = mtbdd_maketerminal(onePtr, lt_classic);

    BDD cube_bdd = mtbdd_cube2(0x0, varNum, variables, leaf1, bdd_false());

    free(variables);
    *c = cube_bdd;
}

LEAF_TYPE qBDD_getTerminalValue(qBDD a) {
    return *(LEAF_TYPE*)mtbdd_getTerminalValue(a);
}

/** TODO */
void deleteCircuit(qBDD *circ) {
    return;
}

void forceGC() {
    bdd_gbc();
}

////////////////////////////////////////////////////
//           SYMB MAP //////////////

void terminal_symb_map_free(void *val) {
    return; // todo
}
CUSTOM_COMPARE_DECLARE(terminal_symb_map_compare)
CUSTOM_COMPARE_DEFINE_START(terminal_symb_map_compare, l_a, l_b)

    if (l_a == NULL && l_b == NULL) return 1;
    else if ((l_a == NULL && l_b != NULL) ||
            (l_a != NULL && l_b == NULL) ) {
                return 0;
            }

    LEAF_TYPE *leaf_a = (LEAF_TYPE*)l_a;
    LEAF_TYPE *leaf_b = (LEAF_TYPE*)l_b;

    sl_map_t *ldata_a = (sl_map_t*) leaf_a->pImpl;
    sl_map_t *ldata_b = (sl_map_t*) leaf_b->pImpl;

    return (ldata_a->va == ldata_b->va) && (ldata_a->vb == ldata_b->vb) && (ldata_a->vc == ldata_b->vc) \
           && (ldata_a->vd == ldata_b->vd);
CUSTOM_COMPARE_DEFINE_END

CUSTOM_HASH_DECLARE(terminal_symb_map_hash)
CUSTOM_HASH_DEFINE_START(terminal_symb_map_hash, l_a)
    LEAF_TYPE *leaf_a = (LEAF_TYPE*) l_a;
    sl_map_t *ldata = (sl_map_t*) leaf_a->pImpl;
    uint64_t val = 1;
    val = MY_HASH_COMB(val, ldata->va);
    val = MY_HASH_COMB(val, ldata->vb);
    val = MY_HASH_COMB(val, ldata->vc);
    val = MY_HASH_COMB(val, ldata->vd);

    return (unsigned)val;
CUSTOM_HASH_DEFINE_END

void init_terminal_symb_map_i() {
    lt_symb_map = mtbdd_new_terminal_type();
    mtbdd_register_compare_function(lt_symb_map, terminal_symb_map_compare);
    mtbdd_register_free_function(lt_symb_map, terminal_symb_map_free);
    mtbdd_register_hash_function(lt_symb_map, terminal_symb_map_hash);
}





///////////////////////////////////////////////////////////////////
//////////////////////// SYMB VAL////////////////////////////////

CUSTOM_COMPARE_DECLARE(terminal_symb_val_compare)
CUSTOM_COMPARE_DEFINE_START(terminal_symb_val_compare, l_a, l_b)
    LEAF_TYPE *leaf_a = (LEAF_TYPE*)l_a;
    LEAF_TYPE *leaf_b = (LEAF_TYPE*)l_b;

    sl_val_t *ldata_a = (sl_val_t *) leaf_a->pImpl;
    sl_val_t *ldata_b = (sl_val_t *) leaf_b->pImpl;

    return symexp_cmp(ldata_a->a, ldata_b->a) && symexp_cmp(ldata_a->b, ldata_b->b) && symexp_cmp(ldata_a->c, ldata_b->c) \
           && symexp_cmp(ldata_a->d, ldata_b->d);
CUSTOM_COMPARE_DEFINE_END

CUSTOM_HASH_DEFINE_START(terminal_symb_val_hash, l_a)
    LEAF_TYPE *leaf_a = (LEAF_TYPE*) l_a;
    sl_val_t *ldata = (sl_val_t*) leaf_a->pImpl;

    uint64_t val = 1;
    val = MY_HASH_COMB(val, ldata->a);
    val = MY_HASH_COMB(val, ldata->b);
    val = MY_HASH_COMB(val, ldata->c);
    val = MY_HASH_COMB(val, ldata->d);

    return val;
CUSTOM_HASH_DEFINE_END

void init_terminal_symb_val_i() {
    lt_symb_val = mtbdd_new_terminal_type();
    mtbdd_register_compare_function(lt_symb_val, terminal_symb_map_compare);
    mtbdd_register_free_function(lt_symb_val, NULL); // TODO ADD FREE
    mtbdd_register_hash_function(lt_symb_val, terminal_symb_map_hash);
}


LEAF_TYPE mtbdd_symb_neg_i(LEAF_TYPE t) {
    // Partial function check
    if (t.pImpl == NULL) return t;

    // Compute -a if mtbdd is a leaf
    sl_val_t *ldata = (sl_val_t*) t.pImpl;

    sl_val_t *res_data = (sl_val_t*)malloc(sizeof(sl_val_t));
    res_data->a = symexp_op(NULL, ldata->a, SYMEXP_SUB);
    res_data->b = symexp_op(NULL, ldata->b, SYMEXP_SUB);
    res_data->c = symexp_op(NULL, ldata->c, SYMEXP_SUB);
    res_data->d = symexp_op(NULL, ldata->d, SYMEXP_SUB);

    LEAF_TYPE res = {.pImpl = (LEAF_TYPE_IMPL*)res_data};
    return res;
}

LEAF_TYPE mtbdd_symb_coef_rot1_i(LEAF_TYPE t)
{
    // Partial function check
    if (t.pImpl == NULL) return t;

    // Compute coeficient rotation if mtbdd is a leaf
    sl_val_t *ldata = (sl_val_t*) t.pImpl;

    sl_val_t *res_data = (sl_val_t *)malloc(sizeof(sl_val_t));
    res_data->a = symexp_op(NULL, ldata->d, SYMEXP_SUB);
    res_data->b = ldata->a;
    res_data->c = ldata->b;
    res_data->d = ldata->c;

    LEAF_TYPE res = {.pImpl = (LEAF_TYPE_IMPL*)res_data};
    return res;
}

LEAF_TYPE mtbdd_symb_coef_rot2_i(LEAF_TYPE t)
{
    // Partial function check
    if (t.pImpl == NULL) return t;

    // Compute coeficient rotation if mtbdd is a leaf

    sl_val_t *ldata = (sl_val_t*) t.pImpl;

    sl_val_t *res_data = (sl_val_t*)malloc(sizeof(sl_val_t));
    res_data->a = symexp_op(NULL, ldata->c, SYMEXP_SUB);
    res_data->b = symexp_op(NULL, ldata->d, SYMEXP_SUB);
    res_data->c = ldata->a;
    res_data->d = ldata->b;

    LEAF_TYPE res = {.pImpl = (LEAF_TYPE_IMPL*)res_data};
    return res;
}

LEAF_TYPE mtbdd_symb_plus_i(LEAF_TYPE a, LEAF_TYPE b) {

    // Partial function check
    if (a.pImpl == NULL) return b;
    if (b.pImpl == NULL) return a;

    sl_val_t *a_data = (sl_val_t*) a.pImpl;
    sl_val_t *b_data = (sl_val_t*) b.pImpl;

    sl_val_t *res_data = (sl_val_t*) malloc(sizeof(sl_val_t));
    res_data->a = symexp_op(a_data->a, b_data->a, SYMEXP_ADD);
    res_data->b = symexp_op(a_data->b, b_data->b, SYMEXP_ADD);
    res_data->c = symexp_op(a_data->c, b_data->c, SYMEXP_ADD);
    res_data->d = symexp_op(a_data->d, b_data->d, SYMEXP_ADD);

    if (!res_data->a && !res_data->b && !res_data->c && !res_data->d) {
        free(res_data);
        LEAF_TYPE zero_leaf = { .pImpl = NULL };
        return zero_leaf;
    }
    
    LEAF_TYPE res = {.pImpl = (LEAF_TYPE_IMPL*) res_data};
    return res;
}

LEAF_TYPE mtbdd_symb_minus_i(LEAF_TYPE a, LEAF_TYPE b) {

    // Partial function check
    if (a.pImpl == NULL){
        LEAF_TYPE b_minus = mtbdd_symb_neg_i(b);
        return b_minus;
    } 
    if (b.pImpl == NULL) return a;

    // Compute a - b if both mtbdds are leaves
    sl_val_t *a_data = (sl_val_t*) a.pImpl;
    sl_val_t *b_data = (sl_val_t*) b.pImpl;
    
    sl_val_t *res_data = (sl_val_t *)malloc(sizeof(sl_val_t));
    res_data->a = symexp_op(a_data->a, b_data->a, SYMEXP_SUB);
    res_data->b = symexp_op(a_data->b, b_data->b, SYMEXP_SUB);
    res_data->c = symexp_op(a_data->c, b_data->c, SYMEXP_SUB);
    res_data->d = symexp_op(a_data->d, b_data->d, SYMEXP_SUB);

    LEAF_TYPE resleaf;

    if (!res_data->a && !res_data->b && !res_data->c && !res_data->d) {
        resleaf.pImpl = NULL;
        return resleaf;
    }
    
    resleaf.pImpl = res_data;
    return resleaf;
}