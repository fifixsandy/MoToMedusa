#include "bdd.h"
#include "mtbdd.h"
#define QBDD_TYPE_DEFINED
typedef BDD qBDD;

#include "../interface.h"

#include <math.h>
#include <stdint.h>
#include <string.h>


struct LEAF_TYPE_IMPL {
    double real;
    double imag;
};


unsigned int globalSquareRootCoeff;

static inline void allocPimpl(LEAF_TYPE *leaf) {
    leaf->pImpl = malloc(sizeof(LEAF_TYPE_IMPL));
    if (leaf->pImpl) {
        ((LEAF_TYPE_IMPL *)leaf->pImpl)->real = 0.0;
        ((LEAF_TYPE_IMPL *)leaf->pImpl)->imag = 0.0;
    }
}

void freePimpl(void *leafpointer) {
    LEAF_TYPE *leaf = (LEAF_TYPE*)leafpointer;
    if(leaf == NULL || leaf->pImpl == 0) return;
    free(leaf->pImpl);
    leaf->pImpl == NULL;
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
    
    // Apply operation on LEAF_TYPE values to get a new LEAF_TYPE (by value)
    LEAF_TYPE result = applyOperationToConvert(leaf_a_val, leaf_b_val);

    if (result.pImpl == NULL) {
        // No internal data, just copy pointer as NULL
        return NULL;
    }

        // Allocate memory for the new LEAF_TYPE pointer to return
    LEAF_TYPE *result_ptr = (LEAF_TYPE*)malloc(sizeof(LEAF_TYPE));
    if (result_ptr == NULL) {
        fprintf(stderr, "ERROR: malloc failed in convertedApplyOperation\n");
        return NULL;
    }

    // Allocate new LEAF_TYPE_IMPL for pImpl
    LEAF_TYPE_IMPL *new_impl = (LEAF_TYPE_IMPL*)malloc(sizeof(LEAF_TYPE_IMPL));
    if (new_impl == NULL) {
        fprintf(stderr, "ERROR: malloc failed for LEAF_TYPE_IMPL\n");
        free(result_ptr);
        return NULL;
    }


    new_impl->imag = result.pImpl->imag;
    new_impl->real = result.pImpl->real;
    // Assign the new_impl to result_ptr
    result_ptr->pImpl = new_impl;

    return result_ptr;
}

void *convertedApplyOperationAdd(void *a, void *b) {
    return convertedApplyOperation(a, b);
}

void *convertedApplyOperationSub(void *a, void *b) {
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
    LEAF_TYPE *resultPointer = (LEAF_TYPE*)malloc(sizeof(LEAF_TYPE)); // this malloc is wrong?
    if(resultPointer == NULL){
        printf("ERROR\n");
        return NULL;
        // some err
    }
    LEAF_TYPE_IMPL *new_impl = (LEAF_TYPE_IMPL*)malloc(sizeof(LEAF_TYPE_IMPL));
    if (new_impl == NULL) {
        fprintf(stderr, "ERROR: malloc failed for LEAF_TYPE_IMPL\n");
        free(resultPointer);
        return NULL;
    }


    new_impl->imag = result.pImpl->imag;
    new_impl->real = result.pImpl->real;
    // Assign the new_impl to result_ptr
    resultPointer->pImpl = new_impl;
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

qBDD bdd_operation(qBDD operand, size_t* targets, size_t controlNum, qBDD(*op)(size_t tgt, qBDD, qBDD)){
    BDD returned = bdd_addref(mtbdd_operation(operand, targets, controlNum, op));
    return returned;
}
qBDD binary_apply(qBDD l, qBDD r, LEAF_TYPE(*op)(LEAF_TYPE, LEAF_TYPE)){
    applyOperationToConvert = op;
    if (op == addLeaf) {
        return bdd_addref(mtbdd_apply(l, r, convertedApplyOperationAdd));
        }
    else if (op == subLeaf){
        return bdd_addref(mtbdd_apply(l, r, convertedApplyOperationSub));
    }
    else {
        return bdd_addref(mtbdd_apply(l, r, convertedApplyOperation));
    }
}

qBDD unary_apply(qBDD l, LEAF_TYPE(*op)(LEAF_TYPE)){
    applyOperationToConvertUnary = op;
    if (op == rotateCoef1) {
        return bdd_addref(mtbdd_apply_unary(l,convertedApplyOperationRot1));
    } else if (op == rotateCoef2) {
        return bdd_addref(mtbdd_apply_unary(l,convertedApplyOperationRot2));
    } else if (op == times2Leaf) {
        return bdd_addref(mtbdd_apply_unary(l,convertedApplyOperationTimes2));
    } else if (op == invertLeaf){
        return bdd_addref(mtbdd_apply_unary(l, convertedApplyOperationInv));
    }
    return bdd_addref(mtbdd_apply_unary(l, convertedApplyOperationUnary));
}

qBDD newqBDD(unsigned int target, qBDD lhs, qBDD rhs){
    return (bdd_makenode(target, lhs, rhs)); 
}

LEAF_TYPE invertLeaf(LEAF_TYPE a){
    if (a.pImpl == NULL) return a;
    LEAF_TYPE result;
    allocPimpl(&result);
    if (result.pImpl == NULL) {
        LEAF_TYPE empty = { .pImpl = NULL };
        return empty;
    }

    LEAF_TYPE_IMPL *a_data = (LEAF_TYPE_IMPL *) a.pImpl;
    LEAF_TYPE_IMPL *res_data = (LEAF_TYPE_IMPL *) result.pImpl;

    // Negate the real and imaginary parts
    res_data->real = -a_data->real;
    res_data->imag = -a_data->imag;

    return result;
};


LEAF_TYPE addLeaf(LEAF_TYPE a, LEAF_TYPE b){
    if (a.pImpl == NULL) return b; // zero
    if (b.pImpl == NULL) return a; // zero
    LEAF_TYPE_IMPL *a_data = (LEAF_TYPE_IMPL*) a.pImpl;
    LEAF_TYPE_IMPL *b_data = (LEAF_TYPE_IMPL*) b.pImpl;

    LEAF_TYPE res;
    allocPimpl(&res);
    if (res.pImpl == NULL) {
        // Allocation failed, return empty leaf
        LEAF_TYPE zero_leaf = { .pImpl = NULL };
        return zero_leaf;
    }

    LEAF_TYPE_IMPL *res_data = (LEAF_TYPE_IMPL *) res.pImpl;

    // Add real and imaginary parts
    res_data->real = a_data->real + b_data->real;
    res_data->imag = a_data->imag + b_data->imag;

    // If result is zero in both real and imag, free and return empty
    if (res_data->real == 0.0 && res_data->imag == 0.0) {
        free(res.pImpl);
        LEAF_TYPE zero_leaf = { .pImpl = NULL };
        return zero_leaf;
    }

    return res;
}

LEAF_TYPE subLeaf(LEAF_TYPE a, LEAF_TYPE b){
    if (a.pImpl == NULL) return invertLeaf(b);
    if (b.pImpl == NULL) return a;
    LEAF_TYPE_IMPL *a_data = (LEAF_TYPE_IMPL*) a.pImpl;
    LEAF_TYPE_IMPL *b_data = (LEAF_TYPE_IMPL*) b.pImpl;

    LEAF_TYPE res;
    allocPimpl(&res);
    if (res.pImpl == NULL) {
        LEAF_TYPE zero = { .pImpl = NULL };
        return zero;
    }

    LEAF_TYPE_IMPL *res_data = (LEAF_TYPE_IMPL*) res.pImpl;

    res_data->real = a_data->real - b_data->real;
    res_data->imag = a_data->imag - b_data->imag;

    // Check if result is zero (both real and imag are zero)
    if (res_data->real == 0.0 && res_data->imag == 0.0) {
        free(res.pImpl);
        LEAF_TYPE zero = { .pImpl = NULL };
        return zero;
    }

    return res;
}

LEAF_TYPE mulLeaf(LEAF_TYPE a, LEAF_TYPE b){
    // LEAF_TYPE result;
    // allocPimpl(&result);

    // mpz_t *aa = leafA(&a);
    // mpz_t *ab = leafB(&a);
    // mpz_t *ac = leafC(&a);
    // mpz_t *ad = leafD(&a);

    // mpz_t *ba = leafA(&b);
    // mpz_t *bb = leafB(&b);
    // mpz_t *bc = leafC(&b);
    // mpz_t *bd = leafD(&b);

    // mpz_mul(result.pImpl->a, *aa, *ba);
    // mpz_mul(result.pImpl->b, *ab, *bb);
    // mpz_mul(result.pImpl->c, *ac, *bc);
    // mpz_mul(result.pImpl->d, *ad, *bd);

    // return result;
}
LEAF_TYPE divLeaf(LEAF_TYPE a, LEAF_TYPE b){
    // LEAF_TYPE result;
    // allocPimpl(&result);

    // mpz_t *aa = leafA(&a);
    // mpz_t *ab = leafB(&a);
    // mpz_t *ac = leafC(&a);
    // mpz_t *ad = leafD(&a);

    // mpz_t *ba = leafA(&b);
    // mpz_t *bb = leafB(&b);
    // mpz_t *bc = leafC(&b);
    // mpz_t *bd = leafD(&b);

    // mpz_div(result.pImpl->a, *aa, *ba);
    // mpz_div(result.pImpl->b, *ab, *bb);
    // mpz_div(result.pImpl->c, *ac, *bc);
    // mpz_div(result.pImpl->d, *ad, *bd);

    //return result;
}

LEAF_TYPE sqrtLeaf(LEAF_TYPE a){
    // TODO
    return a;
}

LEAF_TYPE rotateCoef1(LEAF_TYPE l) {
    if (l.pImpl == NULL) return l;
    LEAF_TYPE result;
    allocPimpl(&result);
    if (result.pImpl == NULL) {
        LEAF_TYPE empty = { .pImpl = NULL };
        return empty;
    }

    LEAF_TYPE_IMPL *l_data = (LEAF_TYPE_IMPL *) l.pImpl;
    LEAF_TYPE_IMPL *res_data = (LEAF_TYPE_IMPL *) result.pImpl;

    double x = l_data->real;
    double y = l_data->imag;

    res_data->real = x - y;
    res_data->imag = x + y;

    return result;
}

LEAF_TYPE rotateCoef2(LEAF_TYPE l) {
    if (l.pImpl == NULL) return l;
    if (l.pImpl == NULL) return l;
    LEAF_TYPE result;
    allocPimpl(&result);
    if (result.pImpl == NULL) {
        LEAF_TYPE empty = { .pImpl = NULL };
        return empty;
    }

    LEAF_TYPE_IMPL *l_data = (LEAF_TYPE_IMPL *) l.pImpl;
    LEAF_TYPE_IMPL *res_data = (LEAF_TYPE_IMPL *) result.pImpl;

    double x = l_data->real;
    double y = l_data->imag;

    res_data->real = -y;
    res_data->imag = x;

    return result;
}
LEAF_TYPE times2Leaf(LEAF_TYPE l) {
    if (l.pImpl == NULL) {return l;}
    LEAF_TYPE result;
    allocPimpl(&result);
    if (result.pImpl == NULL) {
        LEAF_TYPE zero = { .pImpl = NULL };
        return zero;
    }

    LEAF_TYPE_IMPL *l_data = (LEAF_TYPE_IMPL*) l.pImpl;
    LEAF_TYPE_IMPL *res_data = (LEAF_TYPE_IMPL*) result.pImpl;

    res_data->real = 2.0 * l_data->real;
    res_data->imag = 2.0 * l_data->imag;

    return result;
}


void incInvSqrtCoeff() {
    globalSquareRootCoeff++;
}


CUSTOM_COMPARE_DECLARE(myCompare);
CUSTOM_COMPARE_DEFINE_START(myCompare, a, b)

    const double EPSILON = 1e-9;

    LEAF_TYPE *ldata_a = (LEAF_TYPE *)a;
    LEAF_TYPE *ldata_b = (LEAF_TYPE *)b;

    // Handle null pointers
    if (ldata_a == NULL || ldata_a->pImpl == NULL) {
        return 0;
    }
    if (ldata_b == NULL || ldata_b->pImpl == NULL) {
        return 0;
    }

    LEAF_TYPE_IMPL *a_data = (LEAF_TYPE_IMPL *)ldata_a->pImpl;
    LEAF_TYPE_IMPL *b_data = (LEAF_TYPE_IMPL *)ldata_b->pImpl;

    // Compare real and imag parts with tolerance
    return (fabs(a_data->real - b_data->real) < EPSILON) &&
           (fabs(a_data->imag - b_data->imag) < EPSILON);
CUSTOM_COMPARE_DEFINE_END

static inline uint64_t hash_combine(uint64_t seed, uint64_t value) {
    return seed ^ (value + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}
CUSTOM_HASH_DECLARE(myHash);
CUSTOM_HASH_DEFINE_START(myHash, q)
    uint64_t seed = 1; // todo?
    double eps = 1e-9;
    LEAF_TYPE *ldata = (LEAF_TYPE*)q;
    
    if (ldata == NULL || ldata->pImpl == NULL) return seed;
    LEAF_TYPE_IMPL *data = (LEAF_TYPE_IMPL*)ldata->pImpl;

    // Quantize real and imag parts based on eps

    double q_real = round(data->real / eps) * eps;
    double q_imag = round(data->imag / eps) * eps;

    uint64_t real_bits, imag_bits;
    memcpy(&real_bits, &q_real, sizeof(double));
    memcpy(&imag_bits, &q_imag, sizeof(double));

    uint64_t val = seed;
    val = hash_combine(val, real_bits);
    val = hash_combine(val, imag_bits);

    return val;
CUSTOM_COMPARE_DEFINE_END

void myPrint(qBDD a, FILE *f)
{

    LEAF_TYPE* b = (LEAF_TYPE*)mtbdd_getTerminalValue(a);
    if (!b || !b->pImpl) {
        fprintf(f, "Empty or invalid leaf\n");
        return;
    }

    fprintf(
        f,
        "(1/√2)^(%d) * (%.6f %+.6fi)",
        globalSquareRootCoeff,
        b->pImpl->real,
        b->pImpl->imag
    );
}

// void q_printdot(qBDD a){
//     bdd_printdot(a, myPrint);
// }

void q_fprintdot(FILE* o` {
    bdd_fprintdot(out, a, myPrint);
}

void initPackage(unsigned cacheSize, unsigned nodeSize, unsigned varNum){
    mtbdd = 1;
    bdd_init(10000, 10000);
    bdd_setvarnum(1);
    customCompare = myCompare;
    customHash = myHash;
    customFree = freePimpl;
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
    
    globalSquareRootCoeff = 0;
    LEAF_TYPE* onePtr = malloc(sizeof(LEAF_TYPE));
    if(onePtr == NULL){bdd_error(BDD_MEMORY);}
    onePtr->pImpl = NULL;
    allocPimpl(onePtr);
    onePtr->pImpl->imag = 0.0;
    onePtr->pImpl->real = 1.0;
    
    BDD* variables = malloc(varNum * sizeof(BDD));
    if (variables == NULL) { bdd_error(BDD_MEMORY); }

    for (unsigned i = 0; i < varNum; i++) {
        variables[i] = bdd_ithvar(i); 
    }
    
    
    BDD leaf1 = bdd_addref(mtbdd_maketerminal(onePtr));

    BDD cube_bdd = mtbdd_cube2(0x0, varNum, variables, leaf1, bdd_false());
    bdd_addref(cube_bdd);
    bdd_delref(leaf1);
    free(variables);
    *c = cube_bdd;
}

/** TODO */
void deleteCircuit(qBDD *circ) {
    bdd_delref(*circ);
    return;
}

