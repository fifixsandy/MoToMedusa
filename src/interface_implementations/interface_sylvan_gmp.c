#include <sylvan.h>
#include <gmp.h>
#include "../gates.h"
#include "sylvan_int.h"
#include "../error.h"
#include "../hash.h"
#include "../mtbdd_out.h"
#define MAX_LEAF_STR_LEN 250
#define QBDD_TYPE_DEFINED
typedef MTBDD qBDD;
#define GET_MIN(a, b) ((a) < (b))? (a) : (b)

#include "../interface.h"
#include <string.h>


extern uint32_t ltype_id; // leafType for sylvan
typedef mpz_t coef_t;
struct LEAF_TYPE_IMPL {
        /// a * 1
        coef_t a;
        /// b * ω
        coef_t b;
        /// c * ω²
        coef_t c;
        /// d * ω³
        coef_t d;
};

coef_t globalInvSqrtCoeff; // coefficient denoting leaf*(1/sqrt(2))^coeff

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


void incInvSqrtCoeff() {
    mpz_add_ui(globalInvSqrtCoeff,globalInvSqrtCoeff,1);
}


LEAF_TYPE (*global_binary_op_to_transform)(LEAF_TYPE, LEAF_TYPE);
LEAF_TYPE (*global_unary_op_to_transform)(LEAF_TYPE);
qBDD (*global_operation_to_transform)(size_t, qBDD, qBDD);



qBDD mtbdd_binary_apply_op_transformed(qBDD *p_a, qBDD *p_b) {
    qBDD a = *p_a;
    qBDD b = *p_b;

    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {

        LEAF_TYPE aVal;
        LEAF_TYPE bVal;
        
        if (a == mtbdd_false) {
            aVal.pImpl = 0;
        } else {
            aVal = *(LEAF_TYPE*)mtbdd_getvalue(a);
        }
        if (b == mtbdd_false) {
            bVal.pImpl = 0;
        } else {
            bVal = *(LEAF_TYPE*)mtbdd_getvalue(b);
        }
        
        LEAF_TYPE result = global_binary_op_to_transform(aVal, bVal);

        if (result.pImpl == 0) {
            return mtbdd_false;
        } else {
            MTBDD l = mtbdd_makeleaf(ltype_id, (uint64_t) &result);
            mtbdd_protect(&l);
            return l;
        }

    }

    // only for commutative
    // if (a < b) { 
    //     *p_a = b;
    //     *p_b = a;
    // }

    return mtbdd_invalid; // Recurse deeper
}

TASK_DECL_2(qBDD, mtbdd_binary_apply_add_op, qBDD*, qBDD*);
TASK_IMPL_2(qBDD, mtbdd_binary_apply_add_op, qBDD*, p_a, qBDD*, p_b) {
    global_binary_op_to_transform = addLeaf;
    return mtbdd_binary_apply_op_transformed(p_a, p_b);
}

TASK_DECL_2(qBDD, mtbdd_binary_apply_sub_op, qBDD*, qBDD*);
TASK_IMPL_2(qBDD, mtbdd_binary_apply_sub_op, qBDD*, p_a, qBDD*, p_b) {
    global_binary_op_to_transform = subLeaf;
    return mtbdd_binary_apply_op_transformed(p_a, p_b);
}

TASK_DECL_2(qBDD, mtbdd_unary_apply_op_transformed, qBDD, size_t)
TASK_IMPL_2(qBDD, mtbdd_unary_apply_op_transformed, qBDD, a, size_t, x)
{
    (void)x;

    if (mtbdd_isleaf(a)) {

        LEAF_TYPE aVal;

        if (a == mtbdd_false) {
            aVal.pImpl = 0;
        } else {
            aVal = *(LEAF_TYPE*)mtbdd_getvalue(a);
        }
        
        LEAF_TYPE result = global_unary_op_to_transform(aVal);
        if (result.pImpl == 0) {
            return mtbdd_false;
        } else {
            MTBDD l = mtbdd_makeleaf(ltype_id, (uint64_t) &result);
            return l;
        }
    }


    return mtbdd_invalid; // Recurse deeper
}

TASK_DECL_3(MTBDD, mtbdd_operation_transformed, MTBDD, MTBDD, uint32_t);
TASK_IMPL_3(MTBDD, mtbdd_operation_transformed, MTBDD, low, MTBDD, high, uint32_t, xt)
{
    return global_operation_to_transform(xt, low, high);
}

qBDD newqBDD(unsigned int target, qBDD low, qBDD high) {
    MTBDD new = mtbdd_makenode(target, low, high);
    mtbdd_protect(&new);
    return new;
}

qBDD bdd_operation(qBDD operand, size_t* targets, size_t controlNum, qBDD(*op)(size_t tgt, qBDD, qBDD)) {
    global_operation_to_transform = op;
    if (controlNum == 0) {
        return my_mtbdd_apply_gate(operand, TASK(mtbdd_operation_transformed), targets[0]);
    } else {
        return my_mtbdd_apply_cgate(operand, TASK(mtbdd_operation_transformed), targets[0], targets[1]);
    }
}

qBDD binary_apply(qBDD l, qBDD r, LEAF_TYPE(*op)(LEAF_TYPE, LEAF_TYPE)) {
    if (op == addLeaf) {
        return mtbdd_apply(l, r, TASK(mtbdd_binary_apply_add_op));
    } else if (op == subLeaf) {
        return mtbdd_apply(l, r, TASK(mtbdd_binary_apply_sub_op));
    } else {return mtbdd_false;}
}

qBDD unary_apply(qBDD l, LEAF_TYPE(*op)(LEAF_TYPE)) {
    global_unary_op_to_transform = op;
    return mtbdd_uapply(l, TASK(mtbdd_unary_apply_op_transformed), 0);
}

LEAF_TYPE addLeaf(LEAF_TYPE a, LEAF_TYPE b) {
    if (a.pImpl == 0) return b;
    if (b.pImpl == 0) return a;

    LEAF_TYPE_IMPL *a_data = (LEAF_TYPE_IMPL*) a.pImpl;
    LEAF_TYPE_IMPL *b_data = (LEAF_TYPE_IMPL*) b.pImpl;
    
    LEAF_TYPE res;
    // Allocate new result structures
    allocPimpl(&res);

    // Perform addition
    mpz_add(res.pImpl->a, a_data->a, b_data->a);
    mpz_add(res.pImpl->b, a_data->b, b_data->b);
    mpz_add(res.pImpl->c, a_data->c, b_data->c);
    mpz_add(res.pImpl->d, a_data->d, b_data->d);


    // If result is zero in all fields, clean up and return empty
    if (!mpz_sgn(res.pImpl->a) && !mpz_sgn(res.pImpl->b) &&
        !mpz_sgn(res.pImpl->c) && !mpz_sgn(res.pImpl->d)) {
        mpz_clears(res.pImpl->a, res.pImpl->b, res.pImpl->c, res.pImpl->d, NULL);
        free(res.pImpl);
        LEAF_TYPE zero_leaf = { .pImpl = NULL };
        return zero_leaf;
    }
    

    return res;
}

LEAF_TYPE invertLeaf(LEAF_TYPE a) {

    if (a.pImpl == 0) return a;

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

LEAF_TYPE rotateCoef1(LEAF_TYPE l) {
    if (l.pImpl == NULL) return l;  // Return unchanged if input is null

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

    if (l.pImpl == 0) return l; // l is zero, return 0

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
    if (l.pImpl == 0) return l;

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
LEAF_TYPE subLeaf(LEAF_TYPE a, LEAF_TYPE b) {
    if (a.pImpl == 0) return invertLeaf(b);
    if (b.pImpl == 0) return a;

    LEAF_TYPE_IMPL *a_data = (LEAF_TYPE_IMPL*) a.pImpl;
    LEAF_TYPE_IMPL *b_data = (LEAF_TYPE_IMPL*) b.pImpl;

    // Allocate result
    LEAF_TYPE res;
    allocPimpl(&res);

    // Perform subtraction
    mpz_sub(res.pImpl->a, a_data->a, b_data->a);
    mpz_sub(res.pImpl->b, a_data->b, b_data->b);
    mpz_sub(res.pImpl->c, a_data->c, b_data->c);
    mpz_sub(res.pImpl->d, a_data->d, b_data->d);


    // Return zero leaf if all values are zero
    if (!mpz_sgn(res.pImpl->a) && !mpz_sgn(res.pImpl->b) &&
        !mpz_sgn(res.pImpl->c) && !mpz_sgn(res.pImpl->d)) {
        mpz_clears(res.pImpl->a, res.pImpl->b, res.pImpl->c, res.pImpl->d, NULL);
        free(res.pImpl);
        LEAF_TYPE zero_leaf = { .pImpl = NULL };
        return zero_leaf;
    }


    return res;
}



void init_sylvan_interface() {
    lace_start(1, 0); // 1 thread, default task queue size
    sylvan_set_limits(2000LL*1024*1024, 3, 5); // Allocate 2GB
    sylvan_init_package();
    sylvan_init_mtbdd();
}



void my_leaf_destroy_interface(uint64_t ldata)
{
    if (!ldata) return;
    LEAF_TYPE *data_p = (LEAF_TYPE*) ldata; // Data in leaf = pointer to my data
    mpz_clears(*leafA(data_p),*leafB(data_p), *leafC(data_p), *leafD(data_p), NULL);
    free(data_p->pImpl);
    free(data_p);
}

int my_leaf_equals_interface(const uint64_t ldata_a_raw, const uint64_t ldata_b_raw)
{
    LEAF_TYPE *ldata_a = (LEAF_TYPE *) ldata_a_raw;
    LEAF_TYPE *ldata_b = (LEAF_TYPE *) ldata_b_raw;

    return !mpz_cmp(*leafA(ldata_a), *leafA(ldata_b)) &&
        !mpz_cmp(*leafB(ldata_a), *leafB(ldata_b)) &&
        !mpz_cmp(*leafC(ldata_a), *leafC(ldata_b)) &&
        !mpz_cmp(*leafD(ldata_a), *leafD(ldata_b));
}

void my_leaf_create_interface(uint64_t *ldata_p_raw)
{
    // Step 1: Interpret the raw MTBDD leaf value as a LEAF_TYPE*
    
    LEAF_TYPE *orig_leaf = (LEAF_TYPE *)(*ldata_p_raw);
    if (!orig_leaf || !orig_leaf->pImpl) {
        fprintf(stderr, "Error: invalid leaf or pImpl in my_leaf_create_interface\n");
        exit(EXIT_FAILURE); // Or handle gracefully
    }
    // Step 2: Copy the internal LEAF_TYPE_IMPL structure
    LEAF_TYPE_IMPL *new_impl = (LEAF_TYPE_IMPL *)my_malloc(sizeof(LEAF_TYPE_IMPL));
    mpz_init_set(new_impl->a, orig_leaf->pImpl->a);
    mpz_init_set(new_impl->b, orig_leaf->pImpl->b);
    mpz_init_set(new_impl->c, orig_leaf->pImpl->c);
    mpz_init_set(new_impl->d, orig_leaf->pImpl->d);

    // Step 3: Allocate a new LEAF_TYPE and attach the new pImpl
    LEAF_TYPE *new_leaf = (LEAF_TYPE *)my_malloc(sizeof(LEAF_TYPE));
    new_leaf->pImpl = new_impl;

    // Step 4: Store the pointer as the MTBDD leaf value
    *ldata_p_raw = (uint64_t)new_leaf;
}

uint64_t my_leaf_hash_interface(const uint64_t ldata_raw, const uint64_t seed)
{
    LEAF_TYPE *ldata = (LEAF_TYPE*) ldata_raw;

    uint64_t val = seed;
    
    val = MY_HASH_COMB_GMP(val, *leafD(ldata));
    val = MY_HASH_COMB_GMP(val, *leafC(ldata));
    val = MY_HASH_COMB_GMP(val, *leafB(ldata));
    val = MY_HASH_COMB_GMP(val, *leafA(ldata));

    return val;
}
#define VAR_NAME_FMT "large-number[%ld]"
#define MAX_NUM_LEN 50
static int _leaf_to_str_output(char *buf, LEAF_TYPE_IMPL *ldata, mpz_t a, mpz_t b, mpz_t c, mpz_t d, mpz_t k, mp_bitcnt_t shift_cnt)
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

char* my_leaf_to_str_interface(int complemented, uint64_t ldata_raw, char *sylvan_buf, size_t sylvan_bufsize)
{
    (void) complemented;
    LEAF_TYPE *leaf = (LEAF_TYPE*) ldata_raw;
    char ldata_string[MAX_LEAF_STR_LEN] = {0};
    int chars_written;
    // Result in leaf will be divided by GCD for better clarity
    if (mpz_even_p(*leafA(leaf)) && mpz_even_p(*leafB(leaf)) && mpz_even_p(*leafC(leaf)) && mpz_even_p(*leafD(leaf))) {
        mp_bitcnt_t greatest_pow2_a = mpz_scan1(*leafA(leaf), 0);
        mp_bitcnt_t greatest_pow2_b = mpz_scan1(*leafB(leaf), 0);
        mp_bitcnt_t greatest_pow2_c = mpz_scan1(*leafC(leaf), 0);
        mp_bitcnt_t greatest_pow2_d = mpz_scan1(*leafD(leaf), 0);

        mp_bitcnt_t greatest_pow2 = GET_MIN(GET_MIN(GET_MIN(greatest_pow2_a, greatest_pow2_b), greatest_pow2_c), greatest_pow2_d);

        // Get the power of 2, by which the result will be divided
        // Will be k/2 for even k and (k-1)/2 for odd k
        assert(mpz_fits_uint_p(globalInvSqrtCoeff));
        mp_bitcnt_t pow2 = mpz_get_ui(globalInvSqrtCoeff) >> 1;

        mp_bitcnt_t shift_cnt = GET_MIN(greatest_pow2, pow2);

        mpz_t a,b,c,d,k;
        mpz_inits(a, b, c, d, k, NULL);
        mpz_fdiv_q_2exp(a, *leafA(leaf), shift_cnt);
        mpz_fdiv_q_2exp(b, *leafB(leaf), shift_cnt);
        mpz_fdiv_q_2exp(c, *leafC(leaf), shift_cnt);
        mpz_fdiv_q_2exp(d, *leafD(leaf), shift_cnt);
        unsigned long k_decrement = shift_cnt << 1; // need to subtract 2*shift_cnt from globalInvSqrtCoeff
        mpz_sub_ui(k, globalInvSqrtCoeff, k_decrement);

        chars_written = _leaf_to_str_output(ldata_string, leaf->pImpl, a, b, c, d,k, shift_cnt);
        mpz_clears(a, b, c, d, k, NULL);
    }
    else {
        chars_written = _leaf_to_str_output(ldata_string, leaf->pImpl, *leafA(leaf), *leafB(leaf), *leafC(leaf), *leafD(leaf), globalInvSqrtCoeff, 0);
    }

    // Was string truncated?
    if (chars_written >= MAX_LEAF_STR_LEN) {
        error_exit("Allocated string length for leaf value output has not been sufficient.\n");
    }
    else if (chars_written < 0) {
        error_exit("An encoding error has occured when producing leaf value output.\n");
    }

    // Is buffer large enough?
    if (chars_written < sylvan_bufsize) {
        memcpy(sylvan_buf, ldata_string, chars_written * sizeof(char));
        sylvan_buf[chars_written] = '\0';
        return sylvan_buf;
    }
    
    // Else return newly allocated string
    char *new_buf = (char*)my_malloc((chars_written + 1) * sizeof(char));
    memcpy(new_buf, ldata_string, chars_written * sizeof(char));
    new_buf[chars_written] = '\0';
    return new_buf;
}

void init_my_leaf_interface(bool is_prob)
{
    ltype_id = sylvan_mt_create_type();
    mtbdd_apply_gate_id = cache_next_opid();
    mtbdd_apply_cgate_id = cache_next_opid();

    sylvan_mt_set_create(ltype_id, my_leaf_create_interface);
    sylvan_mt_set_destroy(ltype_id, my_leaf_destroy_interface);
    sylvan_mt_set_equals(ltype_id, my_leaf_equals_interface);
    if (is_prob) {
        sylvan_mt_set_to_str(ltype_id, my_leaf_to_str_prob);
    }
    else {
        sylvan_mt_set_to_str(ltype_id, my_leaf_to_str_interface);
    }
    sylvan_mt_set_hash(ltype_id, my_leaf_hash_interface);
}

void circuit_init_interface(qBDD *c, const uint32_t n)
{
    BDDSET variables = mtbdd_set_empty();
    for (uint32_t i = 0; i < n; i++) {
        variables = mtbdd_set_add(variables, i);
    }

    LEAF_TYPE_IMPL point; // can be local, mtbdd_makeleaf makes realloc
    LEAF_TYPE type;
    mpz_init_set_ui(point.a, 1);
    mpz_inits(point.b, point.c, point.d, globalInvSqrtCoeff, NULL);
    uint8_t point_symbol[n];
    memset(point_symbol, 0, n*sizeof(uint8_t));
    type.pImpl = &point;
    qBDD leaf  = mtbdd_makeleaf(ltype_id, (uint64_t) &type);
    *c = mtbdd_cube(variables, point_symbol, leaf);

    mpz_clears(point.a, point.b, point.c, point.d, NULL);
    mtbdd_protect(c);
}

void q_fprintdot(FILE* out, qBDD a) {
    mtbdd_fprintdot(out, a);
}


void initPackage(unsigned cacheSize, unsigned nodeSize, unsigned varNum) {
    init_sylvan_interface();
    init_my_leaf_interface(false);
}

void deleteCircuit(qBDD *circ) {
    circuit_delete(circ);
}

void freePackage() {
    stop_sylvan();
}






