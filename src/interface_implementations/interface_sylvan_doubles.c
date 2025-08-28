#include <sylvan.h>
#include "../gates.h"
#include "sylvan_int.h"
#include "../error.h"
#include "../hash.h"
#include "../mtbdd_out.h"
#include <math.h>
#define MAX_LEAF_STR_LEN 250
#define QBDD_TYPE_DEFINED
typedef MTBDD qBDD;
#define GET_MIN(a, b) ((a) < (b))? (a) : (b)

#include "../interface.h"
#include <string.h>


extern uint32_t ltype_id; // leafType for sylvan
typedef mpz_t coef_t;
struct LEAF_TYPE_IMPL {
    double real;
    double imag;
};

unsigned int globalInvSqrtCoeff; // coefficient denoting leaf*(1/sqrt(2))^coeff


static inline void allocPimpl(LEAF_TYPE *leaf) {
    leaf->pImpl = malloc(sizeof(LEAF_TYPE_IMPL));
    if (leaf->pImpl) {
        ((LEAF_TYPE_IMPL *)leaf->pImpl)->real = 0.0;
        ((LEAF_TYPE_IMPL *)leaf->pImpl)->imag = 0.0;
    }
}


void incInvSqrtCoeff() {
    globalInvSqrtCoeff++;
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
        qBDD returned = my_mtbdd_apply_gate(operand, TASK(mtbdd_operation_transformed), targets[0]);
        return returned;
    } else {
        return my_mtbdd_apply_cgate(operand, TASK(mtbdd_operation_transformed), targets[0], targets[1]);
    }
}

qBDD binary_apply(qBDD l, qBDD r, LEAF_TYPE(*op)(LEAF_TYPE, LEAF_TYPE)) {
    if (op == addLeaf) {
        return mtbdd_apply(l, r, TASK(mtbdd_binary_apply_add_op));
    } else if (op == subLeaf) {
        return mtbdd_apply(l, r, TASK(mtbdd_binary_apply_sub_op));
    } else {printf("here"); fflush(stdout);return mtbdd_false;}
}

qBDD unary_apply(qBDD l, LEAF_TYPE(*op)(LEAF_TYPE)) {
    global_unary_op_to_transform = op;
    return mtbdd_uapply(l, TASK(mtbdd_unary_apply_op_transformed), 0);
}

LEAF_TYPE addLeaf(LEAF_TYPE a, LEAF_TYPE b) {
    if (a.pImpl == NULL) return b;
    if (b.pImpl == NULL) return a;

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

    
    if (isnan(res_data->real) || isnan(res_data->imag)) {
        printf("Detected NaN!\n");
        exit(1);
    }
    if (isinf(res_data->real) || isinf(res_data->imag)) {
        printf("Detected Inf!\n");
        exit(1);
    }

    // If result is zero in both real and imag, free and return empty
    if (res_data->real == 0.0 && res_data->imag == 0.0) {
        free(res.pImpl);
        LEAF_TYPE zero_leaf = { .pImpl = NULL };
        return zero_leaf;
    }

    return res;
}

LEAF_TYPE invertLeaf(LEAF_TYPE a) {
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
}

LEAF_TYPE rotateCoef1(LEAF_TYPE l) {
    if (l.pImpl == NULL) return l;  // Return unchanged if input is null

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
    if (l.pImpl == NULL) return l;

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

    if (isnan(res_data->real) || isnan(res_data->imag)) {
        printf("Detected NaN!\n");
        exit(1);
    }
    if (isinf(res_data->real) || isinf(res_data->imag)) {
        printf("Detected Inf!\n");
        exit(1);
    }

    return result;
}
LEAF_TYPE subLeaf(LEAF_TYPE a, LEAF_TYPE b) {
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

    
    if (isnan(res_data->real) || isnan(res_data->imag)) {
        printf("Detected NaN!\n");
        exit(1);
    }
    if (isinf(res_data->real) || isinf(res_data->imag)) {
        printf("Detected Inf!\n");
        exit(1);
    }

    // Check if result is zero (both real and imag are zero)
    if (res_data->real == 0.0 && res_data->imag == 0.0) {
        free(res.pImpl);
        LEAF_TYPE zero = { .pImpl = NULL };
        return zero;
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
    free(data_p->pImpl);
    free(data_p);
}

int my_leaf_equals_interface(const uint64_t ldata_a_raw, const uint64_t ldata_b_raw)
{
    const double EPSILON = 1e-9;

    LEAF_TYPE *ldata_a = (LEAF_TYPE *)ldata_a_raw;
    LEAF_TYPE *ldata_b = (LEAF_TYPE *)ldata_b_raw;

    // Handle null pointers
    if (ldata_a == NULL || ldata_a->pImpl == NULL) {
        return (ldata_b == NULL || ldata_b->pImpl == NULL);
    }
    if (ldata_b == NULL || ldata_b->pImpl == NULL) {
        return 0;
    }

    LEAF_TYPE_IMPL *a_data = (LEAF_TYPE_IMPL *)ldata_a->pImpl;
    LEAF_TYPE_IMPL *b_data = (LEAF_TYPE_IMPL *)ldata_b->pImpl;

    // Compare real and imag parts with tolerance
    return (fabs(a_data->real - b_data->real) < EPSILON) &&
           (fabs(a_data->imag - b_data->imag) < EPSILON);
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
    new_impl->imag = orig_leaf->pImpl->imag;
    new_impl->real = orig_leaf->pImpl->real;

    // Step 3: Allocate a new LEAF_TYPE and attach the new pImpl
    LEAF_TYPE *new_leaf = (LEAF_TYPE *)my_malloc(sizeof(LEAF_TYPE));
    new_leaf->pImpl = new_impl;

    // Step 4: Store the pointer as the MTBDD leaf value
    *ldata_p_raw = (uint64_t)new_leaf;
}
static inline uint64_t hash_combine(uint64_t seed, uint64_t value) {
    return seed ^ (value + 0x9e3779b9 + (seed << 6) + (seed >> 2));
}
uint64_t my_leaf_hash_interface(const uint64_t ldata_raw, const uint64_t seed)
{
    double eps = 1e-9;
    LEAF_TYPE *ldata = (LEAF_TYPE*)ldata_raw;
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
}
#define VAR_NAME_FMT "large-number[%ld]"
#define MAX_NUM_LEN 50
static int leaf_to_str_output(char *buf, LEAF_TYPE_IMPL *ldata, double k, int shift_cnt) {
    // Compose string representation:
    // Format: "(1/√2)^(k) * (real + imag i)"
    // shift_cnt can be appended as a power or just noted, here simplified.

    // For the sake of example, ignore lnum_map_add and assume k is just a double exponent.

    int chars_written = snprintf(
        buf, MAX_LEAF_STR_LEN,
        "(1/√2)^(%.2f) * (%.6f %+.6fi)",
        k, 
        ldata->real,
        ldata->imag
    );

    return chars_written;
}
char* my_leaf_to_str_interface(int complemented, uint64_t ldata_raw, char *sylvan_buf, size_t sylvan_bufsize) {
    (void) complemented;

    LEAF_TYPE *leaf = (LEAF_TYPE *)ldata_raw;  // First cast to LEAF_TYPE*
    if (!leaf || !leaf->pImpl) {
        // Empty leaf or invalid data
        snprintf(sylvan_buf, sylvan_bufsize, "(null leaf)");
        return sylvan_buf;
    }

    LEAF_TYPE_IMPL *ldata = (LEAF_TYPE_IMPL *)leaf->pImpl;  // Now get the internal data pointer
    char ldata_string[MAX_LEAF_STR_LEN] = {0};

    // We assume shift count and gcd logic is not needed for doubles,
    // so just format directly with the exponent globalInvSqrtCoeff.

    int chars_written = snprintf(
        ldata_string, MAX_LEAF_STR_LEN,
        "(1/√2)^(%d) * (%.6f %+ .6fi)",
        globalInvSqrtCoeff,
        ldata->real,
        ldata->imag
    );

    if (chars_written >= MAX_LEAF_STR_LEN) {
        fprintf(stderr, "Allocated string length for leaf value output has not been sufficient.\n");
        exit(EXIT_FAILURE);
    } else if (chars_written < 0) {
        fprintf(stderr, "An encoding error has occurred when producing leaf value output.\n");
        exit(EXIT_FAILURE);
    }

    if ((size_t)chars_written < sylvan_bufsize) {
        memcpy(sylvan_buf, ldata_string, chars_written * sizeof(char));
        sylvan_buf[chars_written] = '\0';
        return sylvan_buf;
    }

    char *new_buf = (char *)malloc((chars_written + 1) * sizeof(char));
    if (!new_buf) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
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
    point.imag = 0.0;
    point.real = 1.0;
    globalInvSqrtCoeff = 0;
    uint8_t point_symbol[n];
    memset(point_symbol, 0, n*sizeof(uint8_t));
    type.pImpl = &point;
    qBDD leaf  = mtbdd_makeleaf(ltype_id, (uint64_t) &type);
    *c = mtbdd_cube(variables, point_symbol, leaf);

    mtbdd_protect(c);

}

void q_fprintdot(FILE* out, qBDD a) {
    mtbdd_fprintdot(out, a);
}


void initPackage(unsigned cacheSize, unsigned nodeSize, unsigned varNum) {
    //feraiseexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    init_sylvan_interface();
    init_my_leaf_interface(false);
}

void deleteCircuit(qBDD *circ) {
    circuit_delete(circ);
}

void freePackage() {
    stop_sylvan();
}






