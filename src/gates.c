#include "gates.h"
#include "sylvan_int.h" // Because of 'mtbddnode_t'
#include "interface.h"
prob_t measure(qBDD *a, uint32_t xt, char *curr_state, int n)
{
    uint32_t var_t;
    if (qBDD_isTerminal(*a)) {
        var_t = n;
    }
    else {
        var_t = qBDD_getVar(*a);
    }
    long skip_coef;
    if (var_t == 0) {
        skip_coef = 1;
    }
    else {
        skip_coef = get_coef(0, var_t, xt, curr_state);
    }
    prob_t prob = interface_prob_sum(*a, xt, curr_state, n);
    prob *= skip_coef;
    
    return prob;
}

qBDD interface_gate_x(size_t xt, qBDD low, qBDD high) {
    return newqBDD(xt, high, low);
}

qBDD interface_gate_y(size_t xt, qBDD low, qBDD high) {
    qBDD nhigh = unary_apply(high, invertLeaf);
    qBDD_protect(nhigh);
    qBDD updated = newqBDD(xt, nhigh, low);
    qBDD_protect(updated);
    qBDD_unprotect(nhigh);
    qBDD res = unary_apply(updated, rotateCoef2);
    qBDD_unprotect(updated);
    return res;
}

qBDD interface_gate_z(size_t xt, qBDD low, qBDD high) {
    qBDD nhigh = unary_apply(high, invertLeaf);
    qBDD_protect(nhigh);
    qBDD res = newqBDD(xt, low, nhigh);
    qBDD_unprotect(nhigh);
    return res;
}

qBDD interface_gate_s(size_t xt, qBDD low, qBDD high) {
    qBDD newH = unary_apply(high, rotateCoef2);
    qBDD_protect(newH);
    qBDD res = newqBDD(xt, low, newH);
    qBDD_unprotect(newH);
    return res;
}

qBDD interface_gate_t(size_t xt, qBDD low, qBDD high) {
    qBDD newH = unary_apply(high, rotateCoef1);
    qBDD_protect(newH);
    qBDD res = newqBDD(xt, low, newH);
    qBDD_unprotect(newH);
    return res;
}

qBDD interface_gate_h(size_t xt, qBDD low, qBDD high) {
    if (low == high) {
        qBDD nlow = unary_apply(low , times2Leaf);
        qBDD_protect(nlow);
        qBDD res = newqBDD(xt, nlow, qBDD_false());
        qBDD_unprotect(nlow);
        return res;
    } else {
        qBDD nh = qBDD_protect(binary_apply(low, high, subLeaf));
        qBDD nl = qBDD_protect(binary_apply(low, high, addLeaf));
        qBDD new = newqBDD(xt, nl, nh);
        qBDD_unprotect(nh);
        qBDD_unprotect(nl);

        return new;
    }
}

qBDD interface_gate_rx_pihalf(size_t xt, qBDD low, qBDD high) {
    qBDD rot_low, rot_high;

    rot_low = unary_apply(low, rotateCoef2);
    qBDD_protect(rot_low);
    rot_high = unary_apply(high, rotateCoef2);
    qBDD_protect(rot_high);
    qBDD newLow = binary_apply(low, rot_high, subLeaf);
    qBDD_protect(newLow);
    qBDD newHigh = binary_apply(high, rot_low, subLeaf);
    qBDD_protect(newHigh);
    qBDD_unprotect(rot_high);
    qBDD_unprotect(rot_low);
    qBDD res = newqBDD(xt, newLow, newHigh);
    qBDD_unprotect(newLow);
    qBDD_unprotect(newHigh);
    return res;
}

qBDD interface_gate_ry_pihalf(size_t xt, qBDD low, qBDD high) {
    qBDD newLow = binary_apply(low, high, subLeaf);
    qBDD_protect(newLow);
    qBDD newHigh = binary_apply(low, high, addLeaf);
    qBDD_protect(newHigh);
    qBDD res = newqBDD(xt, newLow, newHigh);
    qBDD_unprotect(newLow);
    qBDD_unprotect(newHigh);
    return res;
}

void gate_x(qBDD *p_t, uint32_t xt)
{
    // Check if xt is a missing root not needed -> swap is identical with the original MTBDD
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_x));
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_y(qBDD *p_t, uint32_t xt)
{
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_y));
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_z(qBDD *p_t, uint32_t xt)
{
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_z));
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_s(qBDD *p_t, uint32_t xt)
{
    size_t newxt = (size_t) xt;
    //printf("protect 7\n");
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_s));
    //printf("unprotect number 4\n"); 
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_t(qBDD *p_t, uint32_t xt)
{   
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_t));
    //incInvSqrtCoeff(); // TODO CORRECT?
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_h(qBDD *p_t, uint32_t xt)
{   
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_h));
    incInvSqrtCoeff();
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_rx_pihalf(qBDD *p_t, uint32_t xt)
{
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_rx_pihalf));
    incInvSqrtCoeff();
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_ry_pihalf(qBDD *p_t, uint32_t xt)
{
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_ry_pihalf));
    incInvSqrtCoeff();
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_cnot(qBDD *p_t, uint32_t xt, uint32_t xc)
{
    qBDD t = *p_t;
    qBDD_protect(t);
    qBDD res;

    res = my_mtbdd_b_xt_comp_mul_i(t, xc); // Bxc_c * T
    qBDD_protect(res);

    qBDD t_xt = my_mtbdd_t_xt_i(t, xt);
    qBDD_protect(t_xt);
    qBDD bracket_left = my_mtbdd_b_xt_comp_mul_i(t_xt, xt); // Bxt_c * Txt
    qBDD_protect(bracket_left);
    qBDD_unprotect(t_xt);

    qBDD t_xt_comp = my_mtbdd_t_xt_comp_i(t, xt);
    qBDD_protect(t_xt_comp);
    qBDD_unprotect(t);
    qBDD bracket_right =  my_mtbdd_b_xt_mul_i(t_xt_comp, xt); // Bxt * Txt_c
    qBDD_protect(bracket_right);
    qBDD_unprotect(t_xt_comp);

    qBDD inter_res = binary_apply(bracket_left, bracket_right, addLeaf); // (Bxt_c * Txt) + (Bxt * Txt_c)
    qBDD_protect(inter_res);
    qBDD_unprotect(bracket_left);
    qBDD_unprotect(bracket_right);
    qBDD inter_res2 = my_mtbdd_b_xt_mul_i(inter_res, xc); // Bxc * (Bxt_c * Txt + Bxt * Txt_c)
    qBDD_unprotect(inter_res);
    qBDD_protect(inter_res2);
    qBDD res2 = qBDD_protect(binary_apply(res, inter_res2, addLeaf)); // (Bxc_c * T) + (Bxc * (Bxt_c * Txt + Bxt * Txt_c))
    qBDD_unprotect(inter_res2);
    qBDD_unprotect(*p_t);

    *p_t = res2;
}

void gate_cz(qBDD *p_t, uint32_t xt, uint32_t xc)
{
    size_t newxt = (size_t) xt;
    size_t newxc = (size_t) xc;
    size_t ctrls[2] = {newxc, newxt};
    //printf("protect 21\n");
    qBDD res = qBDD_protect(bdd_operation(*p_t, ctrls, 1, interface_gate_z));

    //printf("unprotect number 16\n"); 
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_toffoli(qBDD *p_t, uint32_t xt, uint32_t xc1, uint32_t xc2) {
    
    qBDD t = *p_t;
    qBDD_protect(t);
    qBDD res = my_mtbdd_b_xt_comp_mul_i(t, xc1); // Bxc_c * T
    qBDD_protect(res);
    qBDD t_xt = my_mtbdd_t_xt_i(t, xt);
    qBDD_protect(t_xt);
    qBDD bracket_left = my_mtbdd_b_xt_comp_mul_i(t_xt, xt); // Bxt_c * Txt
    qBDD_protect(bracket_left);
    qBDD_unprotect(t_xt);
    qBDD t_xt_comp = my_mtbdd_t_xt_comp_i(t, xt);
    qBDD_protect(t_xt_comp);
    qBDD bracket_right = my_mtbdd_b_xt_mul_i(t_xt_comp, xt); // Bxt * Txt_c
    qBDD_protect(bracket_right);
    qBDD_unprotect(t_xt_comp);
    qBDD inter_res = binary_apply(bracket_left, bracket_right, addLeaf); // (Bxt_c * Txt) + (Bxt * Txt_c)
    qBDD_protect(inter_res);
    qBDD_unprotect(bracket_left);
    qBDD_unprotect(bracket_right);
    qBDD bracket_right2 = my_mtbdd_b_xt_mul_i(inter_res, xc2); // Bxc' * (Bxt_c * Txt + Bxt * Txt_c)
    qBDD_protect(bracket_right2);
    qBDD_unprotect(inter_res);
    qBDD bracket_left2 = my_mtbdd_b_xt_comp_mul_i(t, xc2); // Bxc'_c * T
    qBDD_protect(bracket_left2);
    qBDD_unprotect(t);
    qBDD inter_res2 = binary_apply(bracket_left2, bracket_right2, addLeaf); // (Bxc'_c * T) + (Bxc' * (Bxt_c * Txt + Bxt * Txt_c))
    qBDD_protect(inter_res2);
    qBDD_unprotect(bracket_left2);
    qBDD_unprotect(bracket_right2);
    qBDD inter_res3 = my_mtbdd_b_xt_mul_i(inter_res2, xc1); // Bxc * (Bxc'_c * T + Bxc' * (Bxt_c * Txt + Bxt * Txt_c))
    qBDD_protect(inter_res3);
    qBDD_unprotect(inter_res2);
    qBDD res2 = binary_apply(res, inter_res3, addLeaf); // (Bxc_c * T) + (Bxc * (Bxc'_c * T + Bxc' * (Bxt_c * Txt + Bxt * Txt_c)))
    qBDD_protect(res2);
    qBDD_unprotect(res);
    qBDD_unprotect(inter_res3);
    *p_t = res2;
}

void gate_mcx(qBDD *p_t, qparam_list_t *qparams) {
// Assumes indices in qparams are in the reverse order of the input file! (target qubit first, qc1 last)

    qBDD t = *p_t;
    qBDD_protect(t);
    qBDD res;

    qparam_list_first(qparams);

    // Target qubit part
    qBDD t_xt = my_mtbdd_t_xt_i(t, qparams->active->q_index);
    qBDD_protect(t_xt);
    qBDD bracket_left = my_mtbdd_b_xt_comp_mul_i(t_xt, qparams->active->q_index); // Bxt_c * Txt
    qBDD_protect(bracket_left);
    qBDD_unprotect(t_xt);

    qBDD t_xt_comp = my_mtbdd_t_xt_comp_i(t, qparams->active->q_index);
    qBDD_protect(t_xt_comp);
    qBDD bracket_right = my_mtbdd_b_xt_mul_i(t_xt_comp, qparams->active->q_index); // Bxt * Txt_c
    qBDD_protect(bracket_right);
    qBDD_unprotect(t_xt_comp);

    res = binary_apply(bracket_left, bracket_right, addLeaf); // (Bxt_c * Txt) + (Bxt * Txt_c)
    qBDD_protect(res);

    // Get the last control qubit
    qparam_list_next(qparams);

    // Handling the control qubits
    while (qparams->active) {
        qBDD bracket_right1 = my_mtbdd_b_xt_mul_i(res, qparams->active->q_index);   // BxcN * (result so far)
        qBDD_protect(bracket_right1);
        qBDD bracket_left1 = my_mtbdd_b_xt_comp_mul_i(t, qparams->active->q_index); // BxcN_c * T
        qBDD_protect(bracket_left1);
        qBDD_unprotect(res);
        res = binary_apply(bracket_left1, bracket_right1, addLeaf);
        qBDD_protect(res);                   // (BxcN_c * T) + (BxcN * (result so far))
        qBDD_unprotect(bracket_left1);
        qBDD_unprotect(bracket_right1);
        qparam_list_next(qparams);  // move to the next control qubit
    }
    qBDD_unprotect(t);
    qBDD_unprotect(bracket_left);
    qBDD_unprotect(bracket_right);
    qBDD_unprotect(*p_t);
    qBDD_protect(res);
    *p_t = res;
}
/* end of "gates.c" */