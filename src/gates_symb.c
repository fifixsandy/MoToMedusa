#include "gates_symb.h"
#include "gates.h"

void gate_symb_x(qBDD *p_t, uint32_t xt) {
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_x));
    qBDD_unprotect(*p_t);
    *p_t = res;
}

qBDD interface_gate_symb_y(size_t xt, qBDD low, qBDD high) {
    qBDD nhigh = my_mtbdd_symb_neg_i(high);
    qBDD_protect(nhigh);
    qBDD updated = newqBDD(xt, nhigh, low);
    qBDD_protect(updated);
    qBDD_unprotect(nhigh);
    qBDD res = my_mtbdd_symb_coef_rot2_i(updated);
    qBDD_unprotect(updated);
    return res;
}

void gate_symb_y(qBDD *p_t, uint32_t xt) {
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_symb_y));
    qBDD_unprotect(*p_t);
    *p_t = res;
}

qBDD interface_gate_symb_z(size_t xt, qBDD low, qBDD high) {
    qBDD nhigh = my_mtbdd_symb_neg_i(high);
    qBDD_protect(nhigh);
    qBDD res = newqBDD(xt, low, nhigh);
    qBDD_unprotect(nhigh);
    return res;
}

void gate_symb_z(qBDD *p_t, uint32_t xt) {
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_symb_z));
    qBDD_unprotect(*p_t);
    *p_t = res;
}

qBDD interface_gate_symb_s(size_t xt, qBDD low, qBDD high) {
    qBDD newH = my_mtbdd_symb_coef_rot2_i(high);
    qBDD_protect(newH);
    qBDD res = newqBDD(xt, low, newH);
    qBDD_unprotect(newH);
    return res;
}

void gate_symb_s(qBDD *p_t, uint32_t xt) {
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_symb_s));
    qBDD_unprotect(*p_t);
    *p_t = res;
}

qBDD interface_gate_symb_t(size_t xt, qBDD low, qBDD high) {
    qBDD newH = my_mtbdd_symb_coef_rot1_i(high);
    qBDD_protect(newH);
    qBDD res = newqBDD(xt, low, newH);
    qBDD_unprotect(newH);
    return res;
}

void gate_symb_t(qBDD *p_t, uint32_t xt) {
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_symb_t));
    //incInvSqrtCoeff(); // TODO CORRECT?
    qBDD_unprotect(*p_t);
    *p_t = res;
}

qBDD interface_gate_symb_h(size_t xt, qBDD low, qBDD high) {
    // if (low == high) {
    //     qBDD nlow = my_mtbdd_symb_times_c_i(low, 2);
    //     qBDD_protect(nlow);
    //     qBDD res = newqBDD(xt, nlow, qBDD_false());
    //     qBDD_unprotect(nlow);
    //     return res;
    // } else {
        qBDD nh = qBDD_protect(my_mtbdd_symb_minus_i(low, high));
        qBDD nl = qBDD_protect(my_mtbdd_symb_plus_i(low, high));
        qBDD new = newqBDD(xt, nl, nh);
        qBDD_unprotect(nh);
        qBDD_unprotect(nl);
        return new;
    //}
}

void gate_symb_h(qBDD *p_t, uint32_t xt)
{   
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_symb_h));
    incInvSqrtCoeffSymb();
    qBDD_unprotect(*p_t);
    *p_t = res;
}

qBDD interface_gate_symb_rx_pihalf(size_t xt, qBDD low, qBDD high) {
    qBDD rot_low, rot_high;

    rot_low = my_mtbdd_symb_coef_rot2_i(low);
    qBDD_protect(rot_low);
    rot_high = my_mtbdd_symb_coef_rot2_i(high);
    qBDD_protect(rot_high);
    qBDD newLow = my_mtbdd_symb_minus_i(low, rot_high);
    qBDD_protect(newLow);
    qBDD newHigh = my_mtbdd_symb_minus_i(high, rot_low);
    qBDD_protect(newHigh);
    qBDD_unprotect(rot_high);
    qBDD_unprotect(rot_low);
    qBDD res = newqBDD(xt, newLow, newHigh);
    qBDD_unprotect(newLow);
    qBDD_unprotect(newHigh);
    return res;
}

void gate_symb_rx_pihalf(qBDD *p_t, uint32_t xt) {
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_symb_rx_pihalf));
    incInvSqrtCoeffSymb();
    qBDD_unprotect(*p_t);
    *p_t = res;
}

qBDD interface_gate_symb_ry_pihalf(size_t xt, qBDD low, qBDD high) {
    qBDD newLow = my_mtbdd_symb_minus_i(low, high);
    qBDD_protect(newLow);
    qBDD newHigh = my_mtbdd_symb_plus_i(low, high);
    qBDD_protect(newHigh);
    qBDD res = newqBDD(xt, newLow, newHigh);
    qBDD_unprotect(newLow);
    qBDD_unprotect(newHigh);
    return res;
}

void gate_symb_ry_pihalf(qBDD *p_t, uint32_t xt) {
    size_t newxt = (size_t) xt;
    qBDD res = qBDD_protect(bdd_operation(*p_t, &newxt, 0, interface_gate_symb_ry_pihalf));
    incInvSqrtCoeffSymb();
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_symb_cnot(qBDD *p_t, uint32_t xt, uint32_t xc)
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

    qBDD inter_res = my_mtbdd_symb_plus_i(bracket_left, bracket_right); // (Bxt_c * Txt) + (Bxt * Txt_c)
    qBDD_protect(inter_res);
    qBDD_unprotect(bracket_left);
    qBDD_unprotect(bracket_right);
    qBDD inter_res2 = my_mtbdd_b_xt_mul_i(inter_res, xc); // Bxc * (Bxt_c * Txt + Bxt * Txt_c)
    qBDD_unprotect(inter_res);
    qBDD_protect(inter_res2);
    qBDD res2 = my_mtbdd_symb_plus_i(res, inter_res2); // (Bxc_c * T) + (Bxc * (Bxt_c * Txt + Bxt * Txt_c))
    qBDD_unprotect(inter_res2);
    qBDD_unprotect(*p_t);

    *p_t = res2;
}

void gate_symb_cz(qBDD *p_t, uint32_t xt, uint32_t xc)
{
    size_t newxt = (size_t) xt;
    size_t newxc = (size_t) xc;
    size_t ctrls[2] = {newxc, newxt};
    qBDD res = qBDD_protect(bdd_operation(*p_t, ctrls, 1, interface_gate_symb_z));
    qBDD_unprotect(*p_t);
    *p_t = res;
}

void gate_symb_toffoli(qBDD *p_t, uint32_t xt, uint32_t xc1, uint32_t xc2)
{
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
    qBDD inter_res = my_mtbdd_symb_plus_i(bracket_left, bracket_right); // (Bxt_c * Txt) + (Bxt * Txt_c)
    qBDD_protect(inter_res);
    qBDD_unprotect(bracket_left);
    qBDD_unprotect(bracket_right);
    qBDD bracket_right2 = my_mtbdd_b_xt_mul_i(inter_res, xc2); // Bxc' * (Bxt_c * Txt + Bxt * Txt_c)
    qBDD_protect(bracket_right2);
    qBDD_unprotect(inter_res);
    qBDD bracket_left2 = my_mtbdd_b_xt_comp_mul_i(t, xc2); // Bxc'_c * T
    qBDD_protect(bracket_left2);
    qBDD_unprotect(t);
    qBDD inter_res2 = my_mtbdd_symb_plus_i(bracket_left2, bracket_right2); // (Bxc'_c * T) + (Bxc' * (Bxt_c * Txt + Bxt * Txt_c))
    qBDD_protect(inter_res2);
    qBDD_unprotect(bracket_left2);
    qBDD_unprotect(bracket_right2);
    qBDD inter_res3 = my_mtbdd_b_xt_mul_i(inter_res2, xc1); // Bxc * (Bxc'_c * T + Bxc' * (Bxt_c * Txt + Bxt * Txt_c))
    qBDD_protect(inter_res3);
    qBDD_unprotect(inter_res2);
    qBDD res2 = my_mtbdd_symb_plus_i(res, inter_res3); // (Bxc_c * T) + (Bxc * (Bxc'_c * T + Bxc' * (Bxt_c * Txt + Bxt * Txt_c)))
    qBDD_protect(res2);
    qBDD_unprotect(res);
    qBDD_unprotect(inter_res3);
    *p_t = res2;
}

void gate_symb_mcx(qBDD *p_t, qparam_list_t *qparams)
{
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

    res = my_mtbdd_symb_plus_i(bracket_left, bracket_right); // (Bxt_c * Txt) + (Bxt * Txt_c)
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
        res = my_mtbdd_symb_plus_i(bracket_left1, bracket_right1);                   // (BxcN_c * T) + (BxcN * (result so far))
        qBDD_protect(res);
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

/* end of "gates_symb.c" */