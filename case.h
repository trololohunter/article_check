//
// Created by vover on 4/4/18.
//

#ifndef UNTITLED2_CASE_H
#define UNTITLED2_CASE_H

#include "functions.h"
#include <stdio.h>
#include "gas_two.h"
#include "laspack/itersolv.h"
#include "functions.h"
#include <math.h>
#include <unistd.h>
#include "residuals.h"

#define DEBUG 3

void first_fill (double *V1, double *V2, double *P, P_she p_s,
                 double w, func u1, func u2, func ro);
void first_fill___ (double *V1, double *V2, double *P, P_she p_s,
                 double w);

void first_fill_check (double *V1, double *V2, double *P, P_she p_s,
                       double w);
void first_fill_check_vjump (double *V1, double *V2, double *P, P_she p_s);
void first_fill_check_pjump (double *V1, double *V2, double *P, P_she p_s);

size_t case0 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case1 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case2 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case3 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case4 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case5 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case6 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case7 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case8 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case9 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case10 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case11 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);
size_t case12 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm);

#endif //UNTITLED2_CASE_H
