//
// Created by vover on 3/26/18.
//

#include <stdio.h>
#include "gas_two.h"
#include "laspack/itersolv.h"
#include "functions.h"
#include "case.h"

typedef struct
{
    double p00;
    double pL0;
    double pR0;
    double p0L;
    double p0R;
    double pLL;
    double pRL;
    double pLR;
    double pRR;
    double v100;
    double v1L0;
    double v1R0;
    double v10L;
    double v10R;
    double v1LL;
    double v1RL;
    double v1LR;
    double v1RR;
    double v200;
    double v2L0;
    double v2R0;
    double v20L;
    double v20R;
    double v2LL;
    double v2RL;
    double v2LR;
    double v2RR;

} gv1v2;

void printl_m_s (MM_step m_s, int m) {
    printf("m = %d \n", m);
    printf("mmp00 = %d \n", m_s.mmp00);
    printf("mmpL0 = %d \n", m_s.mmpL0);
    printf("mmp0L = %d \n", m_s.mmp0L);
    printf("mmv100 = %d \n", m_s.mmv100);
    printf("mmv1R0 = %d \n", m_s.mmv1R0);
    printf("mmv1L0 = %d \n", m_s.mmv1L0);
    printf("mmv10R = %d \n", m_s.mmv10R);
    printf("mmv10L = %d \n", m_s.mmv10L);
    printf("mmv1RL = %d \n", m_s.mmv1RL);
    printf("mmv200 = %d \n", m_s.mmv200);
    printf("mmv2R0 = %d \n", m_s.mmv2R0);
    printf("mmv2L0 = %d \n", m_s.mmv2L0);
    printf("mmv20R = %d \n", m_s.mmv20R);
    printf("mmv20L = %d \n", m_s.mmv20L);
    printf("mmv2LR = %d \n", m_s.mmv2LR);
    printf("m = %d \n", m);
    printf("\n");
}



void first_fill(double *V1, double *V2, double *G, P_she p_s, double w, func u1, func u2, func ro) {
    int i, j;
    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1[i * (p_s.M_x + 1) + j] = u1(0, j * p_s.h_x, i * p_s.h_y);
            V2[i * (p_s.M_x + 1) + j] = u2(0, j * p_s.h_x, i * p_s.h_y);
            G[i * (p_s.M_x + 1) + j] = ro(0, j * p_s.h_x, i * p_s.h_y);
            /*if (i > 0 && i < (p_s.M_y/2) && j == p_s.M_x/2)
            {
                V1[i * (p_s.M_x + 1) + j] = 0;
                V2[i * (p_s.M_x + 1) + j] = 0;
            }*/
            if (i == p_s.M_y / 2 + 2 && j > 0 && j < p_s.M_x) {
                V1[i * (p_s.M_x + 1) + j] = 0;
                V2[i * (p_s.M_x + 1) + j] = sin(j * p_s.h_x / 2);
            }

        }
    //for (j = 1; j < p_s.M_x/2; ++j)
    //V2[j] = w;
}

void first_fill___ (double *V1, double *V2, double *G, P_she p_s,
                    double w)
{
    int i, j;
    /* for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1[i * (p_s.M_x + 1) + j] = 0;
            V2[i * (p_s.M_x + 1) + j] = 0;
            G[i * (p_s.M_x + 1) + j] = 0;
        }
    for (j = 1; j < p_s.M_x/2; ++j)
    V2[j] = w;
    */
    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1[i * (p_s.M_x + 1) + j] = 0;
            V2[i * (p_s.M_x + 1) + j] = 0;
            G[i * (p_s.M_x + 1) + j] = 1;
            /*if (i > 0 && i < (p_s.M_y/2) && j == p_s.M_x/2)
            {
                V1[i * (p_s.M_x + 1) + j] = 0;
                V2[i * (p_s.M_x + 1) + j] = 0;
            }*/
            if (i > 2 * p_s.M_y / 5 && i < 3 * p_s.M_y / 5
            &&
            j > 2 * p_s.M_x / 5 && j < 3 * p_s.M_x / 5) {
                V1[i * (p_s.M_x + 1) + j] = 1;
                V2[i * (p_s.M_x + 1) + j] = 1;
                //G[i * (p_s.M_x + 1) + j] = 2;
            }

        }

}

void first_fill_check (double *V1, double *V2, double *P, P_she p_s,
                       double w) {
    int i, j;
    for (i = 0; i < p_s.M_y; ++i)
        for (j = 0; j < p_s.M_x; ++j)
        {
            V1[i * (p_s.M_x) + j] = 0;
            V2[i * (p_s.M_x) + j] = 0;
            P[i * (p_s.M_x) + j] = 0;
        }

    V1[3] = 1;
    V1[9] = 1;
    V2[5] = 1;
    V2[14] = 1;
    P[0] = 1;
    P[7] = 1;

    for (i = p_s.M_x * p_s.M_y; i < p_s.M_x * p_s.M_y + p_s.M_y; ++i)
        V1[i] = 0;
    for (i = p_s.M_x * p_s.M_y; i < p_s.M_x * p_s.M_y + p_s.M_x; ++i)
        V2[i] = 0;

    return;
}

size_t case0 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("case0 P[m] = %e \n V1[m] = %e \n V2[m] = %e \n\n\n", P[m], V1[m], V2[m]);
    }
    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case0 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0 + 2 * t_c.thxx43 + 2 * t_c.thyymu;
    if (DEBUG == 1) {
        printf("case0 Av1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, m_s.mmv1R0, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10R, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thyymu);
    Q_SetEntry(A, mm, 10, m_s.mmv2L0, -t_c.thxy);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case0 Bv1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0+ 2 * t_c.thyy43 + 2 * t_c.thxxmu;
    if (DEBUG == 1) {
        printf("vase0 Av2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thyy43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thyy43);
    Q_SetEntry(A, mm, 8, m_s.mmv2R0, -t_c.thxxmu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thxxmu);
    Q_SetEntry(A, mm, 10, m_s.mmv10L, -t_c.thxy);

    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case0 Bv2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}

size_t case1 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("case1 P[m] = %e \n V1[m] = %e \n V2[m] = %e \n\n\n", P[m], V1[m], V2[m]);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case1 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, 0.);
    if (DEBUG == 1) {
        printf("case1 Bv1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 10);

    tmp = ro0+ 2 * t_c.thyy43 + t_c.thxxmu;
    if (DEBUG == 1) {
        printf("case1 Av2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy); //0
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thyy43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thyy43);
    Q_SetEntry(A, mm, 8, m_s.mmv2R0, -t_c.thxxmu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thxy);


    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case1 Bv2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}

size_t case2 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp = 0;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
    }

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);
    V_SetCmp(B, mm, 0.);

    if (DEBUG == 1) {
        printf("case2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}


size_t case3 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("csae3 P[m] = %e \n V1[m] = %e \n\n\n", P[m], V1[m]);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case3 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 10);

    tmp = ro0 + 2 * t_c.thxx43 + t_c.thyymu;
    if (DEBUG == 1) {
        printf("case3 Av1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, m_s.mmv1R0, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10R, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thxy);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case3 Bv1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.); //mm=mmv200
    V_SetCmp(B, mm, 0.);
    if (DEBUG == 1) {
        tmp = 0;
        printf("case3 Bv2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}

size_t case4 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp = 0;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
    }

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);
    V_SetCmp(B, mm, 0.);
    if (DEBUG == 1) {
        printf("case4 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}


size_t case5 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp = 0;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("P[m] = %e \n\n\n", P[m]);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case5 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 1); //mm=mmv100
    Q_SetEntry(A, mm, 0, mm, 1.);
    V_SetCmp(B, mm, 0.);

    mm++;

    Q_SetLen(A, mm, 1); //mm=mmv200
    Q_SetEntry(A, mm, 0, mm, 1.);
    V_SetCmp(B, mm, 0.);

    return mm;
}


size_t case6 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    m_s.mmv1R0 = (size_t) 3*p_s.M_x*p_s.M_y + 1;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("case6 P[m] = %e \n V1[m] = %e \n\n\n", P[m], V1[m]);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case6 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 10);

    tmp = ro0 + 2 * t_c.thxx43 + t_c.thyymu;
    if (DEBUG == 1) {
        printf("case6 Av1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, m_s.mmv1R0, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10R, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thxy);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case6 Bv1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.); //mm=mmv200
    V_SetCmp(B, mm, 0.);

    return mm;
}


size_t case7 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    m_s.mmv20R = (size_t) 3*p_s.M_x*p_s.M_y + p_s.M_y;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("case7 P[m] = %e \n V1[m] = %e \n V2[m] = %e \n\n\n", P[m], V1[m], V2[m]);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case7 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, 0.);

    mm++;

    Q_SetLen(A, mm, 10);

    tmp = ro0+ 2 * t_c.thyy43 + t_c.thxxmu;
    if (DEBUG == 1) {
        printf("case7 Av2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy); //0
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thyy43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thyy43);
    Q_SetEntry(A, mm, 8, m_s.mmv2R0, -t_c.thxxmu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thxy);


    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case7 Bv2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}


size_t case8 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    m_s.mmv1R0 = (size_t) 3*p_s.M_x*p_s.M_y + m/(p_s.M_x) + 1;
    m_s.mmv20R = (size_t) 3*p_s.M_x*p_s.M_y + m - p_s.M_x*(p_s.M_y-1) + p_s.M_y +1;
    m_s.mmv2LR = (size_t) 3*p_s.M_x*p_s.M_y + m - p_s.M_x*(p_s.M_y-1) + p_s.M_y +2;
    m_s.mmv1RL = (size_t) 3*p_s.M_x*p_s.M_y + m/(p_s.M_x);

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("case8 P[m] = %e \n V1[m] = %e \n V2[m] = %e \n\n\n", P[m], V1[m], V2[m]);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case8 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 10);

    tmp = ro0 + 2 * t_c.thxx43 + t_c.thyymu;
    if (DEBUG == 1) {
        printf("case8 Av1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, m_s.mmv1R0, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10L, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thxy);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case8 Bv1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0+ 2 * t_c.thyy43 + t_c.thxxmu;
    if (DEBUG == 1) {
        printf("case8 Av2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv2L0, -t_c.thxxmu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thxy);

    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case8 Bv2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}

size_t case9 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0 + 2 * t_c.thxx43 + 2 * t_c.thyymu;
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, m_s.mmv1R0, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10R, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thyymu);
    Q_SetEntry(A, mm, 10, m_s.mmv2L0, -t_c.thxy);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0+ 2 * t_c.thyy43 + 2 * t_c.thxxmu;
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv2R0, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thyymu);
    Q_SetEntry(A, mm, 10, m_s.mmv10L, -t_c.thxy);

    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);

    return mm;
}

size_t case10 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0 + 2 * t_c.thxx43 + 2 * t_c.thyymu;
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, m_s.mmv1R0, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10R, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thyymu);
    Q_SetEntry(A, mm, 10, m_s.mmv2L0, -t_c.thxy);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0+ 2 * t_c.thyy43 + 2 * t_c.thxxmu;
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv2R0, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thyymu);
    Q_SetEntry(A, mm, 10, m_s.mmv10L, -t_c.thxy);

    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);

    return mm;
}


size_t case11 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    m_s.mmv20R = (size_t) 3*p_s.M_x*p_s.M_y /*+ p_s.M_y */ + m - p_s.M_x*(p_s.M_y-1) + p_s.M_y + 1;
    m_s.mmv2LR = (size_t) 3*p_s.M_x*p_s.M_y /*+ p_s.M_y */ + m - p_s.M_x*(p_s.M_y-1) + p_s.M_y + 2;

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("case11 P[m] = %e \n V1[m] = %e \n V2[m] = %e \n\n\n", P[m], V1[m], V2[m]);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case11 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0 + 2 * t_c.thxx43 + t_c.thyymu;
    if (DEBUG == 1) {
        printf("case11 Av1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, m_s.mmv1R0, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10L, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thxy);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case11 Bv1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0+ 2 * t_c.thyy43 + 2 * t_c.thxxmu;
    if (DEBUG == 1) {
        printf("case11 Av2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thyy43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thyy43);
    Q_SetEntry(A, mm, 8, m_s.mmv2R0, -t_c.thxxmu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thxxmu);
    Q_SetEntry(A, mm, 10, m_s.mmv10L, -t_c.thxy);

    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case11 Bv2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}


size_t case12 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    m_s.mmv1R0 = (size_t) 3*p_s.M_x*p_s.M_y + m/(p_s.M_x) + 1;
    m_s.mmv1RL = (size_t) 3*p_s.M_x*p_s.M_y + m/(p_s.M_x);

    if (DEBUG == 1) {
        printl_m_s(m_s, m);
        printf("case12 P[m] = %e \n V1[m] = %e \n V2[m] = %e \n\n\n", P[m], V1[m], V2[m]);
    }

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case12 Bp mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0 + 2 * t_c.thxx43 + 2 * t_c.thyymu;
    if (DEBUG == 1) {
        printf("case12 Av1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, m_s.mmv1R0, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10R, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thyymu);
    Q_SetEntry(A, mm, 10, m_s.mmv2L0, -t_c.thxy);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case12 Bv1 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    mm++;

    Q_SetLen(A, mm, 11);

    tmp = ro0+ 2 * t_c.thyy43 + t_c.thxxmu;
    if (DEBUG == 1) {
        printf("case12 Av2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thyy43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thyy43);
    Q_SetEntry(A, mm, 8, m_s.mmv2L0, -t_c.thxxmu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thxy);

    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);
    if (DEBUG == 1) {
        printf("case12 Bv2 mm = %d \t tmp = %e \n\n\n", mm, tmp);
    }

    return mm;
}