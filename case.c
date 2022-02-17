//
// Created by vover on 3/26/18.
//

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

void first_fill_check (double *V1, double *V2, double *G, P_she p_s,
                       double w) {
    int i, j;
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

    return;
}

size_t case0 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

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

size_t case1 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;


    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);

    tmp = ro0 * V1[m];
    V_SetCmp(B, mm, 0.);

    mm++;

    Q_SetLen(A, mm, 10);

    tmp = ro0+ 2 * t_c.thyy43 + t_c.thxxmu;
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv200
    Q_SetEntry(A, mm, 1, m_s.mmv100, t_c.thxy); //0
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohy);
    Q_SetEntry(A, mm, 3, m_s.mmp0L, -t_c.crohy);
    Q_SetEntry(A, mm, 4, m_s.mmv20R, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv1R0, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv1RL, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv20L, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv2R0, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv10L, -t_c.thxy);


    tmp = ro0 * V2[m];
    V_SetCmp(B, mm, tmp);


    return mm;
}

size_t case2 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);
    V_SetCmp(B, mm, 0.);

    return mm;
}


size_t case3 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

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

    tmp = ro0 + 2 * t_c.thxx43 + t_c.thyymu;
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

    mm++;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.); //mm=mmv100
    V_SetCmp(B, mm, 0.);


    return mm;
}

size_t case4 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);
    V_SetCmp(B, mm, 0.);

    return mm;
}


size_t case5 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, m_s.mmv1R0, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, m_s.mmv20R, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);

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

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);
    V_SetCmp(B, mm, 0.);

    return mm;
}


size_t case7 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    Q_SetLen(A, mm, 1);
    Q_SetEntry(A, mm, 0, mm, 1.);
    V_SetCmp(B, mm, 0.);

    return mm;
}


size_t case8 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

    Q_SetLen(A, mm, 5);
    Q_SetEntry(A, mm, 0, mm, 1.);
    Q_SetEntry(A, mm, 1, (size_t) 3*p_s.M_x*p_s.M_y + p_s.M_x  /* mmv1R0 */, t_c.ro0hx );
    Q_SetEntry(A, mm, 2, m_s.mmv100, -t_c.ro0hx);
    Q_SetEntry(A, mm, 3, (size_t) 3*p_s.M_x*p_s.M_y + p_s.M_y + p_s.M_x  /* mmv20R */, t_c.ro0hy);
    Q_SetEntry(A, mm, 4, m_s.mmv200, -t_c.ro0hy);

    tmp = P[m];
    V_SetCmp(B, mm, tmp);

    mm++;

    Q_SetLen(A, mm, 10);

    tmp = ro0 + 2 * t_c.thxx43 + t_c.thyymu;
    Q_SetEntry(A, mm, 0, mm, tmp); //mm=mmv100
    Q_SetEntry(A, mm, 1, m_s.mmv200, t_c.thxy);
    Q_SetEntry(A, mm, 2, m_s.mmp00, t_c.crohx);
    Q_SetEntry(A, mm, 3, m_s.mmpL0, -t_c.crohx);
    Q_SetEntry(A, mm, 4, (size_t) 3*p_s.M_x*p_s.M_y + p_s.M_x  /* mmv1R0 */, -t_c.thxx43);
    Q_SetEntry(A, mm, 5, m_s.mmv20R, -t_c.thxy);
    Q_SetEntry(A, mm, 6, m_s.mmv2LR, t_c.thxy);
    Q_SetEntry(A, mm, 7, m_s.mmv1L0, -t_c.thxx43);
    Q_SetEntry(A, mm, 8, m_s.mmv10L, -t_c.thyymu);
    Q_SetEntry(A, mm, 9, m_s.mmv2L0, -t_c.thxy);

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

size_t case9 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
              double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

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


size_t case12 (QMatrix_L *A, Vector *B, T_const t_c, MM_step m_s, int k,
               double *V1, double *V2, double *P, int m, P_she p_s, double ro0, size_t mm)
{
    double tmp;

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