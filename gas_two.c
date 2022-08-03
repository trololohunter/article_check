//
// Created by vover on 3/5/18.
//

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "laspack/itersolv.h"
#include "gas_two.h"
#include "case.h"
#include "residuals.h"
#include "laspack/rtc.h"
#include "functions.h"
#include "gnuploting.h"


#define EPS 1e-8
#define EEPPSS 1e-15
#define MAX_ITER 2000


void param_dif (P_gas *p_d)
{
    p_d->Segm_X = 2*M_PI;
    p_d->Segm_Y = 2*M_PI;
    p_d->Segm_T = 10;
    p_d->mu = 1;
    p_d->p_ro = 1;
    p_d->omega = 1;
    p_d->p_ro_0 = 1;
    p_d->C_ro = 1;

    return;
}

void param_she_step(P_she *p_s, P_gas p_d, int it_t, int it_sp)
{
    p_s->M_x = zero_spl_x * it_sp;
    p_s->M_y = zero_spl_y * it_sp;
    p_s->N = zero_spl_t * it_t;
    p_s->eta = 1;
    p_s->h_x = p_d.Segm_X / p_s->M_x;
    p_s->h_y = p_d.Segm_Y / p_s->M_y;
    p_s->tau = p_d.Segm_T / p_s->N;
    p_s->Dim = (p_s->M_x + 1) * (p_s->M_y + 1) - 1;
    p_s->vecDim = 3 * p_s->M_x * p_s->M_y + p_s->M_x + p_s->M_y;

    return;
}

void param_t_const (T_const *t_c, P_she p_s, P_gas p_d)
{
    double thx_ = p_s.tau/p_s.h_x;
    double thy_ = p_s.tau/p_s.h_y;
    t_c->thx = thx_;
    t_c->thy = thy_;
    t_c->thx2 = 2.*thx_;
    t_c->thy2 = 2.*thy_;
    t_c->thxx43 = 4.*p_d.mu*thx_/(3. * p_s.h_x);
    t_c->thyy43 = 4.*p_d.mu*thy_/(3. * p_s.h_y);
    t_c->thxxmu = p_d.mu * thx_/p_s.h_x;
    t_c->thyymu = p_d.mu * thy_/p_s.h_y;
    t_c->thxy = p_d.mu * p_s.tau / (3. * p_s.h_x * p_s.h_y);
    //t_c->Max = thx_ * p_d.p_ro;
    //t_c->May = thy_ * p_d.p_ro;
    t_c->crohx = p_d.C_ro * thx_;
    t_c->crohy = p_d.C_ro * thy_;
    t_c->ro0hx = p_d.p_ro_0 * thx_;
    t_c->ro0hy = p_d.p_ro_0 * thy_;

    return;
}

void param_MM_step (MM_step *MM_s, size_t mm, int n, int m)
{
    MM_s->mmp00 = mm;
    MM_s->mmv100 = mm + 1;
    MM_s->mmv200 = mm + 2;
    MM_s->mmpL0 = mm - 3;
    MM_s->mmv1L0 = mm - 2;
    MM_s->mmv2L0 = mm - 1;
    MM_s->mmpR0 = mm + 3;
    MM_s->mmv1R0 = mm + 4;
    MM_s->mmv2R0 = mm + 5;
    MM_s->mmp0L = (size_t) 3 * (m - n) + 1;
    MM_s->mmv10L = MM_s->mmp0L + 1;
    MM_s->mmv20L = MM_s->mmv10L + 1;
    MM_s->mmp0R = (size_t) 3 * (m + n) + 1;
    MM_s->mmv10R = MM_s->mmp0R + 1;
    MM_s->mmv20R = MM_s->mmv10R + 1;
    MM_s->mmpLL = MM_s->mmp0L - 3;
    MM_s->mmv1LL = MM_s->mmp0L - 2;
    MM_s->mmv2LL = MM_s->mmp0L - 1;
    MM_s->mmpLR = MM_s->mmp0R - 3;
    MM_s->mmv1LR = MM_s->mmp0R - 2;
    MM_s->mmv2LR = MM_s->mmp0R - 1;
    MM_s->mmpRL = MM_s->mmp0L + 3;
    MM_s->mmv1RL = MM_s->mmp0L + 4;
    MM_s->mmv2RL = MM_s->mmp0L + 5;
    MM_s->mmpRR = MM_s->mmp0R + 3;
    MM_s->mmv1RR = MM_s->mmp0R + 4;
    MM_s->mmv2RR = MM_s->mmp0R + 5;

    return;
}

int son (int i, int j, P_she *p_s) //  the status of the node
{
    if (i == 0 && j == 0) // x = 0; y = 0
        return 5;
    if (i == p_s->M_y-1 && j == 0) // x = 0; y = b
        return 7;
    if (i == 0 && j == p_s->M_x-1) // x = a; y = 0;
        return 6;
    if (i == p_s->M_y-1 && j == p_s->M_x-1)
        return 8;
        // x
        //if (i > 0 && i < (p_s->M_y/2) && j == p_s->M_x/2) // Interpretation
        //    return 9; = a; y = b
    if (i == 0 && j > 0 && j < p_s->M_x-1) // x in; y = 0
        return 3;
    if (i == p_s->M_y && j >= 0 && j < p_s->M_x) // x in; y = b
        return 4;
    if (i > 0 && i < p_s->M_y && j == 0) // x = 0; y in
        return 1;
    if (i >= 0 && i < p_s->M_y && j == p_s->M_x) // x = a; y in
        return 2;
   // if (i == 0 && j > p_s->M_x/2-1 && j < p_s->M_x)
    //    return 10;
    if (i == p_s->M_y-1)
        return 11;
    if (j == p_s->M_x-1)
        return 12;

    return 0; // internal mesh node
}

void print_vector_int(int* a, int n, int start)
{
    int i;
    printf("\n\n start print vector\n");
    for (i = 0; i < n; ++i)
        printf ("%d \n", a[start + i]);
    printf("let's going on\n\n\n");

    return;
}

void print_vector(double* a, int n, int start)
{
    int i;
    printf("\n\n start print vector\n");
    for (i = 0; i < n; ++i)
        printf ("%e \n", a[start + i]);
    printf("let's going on\n\n\n");

    return;
}

void Setka (int *st, P_she *p_s)
{
    int k = 0;
    int i = 0;
    int j = 0;

    for (i = 0; i < p_s->M_y; ++i) {
        for (j = 0; j < p_s->M_x; ++j) {
            st[k] = son(i, j, p_s);
          printf("k=%d  \t st[k] = %d \n", k, st[k]);
            ++k;
        }
//        printf("\n");
    }

    for (i = 0; i < p_s->M_y; ++i) {
        st[k] = son(i, p_s->M_x, p_s);
        printf("k=%d  \t st[k] = %d \n", k, st[k]);
        ++k;
    }

    for (j = 0; j < p_s->M_x; ++j) {
        st[k] = son(p_s->M_y, j, p_s);
        printf("k=%d  \t st[k] = %d \n", k, st[k]);
        ++k;
    }
//    sleep(10);
    return;

}

void printl_t_c (T_const t_c) {
    printf("t_c.thx     = %e \n", t_c.thx);
    printf("t_c.thy     = %e \n",    t_c.thy);
    printf("t_c.thx2    = %e \n",   t_c.thx2);
    printf("t_c.thy2    = %e \n",   t_c.thy2);
    printf("t_c.thxx43  = %e \n",   t_c.thxx43);
    printf("t_c.thyy43  = %e \n",    t_c.thyy43);
    printf("t_c.thxxmu  = %e \n",   t_c.thxxmu);
    printf("t_c.thyymu  = %e \n",    t_c.thyymu);
    printf("t_c.thxy    = %e \n",    t_c.thxy);
    printf("t_c.crohx   = %e \n",  t_c.crohx);
    printf("t_c.crohy   = %e \n",    t_c.crohy);
    printf("t_c.ro0hx   = %e \n",    t_c.ro0hx);
    printf("t_c.ro0hy   = %e \n",    t_c.ro0hy);
    printf("\n\n");
    //t_c->Max = thx
}

void Sxema (double *P, double *V1, double *V2, int *st, P_she p_s, P_gas p_d)
{
    int k;

    QMatrix_L A;
    Vector b, x;

    T_const t_c;
    MM_step m_s;
    Norm_Step n_s;

    double GG = 0;
    double tmp;
    size_t mm;
    int m;
    double start_norm;
    double now_norm;
    FILE *fp;

    fp = fopen ("norma.txt","w");

    first_fill_check(V1, V2, P, p_s, p_d.omega);

    print_vector(V1, p_s.M_x * p_s.M_y + p_s.M_y, 0);
    print_vector(V2, p_s.M_x * p_s.M_y + p_s.M_x, 0);
    print_vector(P, p_s.M_x * p_s.M_y, 0);

    param_t_const(&t_c, p_s, p_d);
    printl_t_c (t_c);
    SetRTCAccuracy(EPS);

    if (DEBUG == 1) {
        print_vector_int(st, p_s.Dim, 0);
        printf("Dim = %d \t vecDim = %d \n\n", p_s.Dim, p_s.vecDim);
    }

    //run_gnuplot(p_s, V1, V2, P, 0);
    //residual_Ch(V1, V2, P, p_s, &n_s, u1, u2, g);
    //start_norm = sqrt(n_s.V1norm * n_s.V1norm + n_s.V2norm * n_s.V2norm);
    //now_norm = start_norm;

    for (k = 1; k < p_s.N + 1; ++k) {
    /*k = 0;
    while (now_norm > start_norm * 0.001)
    {
        ++k;*/
        mm = 1;

        Q_Constr(&A, "A", (size_t) p_s.vecDim, False, Rowws, Normal, True);
        V_Constr(&b, "b", (size_t) p_s.vecDim, Normal, True);
        V_Constr(&x, "x", (size_t) p_s.vecDim, Normal, True);

 /*       GG = -1000000;
        for (m = 0; m < p_s.Dim; ++m)
        {
//            printf("%lf \n", GG);
            if (exp(-G[m]) > GG) GG = exp(-G[m]);
        }
        //printf("%lf \n", GG);
        param_MUM_const(&m_c, p_s, GG, p_d);
*/

        for (size_t i = 0; i < p_s.M_y * p_s.M_x; ++i)
        {
            V_SetCmp(&x, 3 * i + 1, P[i]);
            V_SetCmp(&x, 3 * i + 2, V1[i]);
            V_SetCmp(&x, 3 * i + 3, V2[i]);
        }
        for (size_t i = p_s.M_x * p_s.M_y; i < p_s.M_x * p_s.M_y + p_s.M_y; ++i)
            V_SetCmp(&x, 2*p_s.M_x * p_s.M_y+i+1, V1[i]);
        for (size_t i = p_s.M_x * p_s.M_y + p_s.M_y; i < p_s.Dim; ++i)
            V_SetCmp(&x, 2*p_s.M_x * p_s.M_y+i+1, V2[i-p_s.M_y]);

        if (DEBUG == 1) {
            printf("x before begin \n");
            for (size_t i = 1; i < p_s.vecDim + 1; ++i)
                printf("%e \n", V_GetCmp(&x, i));
            printf("x before end \n");
        }

        for (m = 0; m < p_s.Dim; ++m)
        {
            param_MM_step(&m_s,mm,p_s.M_x, m);
            switch (st[m])
            {
                case 0:
                    mm = case0(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 1:
                    mm = case1(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 2:
                    mm = case2(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 3:
                    mm = case3(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 4:
                    mm = case4(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 5:
                    mm = case5(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 6:
                    mm = case6(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 7:
                    mm = case7(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 8:
                    mm = case8(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                /*case 9:
                    mm = case0(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 10:
                    mm = case10(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break; */
                case 11:
                    mm = case11(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                case 12:
                    mm = case12(&A, &b, t_c, m_s, k, V1, V2, P, m, p_s, p_d.p_ro_0, mm);
                    break;
                default:
                {
                    printf("FATAL ERROR: unknown type of node");
                    exit(1);
                }

            }
           // printf("%d \n", m);
            ++mm;
        }

        if (DEBUG == 1) {
            printf("b before begin \n");
            for (size_t i = 1; i < p_s.vecDim + 1; ++i)
                printf("%e \n", V_GetCmp(&b, i));
            printf("b before end \n");
        }

        CGSIter(&A, &x, &b, MAX_ITER, SSORPrecond, 1);
        //CGSIter(&A, &x, &b, MAX_ITER, NULL, 1);
        //JacobiIter(&A, &x, &b, MAX_ITER, NULL, 1);

        for (size_t i = 0; i < p_s.M_y * p_s.M_x; ++i)
        {
            P[i] = V_GetCmp(&x, 3 * i + 1);
            V1[i] = V_GetCmp(&x, 3 * i + 2);
            V2[i] = V_GetCmp(&x, 3 * i + 3);
        }
        for (size_t i = p_s.M_x * p_s.M_y; i < p_s.M_x * p_s.M_y + p_s.M_y; ++i)
            V1[i] = V_GetCmp(&x, 2*p_s.M_x * p_s.M_y+i+1);
        for (size_t i = p_s.M_x * p_s.M_y + p_s.M_y; i < p_s.Dim; ++i)
            V2[i-p_s.M_y] = V_GetCmp(&x, 2*p_s.M_x * p_s.M_y+i+1);

if (DEBUG == 1) {
    printf("x after begin \n");
    for (size_t i = 1; i < p_s.vecDim+1; ++i )
        printf("%e \n", V_GetCmp(&x, i));
    printf("x after end \n");
    print_vector(V1, p_s.M_x * p_s.M_y + p_s.M_y, 0);
    print_vector(V2, p_s.M_x * p_s.M_y + p_s.M_x, 0);
    print_vector(P, p_s.M_x * p_s.M_y, 0);
}
        for (m = 0; m < p_s.Dim; ++m)
        {
            if (fabs(P[m])  < EEPPSS ) P[m] = 0;
            if (fabs(V1[m]) < EEPPSS ) V1[m] = 0;
            if (fabs(V2[m]) < EEPPSS ) V2[m] = 0;
        }
        //usleep(10);


        Q_Destr(&A);
        V_Destr(&b);
        V_Destr(&x);

        //residual_Ch(V1, V2, P, p_s, &n_s, u1, u2, g);
        //fprintf(fp, "%lf \t %lf \n", k*p_s.tau, sqrt(n_s.V1norm * n_s.V1norm + n_s.V2norm * n_s.V2norm));
        //if (SMOOTH_SOLUTION != 1)
        now_norm = sqrt(n_s.V1norm * n_s.V1norm + n_s.V2norm * n_s.V2norm);
        //if (k % 20 == 1)    run_gnuplot(p_s, V1, V2, P, k);
    }

    printf("\n\n time: %lf \n\n", k*p_s.tau);

    fclose(fp);

    return;
}


