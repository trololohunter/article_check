//
// Created by vover on 4/3/18.
//

#include <stdio.h>
#include <unistd.h>
#include "residuals.h"
#include "gas_two.h"

double residual_Ch_step (double *V1, double *V2, double *G, P_she p_s, int k,
                         func u1, func u2, func ro)
{
    double V1max = 0, V2max = 0, Gmax = 0;
    int i, j;

    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1max = (fabs(u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) > V1max) ? fabs(u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) : V1max;
            V2max = (fabs(u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) > V2max) ? fabs(u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) : V2max;
            Gmax = (fabs(ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]) > Gmax) ? fabs(ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]) : Gmax;
 /*           printf ("i = %d \t j = %d \n V1max = %e \t V2max = %e \t Gmax = %e \n", i, j,
                    fabs(u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]),
                    fabs(u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]),
                    fabs(ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]));
  */      }

//    printf ("V1max = %e \t V2max = %e \t Gmax = %e \n", V1max, V2max, Gmax);
    if (V1max > V2max && V1max > Gmax) return V1max;
    if (V2max > V1max && V2max > Gmax) return V2max;
    return Gmax;

}

double residual_Ch (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s,
                    func u1, func u2, func ro)
{
    double V1max = 0, V2max = 0, Gmax = 0;
    int i, j;

    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1max = (fabs(u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) > V1max) ? fabs(u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) : V1max;
            V2max = (fabs(u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) > V2max) ? fabs(u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) : V2max;
            Gmax = (fabs(ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]) > Gmax) ? fabs(ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]) : Gmax;
        }

    printf ("armadas V1max = %e \t V2max = %e \t Gmax = %e \n", V1max, V2max, Gmax);

    n_s->Gnorm = Gmax;
    n_s->V1norm = V1max;
    n_s->V2norm = V2max;

    if (V1max > V2max && V1max > Gmax) return V1max;
    if (V2max > V1max && V2max > Gmax) return V2max;
    return Gmax;
}

double residual_L2h_step (double *V1, double *V2, double *G, P_she p_s, int k,
                          func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(k * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);

}

double residual_L2h (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s,
                     func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += (ro(p_s.N * p_s.tau, j * p_s.h_x , i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; i < p_s.M_y; ++i)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; i < p_s.M_y; ++i)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    n_s->V1norm = sqrt(V1res);
    n_s->V2norm = sqrt(V2res);
    n_s->Gnorm  = sqrt(Gres);

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);
}



double residual_W12 (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += 2 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += 2 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += 2 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    n_s->V1norm = sqrt(V1res);
    n_s->V2norm = sqrt(V2res);
    n_s->Gnorm  = sqrt(Gres);

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);
}

double residual_Ch_Leb (double *V1, double *V2, double *P, P_she p_s, Norm_Step *n_s,
                    func u1, func u2, func ro)
{
    double V1max = 0, V2max = 0, Gmax = 0;
    int i, j;
/*
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
 */

    for (i = 0; i < p_s.M_y + 1; ++i)
        for (j = 0; j < p_s.M_x + 1; ++j)
        {
            V1max = (fabs(u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) > V1max) ? fabs(u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) : V1max;
            V2max = (fabs(u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) > V2max) ? fabs(u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) : V2max;
            if ((j < p_s.M_x) && (i < p_s.M_y))
            Gmax = (fabs(ro(p_s.N * p_s.tau, j * p_s.h_x + p_s.h_x/2., i * p_s.h_y + p_s.h_y/2.) - P[i * (p_s.M_x + 1) + j]) > Gmax) ? fabs(ro(p_s.N * p_s.tau, j * p_s.h_x + p_s.h_x/2., i * p_s.h_y + p_s.h_y/2.) - P[i * (p_s.M_x + 1) + j]) : Gmax;
        }

    printf ("armadas V1max = %e \t V2max = %e \t Gmax = %e \n", V1max, V2max, Gmax);

    n_s->Gnorm = Gmax;
    n_s->V1norm = V1max;
    n_s->V2norm = V2max;

    if (V1max > V2max && V1max > Gmax) return V1max;
    if (V2max > V1max && V2max > Gmax) return V2max;
    return Gmax;
}

double residual_L2h_Sokol (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s,
                     func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += (ro(p_s.N * p_s.tau, (j-1) * p_s.h_x + p_s.h_x/2., (i-1) * p_s.h_y + p_s.h_y/2.) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, (j-1) * p_s.h_x + p_s.h_x/2., (i-1) * p_s.h_y + p_s.h_y/2.) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; i < p_s.M_x; ++i)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; i < p_s.M_x; ++i)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    n_s->V1norm = sqrt(V1res);
    n_s->V2norm = sqrt(V2res);
    n_s->Gnorm  = sqrt(Gres);

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);
}



double residual_W12_Sokol (double *V1, double *V2, double *G, P_she p_s, Norm_Step *n_s, func u1, func u2, func ro)
{
    double V1res = 0, V2res = 0, Gres = 0;
    int i,j;

    for (i = 1; i < p_s.M_y; ++i)
        for (j = 1; j < p_s.M_x; ++j)
        {
            V1res += 2 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
            V2res += 2 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
            Gres  += 2 * (ro(p_s.N * p_s.tau, (j-1) * p_s.h_x + p_s.h_x/2., (i-1) * p_s.h_y + p_s.h_y/2.) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, (j-1) * p_s.h_x + p_s.h_x/2., (i-1) * p_s.h_y + p_s.h_y/2.) - G[i * (p_s.M_x + 1) + j]);
        }

    i = 0;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }
    i = p_s.M_y;
    for (j = 0; j < p_s.M_x + 1; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = 0;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    j = p_s.M_x;
    for (i = 1; j < p_s.M_x; ++j)
    {
        V1res += 0.5 * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]) * (u1(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V1[i * (p_s.M_x + 1) + j]);
        V2res += 0.5 * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]) * (u2(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - V2[i * (p_s.M_x + 1) + j]);
//        Gres  += 0.5 * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) -  G[i * (p_s.M_x + 1) + j]) * (ro(p_s.N * p_s.tau, j * p_s.h_x, i * p_s.h_y) - G[i * (p_s.M_x + 1) + j]);
    }

    V1res *= p_s.h_x * p_s.h_y;
    V2res *= p_s.h_x * p_s.h_y;
    Gres  *= p_s.h_x * p_s.h_y;

    n_s->V1norm = sqrt(V1res);
    n_s->V2norm = sqrt(V2res);
    n_s->Gnorm  = sqrt(Gres);

    if (V1res > V2res && V1res > Gres) return sqrt(V1res);
    if (V2res > V1res && V2res > Gres) return sqrt(V2res);
    return sqrt(Gres);
}

double residial_Ch_ (double *u, double *v, double *p, int u_size, int v_size, int p_size) {
    double Vmax = 0, Umax = 0, Pmax = 0;
    int i, j;

    for (i = 0; i < u_size; ++i)
        Umax = (u[i] > Umax) ? u[i] : Umax;

    for (j = 0; j < u_size; ++j)
        Vmax = (v[j] > Vmax) ? v[j] : Vmax;

    return (Vmax > Umax) ? Vmax : Umax;

}


void residual_matrix_schema (double *u_, double *v_, double *p_, //from previous time layer
                             double *u, double *v, double *p,   //from current time layer
                             P_she p_s, P_gas p_d) {

    int u_i = 0, v_i = 0, p_i = 0;
    double tmp, tmp_u, tmp_v, tmp_p;

    FILE *fp;
    fp = fopen ("residual_matrix_schema.txt","w");

    for (p_i = 0; p_i < p_s.M_x * p_s.M_y; ++p_i)
    {
        tmp_p = p[p_i] - p_[p_i];
        tmp_u = (p_i+1) % p_s.M_x != 0 ? u[p_i+1] - u[p_i] : -u[p_i];
        tmp_v = v[p_i + p_s.M_x] - v[p_i];
        tmp = tmp_p / p_s.tau
                + p_d.p_ro_0 * tmp_u / p_s.h_x
                + p_d.p_ro_0 * tmp_v / p_s.h_y;
        fprintf (fp, "%e \t", tmp);
        if ((p_i+1) % p_s.M_x == 0)
            fprintf(fp, "\n");
    }

    fprintf (fp, "\n\n\n");

    fprintf (fp, "%e \t", u[0]);
    for (u_i = 1; u_i < p_s.M_x; ++u_i)
    {
        tmp = 0;

    }

    for (u_i = 1; u_i < p_s.M_x * p_s.M_y; ++u_i)
    {
        tmp = p_d.p_ro_0 * (u[u_i] - u_[u_i]) / p_s.tau
              + p_d.C_ro * (p[u_i] - p[u_i-1]) / p_s.h_x
              - p_d.mu * (
                        4./3. * (u[u_i + 1] - 2 * u[u_i] + u[u_i - 1]) / (p_s.h_x*p_s.h_x)
                        + 1
                      )
                ;
        if (u_i % p_s.M_x == 0) fprintf (fp, "%e \t", u[u_i]);
        if ((u_i+1) % p_s.M_x == 0)
            fprintf(fp, "\n");
    }
    for (u_i = p_s.M_x * p_s.M_y; u_i < p_s.M_x * p_s.M_y + p_s.M_y; ++u_i)
    {


    }

    fclose(fp);

    return;
}