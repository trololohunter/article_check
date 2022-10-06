//
// Created by vover on 3/26/18.
//

#ifndef UNTITLED2_GAS_TWO_H
#define UNTITLED2_GAS_TWO_H

#include <math.h>
#include <glob.h>

#define zero_spl_x  32 //zero splitting x
#define zero_spl_y  32
#define zero_spl_t  32

#define SMOOTH_SOLUTION 1
#define SOKOL 0
#define GAMMA 1.4
#define WALL 0

typedef struct
{

    double Segm_T;
    double Segm_X;
    double Segm_Y;
    double p_ro;
    double mu;
    double omega;
    double p_ro_0;
    double C_ro;

} P_gas;

typedef struct
{

    int M_x;
    int M_y;
    int N;
    int Dim;
    int vecDim;
    double h_x;
    double h_y;
    double tau;
    double eta;

} P_she;

typedef struct
{

    double thx;
    double thy;
    double thx2;
    double thy2;
    double thx4;
    double thy4;
    double thxx43;
    double thyy43;
    double thxxmu;
    double thyymu;
    double thxy;
    double crohx;
    double crohy;
    double ro0hx;
    double ro0hy;


} T_const;


typedef struct
{

    size_t mmp00;
    size_t mmv100;
    size_t mmv200;
    size_t mmpL0;
    size_t mmv1L0;
    size_t mmv2L0;
    size_t mmpR0;
    size_t mmv1R0;
    size_t mmv2R0;
    size_t mmp0L;
    size_t mmv10L;
    size_t mmv20L;
    size_t mmp0R;
    size_t mmv10R;
    size_t mmv20R;
    size_t mmpLL;
    size_t mmv1LL;
    size_t mmv2LL;
    size_t mmpLR;
    size_t mmv1LR;
    size_t mmv2LR;
    size_t mmpRL;
    size_t mmv1RL;
    size_t mmv2RL;
    size_t mmpRR;
    size_t mmv1RR;
    size_t mmv2RR;


} MM_step;

void param_dif (P_gas *p_d, double mu);
void param_she_step(P_she *p_s, P_gas p_d, int it_t, int it_sp);
void param_t_const (T_const *t_c, P_she p_s, P_gas p_d);
void param_MM_step (MM_step *MM_s, size_t mm, int n, int m);
void Setka (int *st, P_she *p_s);
void Sxema (double *P, double *V1, double *V2, int *st, P_she p_s, P_gas p_d);

#endif //UNTITLED2_GAS_TWO_H
