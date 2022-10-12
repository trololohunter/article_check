//
// Created by vover on 06.10.22.
//

#include "exact_fading.h"

#define jump_type 1 //1 - vjump 2 - pjump

void print_norms(FILE *f, double *norm, int it_max);

int deg2inx (int x);

void article_check() {

    P_gas p_d;
    P_she p_s;
    double *P;
    double *V1;
    double *V2;
    int *st;

    FILE *f;

    param_dif(&p_d, 1);


    param_she_step(&p_s, p_d, 1, 1);

    st = (int*) malloc(p_s.Dim * sizeof(int));
    P = (double*) malloc((p_s.M_x * p_s.M_y) * sizeof(double));
    V1 = (double*) malloc((p_s.M_x * p_s.M_y + p_s.M_x) * sizeof(double));
    V2 = (double*) malloc((p_s.M_x * p_s.M_y + p_s.M_y) * sizeof(double));

    Setka(st, &p_s);

    if (jump_type == 1)
        first_fill_check_vjump(V1, V2, P, p_s);
    if (jump_type == 2)
        first_fill_check_pjump(V1, V2, P, p_s);

    Sxema_particularT (P, V1, V2, st, p_s, p_d);

    printf("asdfsdfa \n M_x = %d \n M_y = %d \n N   = %d \n", p_s.M_x, p_s.M_y, p_s.N);

    free(V1);
    free(V2);
    free(P);


    return;
}


void T_mu_jump() {

    P_gas p_d;
    P_she p_s;
    double *P;
    double *V1;
    double *V2;
    int *st;
    double T = 0;
    double mu;

    FILE *f;
    f = fopen("mu_T.txt", "w");

    for (mu = 0.01; mu < 1.01; mu += 0.01) {

        param_dif(&p_d, mu);
        param_she_step(&p_s, p_d, 1, 1);

        st = (int *) malloc(p_s.Dim * sizeof(int));
        P = (double *) malloc((p_s.M_x * p_s.M_y) * sizeof(double));
        V1 = (double *) malloc((p_s.M_x * p_s.M_y + p_s.M_x) * sizeof(double));
        V2 = (double *) malloc((p_s.M_x * p_s.M_y + p_s.M_y) * sizeof(double));

        Setka(st, &p_s);

        if (jump_type == 1)
            first_fill_check_vjump(V1, V2, P, p_s);
        if (jump_type == 2)
            first_fill_check_pjump(V1, V2, P, p_s);

        Sxema_fading_in_q_times(P, V1, V2, st, p_s, p_d, &T);

        printf("asdfsdfa \n M_x = %d \n M_y = %d \n N   = %d \n T = %lf", p_s.M_x, p_s.M_y, p_s.N, T);
        fprintf(f, "%lf \t %lf \n", mu, T);

        free(V1);
        free(V2);
        free(P);

    }

    fclose(f);

    return;
}

int deg2inx (int x)
{
    int i;
    int ans = 1;

    if (x < 0) return 0;
    if (x == 0) return 1;
    if (x == 1) return 2;
    for (i = 1; i < x+1; ++i)
        ans *= 2;
    return ans;
}

void print_norms(FILE *f, double *norm, int it_max)
{
    int k, it_t, it_sp;

    printf("start print norm \n");
    fprintf(f, "\n\n");
    fprintf(f, "\\begin{tabular}{c r r r r}\n");
    fprintf(f, "\\hline \n");
    fprintf(f, "N \\texttt{\\char`\\\\} M ");
    for (it_sp = 1; it_sp < deg2inx(it_max) + 1; it_sp *= 2)
        fprintf(f, "& %d", zero_spl_x * it_sp);
    fprintf(f, "\\\\ \n\\hline \n");

    k = 0;
    for (it_t = 1; it_t < deg2inx(it_max) + 1; it_t *= 2)
    {
        fprintf(f,"%d ", zero_spl_t * it_t);
        for (it_sp = 1; it_sp < deg2inx(it_max) + 1; it_sp *= 2)
        {
            fprintf (f, "& %e", norm[k]);
            ++k;
        }
        fprintf(f, "\\\\ \n");
    }
    fprintf(f, "\\hline \n\\end{tabular}");
}