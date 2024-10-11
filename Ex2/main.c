/*
i_tel
11

i_le
25

i_teu
39

j_tel
0

j_le
0

j_teu
0

mach
0.3

angle_of_attack_deg
0

density
1.225

environment_pressure
101325

delta_t
0.0000001

gamma
1.4

epse
0.06

iterations
2
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>

#define MAXDIR 1000
#define MAXWORD MAXDIR
#define PI M_PI
#define dprintSTRING(expr) do{printf("%d: ", __LINE__); printf(#expr " = %s\n", expr);} while(0)   /* macro for easy debuging*/
#define dprintINT(expr) do{printf("%d: ", __LINE__); printf(#expr " = %d\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintF(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintD(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/

void read_input(char *input_dir, char *mesh_dir, double *x_vals_mat,
                double *y_vals_mat);
void read_mesh_file(FILE *mesh_fp, double *x_vals_mat,
                   double *y_vals_mat);
void read_input_file(FILE *fp);
void read_mat_from_file(FILE *fp, double *des);
void output_solution(double *current_Q, double *U_mat, double *V_mat,
                     double *x_vals_mat, double *y_vals_mat,
                     double *dxi_dx_mat, double *dxi_dy_mat,
                     double *deta_dx_mat, double *deta_dy_mat);
int offset2d(int i, int j, int ni);
int offset3d(int i, int j, int k, int ni, int nj);
void print_mat2D(double *data);
void print_layer_of_mat3D(double *data, int layer);
double first_deriv(double *mat, char diraction, int i, int j);
double second_deriv(double *mat, char diraction, int i, int j);
double calculate_one_over_jacobian_at_a_point(double *x_vals_mat,
                                              double *y_vals_mat, int i,
                                              int j);
void contravariant_velocities(double *U, double *V,
                              double *dxi_dx_mat, double *dxi_dy_mat,
                              double *deta_dx_mat, double *deta_dy_mat,
                              double *Q, int i, int j);
void calculate_u_and_v(double *u, double *v, double *Q, int i, int j);
double calculate_p(double energy, double rho, double u, double v);
double calculate_energy(double p, double u, double v, double rho);
void calculate_E_hat_at_a_point(double *E0, double *E1, double *E2,
                                double *E3, double *J_vals_mat,
                                double *dxi_dx_mat, double *dxi_dy_mat,
                                double *deta_dx_mat, double *deta_dy_mat,
                                double *Q, int i, int j);
void calculate_F_hat_at_a_point(double *F0, double *F1, double *F2,
                                double *F3, double *J_vals_mat,
                                double *dxi_dx_mat, double *dxi_dy_mat,
                                double *deta_dx_mat, double *deta_dy_mat,
                                double *Q, int i, int j);                               
void initialize_flow_field(double *Q);
void RHS(double *S, double *W, double *Q, double *J_vals_mat,
         double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat,
         double *deta_dy_mat, double *s2, double *rspec, double *qv, double *dd);
void advance_Q(double *next_Q, double *current_Q ,double *S, double *J_vals_mat);
void copy_3Dmat_to_3Dmat(double *dst, double *src);
int smooth(double *q, double *s, double *jac, double *xx, double *xy, 
           double *yx, double *yy, int id, int jd, double *s2, 
           double *rspec, double *qv, double *dd,
           double epse, double gamma, double fsmach, double dt);
void apply_BC(double *Q, double *J_vals_mat, double *dxi_dx_mat,
              double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat);
void output_mat2D_to_file(FILE *fp, double *data);
void output_layer_of_mat3D_to_file(FILE *fp, double *data, int layer);
void calculate_A_hat_j_const(double *dst, double *Q, double *dxi_dx_mat,
                             double *dxi_dy_mat, int i, int j);
void calculate_B_hat_i_const(double *dst, double *Q, double *deta_dx_mat,
                             double *deta_dy_mat, int i, int j);
int smoothx(double *q, double *xx, double *xy, int id, int jd, double *a,
	        double *b, double *c, int j,double *jac, double *drr, double *drp, 
            double *rspec, double *qv, double *dd,
            double epsi, double gamma, double fsmach, double dt);
int smoothy(double *q, double *yx, double *yy, int id, int jd, double *a,
	        double *b, double *c, int i,double *jac, double *drr, double *drp, 
            double *rspec, double *qv, double *dd,
            double epsi, double gamma, double fsmach, double dt);
void LHSX(double *A, double *B, double *C, double *Q,
          double *dxi_dx_mat, double *dxi_dy_mat, int j);
void LHSY(double *A, double *B, double *C, double *Q,
          double *deta_dx_mat, double *deta_dy_mat, int i);
int btri4s(double *a, double *b, double *c, double *f, int kd, int ks, int ke);
double calculate_S_norm(double *S);
double step(double *A, double *B, double *C, double *D, double *current_Q,
          double *S, double *W, double *J_vals_mat, double *dxi_dx_mat,
          double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat,
          double *s2, double *drr, double *drp, double *rspec, double *qv, double *dd);
void matrices_coeffic_and_Jacobian(double *J_vals_mat, double *dxi_dx_mat,
                                   double *dxi_dy_mat, double *deta_dx_mat,
                                   double *deta_dy_mat, double *x_vals_mat,
                                   double *y_vals_mat);
void initialize(double *Q, double *J_vals_mat, double *dxi_dx_mat,
                double *dxi_dy_mat, double *deta_dx_mat,
                double *deta_dy_mat, double *x_vals_mat, double *y_vals_mat);

/* global variables */
int ni, nj, max_ni_nj, i_TEL, i_LE, i_TEU, j_TEL = 0, j_LE = 0,
j_TEU = 0, max_iterations = 1e4;

double Mach_inf, angle_of_attack_deg = 0, angle_of_attack_rad, density = 1.225,     /* default values */
environment_pressure = 101325, delta_t = 1e-6, Gamma = 1.4, epse = 0.06, epsi;

int main(int argc, char const *argv[])
{
/* declerations */
    char input_dir[MAXDIR], mesh_dir[MAXDIR], current_word[MAXWORD];

    double *x_vals_mat, *y_vals_mat, *J_vals_mat, *first_Q,
    *current_Q, *next_Q, *S, *W, *dxi_dx_mat, *dxi_dy_mat, *deta_dx_mat,
    *deta_dy_mat, *s2, *rspec, *qv, *dd, *U_mat, *V_mat, *A, *B, *C, *D,
    *drr, *drp, max_S_norm = 0, current_S_norm;
    
    int i_index, j_index, k_index;

/* getting the input directory and mesh directory*/
    if (--argc != 2) {
        fprintf(stderr, "%s:%d: [Error] not right usage... Usage: main 'input dir' 'mesh dir'\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(input_dir, (*(++argv)), MAXDIR);

    if (input_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [Error] input too long\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(mesh_dir, (*(++argv)), MAXDIR);

    if (mesh_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [Error] input too long\n", __FILE__, __LINE__);
        return -1;
    }

    // dprintSTRING(input_dir);
    // dprintSTRING(mesh_dir);

/*------------------------------------------------------------*/

/* getting ni, nj*/
    FILE *mesh_fp = fopen(mesh_dir, "rt");
    if (!mesh_fp) {
        fprintf(stderr, "%s:%d: [Error] opening input file: %s\n",__FILE__, __LINE__, strerror(errno));
        exit(1);
    }
    
    while(fscanf(mesh_fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "ni")) {
            fscanf(mesh_fp, "%d", &ni);
        } else if (!strcmp(current_word, "nj")) {
            fscanf(mesh_fp, "%d", &nj);
        }
    }
    fclose(mesh_fp);
    max_ni_nj = (int)fmax((double)ni, (double)nj);

/*------------------------------------------------------------*/
/* allocating the matrices */
    x_vals_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            x_vals_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    y_vals_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            y_vals_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    J_vals_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            J_vals_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    dxi_dx_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            dxi_dx_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    dxi_dy_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            dxi_dy_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    deta_dx_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            deta_dx_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    deta_dy_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            deta_dy_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    U_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            U_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    V_mat = (double *)malloc(sizeof(double) * ni * nj);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            V_mat[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    s2 = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        s2[i_index] = 0;
    }
    rspec = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        rspec[i_index] = 0;
    }
    qv = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        qv[i_index] = 0;
    }
    dd = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        dd[i_index] = 0;
    }
    drr = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        drr[i_index] = 0;
    }
    drp = (double *)malloc(sizeof(double) * max_ni_nj);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        drp[i_index] = 0;
    }
    W = (double *)malloc(sizeof(double) * max_ni_nj * 4);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            W[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    D = (double *)malloc(sizeof(double) * max_ni_nj * 4);
    for (i_index = 0; i_index < max_ni_nj; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            D[offset2d(i_index, j_index, ni)] = 0;
        }
    }
    A = (double *)malloc(sizeof(double) * 4 * 4 * max_ni_nj);
    for (i_index = 0; i_index < 4; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            for (k_index = 0; k_index < max_ni_nj; k_index++) {
                A[offset2d(i_index, j_index, ni)] = 0;
            }
        }
    }
    B = (double *)malloc(sizeof(double) * 4 * 4 * max_ni_nj);
    for (i_index = 0; i_index < 4; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            for (k_index = 0; k_index < max_ni_nj; k_index++) {
                B[offset2d(i_index, j_index, ni)] = 0;
            }
        }
    }
    C = (double *)malloc(sizeof(double) * 4 * 4 * max_ni_nj);
    for (i_index = 0; i_index < 4; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < 4; j_index++) {
            for (k_index = 0; k_index < max_ni_nj; k_index++) {
                C[offset2d(i_index, j_index, ni)] = 0;
            }
        }
    }
    first_Q = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            for (k_index = 0; k_index < 4; k_index++) {
                first_Q[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }
    current_Q = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            for (k_index = 0; k_index < 4; k_index++) {
                current_Q[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }
    next_Q = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            for (k_index = 0; k_index < 4; k_index++) {
                next_Q[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }
    S = (double *)malloc(sizeof(double) * ni * nj * 4);
    for (i_index = 0; i_index < ni; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < nj; j_index++) {
            for (k_index = 0; k_index < 4; k_index++) {
                S[offset3d(i_index, j_index, k_index, ni, nj)] = 0;
            }
        }
    }

/*------------------------------------------------------------*/

    read_input(input_dir, mesh_dir, x_vals_mat, y_vals_mat);

/* Checking the input */
    printf("--------------------\n");
    dprintINT(ni);
    dprintINT(nj);
    dprintINT(i_TEL);
    dprintINT(i_LE);
    dprintINT(i_TEU);
    dprintINT(j_TEL);
    dprintINT(j_LE);
    dprintINT(j_TEU);
    dprintD(Mach_inf);
    dprintD(angle_of_attack_deg);
    dprintD(density);
    dprintD(environment_pressure);
    dprintD(delta_t);
    dprintD(Gamma);
    dprintINT(max_ni_nj);
    dprintD(epse);
    dprintINT(max_iterations);
    printf("--------------------\n");

/*------------------------------------------------------------*/

    initialize(current_Q, J_vals_mat, dxi_dx_mat, dxi_dy_mat, deta_dx_mat,
               deta_dy_mat, x_vals_mat, y_vals_mat);
    copy_3Dmat_to_3Dmat(first_Q, current_Q);
    
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        current_S_norm = step(A, B, C, D, current_Q, S, W, J_vals_mat,
                              dxi_dx_mat, dxi_dy_mat, deta_dx_mat,
                              deta_dy_mat, s2, drr, drp, rspec, qv, dd);
        if (max_S_norm < fabs(current_S_norm)) {
            max_S_norm = fabs(current_S_norm);
        }
        advance_Q(next_Q, current_Q, S, J_vals_mat);
        copy_3Dmat_to_3Dmat(current_Q, next_Q);
        apply_BC(current_Q, J_vals_mat, dxi_dx_mat, dxi_dy_mat, deta_dx_mat, deta_dy_mat);
        
        printf("%d: %f\n", iteration, current_S_norm);

        if (fabs(current_S_norm) / max_S_norm < 1e-5 || current_S_norm == 0 || isnan(current_S_norm)) {
            break;
        }
    }
    dprintD(max_S_norm);

        /*test*/
        for (int j = nj-1; j >=0; j--) {
            for (int i = 0; i < ni; i++) {
                if (i == i_LE) {
                    printf("  ");
                }
                double e = current_Q[offset3d(i, j, 3, ni, nj)];
                double rho = current_Q[offset3d(i, j, 0, ni, nj)]; 
                double u, v;
                calculate_u_and_v(&u, &v, current_Q, i, j);
                double p = calculate_p(e, rho, u, v);
                printf("%g ", p);
            }
            printf("\n");
        }

        // for (int i = 0; i < ni; i++) {
        //     for (int j = 0; j < nj; j++) {
        //         double U, V;
        //         contravariant_velocities(&U, &V, dxi_dx_mat, dxi_dy_mat, deta_dx_mat, deta_dy_mat, current_Q, i, j);
        //         int index = offset2d(i, j, ni);
        //         U_mat[index] = U;
        //         V_mat[index] = V;
        //     }
        // }
        // print_mat2D(U_mat);

        // print_mat2D(J_vals_mat);
        /*test*/

    // int layer = 2;
    // print_layer_of_mat3D(first_Q, layer);
    // print_layer_of_mat3D(current_Q, layer);

    // output_solution(current_Q, U_mat, V_mat, x_vals_mat, y_vals_mat);
    
/*------------------------------------------------------------*/

    free(x_vals_mat);
    free(y_vals_mat);
    free(J_vals_mat);
    free(current_Q);
    free(next_Q);
    free(S);  
    free(W);  
    free(dxi_dx_mat);
    free(dxi_dy_mat);
    free(deta_dx_mat);
    free(deta_dy_mat);
    free(s2); 
    free(rspec);
    free(qv); 
    free(dd); 

    return 0;
}

/* sets 'flags' and variables according to the input file
argument list:
input_dir - the directory of the input file 
mesh_dir - the directory of the mesh file
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void read_input(char *input_dir, char *mesh_dir, double *x_vals_mat,
                double *y_vals_mat)
{
    FILE *input_fp = fopen(input_dir, "rt");
    FILE *mesh_fp = fopen(mesh_dir, "rt");

    if (!input_fp) {
        fprintf(stderr, "%s:%d: [Error] failed opening input file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }
    if (!mesh_fp) {
        fprintf(stderr, "%s:%d: [Error] failed opening mesh file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    read_mesh_file(mesh_fp, x_vals_mat, y_vals_mat);
    read_input_file(input_fp);

    fclose(input_fp);
    fclose(mesh_fp);
}

/* reading the mesh file
argument list:
mesh_dir - the directory of the mesh file
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void read_mesh_file(FILE *mesh_fp, double *x_vals_mat,
                   double *y_vals_mat)
{
    char current_word[MAXWORD];

    /* Seting the input varibles according to the mesh file */
    while(fscanf(mesh_fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "x_vals")) {
            read_mat_from_file(mesh_fp, x_vals_mat);
        } else if (!strcmp(current_word, "y_vals")) {
            read_mat_from_file(mesh_fp, y_vals_mat);
        }
    }
}

/* read input parameters from input file
fp - file pointer */
void read_input_file(FILE *fp)
{
    char current_word[MAXWORD];
    float temp;

    while(fscanf(fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "Mach_inf")) {
            fscanf(fp, "%g", &temp);
            Mach_inf = (double)temp;
        } else if (!strcmp(current_word, "angle_of_attack_deg")) {
            fscanf(fp, "%g", &temp);
            angle_of_attack_deg = (double)temp;
            angle_of_attack_rad = PI*angle_of_attack_deg/180.0;
        } else if (!strcmp(current_word, "density")) {
            fscanf(fp, "%g", &temp);
            density = (double)temp;
        } else if (!strcmp(current_word, "environment_pressure")) {
            fscanf(fp, "%g", &temp);
            environment_pressure = (double)temp;
        } else if (!strcmp(current_word, "delta_t")) {
            fscanf(fp, "%g", &temp);
            delta_t = (double)temp;
        } else if (!strcmp(current_word, "Gamma")) {
            fscanf(fp, "%g", &temp);
            Gamma = (double)temp;
        } else if (!strcmp(current_word, "epse")) {
            fscanf(fp, "%g", &temp);
            epse = (double)temp;
            epsi = epse * 2;
        } else if (!strcmp(current_word, "i_TEL")) {
            fscanf(fp, "%d", &i_TEL);
        } else if (!strcmp(current_word, "i_LE")) {
            fscanf(fp, "%d", &i_LE);
        } else if (!strcmp(current_word, "i_TEU")) {
            fscanf(fp, "%d", &i_TEU);
        } else if (!strcmp(current_word, "j_TEL")) {
            fscanf(fp, "%d", &j_TEL);
        } else if (!strcmp(current_word, "j_LE")) {
            fscanf(fp, "%d", &j_LE);
        } else if (!strcmp(current_word, "j_TEU")) {
            fscanf(fp, "%d", &j_TEU);
        } else if (!strcmp(current_word, "iterations")) {
            fscanf(fp, "%d", &max_iterations);
        }
    }
}

/* read matrix from file into a 2D matrix (1D array)
fp - file pointer
des - 1D array */
void read_mat_from_file(FILE *fp, double *des)
{
    float temp;

    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            fscanf(fp, "%g", &temp);
            des[offset2d(i, j, ni)] = (double)temp;
        }
    }
}

void output_solution(double *current_Q, double *U_mat, double *V_mat,
                     double *x_vals_mat, double *y_vals_mat,
                     double *dxi_dx_mat, double *dxi_dy_mat,
                     double *deta_dx_mat, double *deta_dy_mat)
{
    int i, j, index;
    double U, V;

    for (i = 0; i < ni; i++) {
        for (j = 0; j < nj; j++) {
            contravariant_velocities(&U, &V, dxi_dx_mat,
                                     dxi_dy_mat, deta_dx_mat,
                                     deta_dy_mat, current_Q, i, j);
            index = offset2d(i, j, ni);
            U_mat[index] = U;
            V_mat[index] = V;
        }
    }
    FILE *rho_u_fp = fopen("./results/rho_u.txt", "wt");
    FILE *rho_v_fp = fopen("./results/rho_v.txt", "wt");
    FILE *x_fp = fopen("./results/x_mat.txt", "wt");
    FILE *y_fp = fopen("./results/y_mat.txt", "wt");
    FILE *U_fp = fopen("./results/U_mat.txt", "wt");
    FILE *V_fp = fopen("./results/V_mat.txt", "wt");
    output_layer_of_mat3D_to_file(rho_u_fp, current_Q, 1);
    output_layer_of_mat3D_to_file(rho_v_fp, current_Q, 2);
    output_mat2D_to_file(x_fp, x_vals_mat);
    output_mat2D_to_file(y_fp, y_vals_mat);
    output_mat2D_to_file(U_fp, U_mat);
    output_mat2D_to_file(V_fp, V_mat);
}

/* converts a 2D index into 1D index
argument list:
i - first direction
j - second direction
ni - first direction size */
int offset2d(int i, int j, int ni)
{
    return j * ni + i;
}

/* converts a 3D index into 1D index
argument list:
i - first direction
j - second direction
k - third direction
ni - first direction size
nj - second direction size */
int offset3d(int i, int j, int k, int ni, int nj)
{
    return (k * nj + j) * ni + i;
}

void print_mat2D(double *data)
{
    int j_index, i_index;

    for (j_index = nj - 1; j_index >= 0; j_index--) {
        for (i_index = 0; i_index < ni; i_index++) {
                if (i_index == i_LE) {
                    printf("  ");
                }
            printf("%g ", data[offset2d(i_index, j_index, ni)]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_layer_of_mat3D(double *data, int layer)
{
    int j_index, i_index;

    for (j_index = nj - 1; j_index >= 0; j_index--) {
        for (i_index = 0; i_index < ni; i_index++) {
            printf("%g ", data[offset3d(i_index, j_index, layer, ni, nj)]);
        }
        printf("\n");
    }
    printf("\n");
}

/* return the second order first derivitive from the valuse in the mat matrix
argument list:
mat - 1D array of valuse
diraction - i or j
i, j - the points coordinates */
double first_deriv(double *mat, char diraction, int i, int j)
{
    int j_min = 0, j_max = nj-1, i_min = 0, i_max = ni-1;

    if (diraction == 'j') {
        if (j == j_min) {
            return mat[offset2d(i, j+1, ni)] - mat[offset2d(i, j, ni)]; /* (forward) first order first derivitive */
        } else if (j == j_max) {
            return mat[offset2d(i, j, ni)] - mat[offset2d(i, j-1, ni)]; /* (backward) first order first derivitive */
        }
        return (mat[offset2d(i, j+1, ni)] - mat[offset2d(i, j-1, ni)]) / (2); /* (central) second order first derivitive */
    }
    if (diraction == 'i') {
        if (i == i_min) {
            return mat[offset2d(i+1, j, ni)] - mat[offset2d(i, j, ni)]; /* (forward) first order first derivitive */
        } else if (i == i_max) {
            return mat[offset2d(i, j, ni)] - mat[offset2d(i-1, j, ni)]; /* (backward) first order first derivitive */
        } else {
            return (mat[offset2d(i+1, j, ni)] - mat[offset2d(i-1, j, ni)]) / (2); /* (central) second order first derivitive */
        }
    }
    return NAN;
}

/* calculating the jacobian in a single point
argument list:
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse 
i, j - the points coordinates */
double calculate_one_over_jacobian_at_a_point(double *x_vals_mat,
                                              double *y_vals_mat, int i,
                                              int j)
{
    double dx_dxi, dx_deta, dy_dxi, dy_deta;

    dx_dxi = first_deriv(x_vals_mat, 'i', i, j);
    dx_deta = first_deriv(x_vals_mat, 'j', i, j);
    dy_dxi = first_deriv(y_vals_mat, 'i', i, j);
    dy_deta = first_deriv(y_vals_mat, 'j', i, j);

    return (dx_dxi*dy_deta - dy_dxi*dx_deta);
}

/* calculating the contravariant velocities in a single point
argument list:
i, j - the points coordinates */
void contravariant_velocities(double *U, double *V,
                              double *dxi_dx_mat, double *dxi_dy_mat,
                              double *deta_dx_mat, double *deta_dy_mat,
                              double *Q, int i, int j)
{
    double dxi_dx, dxi_dy, deta_dx, deta_dy, u, v;
    int index = offset2d(i, j, ni);

    dxi_dx = dxi_dx_mat[index];
    dxi_dy = dxi_dy_mat[index];
    deta_dx = deta_dx_mat[index];
    deta_dy = deta_dy_mat[index];

    calculate_u_and_v(&u, &v, Q, i, j);

    *U = dxi_dx  * u + dxi_dy  * v;
    *V = deta_dx * u + deta_dy * v;
}

void calculate_u_and_v(double *u, double *v, double *Q, int i, int j)
{
    *u = Q[offset3d(i, j, 1, ni, nj)] / Q[offset3d(i, j, 0, ni, nj)]; /* rho*u / rho */
    *v = Q[offset3d(i, j, 2, ni, nj)] / Q[offset3d(i, j, 0, ni, nj)]; /* rho*v / rho */
}

double calculate_p(double energy, double rho, double u, double v)
{
    return (Gamma - 1) * (energy - 0.5 * rho * (u * u + v * v));
}

double calculate_energy(double p, double u, double v, double rho)
{
    double internal_energy, energy;

    internal_energy = p / (Gamma - 1) / rho;
    energy = rho * internal_energy + rho * (u * u + v * v) / 2;
    return energy;
}

void calculate_E_hat_at_a_point(double *E0, double *E1, double *E2,
                                double *E3, double *J_vals_mat,
                                double *dxi_dx_mat, double *dxi_dy_mat,
                                double *deta_dx_mat, double *deta_dy_mat,
                                double *Q, int i, int j)
{
    double u, v, U, V, dxi_dx, dxi_dy, energy, p, rho, J;
    int index = offset2d(i, j, ni);

    calculate_u_and_v(&u, &v, Q, i, j);
    contravariant_velocities(&U, &V, dxi_dx_mat, dxi_dy_mat, deta_dx_mat,
                             deta_dy_mat, Q, i, j);
    energy = Q[offset3d(i, j, 3, ni, nj)];
    rho = Q[offset3d(i, j, 0, ni, nj)];
    p = calculate_p(energy, rho, u, v);

    dxi_dx = dxi_dx_mat[index];
    dxi_dy = dxi_dy_mat[index];
    J = J_vals_mat[index];

    // if (!one_over_J) {
    //     *E0 = 0;
    //     *E1 = 0;
    //     *E2 = 0;
    //     *E3 = 0;
    // } else {
        *E0 = (1.0 / J) * rho * U;
        *E1 = (1.0 / J) * (rho * u * U + dxi_dx * p); 
        *E2 = (1.0 / J) * (rho * v * U + dxi_dy * p);
        *E3 = (1.0 / J) * (energy + p) * U;
    // }

}

void calculate_F_hat_at_a_point(double *F0, double *F1, double *F2,
                                double *F3, double *J_vals_mat,
                                double *dxi_dx_mat, double *dxi_dy_mat,
                                double *deta_dx_mat, double *deta_dy_mat,
                                double *Q, int i, int j)                               
{
    double u, v, U, V, J, deta_dx, deta_dy,
    energy, p, rho;
    int index = offset2d(i, j, ni);

    calculate_u_and_v(&u, &v, Q, i, j);
    contravariant_velocities(&U, &V, dxi_dx_mat, dxi_dy_mat, deta_dx_mat,
                             deta_dy_mat, Q, i, j);  
    energy = Q[offset3d(i, j, 3, ni, nj)];
    rho = Q[offset3d(i, j, 0, ni, nj)];
    p = calculate_p(energy, rho, u, v);

    J = J_vals_mat[index];
    deta_dx = deta_dx_mat[index];
    deta_dy = deta_dy_mat[index];

    // if (!one_over_J) {
    //     *F0 = 0;
    //     *F1 = 0;
    //     *F2 = 0;
    //     *F3 = 0;
    // } else {
        *F0 = (1.0 / J) * rho * V;
        *F1 = (1.0 / J) * (rho * u * V + deta_dx * p);
        *F2 = (1.0 / J) * (rho * v * V + deta_dy * p);
        *F3 = (1.0 / J) * (energy + p) * V;
    // }

}

void initialize_flow_field(double *Q)
{
    double u, v, energy, p, rho, speed_of_sound, velocity;

    p = environment_pressure;
    rho = density;
    speed_of_sound = sqrt(Gamma * p / rho);

    velocity = Mach_inf * speed_of_sound;
    u = velocity * cos(angle_of_attack_rad);
    v = velocity * sin(angle_of_attack_rad);

    energy = calculate_energy(p, u, v, rho);

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            Q[offset3d(i, j, 0, ni, nj)] = rho;
            Q[offset3d(i, j, 1, ni, nj)] = rho * u;
            Q[offset3d(i, j, 2, ni, nj)] = rho * v;
            Q[offset3d(i, j, 3, ni, nj)] = energy;
        }
    }
}

void RHS(double *S, double *W, double *Q, double *J_vals_mat,
         double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat,
         double *deta_dy_mat, double *s2, double *rspec, double *qv, double *dd)
{
int i, j, k;

    /* zeroing S and W*/
    for (i = 0; i < ni; i++) {   
        for (j = 0; j < nj; j++) {
            for (k = 0; k < 4; k++) {
                S[offset3d(i, j, k, ni, nj)] = 0;
            }
        }
    }
    for (i = 0; i < max_ni_nj; i++) {
        for (j = 0; j < 4; j++) {
            W[offset2d(i, j, max_ni_nj)] = 0;
        }
    }

    /* xi direction (constant j) */
    for (j = 1; j < nj - 1; j++) {
        for (i = 0; i < ni; i++) {
            calculate_E_hat_at_a_point(&W[offset2d(i, 0, max_ni_nj)],
                                       &W[offset2d(i, 1, max_ni_nj)],
                                       &W[offset2d(i, 2, max_ni_nj)],
                                       &W[offset2d(i, 3, max_ni_nj)], J_vals_mat,
                                       dxi_dx_mat, dxi_dy_mat, deta_dx_mat,
                                       deta_dy_mat, Q, i, j);
        }
        for (i = 1; i < ni - 1; i++) {
            for (k = 0; k < 4; k++) {
                S[offset3d(i, j, k, ni, nj)] += -delta_t * 0.5 * (W[offset2d(i+1, k, max_ni_nj)] - W[offset2d(i-1, k, max_ni_nj)]);
            }
        }
    } 

    /* eta direction (constant i) */
    for (i = 1; i < ni - 1; i++) {
        for (j = 0; j < nj; j++) {
            calculate_F_hat_at_a_point(&W[offset2d(j, 0, max_ni_nj)],
                                       &W[offset2d(j, 1, max_ni_nj)],
                                       &W[offset2d(j, 2, max_ni_nj)],
                                       &W[offset2d(j, 3, max_ni_nj)], J_vals_mat,
                                       dxi_dx_mat, dxi_dy_mat, deta_dx_mat,
                                       deta_dy_mat, Q, i, j);
        }
        for (j = 1; j < nj - 1; j++) {
            for (k = 0; k < 4; k++) {
                S[offset3d(i, j, k, ni, nj)] += -delta_t * 0.5 * (W[offset2d(j+1, k, max_ni_nj)] - W[offset2d(j-1, k, max_ni_nj)]);
            }
        }
    }

     if (smooth(Q, S, J_vals_mat, dxi_dx_mat, dxi_dy_mat, deta_dx_mat,
               deta_dy_mat, ni, nj, s2, rspec, qv, dd, epse, Gamma,
               Mach_inf, delta_t)) {
                /* returns zero on success */
                fprintf(stderr, "%s:%d: [Erorr] problem with smooth in RHS\n", __FILE__, __LINE__);
                exit(1);
               }
}

void advance_Q(double *next_Q, double *current_Q ,double *S, double *J_vals_mat)
{
    int index;
    double J;

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            J = J_vals_mat[offset2d(i, j, ni)];
            for (int k = 0; k < 4; k++ ) {
                index = offset3d(i, j, k, ni, nj);
                if (S[index]) {
                    next_Q[index] = current_Q[index] + S[index] * J;
                } else {
                    next_Q[index] = current_Q[index];
                }
            }
        }
    }
}

void copy_3Dmat_to_3Dmat(double *dst, double *src)
{
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < 4; k++) {
                int index = offset3d(i, j, k, ni, nj);
                dst[index] = src[index];
            }
        }
    }
}

int smooth(double *q, double *s, double *jac, double *xx, double *xy, 
           double *yx, double *yy, int id, int jd, double *s2, 
           double *rspec, double *qv, double *dd,
           double epse, double gamma, double fsmach, double dt)
{
    double *rho, *u_vel, *v_vel, *t_e;
    double eratio, smool, gm1, ggm1, cx, cy, eps, ra, u, v, qq, ss, st, 
          qav, qxx, ssfs, qyy;
    int ib, ie, jb, je, i, j, offset, offsetp1, offsetm1, ip, ir, n,
        jp, jr;

    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    smool = 1.0;
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    ib = 1;
    ie = id - 1;
    jb = 1;
    je = jd - 1;

    cx = 2.;
    cy = 1.;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*     smoothing in xi direction */

    for (j = jb; j < je; j++) {
	for (i = 0; i < id; i++) {
	    offset = id * j + i;
	    eps = epse / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq);
            rspec[i] = eps * (fabs(xx[offset] * u + xy[offset] * v) + sqrt((xx[offset] * xx[offset] + xy[offset] * xy[offset]) * ss + 0.01));
	    qv[i] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (i = ib; i < ie; i++) {
	    ip = i + 1;
	    ir = i - 1;
	    qxx = qv[ip] - qv[i] * 2. + qv[ir];
	    qav = (qv[ip] + qv[i] * 2. + qv[ir]) * .25;
	    dd[i] = eratio * fabs(qxx / qav);
	}

	dd[0] = dd[1];
	dd[id - 1] = dd[id - 2];

	for (n = 0; n < 4; n++) {
	    for (i = ib; i < ie; i++) {
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j) * id + i + 1;
		offsetm1 = (jd * n + j) * id + i - 1;
		s2[i] = q[offsetp1] - 2.0 * q[offset] + q[offsetm1];
	    }

	    s2[0] = s2[1] * -1.;
	    s2[id - 1] = s2[id - 2] * -1.;

	    for (i = ib; i < ie; i++) {
		ip = i + 1;
		ir = i - 1;
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j) * id + i + 1;
		offsetm1 = (jd * n + j) * id + i - 1;
		st = ((dd[ip] + dd[i]) * .5 * (q[offsetp1] - q[offset]) - cx / (cx + dd[ip] + dd[i]) * (s2[ip] - s2[i])) * (rspec[ip] + rspec[i]) + ((dd[ir] + dd[i]) * .5 * (q[offsetm1] - q[offset]) - cx / (cx + dd[ir] + dd[i]) * (s2[ir] - s2[i])) * (rspec[ir] + rspec[i]);
		s[offset] += st * .5 * dt;
	    }
	}
    }

/*     smoothing in eta direction */

    ssfs = 1. / (0.001 + fsmach * fsmach);
    for (i = ib; i < ie; i++) {
	for (j = 0; j < jd; j++) {
	    offset = id * j + i;
	    eps = epse / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq) * (1.0 - smool) + smool * qq * ssfs;
            rspec[j] = eps * (fabs(yx[offset] * u + yy[offset] * v) + sqrt((yx[offset] * yx[offset] + yy[offset] * yy[offset]) * ss + 0.01));
	    qv[j] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (j = jb; j < je; j++) {
	    jp = j + 1;
	    jr = j - 1;
	    qyy = qv[jp] - qv[j] * 2. + qv[jr];
	    qav = (qv[jp] + qv[j] * 2. + qv[jr]) * .25;
	    dd[j] = eratio * fabs(qyy / qav);
	}

	dd[0] = dd[1];
	dd[jd - 1] = dd[jd - 2];

	for (n = 0; n < 4; n++) {
	    for (j = jb; j < je; j++) {
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j + 1) * id + i;
		offsetm1 = (jd * n + j - 1) * id + i;
		s2[j] = q[offsetp1] - 2.0 * q[offset] + q[offsetm1];
	    }

	    s2[0] = s2[1] * -1.;
	    s2[jd - 1] = s2[jd - 2] * -1.;

	    for (j = jb; j < je; j++) {
		jp = j + 1;
		jr = j - 1;
		offset = (jd * n + j) * id + i;
		offsetp1 = (jd * n + j + 1) * id + i;
		offsetm1 = (jd * n + j - 1) * id + i;
		st = ((dd[jp] + dd[j]) * .5 * (q[offsetp1] - q[offset]) - cy / (cy + dd[jp] + dd[j]) * (s2[jp] - s2[j])) * (rspec[jp] + rspec[j]) + ((dd[jr] + dd[j]) * .5 * (q[offsetm1] - q[offset]) - cy / (cy + dd[jr] + dd[j]) * (s2[jr] - s2[j])) * (rspec[jr] + rspec[j]);
		s[offset] += st * .5 * dt;
	    }
	}
    }

    return 0;
}

void apply_BC(double *Q, double *J_vals_mat, double *dxi_dx_mat,
              double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat)
{
    int i, k;
    double J_j0, u_j1, v_j1, p_j1, e_j1, rho_j0, rho_j1, p_j0, e_j0, u_j0, v_j0,
    dxi_dx_j0, dxi_dx_j1, dxi_dy_j0, dxi_dy_j1, deta_dx_j0, deta_dy_j0,
    U1, V1,
    e_iTEU_jTEU, rho_iTEU_jTEU, u_iTEU_jTEU, v_iTEU_jTEU, p_iTEU_jTEU,
    e_iTEL_jTEL, rho_iTEL_jTEL, u_iTEL_jTEL, v_iTEL_jTEL, p_iTEL_jTEL,
    e_iTE_jTE, rho_iTE_jTE, u_iTE_jTE, v_iTE_jTE, p_iTE_jTE;
/* wall BC */
    for (i = i_TEL; i <= i_TEU; i++) {
        calculate_u_and_v(&u_j1, &v_j1, Q, i, 1);        
        contravariant_velocities(&U1, &V1, dxi_dx_mat, dxi_dy_mat, deta_dx_mat,
                                 deta_dy_mat, Q, i, 1);

        /* u_i,0 and v_i,0 */
        dxi_dx_j0 = dxi_dx_mat[offset2d(i, 0, ni)];
        dxi_dx_j1 = dxi_dx_mat[offset2d(i, 1, ni)]; 
        dxi_dy_j0 = dxi_dy_mat[offset2d(i, 0, ni)];
        dxi_dy_j1 = dxi_dy_mat[offset2d(i, 1, ni)];
        deta_dx_j0 = deta_dx_mat[offset2d(i, 0, ni)];
        deta_dy_j0 = deta_dy_mat[offset2d(i, 0, ni)];
        J_j0 = J_vals_mat[offset2d(i, 0, ni)];

        // u_j0 = (dxi_dx_j1 * u_j1 + dxi_dy_j1 * v_j1) / (dxi_dx_j0 - dxi_dy_j0 * deta_dx_j0 / deta_dy_j0);
        // v_j0 = -(dxi_dx_j1 * u_j1 + dxi_dy_j1 * v_j1) / (deta_dy_j0 / deta_dx_j0 * dxi_dx_j0 - dxi_dy_j0);
        /*test*/
        u_j0 = U1 * deta_dy_j0 / J_j0;
        v_j0 = - U1 * deta_dx_j0 / J_j0;
        /*test*/

        /* rho_i,0 */
        rho_j1 = Q[offset3d(i, 1, 0, ni, nj)];
        rho_j0 = rho_j1;
        
        /* p_i,0 */
        e_j1 = Q[offset3d(i, 1, 3, ni, nj)];
        p_j1 = calculate_p(e_j1, rho_j1, u_j1, v_j1);
        p_j0 = p_j1;

        /* e_i,0*/
        e_j0 = p_j0 / (Gamma -1) + 0.5 * rho_j0 * (u_j0 * u_j0 + v_j0 *v_j0);
        p_j0 = calculate_p(e_j0, rho_j0, u_j0, v_j0);
        // printf("%d,%g, %g, %g, %g\n", i,p_j1, p_j0, e_j1, e_j0);
        // printf("%g, %g, %g, %g, %g, %g, %g\n", dxi_dx_j0, dxi_dx_j1, dxi_dy_j0, dxi_dy_j1, deta_dx_j0, deta_dy_j0, J_j0);

        Q[offset3d(i, 0, 0, ni, nj)] = rho_j0;
        Q[offset3d(i, 0, 1, ni, nj)] = rho_j0 * u_j0;
        Q[offset3d(i, 0, 2, ni, nj)] = rho_j0 * v_j0;
        Q[offset3d(i, 0, 3, ni, nj)] = e_j0;
    }

/* trailing edge */
    calculate_u_and_v(&u_iTEU_jTEU, &v_iTEU_jTEU, Q, i_TEU, 0);
    calculate_u_and_v(&u_iTEL_jTEL, &v_iTEL_jTEL, Q, i_TEL, 0);

    rho_iTEU_jTEU = Q[offset3d(i_TEU, 0, 0, ni, nj)];
    rho_iTEL_jTEL = Q[offset3d(i_TEL, 0, 0, ni, nj)];

    e_iTEU_jTEU = Q[offset3d(i_TEU, 0, 3, ni, nj)];
    e_iTEL_jTEL = Q[offset3d(i_TEL, 0, 3, ni, nj)];

    p_iTEU_jTEU = calculate_p(e_iTEU_jTEU, rho_iTEU_jTEU, u_iTEU_jTEU, v_iTEU_jTEU);
    p_iTEL_jTEL = calculate_p(e_iTEL_jTEL, rho_iTEL_jTEL, u_iTEL_jTEL, v_iTEL_jTEL);
    
    p_iTE_jTE   = 0.5 * (p_iTEU_jTEU   + p_iTEL_jTEL);
    u_iTE_jTE   = 0.5 * (u_iTEU_jTEU   + u_iTEL_jTEL);
    v_iTE_jTE   = 0.5 * (v_iTEU_jTEU   + v_iTEL_jTEL);
    rho_iTE_jTE = 0.5 * (rho_iTEU_jTEU + rho_iTEL_jTEL);
    e_iTE_jTE   = calculate_energy(p_iTE_jTE, u_iTE_jTE, v_iTE_jTE, rho_iTE_jTE);

    Q[offset3d(i_TEL, 0, 0, ni, nj)] = rho_iTE_jTE;
    Q[offset3d(i_TEU, 0, 0, ni, nj)] = rho_iTE_jTE;
    Q[offset3d(i_TEU, 0, 1, ni, nj)] = rho_iTE_jTE * u_iTE_jTE;
    Q[offset3d(i_TEL, 0, 1, ni, nj)] = rho_iTE_jTE * u_iTE_jTE;
    Q[offset3d(i_TEL, 0, 2, ni, nj)] = rho_iTE_jTE * v_iTE_jTE;
    Q[offset3d(i_TEU, 0, 2, ni, nj)] = rho_iTE_jTE * v_iTE_jTE;
    Q[offset3d(i_TEL, 0, 3, ni, nj)] = e_iTE_jTE;
    Q[offset3d(i_TEU, 0, 3, ni, nj)] = e_iTE_jTE;

/* wake BC */
    for (i = 1; i < i_TEL; i++) {
        for (k = 0; k < 4; k++) {
            Q[offset3d(i, 0, k, ni, nj)] = 0.5 * (Q[offset3d(i, 1, k, ni, nj)] + Q[offset3d(ni-1-i, 1, k, ni, nj)]);
            Q[offset3d(ni-1-i, 0, k, ni, nj)] =  Q[offset3d(i, 0, k, ni, nj)];
        }
    }

/* outflow BC */
    for (int j = 1; j < nj; j++) {
        for (k = 0; k < 4; k++) {
            Q[offset3d(0, j, k, ni, nj)]    = Q[offset3d(1, j, k, ni, nj)];
            Q[offset3d(ni-1, j, k, ni, nj)] = Q[offset3d(ni-1-1, j, k, ni, nj)];
        }
    }

}

/* output data;
argument list:
fp - file pointer to output file
data - 1D array of the output valuse */
void output_mat2D_to_file(FILE *fp, double *data)
{
    int i, j, j_max = nj - 1, i_max = ni - 1;
    
    for (j = 0; j < j_max+1; j++) {
        for (i = 0; i < i_max+1; i++) {
            fprintf(fp, "%g ", data[offset2d(i, j, i_max+1)]);
        }
        fprintf(fp, "\n");
    }
}

/* output data;
argument list:
fp - file pointer to output file
data - 1D array of 3D matrix
layer - the k value */
void output_layer_of_mat3D_to_file(FILE *fp, double *data, int layer)
{
    int i, j, j_max = nj - 1, i_max = ni - 1;
    
    for (j = 0; j < j_max+1; j++) {
        for (i = 0; i < i_max+1; i++) {
            fprintf(fp, "%g ", data[offset3d(i, j, layer, i_max+1, j_max+1)]);
        }
        fprintf(fp, "\n");
    }
}

void calculate_A_hat_j_const(double *dst, double *Q, double *dxi_dx_mat,
                             double *dxi_dy_mat, int i, int j)
{
    /* xi (i) direction */
    double u, v, phi_square, theta, gamma1, gamma2, beta, energy, rho,
    dxi_dx, dxi_dy;


    calculate_u_and_v(&u, &v, Q, i, j);

    dxi_dx = dxi_dx_mat[offset2d(i, j, ni)];
    dxi_dy = dxi_dy_mat[offset2d(i, j, ni)];

    rho = Q[offset3d(i, j, 0, ni, nj)];
    energy = Q[offset3d(i, j, 3, ni, nj)];

    phi_square = 0.5 * (Gamma - 1) * (u * u + v * v);
    theta = dxi_dx * u + dxi_dy * v;
    gamma1 = Gamma - 1;
    gamma2 = Gamma - 2;
    beta = Gamma * energy / rho - phi_square;

    dst[offset3d(i, 0, 0, ni, 4)] = 0;
    dst[offset3d(i, 0, 1, ni, 4)] = dxi_dx;
    dst[offset3d(i, 0, 2, ni, 4)] = dxi_dy;
    dst[offset3d(i, 0, 3, ni, 4)] = 0;

    dst[offset3d(i, 1, 0, ni, 4)] = dxi_dx * phi_square - u * theta;
    dst[offset3d(i, 1, 1, ni, 4)] = theta - dxi_dx * gamma2 * u;
    dst[offset3d(i, 1, 2, ni, 4)] = dxi_dy * u - gamma1 * dxi_dx * v;
    dst[offset3d(i, 1, 3, ni, 4)] = dxi_dx * gamma1;

    dst[offset3d(i, 2, 0, ni, 4)] = dxi_dy * phi_square - v * theta;
    dst[offset3d(i, 2, 1, ni, 4)] = dxi_dx * v - dxi_dy * gamma1 * u;
    dst[offset3d(i, 2, 2, ni, 4)] = theta - dxi_dy * gamma2 * v;
    dst[offset3d(i, 2, 3, ni, 4)] = dxi_dy * gamma1;

    dst[offset3d(i, 3, 0, ni, 4)] = theta * (2 * phi_square - Gamma * energy / rho);
    dst[offset3d(i, 3, 1, ni, 4)] = dxi_dx * beta - gamma1 * u * theta;
    dst[offset3d(i, 3, 2, ni, 4)] = dxi_dy * beta - gamma1 * v * theta;
    dst[offset3d(i, 3, 3, ni, 4)] = Gamma * theta;

}

void calculate_B_hat_i_const(double *dst, double *Q, double *deta_dx_mat,
                             double *deta_dy_mat, int i, int j)
{
    /* eta (j) direction */
    double u, v, phi_square, theta, gamma1, gamma2, beta, energy, rho,
    deta_dx, deta_dy;


    calculate_u_and_v(&u, &v, Q, i, j);

    deta_dx = deta_dx_mat[offset2d(i, j, ni)];
    deta_dy = deta_dy_mat[offset2d(i, j, ni)];

    rho = Q[offset3d(i, j, 0, ni, nj)];
    energy = Q[offset3d(i, j, 3, ni, nj)];

    phi_square = 0.5 * (Gamma - 1) * (u * u + v * v);
    theta = deta_dx * u + deta_dy * v;
    gamma1 = Gamma - 1;
    gamma2 = Gamma - 2;
    beta = Gamma * energy / rho - phi_square;

    dst[offset3d(j, 0, 0, nj, 4)] = 0;
    dst[offset3d(j, 0, 1, nj, 4)] = deta_dx;
    dst[offset3d(j, 0, 2, nj, 4)] = deta_dy;
    dst[offset3d(j, 0, 3, nj, 4)] = 0;

    dst[offset3d(j, 1, 0, nj, 4)] = deta_dx * phi_square - u * theta;
    dst[offset3d(j, 1, 1, nj, 4)] = theta - deta_dx * gamma2 * u;
    dst[offset3d(j, 1, 2, nj, 4)] = deta_dy * u - gamma1 * deta_dx * v;
    dst[offset3d(j, 1, 3, nj, 4)] = deta_dx * gamma1;

    dst[offset3d(j, 2, 0, nj, 4)] = deta_dy * phi_square - v * theta;
    dst[offset3d(j, 2, 1, nj, 4)] = deta_dx * v - deta_dy * gamma1 * u;
    dst[offset3d(j, 2, 2, nj, 4)] = theta - deta_dy * gamma2 * v;
    dst[offset3d(j, 2, 3, nj, 4)] = deta_dy * gamma1;

    dst[offset3d(j, 3, 0, nj, 4)] = theta * (2 * phi_square - Gamma * energy / rho);
    dst[offset3d(j, 3, 1, nj, 4)] = deta_dx * beta - gamma1 * u * theta;
    dst[offset3d(j, 3, 2, nj, 4)] = deta_dy * beta - gamma1 * v * theta;
    dst[offset3d(j, 3, 3, nj, 4)] = Gamma * theta;
}

int smoothx(double *q, double *xx, double *xy, int id, int jd, double *a,
	        double *b, double *c, int j,double *jac, double *drr, double *drp, 
            double *rspec, double *qv, double *dd,
            double epsi, double gamma, double fsmach, double dt)
{

    double *rho, *u_vel, *v_vel, *t_e;

    double eratio, gm1, ggm1, eps, ra, u, v, qq, ss,
          qav, qxx, rr, rp;
    int ib, ie, i, offset, offsetp1, offsetm1, ip, ir, n;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    ib = 1;
    ie = id - 1;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*  smoothing in xi direction */

        for (i = 0; i < id; i++) {
	    offset = id * j + i;
	    eps = epsi / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq);
            rspec[i] = eps * (fabs(xx[offset] * u + xy[offset] * v) + sqrt((xx[offset] * xx[offset] + xy[offset] * xy[offset]) * ss + 0.01));
	    qv[i] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (i = ib; i < ie; i++) {
	    ip = i + 1;
	    ir = i - 1;
	    qxx = qv[ip] - qv[i] * 2. + qv[ir];
	    qav = (qv[ip] + qv[i] * 2. + qv[ir]) * .25;
	    dd[i] = eratio * fabs(qxx / qav);
	}

	dd[0] = dd[1];
	dd[id - 1] = dd[id - 2];

        for (i = ib; i < ie; i++) {
	    ip = i + 1;
	    ir = i - 1;
	    offset = j * id + i;
	    offsetp1 = j * id + i + 1;
	    offsetm1 = j * id + i - 1;
	    rp = (0.5 * (dd[ip] + dd[i]) + 2.5) * dt * 0.5 * (rspec[ip] + rspec[i]);
	    rr = (0.5 * (dd[ir] + dd[i]) + 2.5) * dt * 0.5 * (rspec[ir] + rspec[i]);
	    qv[i] = (rr + rp) * jac[offset];
	    drr[i] = rr * jac[offsetm1];
	    drp[i] = rp * jac[offsetp1];
	}

	for (n = 0; n < 4; n++) {
	    for (i = ib; i < ie; i++) {
	        offset = (n * 4 + n) * id + i;
		a[offset] -= drr[i];
		b[offset] += qv[i];
		c[offset] -= drp[i];
	    }
        }
    return 0;
}

int smoothy(double *q, double *yx, double *yy, int id, int jd, double *a,
	        double *b, double *c, int i,double *jac, double *drr, double *drp, 
            double *rspec, double *qv, double *dd,
            double epsi, double gamma, double fsmach, double dt)
{

    double *rho, *u_vel, *v_vel, *t_e;

    double eratio, smool, gm1, ggm1, eps, ra, u, v, qq, ss, 
          qav, ssfs, qyy, rp, rr;
    int jb, je, j, offset, offsetp1, offsetm1, n,
        jp, jr;
    eratio = 0.25 + 0.25 * pow(fsmach + 0.0001,gamma);
    smool = 1.0;
    gm1 = gamma - 1.0;
    ggm1 = gamma * gm1;
    jb = 1;
    je = jd - 1;

    rho = q;
    u_vel = &q[id*jd];
    v_vel = &q[2*id*jd];
    t_e = &q[3*id*jd];

/*     smoothing in eta direction */

        ssfs = 1. / (0.001 + fsmach * fsmach);
        for (j = 0; j < jd; j++) {
	    offset = id * j + i;
	    eps = epsi / jac[offset];
	    ra = 1. / rho[offset];
	    u = u_vel[offset] * ra;
	    v = v_vel[offset] * ra;
	    qq = u * u + v * v;
	    ss = ggm1 * (t_e[offset] * ra - 0.5 * qq) * (1.0 - smool) + smool * qq * ssfs;
            rspec[j] = eps * (fabs(yx[offset] * u + yy[offset] * v) + sqrt((yx[offset] * yx[offset] + yy[offset] * yy[offset]) * ss + 0.01));
	    qv[j] = gm1 * (t_e[offset] - 0.5 * qq * rho[offset]);
	}

	for (j = jb; j < je; j++) {
	    jp = j + 1;
	    jr = j - 1;
	    qyy = qv[jp] - qv[j] * 2. + qv[jr];
	    qav = (qv[jp] + qv[j] * 2. + qv[jr]) * .25;
	    dd[j] = eratio * fabs(qyy / qav);
	}

	dd[0] = dd[1];
	dd[jd - 1] = dd[jd - 2];

        for (j = jb; j < je; j++) {
	    jp = j + 1;
	    jr = j - 1;
	    offset = j * id + i;
	    offsetp1 = (j + 1) * id + i;
	    offsetm1 = (j - 1) * id + i;
	    rp = (0.5 * (dd[jp] + dd[j]) + 2.5) * dt * 0.5 * (rspec[jp] + rspec[j]);
	    rr = (0.5 * (dd[jr] + dd[j]) + 2.5) * dt * 0.5 * (rspec[jr] + rspec[j]);
	    qv[j] = (rr + rp) * jac[offset];
	    drr[j] = rr * jac[offsetm1];
	    drp[j] = rp * jac[offsetp1];
	}

	for (n = 0; n < 4; n++) {
	    for (j = jb; j < je; j++) {
	        offset = (n * 4 + n) * jd + j;
		a[offset] -= drr[j];
		b[offset] += qv[j];
		c[offset] -= drp[j];
	    }
        }
    return 0;
}

void LHSX(double *A, double *B, double *C, double *Q,
          double *dxi_dx_mat, double *dxi_dy_mat, int j)
{
    int i, n, m; 

    // for (i = 0; i < max_ni_nj; i++) {
    //     for (n = 0; n < 4; n++) {
    //         for (m = 0; m < 4; m++) {
    //                 A[offset3d(i, m, n, max_ni_nj, 4)] = 0;
    //                 B[offset3d(i, m, n, max_ni_nj, 4)] = 0;
    //                 C[offset3d(i, m, n, max_ni_nj, 4)] = 0;
    //         }
    //     }
    // }

    for (i = 0; i < ni; i++) {
        calculate_A_hat_j_const(B, Q, dxi_dx_mat, dxi_dy_mat, i, j);
    }
    for (i = 1; i < ni - 1; i++) {
        for (n = 0; n < 4; n++) {
            for (m = 0; m < 4; m++) {
                    A[offset3d(i, m, n, ni, 4)] = - delta_t * 0.5 * B[offset3d(i-1, m, n, ni, 4)];
                    C[offset3d(i, m, n, ni, 4)] =   delta_t * 0.5 * B[offset3d(i+1, m, n, ni, 4)];
            }
        }
    }
    /*test*/
    // printf("##########____C____##########\n");
    // for (i = 0; i < ni; i++) {
    //     for (m = 0; m < 4; m++) {
    //         for (n = 0; n < 4; n++) {
    //             printf("%g ", C[offset3d(i, m, n, ni, 4)]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    /*test*/
    for (i = 0; i < ni; i++) {
        for (n = 0; n < 4; n++) {
            for (m = 0; m < 4; m++) {
                if (n == m) {
                    B[offset3d(i, m, n, ni, 4)] = 1;
                } else {
                    B[offset3d(i, m, n, ni, 4)] = 0;
                }
            }
        }
    }



    // /*test*/
    // printf("##########____C____##########\n");
    // for (i = 0; i < ni; i++) {
    //     for (m = 0; m < 4; m++) {
    //         for (n = 0; n < 4; n++) {
    //             printf("%g ", C[offset3d(i, m, n, ni, 4)]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // /*test*/
    
}

void LHSY(double *A, double *B, double *C, double *Q,
          double *deta_dx_mat, double *deta_dy_mat, int i)
{
    int j, n, m; 

    // for (j = 0; j < max_ni_nj; j++) {
    //     for (n = 0; n < 4; n++) {
    //         for (m = 0; m < 4; m++) {
    //                 A[offset3d(j, m, n, max_ni_nj, 4)] = 0;
    //                 B[offset3d(j, m, n, max_ni_nj, 4)] = 0;
    //                 C[offset3d(j, m, n, max_ni_nj, 4)] = 0;
    //         }
    //     }
    // }

    for (j = 0; j < nj; j++) {
        calculate_B_hat_i_const(B, Q, deta_dx_mat, deta_dy_mat, i, j);
    }
    for (j = 1; j < nj - 1; j++) {
        for (n = 0; n < 4; n++) {
            for (m = 0; m < 4; m++) {
                    A[offset3d(j, m, n, nj, 4)] = - delta_t * 0.5 * B[offset3d(j-1, m, n, nj, 4)];
                    C[offset3d(j, m, n, nj, 4)] =   delta_t * 0.5 * B[offset3d(j+1, m, n, nj, 4)];
            }
        }
    }
    for (j = 0; j < nj; j++) {
        for (n = 0; n < 4; n++) {
            for (m = 0; m < 4; m++) {
                if (n == m) {
                    B[offset3d(j, m, n, nj, 4)] = 1;
                } else {
                    B[offset3d(j, m, n, nj, 4)] = 0;
                }
            }
        }
    }
}

int btri4s(double *a, double *b, double *c, double *f, int kd, int ks, int ke)
{
  /* Local variables */
  int k, m, n, nd, md;

  double c1, d1, d2, d3, d4, c2, c3, c4, b11, b21, b22, b31, b32, b33, 
    b41, b42, b43, b44, u12, u13, u14, u23, u24, u34;
  
  
  /*   (A,B,C)F = F, F and B are overloaded, solution in F */

  md = 4;
  nd = 4;

  /*   Part 1. Forward block sweep */
  
  for (k = ks; k <= ke; k++)
    {
      
      /*      Step 1. Construct L in B */
      
      if (k != ks) 
	{
	  for (m = 0; m < md; m++) 
	    {
	      for (n = 0; n < nd; n++) 
		{
		  b[k + kd * (m + md * n)] = b[k + kd * (m + md * n)] 
		    - a[k + kd * (m + md * 0)] * b[k - 1 + kd * (0 + md * n)] 
		    - a[k + kd * (m + md * 1)] * b[k - 1 + kd * (1 + md * n)] 
		    - a[k + kd * (m + md * 2)] * b[k - 1 + kd * (2 + md * n)] 
		    - a[k + kd * (m + md * 3)] * b[k - 1 + kd * (3 + md * n)] ;
		}
	    }
	}
      
      /*      Step 2. Compute L inverse (block matrix) */
      
      /*          A. Decompose L into L and U */
      
      b11 = 1. / b[k + kd * (0 + md * 0)];
      u12 = b[k + kd * (0 + md * 1)] * b11;
      u13 = b[k + kd * (0 + md * 2)] * b11;
      u14 = b[k + kd * (0 + md * 3)] * b11;
      b21 = b[k + kd * (1 + md * 0)];
      b22 = 1. / (b[k + kd * (1 + md * 1)] - b21 * u12);
      u23 = (b[k + kd * (1 + md * 2)] - b21 * u13) * b22;
      u24 = (b[k + kd * (1 + md * 3)] - b21 * u14) * b22;
      b31 = b[k + kd * (2 + md * 0)];
      b32 = b[k + kd * (2 + md * 1)] - b31 * u12;
      b33 = 1. / (b[k + kd * (2 + md * 2)] - b31 * u13 - b32 * u23);
      u34 = (b[k + kd * (2 + md * 3)] - b31 * u14 - b32 * u24) * b33;
      b41 = b[k + kd * (3 + md * 0)];
      b42 = b[k + kd * (3 + md * 1)] - b41 * u12;
      b43 = b[k + kd * (3 + md * 2)] - b41 * u13 - b42 * u23;
      b44 = 1. / (b[k + kd * (3 + md * 3)] - b41 * u14 - b42 * u24 
		  - b43 * u34);
      
      /*      Step 3. Solve for intermediate vector */
      
      /*          A. Construct RHS */
      if (k != ks) 
	{
	  for (m = 0; m < md; m++) 
	    {
	      f[k + kd * m] = f[k + kd * m] 
		- a[k + kd * (m + md * 0)] * f[k - 1 + kd * 0] 
		- a[k + kd * (m + md * 1)] * f[k - 1 + kd * 1] 
		- a[k + kd * (m + md * 2)] * f[k - 1 + kd * 2] 
		- a[k + kd * (m + md * 3)] * f[k - 1 + kd * 3];
	    }
	}
      
      /*          B. Intermediate vector */
      
      /*          Forward substitution */
      
      d1 = f[k + kd * 0] * b11;
      d2 = (f[k + kd * 1] - b21 * d1) * b22;
      d3 = (f[k + kd * 2] - b31 * d1 - b32 * d2) * b33;
      d4 = (f[k + kd * 3] - b41 * d1 - b42 * d2 - b43 * d3) * b44;
      
      /*          Backward substitution */
      
      f[k + kd * 3] = d4;
      f[k + kd * 2] = d3 - u34 * d4;
      f[k + kd * 1] = d2 - u23 * f[k + kd * 2] - u24 * d4;
      f[k + kd * 0] = d1 - u12 * f[k + kd * 1] - u13 * f[k + kd * 2] - u14 * d4;
      
      /*      Step 4. Construct U = L ** (-1) * C */
      /*              by columns and store in B */
      
      if (k != ke) 
	{
	  for (n = 0; n < nd; n++) 
	    {

	      /*          Forward substitution */
	      
	      c1 = c[k + kd * (0 + md * n)] * b11;
	      c2 = (c[k + kd * (1 + md * n)] - b21 * c1) * b22;
	      c3 = (c[k + kd * (2 + md * n)] - b31 * c1 - b32 * c2) * 
		b33;
	      c4 = (c[k + kd * (3 + md * n)] - b41 * c1 - b42 * c2 - 
		    b43 * c3) * b44;
	      
	      /*          Backward substitution */
	      
	      b[k + kd * (3 + md * n)] = c4;
	      b[k + kd * (2 + md * n)] = c3 - u34 * c4;
	      b[k + kd * (1 + md * n)] = c2 - u23 * b[k + kd * (2 + md * n)] - u24 * c4;
	      b[k + kd * (0 + md * n)] = c1 - u12 * b[k + kd * (1 + md * n)] 
		- u13 * b[k + kd * (2 + md * n)] - u14 * c4;
	    }
	}
    }
  
  /*   Part 2. Backward block sweep */
  
  if (ke == ks) 
    {
      return 0;
    }

  for (k = ke - 1; k >= ks; --k) 
    {
      for (m = 0; m < md; m++) 
	{
	  f[k + kd * m] = f[k + kd * m] 
	    - b[k + kd * (m + md * 0)] * f[k + 1 + kd * 0] 
	    - b[k + kd * (m + md * 1)] * f[k + 1 + kd * 1] 
	    - b[k + kd * (m + md * 2)] * f[k + 1 + kd * 2] 
	    - b[k + kd * (m + md * 3)] * f[k + 1 + kd * 3];
	}
    }
  
  return 0;
  
}

double calculate_S_norm(double *S)
{
    double sum = 0, value;

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            for (int k = 0; k < 4; k++) {
                value = S[offset3d(i, j, k, ni, nj)];
                sum += value * value;
                // dprintD(value);
            }
        }
    }
    // dprintD(sum);

    return sqrt(sum);
}

/* returns S_Norm */
double step(double *A, double *B, double *C, double *D, double *current_Q,
          double *S, double *W, double *J_vals_mat, double *dxi_dx_mat,
          double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat,
          double *s2, double *drr, double *drp, double *rspec, double *qv, double *dd)
{
    int i, j, k;

    RHS(S, W, current_Q, J_vals_mat, dxi_dx_mat, dxi_dy_mat, deta_dx_mat, deta_dy_mat,
        s2, rspec, qv, dd);

/* xi inversions */
    for (j = 1; j < nj - 1; j++) {
        LHSX(A, B, C, current_Q, dxi_dx_mat, dxi_dy_mat, j);
        /*test*/
        // if (j == 1) {
        //     dprintINT(j);
        //     printf("##########____C____##########\n");
        //     for (int i = 0; i < ni; i++) {
        //         dprintINT(i);
        //         for (int m = 0; m < 4; m++) {
        //             for (int n = 0; n < 4; n++) {
        //                 printf("%g ", C[offset3d(i, m, n, ni, 4)]);
        //             }
        //             printf("\n");
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }
        /*test*/
        if (smoothx(current_Q, dxi_dx_mat, dxi_dy_mat, ni, nj, A, B, C,
                    j, J_vals_mat, drr, drp, rspec, qv, dd, epsi, Gamma,
                    Mach_inf, delta_t)) {
                    /* returns zero on success */
                    fprintf(stderr, "%s:%d: [Erorr] problem with smoothx in LHSX\n", __FILE__, __LINE__);
                    exit(1);
                }
        for (k = 0; k < 4; k++) {
            for (i = 0; i < ni; i++) {
                D[offset2d(i, k, max_ni_nj)] = S[offset3d(i, j, k, ni, nj)];
            }
        }
        /*test*/
        // if (j == 1) {
        //     dprintINT(j);
        //     printf("##########____B____##########\n");
        //     for (int i = 0; i < ni; i++) {
        //         dprintINT(i);
        //         for (int m = 0; m < 4; m++) {
        //             for (int n = 0; n < 4; n++) {
        //                 printf("%g ", B[offset3d(i, m, n, ni, 4)]);
        //             }
        //             printf("\n");
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }
        /*test*/

        btri4s(A, B, C, D, ni, 1, ni - 2);

        /*test*/
        // if (j == 1) {
        //     dprintINT(j);
        //     for (i = 0; i < ni; i++) {
        //         for (k = 0; k < 4; k++) {
        //             printf("%g ", D[offset2d(i, k, max_ni_nj)]);
        //         }
        //         printf("\n");
        //     }
        // }
        /*test*/

        for (k = 0; k < 4; k++) {
            for (i = 0; i < ni; i++) {
                S[offset3d(i, j, k, ni, nj)] = D[offset2d(i, k, max_ni_nj)];
            }
        }
    }

/* eta inversions */
    for (i = 1; i < ni - 1; i++) {
        LHSY(A, B, C, current_Q, deta_dx_mat, deta_dy_mat, i);
        if (smoothy(current_Q, deta_dx_mat, deta_dy_mat, ni,nj, A, B, C, i, J_vals_mat,
                    drr, drp, rspec, qv, dd, epsi, Gamma, Mach_inf, delta_t)) {
                    /* returns zero on success */
                    fprintf(stderr, "%s:%d: [Erorr] problem with smoothy in LHSY\n", __FILE__, __LINE__);
                    exit(1);
                }
        for (k = 0; k < 4; k++) {
            for (j = 0; j < nj; j++) {
                D[offset2d(j, k, max_ni_nj)] = S[offset3d(i, j, k, ni, nj)];
            }
        }
        btri4s(A, B, C, D, nj, 1, nj - 2);
        for (k = 0; k < 4; k++) {
            for (j = 0; j < nj; j++) {
                S[offset3d(i, j, k, ni, nj)] = D[offset2d(j, k, max_ni_nj)];
            }
        }
    }

    return calculate_S_norm(S);
}

void matrices_coeffic_and_Jacobian(double *J_vals_mat, double *dxi_dx_mat,
                                   double *dxi_dy_mat, double *deta_dx_mat,
                                   double *deta_dy_mat, double *x_vals_mat,
                                   double *y_vals_mat)
{
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            int index = offset2d(i, j, ni);
            double dx_dxi = first_deriv(x_vals_mat, 'i', i, j);
            double dx_deta = first_deriv(x_vals_mat, 'j', i, j);
            double dy_dxi = first_deriv(y_vals_mat, 'i', i, j);
            double dy_deta = first_deriv(y_vals_mat, 'j', i, j);

            double J = 1.0 / (dx_dxi * dy_deta - dy_dxi * dx_deta);
            J_vals_mat[index] = J;

            dxi_dx_mat[index]  =   J * dy_deta;
            dxi_dy_mat[index]  = - J * dx_deta;
            deta_dx_mat[index] = - J * dy_dxi;
            deta_dy_mat[index] =   J * dx_dxi;
        }
    }
}

void initialize(double *Q, double *J_vals_mat, double *dxi_dx_mat,
                double *dxi_dy_mat, double *deta_dx_mat,
                double *deta_dy_mat, double *x_vals_mat, double *y_vals_mat)
{
    initialize_flow_field(Q);
    matrices_coeffic_and_Jacobian(J_vals_mat, dxi_dx_mat, dxi_dy_mat,
                                  deta_dx_mat, deta_dy_mat, x_vals_mat,
                                  y_vals_mat);
}
