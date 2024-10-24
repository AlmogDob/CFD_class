/* example input file:
t
0.12

i_max
50

j_max
25

i_TEL
11

i_LE
25

i_TEU
39

delta_y
0.02

XSF
1.15

YSF
1.15

x_int
1.008930411365

r
0.01

omega
1

phi
-1

psi
-1
*/
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#define MAXDIR 100
#define MAXWORD 100
#define PI 3.14159265359
#define dprintSTRING(expr) printf(#expr " = %s\n", expr)    /* macro for easy debuging*/
#define dprintINT(expr) printf(#expr " = %d\n", expr)
#define dprintF(expr) printf(#expr " = %g\n", expr)
#define dprintD(expr) printf(#expr " = %g\n", expr)

typedef struct {
    double x;
    double y;
} Vec2;

void read_input(char *dir);
void output_solution(FILE *fp, double *x_vals_mat, double *y_vals_mat);
void mat_output_to_file(FILE *fp, double *data);
int offset2d(int i, int j, int ni);
void initialize(double *x_vals_mat, double *y_vals_mat, double *alpha_vals_mat,
                double *beta_vals_mat, double *gama_vals_mat,
                double *psi_vals_mat,double *phi_vals_mat);
void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat);
double airfoil(double x, char side);
void interpulat_mat(double *x_vals_mat, char diraction);
double first_deriv(double *mat, char diraction, int i, int j);
double second_deriv(double *mat, char diraction, int i, int j);
void alpha_beta_gama(double *alpha_vals_mat, double *beta_vals_mat,
                     double *gama_vals_mat, double *x_vals_mat, double *y_vals_mat);
void psi_phi(double *psi_vals_mat, double *phi_vals_mat,
             double *x_vals_mat, double *y_vals_mat);
void copy_mat(double *dst, double *src);
void copy_row_to_mat(double *dst, double *src, int row_num);
void copy_col_to_mat(double *dst, double *src, int col_num);
double L_x(double *x_vals_mat, double *alpha_vals_mat,
           double *phi_vals_mat, double *beta_vals_mat,
           double *gama_vals_mat, double *psi_vals_mat, int i, int j);
double L_y(double *y_vals_mat, double *alpha_vals_mat,
           double *phi_vals_mat, double *beta_vals_mat,
           double *gama_vals_mat, double *psi_vals_mat, int i, int j);
int sweep1(double *fx_vals_mat, double *fy_vals_mat, double *x_vals_mat,
           double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat,
           double *beta_vals_mat, double *gama_vals_mat,
           double *psi_vals_mat, double *A, double *B,
           double *C,  double *D, double *temp_row);
int sweep2(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat,
           double *fy_vals_mat, double *gama_vals_mat, double *A,
           double *B, double *C, double *D, double *temp_row);
void LHS_sweep1(double *A, double *B, double *C, double *alpha_vals_mat, int j);
void LHS_sweep2(double *A, double *B, double *C, double *alpha_vals_mat, int i);
void RHS_sweep1_x(double *D, double *x_vals_mat, double *alpha_vals_mat,
                double *phi_vals_mat, double *beta_vals_mat,
                double *gama_vals_mat, double *psi_vals_mat, int j);
void RHS_sweep2_x(double *D,double *fx_vals_mat, int j);
void RHS_sweep1_y(double *D, double *y_vals_mat, double *alpha_vals_mat,
                double *phi_vals_mat, double *beta_vals_mat,
                double *gama_vals_mat, double *psi_vals_mat, int j);
void RHS_sweep2_y(double *D,double *fy_vals_mat, int i);
void BC_sweep1(double *A, double *B, double *C, double *D);
void BC_sweep2(double *A, double *B, double *C, double *D);
int tridiag(double *a, double *b, double *c, double *d,
            double *u, int is, int ie);
Vec2 step(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat,
         double *fy_vals_mat, double *x_vals_mat_current,
         double *x_vals_mat_next, double *y_vals_mat_current,
         double *y_vals_mat_next, double *alpha_vals_mat,
         double *phi_vals_mat, double *beta_vals_mat,
         double *gama_vals_mat, double *psi_vals_mat,
         double *A, double *B, double *C, double *D, double *temp_row);
void mat_print(double *data);
double calculate_max_L_x(double *x_vals_mat, double *alpha_vals_mat,
           double *phi_vals_mat, double *beta_vals_mat,
           double *gama_vals_mat, double *psi_vals_mat);
double calculate_max_L_y(double *y_vals_mat, double *alpha_vals_mat,
           double *phi_vals_mat, double *beta_vals_mat,
           double *gama_vals_mat, double *psi_vals_mat);

/* Input variables */
double t, delta_x, delta_y, XSF, YSF, x_int = 1, r, omega, psi_valuse = NAN, phi_valuse = NAN;
int i_max, j_max, i_TEL, i_LE, i_TEU, i_min = 0, j_min = 0;


int main()
{
    /* decleraitons */
    char input_dir[MAXDIR];
    int i_index, j_index;
    double *x_vals_mat_init, *y_vals_mat_init, *x_vals_mat_current,
           *y_vals_mat_current, *x_vals_mat_next, *y_vals_mat_next,
           *alpha_vals_mat, *beta_vals_mat, *gama_vals_mat, *psi_vals_mat,
           *phi_vals_mat, *fx_vals_mat, *fy_vals_mat, *Cx_vals_mat,
           *Cy_vals_mat;
    /* matrix diaganosl for different sweeps */
    double *A, *B, *C, *D, *temp_row;
    Vec2 result, first_result;
    FILE *output_fp = fopen("mesh_output.txt", "wt");
    if (output_fp == NULL) {
        fprintf(stderr, "ERROR: problem opening output file | %s\n", strerror(errno));
        return -1;
    }

    /*------------------------------------------------------------*/

    strcpy(input_dir, "mesh_input.txt");

    /*------------------------------------------------------------*/

    read_input(input_dir);

    /*------------------------------------------------------------*/

    /* Checking that I got the right input */
    dprintD(t);
    dprintINT(i_max);
    dprintINT(j_max);
    dprintINT(i_TEL);
    dprintINT(i_LE);
    dprintINT(i_TEU);
    dprintD(delta_x);
    dprintD(delta_y);
    dprintD(XSF);
    dprintD(YSF);
    dprintD(x_int);
    dprintD(r);
    dprintD(omega);
    dprintD(phi_valuse);
    dprintD(psi_valuse);
    printf("--------------------\n");

    /*------------------------------------------------------------*/
    
    /* Memory allocation */
    x_vals_mat_init = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            x_vals_mat_init[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    y_vals_mat_init = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            y_vals_mat_init[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    x_vals_mat_current = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            x_vals_mat_current[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    y_vals_mat_current = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            y_vals_mat_current[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    x_vals_mat_next = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            x_vals_mat_next[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    y_vals_mat_next = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            y_vals_mat_next[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    alpha_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            alpha_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    beta_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            beta_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    gama_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            gama_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    psi_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            psi_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    phi_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            phi_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    fx_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            fx_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    fy_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            fy_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    Cx_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            Cx_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    Cy_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            Cy_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }

    A = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        A[i_index] = 0;
    }
    B = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        B[i_index] = 0;
    }
    C = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        C[i_index] = 0;
    }
    D = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        D[i_index] = 0;
    }
    temp_row = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        temp_row[i_index] = 0;
    }
    
    /*------------------------------------------------------------*/
    
    initialize(x_vals_mat_init, y_vals_mat_init, alpha_vals_mat, beta_vals_mat, gama_vals_mat,
               psi_vals_mat, phi_vals_mat);
    
    copy_mat(x_vals_mat_current, x_vals_mat_init);
    copy_mat(x_vals_mat_next, x_vals_mat_init);
    copy_mat(y_vals_mat_current, y_vals_mat_init);
    copy_mat(y_vals_mat_next, y_vals_mat_init);

    for (i_index = 0; i_index < 5e6; i_index++) {
        result = step(Cx_vals_mat, Cy_vals_mat, fx_vals_mat, fy_vals_mat,
                       x_vals_mat_current, x_vals_mat_next, y_vals_mat_current,
                       y_vals_mat_next, alpha_vals_mat, phi_vals_mat,
                       beta_vals_mat, gama_vals_mat, psi_vals_mat,
                       A, B, C, D, temp_row);
        if (i_index == 0) {
            first_result = result;
        }
        if (result.x == 1 && result.y == 0) {
            fprintf(stderr, "ERROR: Step - sweep 1\n");
            exit(1);
        }
        if (result.x == 2 && result.y == 0) {
            fprintf(stderr, "ERROR: Step - sweep 2\n");
            exit(2);
        }

        copy_mat(x_vals_mat_current, x_vals_mat_next);
        copy_mat(y_vals_mat_current, y_vals_mat_next);

        printf("%5d. Lx_max: %0.10f, Ly_max: %0.10f\n",i_index+1, result.x, result.y);
        /* checking convergenc */
        if (log10(fabs(first_result.x/result.x)) > 6 && log10(fabs(first_result.y/result.y)) > 6) {
            break;
        }
    }

    /* printing Lx and Ly */
    // printf("%4d. Lx_max: %0.10f, Ly_max: %0.10f\n",i_index+1, result.x, result.y);
    
    output_solution(output_fp, x_vals_mat_next, y_vals_mat_next);

    /*------------------------------------------------------------*/

    free(x_vals_mat_init);
    free(y_vals_mat_init);
    free(x_vals_mat_current);
    free(y_vals_mat_current);
    free(x_vals_mat_next);
    free(y_vals_mat_next);
    free(alpha_vals_mat);
    free(beta_vals_mat);
    free(gama_vals_mat);
    free(psi_vals_mat);
    free(phi_vals_mat);
    free(fx_vals_mat);
    free(fy_vals_mat);
    free(Cx_vals_mat);
    free(Cy_vals_mat);
    free(A);
    free(B);
    free(C);
    free(D);
    free(temp_row);

    fclose(output_fp);

    return 0;
}

/* sets 'flags' and variables according to the input file
argument list:
dir - the directory of the input file */
void read_input(char *dir)
{
    FILE *fp = fopen(dir, "rt");
    char current_word[MAXWORD];
    float temp;

    if (!fp) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
        exit(1);
    }

    /* Seting the input varibles according to the input file */
    while(fscanf(fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "t")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            t = (double)temp;
        } else if (!strcmp(current_word, "i_max")) {
            fscanf(fp, "%d", &i_max);
        } else if (!strcmp(current_word, "j_max")) {
            fscanf(fp, "%d ", &j_max);
        } else if (!strcmp(current_word, "i_TEL")) {
            fscanf(fp, "%d", &i_TEL);
        } else if (!strcmp(current_word, "i_LE")) {
            fscanf(fp, "%d", &i_LE);
        } else if (!strcmp(current_word, "i_TEU")) {
            fscanf(fp, "%d", &i_TEU);
        } else if (!strcmp(current_word, "delta_y")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            delta_y = (double)temp;
        } else if (!strcmp(current_word, "XSF")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            XSF = (double)temp;
        } else if (!strcmp(current_word, "YSF")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            YSF = (double)temp;
        } else if (!strcmp(current_word, "x_int")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            x_int = (double)temp;
        } else if (!strcmp(current_word, "r")) {
            fscanf(fp, "%g", &temp);
            r = (double)temp;
        } else if (!strcmp(current_word, "omega")) {
            fscanf(fp, "%g", &temp);
            omega = (double)temp;
        } else if (!strcmp(current_word, "phi")) {
            fscanf(fp, "%g", &temp);
            phi_valuse = (double)temp;
        } else if (!strcmp(current_word, "psi")) {
            fscanf(fp, "%g", &temp);
            psi_valuse = (double)temp;
        }
    }

    delta_x = 1.0/(i_LE - i_TEL);

    fclose(fp);
}

/* output solution;
argument list:
fp - file pointer to output file
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void output_solution(FILE *fp, double *x_vals_mat, double *y_vals_mat)
{
    fprintf(fp, "ni\n%d\n\nnj\n%d\n\n", i_max+1, j_max+1);
    fprintf(fp, "x_vals\n");
    mat_output_to_file(fp, x_vals_mat);
    fprintf(fp, "\n");
    fprintf(fp, "y_vals\n");
    mat_output_to_file(fp, y_vals_mat);
}

/* output data;
argument list:
fp - file pointer to output file
data - 1D array of the y valuse */
void mat_output_to_file(FILE *fp, double *data)
{
    int i, j;
    
    for (j = 0; j < j_max+1; j++) {
        for (i = 0; i < i_max+1; i++) {
            fprintf(fp, "%g ", data[offset2d(i, j, i_max+1)]);
        }
        fprintf(fp, "\n");
    }
}

/* converts a 2D index into 1D index
argument list:
i - x position
j - y position
ni - stride */
int offset2d(int i, int j, int ni)
{
    return j * ni + i;
}

/* set inital valsuse of the mash points 
argument list:
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse
alpha_vals_mat - 1D array for the alpha valus
beta_vals_mat - 1D array for the beta valus
gama_vals_mat - 1D array for the gama valus */
void initialize(double *x_vals_mat, double *y_vals_mat, double *alpha_vals_mat,
                double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat,
                double *phi_vals_mat)
{
    set_grid_boundaries(x_vals_mat, y_vals_mat);
    interpulat_mat(x_vals_mat, 'j');
    interpulat_mat(y_vals_mat, 'j');
    alpha_beta_gama(alpha_vals_mat, beta_vals_mat, gama_vals_mat, x_vals_mat, y_vals_mat);
    psi_phi(psi_vals_mat, phi_vals_mat, x_vals_mat, y_vals_mat);
}

/* set the mash boundaries coorditates 
argument list: x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat)
{
    int i_index, i_min = 0, j_index, j_min = 0, index = 0, num_points_befor_circle,
    num_of_outer_segments, num_of_top_outer_segments;

    double x, x_i_minos_2, x_i_minos_1, y_j_minos_2, y_j_minos_1, y_imax_jmax, x_imax_jmax,
    delta_theta, R, theta = 0, length_top_outer, segment_length, current_x_val = 0;

    /* setting the boundary according to the exercie */
    /* Eq 6 */
    for (i_index = i_TEL, j_index = j_min; i_index < i_TEU+1; i_index++) { 
        x = 1 - cos(0.5*PI*(i_LE-i_index)*delta_x);
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x;
        if (i_index < i_LE) {
            y_vals_mat[offset2d(i_index, j_index, i_max+1)] = airfoil(x, 'l');
        } else if (i_index >= i_LE) {
            y_vals_mat[offset2d(i_index, j_index, i_max+1)] = airfoil(x, 'u');
        }
    }
    /* Eq 7 */
    for (i_index = i_TEU + 1, j_index = j_min; i_index < i_max+1; i_index++) {
        x_i_minos_1 = x_vals_mat[offset2d(i_index-1, j_index, i_max+1)];  
        x_i_minos_2 = x_vals_mat[offset2d(i_index-2, j_index, i_max+1)];  
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x_i_minos_1 + (x_i_minos_1 - x_i_minos_2) * XSF;
    }
    for (i_index = i_min, j_index = j_min; i_index < i_TEL; i_index++) {
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x_vals_mat[offset2d(i_max-i_index, j_index, i_max+1)];
    }
    /* Eq 8 */
    y_vals_mat[offset2d(i_max, j_min+1, i_max+1)] = delta_y;
    for (i_index = i_max, j_index = j_min+2; j_index < j_max+1; j_index++) {
        y_j_minos_1 = y_vals_mat[offset2d(i_index, j_index-1, i_max+1)];  
        y_j_minos_2 = y_vals_mat[offset2d(i_index, j_index-2, i_max+1)];  
        y_vals_mat[offset2d(i_max, j_index, i_max+1)] = y_j_minos_1 + (y_j_minos_1 - y_j_minos_2) * YSF;
    }
    for (i_index = i_max, j_index = j_min+1; j_index < j_max+1; j_index++) {
        x_vals_mat[offset2d(i_max, j_index, i_max+1)] = x_vals_mat[offset2d(i_max, j_min, i_max+1)];
        y_vals_mat[offset2d(i_min, j_index, i_max+1)] = -y_vals_mat[offset2d(i_max, j_index, i_max+1)];
        x_vals_mat[offset2d(i_min, j_index, i_max+1)] = x_vals_mat[offset2d(i_min, j_min, i_max+1)];
    }
    /* Outer boundary */
    y_imax_jmax = y_vals_mat[offset2d(i_max, j_max, i_max+1)];
    x_imax_jmax = x_vals_mat[offset2d(i_max, j_max, i_max+1)];
    R = y_imax_jmax;
    /*test*/
    // dprintD(x_imax_jmax);
    /*test*/

    num_of_outer_segments = i_max;
    num_of_top_outer_segments = num_of_outer_segments/2;
    length_top_outer = x_imax_jmax + 0.5*PI*R; /* length of stright part and quarter of the circle */
    segment_length = length_top_outer/num_of_top_outer_segments;

    /* the stright line part */
    for (num_points_befor_circle = 0;
         num_points_befor_circle < num_of_top_outer_segments + 1;
         num_points_befor_circle++) {
            current_x_val = x_imax_jmax - num_points_befor_circle*segment_length; 
            if (current_x_val < 0) {
                break;
            }
            x_vals_mat[offset2d(i_max-num_points_befor_circle, j_max, i_max+1)] = current_x_val;
            y_vals_mat[offset2d(i_max-num_points_befor_circle, j_max, i_max+1)] = y_imax_jmax;
    }

    theta = PI/2 + atan(x_vals_mat[offset2d(i_max-num_points_befor_circle+1, j_max, i_max+1)] / R);
    delta_theta = theta / (num_of_top_outer_segments - num_points_befor_circle + 1);

    /* the quarter circle part */
    for (index = 0; index < num_of_top_outer_segments - num_points_befor_circle + 1; index++) {
        x_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)] = -R*cos(delta_theta*index);
        y_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)] = R*sin(delta_theta*index);
    }

    /* coping to the bottom side */
    for (index = 1; index < i_max/2 + 1; index++) {
        x_vals_mat[offset2d(i_max/2 - index, j_max, i_max+1)] = x_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)];
        y_vals_mat[offset2d(i_max/2 - index, j_max, i_max+1)] = -y_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)];
    }
}

/* returns the shape of the airfoil as a function of x
argument list: x - x positon 
side - the uppers side are the low side of the airfoil */
double airfoil(double x, char side)
{
    if (side == 'u') {
        return 5 * t * (0.2969*sqrt(x_int*x) - 0.1260*(x_int*x)
               - 0.3516*pow(x_int*x,2) + 0.2843*pow(x_int*x,3)
               - 0.1015*pow(x_int*x,4));
    } else if (side == 'l') {
        return - 5 * t * (0.2969*sqrt(x_int*x) - 0.1260*(x_int*x)
               - 0.3516*pow(x_int*x,2) + 0.2843*pow(x_int*x,3)
               - 0.1015*pow(x_int*x,4));
    }
    fprintf(stderr, "Error with the airfoil usage\n");
    exit(1);
}

/* fill the matrix in a way of interpulation the boundarys
argument list:
mat - 1D array of valsuse
diraction - i direction or j direction */
void interpulat_mat(double *mat, char diraction)
{
    int i, j; 
    double max, min;

    assert(diraction == 'j' || diraction == 'i'); /* different directions are not implemented */

    if (diraction == 'j') {
        for (i = 1; i < i_max; i++) {
            for (j = 1; j < j_max; j++) {
                max = mat[offset2d(i, j_max, i_max+1)];
                min = mat[offset2d(i, j_min, i_max+1)];
                mat[offset2d(i, j, i_max+1)] = (max - min)/(j_max) * j + min;   /* liniar interpulation */
            }
        }
    }
    if (diraction == 'i') {
        for (j = 1; j < j_max; j++) {
            for (i = 1; i < i_max; i++) {
                max = mat[offset2d(i_max, j, i_max+1)];
                min = mat[offset2d(i_min, j, i_max+1)];
                mat[offset2d(i, j, i_max+1)] = (max - min)/(i_max) * i + min;   /* liniar interpulation */
            }
        }
    }
}

/* return the second order first derivitive from the valuse in the mat matrix
argument list:
mat - 1D array of valuse
diraction - i or j
i, j - the points coordinates */
double first_deriv(double *mat, char diraction, int i, int j)
{
    if (diraction == 'j') {
        if (j == j_min || j == j_max) {
            return 0;
        }
        return (mat[offset2d(i, j+1, i_max+1)] - mat[offset2d(i, j-1, i_max+1)]) / (2); /* second order first derivitive */
    }
    if (diraction == 'i') {
        if (i == i_min || i == i_max) {
            return 0;
        }
        return (mat[offset2d(i+1, j, i_max+1)] - mat[offset2d(i-1, j, i_max+1)]) / (2); /* second order first derivitive */
    }
    return NAN;
}

/* return the second order second derivitive from the valuse in the mat matrix
argument list:
mat - 1D array of valuse
diraction - i or j
i, j - the points coordinates */
double second_deriv(double *mat, char diraction, int i, int j)
{
    if (diraction == 'j') {
        if (j == j_min || j == j_max) {
            return 0;
        }
        return (mat[offset2d(i, j+1, i_max+1)] -2*mat[offset2d(i, j, i_max+1)] + mat[offset2d(i, j-1, i_max+1)]) / (1); /* second order second derivitive */
    }
    if (diraction == 'i') {
        if (i == i_min || i == i_max) {
            return 0;
        }
        return (mat[offset2d(i+1, j, i_max+1)] -2*mat[offset2d(i, j, i_max+1)] + mat[offset2d(i-1, j, i_max+1)]) / (1); /* second order second derivitive */
    }
    return NAN;
}

/* fills the alpha and beta and gama matrices
argument list:
alpha_vals_mat - 1D array for the alpha valus
beta_vals_mat - 1D array for the beta valus
gama_vals_mat - 1D array for the gama valus
x_vals_mat - 1D array of the x valus
y_vals_mat - 1D array of the y valus */
void alpha_beta_gama(double *alpha_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *x_vals_mat, double *y_vals_mat)
{
    /* according to equation 3 */
    int i, j;
    double Dx_Deta, Dy_Deta, Dx_Dxai, Dy_Dxai;

    for (i = 0; i < i_max + 1; i++) {
        for (j = 0; j < j_max + 1; j++) {
            Dx_Deta = first_deriv(x_vals_mat, 'j', i, j);
            Dy_Deta = first_deriv(y_vals_mat, 'j', i, j);
            Dx_Dxai = first_deriv(x_vals_mat, 'i', i, j);
            Dy_Dxai = first_deriv(y_vals_mat, 'i', i, j);
            alpha_vals_mat[offset2d(i, j, i_max+1)] = Dx_Deta*Dx_Deta + Dy_Deta*Dy_Deta;
            beta_vals_mat[offset2d(i, j, i_max+1)] = Dx_Dxai*Dx_Deta + Dy_Dxai*Dy_Deta;
            gama_vals_mat[offset2d(i, j, i_max+1)] = Dx_Dxai*Dx_Dxai + Dy_Dxai*Dy_Dxai;
        }
    }
}

/* creat and fill the psi and phi matrices
argument list:
psi_vals_mat - 1D array for the psi valus
phi_vals_mat - 1D array for the phi valus
x_vals_mat - 1D array of the x valus
y_vals_mat - 1D array of the y valus */
void psi_phi(double *psi_vals_mat, double *phi_vals_mat, double *x_vals_mat, double *y_vals_mat)
{
    /* according to equation 4 and 5 */
    int i, j;
    double Dx_Deta_min, Dy_Deta_min, Dx_Dxai_min, Dy_Dxai_min, Dx_Deta_max, Dy_Deta_max, Dx_Dxai_max, Dy_Dxai_max,
    Dx_Deta_Deta_min, Dx_Deta_Deta_max, Dy_Deta_Deta_min, Dy_Deta_Deta_max, Dx_Dxai_Dxai_min, Dx_Dxai_Dxai_max,
    Dy_Dxai_Dxai_min, Dy_Dxai_Dxai_max;

    /*test*/
    // dprintD(second_deriv(x_vals_mat, 'j', i_min, j = 1));
    /*test*/

    /* eq 4 */
    for (j = 0; j < j_max+1; j++) {
        if (psi_valuse == -1) {
            Dx_Deta_min = first_deriv(x_vals_mat, 'j', i_min, j);
            Dy_Deta_min = first_deriv(y_vals_mat, 'j', i_min, j);
            Dx_Deta_max = first_deriv(x_vals_mat, 'j', i_max, j);
            Dy_Deta_max = first_deriv(y_vals_mat, 'j', i_max, j);
            Dx_Deta_Deta_min = second_deriv(x_vals_mat, 'j', i_min, j);
            Dy_Deta_Deta_min = second_deriv(y_vals_mat, 'j', i_min, j);
            Dx_Deta_Deta_max = second_deriv(x_vals_mat, 'j', i_max, j);
            Dy_Deta_Deta_max = second_deriv(y_vals_mat, 'j', i_max, j);

            if (fabs(Dy_Deta_min) > fabs(Dx_Deta_min)) {
                psi_vals_mat[offset2d(i_min, j, i_max+1)] = - Dy_Deta_Deta_min / Dy_Deta_min;
            }
            if (fabs(Dy_Deta_min) < fabs(Dx_Deta_min)) {
                psi_vals_mat[offset2d(i_min, j, i_max+1)] = - Dx_Deta_Deta_min / Dx_Deta_min;
            }
            if (fabs(Dy_Deta_max) > fabs(Dx_Deta_max)) {
                psi_vals_mat[offset2d(i_max, j, i_max+1)] = - Dy_Deta_Deta_max / Dy_Deta_max;
            }
            if (fabs(Dy_Deta_max) < fabs(Dx_Deta_max)) {
                psi_vals_mat[offset2d(i_max, j, i_max+1)] = - Dx_Deta_Deta_max / Dx_Deta_max;
            }
        } else {
            psi_vals_mat[offset2d(i_min, j, i_max+1)] = psi_valuse;
            psi_vals_mat[offset2d(i_max, j, i_max+1)] = psi_valuse;
        }
    }

    /* eq 5 */
    for (i = 0; i < i_max+1; i++) {
        if (phi_valuse == -1) {
            Dx_Dxai_min = first_deriv(x_vals_mat, 'i', i, j_min);
            Dy_Dxai_min = first_deriv(y_vals_mat, 'i', i, j_min);
            Dx_Dxai_max = first_deriv(x_vals_mat, 'i', i, j_max);
            Dy_Dxai_max = first_deriv(y_vals_mat, 'i', i, j_max);
            Dx_Dxai_Dxai_min = second_deriv(x_vals_mat, 'i', i, j_min);
            Dy_Dxai_Dxai_min = second_deriv(y_vals_mat, 'i', i, j_min);
            Dx_Dxai_Dxai_max = second_deriv(x_vals_mat, 'i', i, j_max);
            Dy_Dxai_Dxai_max = second_deriv(y_vals_mat, 'i', i, j_max);

            if (fabs(Dx_Dxai_min) > fabs(Dy_Dxai_min)) {
                phi_vals_mat[offset2d(i, j_min, i_max+1)] = - Dx_Dxai_Dxai_min / Dx_Dxai_min;
            }
            if (fabs(Dx_Dxai_min) < fabs(Dy_Dxai_min)) {
                phi_vals_mat[offset2d(i, j_min, i_max+1)] = - Dy_Dxai_Dxai_min / Dy_Dxai_min;
            }
            if (fabs(Dx_Dxai_max) > fabs(Dy_Dxai_max)) {
                phi_vals_mat[offset2d(i, j_max, i_max+1)] = - Dx_Dxai_Dxai_max / Dx_Dxai_max;
            }
            if (fabs(Dx_Dxai_max) < fabs(Dy_Dxai_max)) {
                phi_vals_mat[offset2d(i, j_max, i_max+1)] = - Dy_Dxai_Dxai_max / Dy_Dxai_max;
            }
        } else {
            phi_vals_mat[offset2d(i_min, j, i_max+1)] = phi_valuse;
            phi_vals_mat[offset2d(i_max, j, i_max+1)] = phi_valuse;
        }
    }

    interpulat_mat(psi_vals_mat, 'i');
    interpulat_mat(phi_vals_mat, 'j');
}

/* copys matrix src to matrix src 
argument list:
dst - 1D array for the destination valus 
src - 1D array of the source valus */
void copy_mat(double *dst, double *src)
{
    int i, j;
    
    for (i = 0; i < i_max+1; i++) {
        for (j = 0; j < j_max+1; j++) {
            dst[offset2d(i, j, i_max+1)] = src[offset2d(i, j, i_max+1)];
        }
    }
}

/* copys the row 'src' into the j row in the destination matrix */
void copy_row_to_mat(double *dst, double *src, int row_num)
{
    int i;
    
    for (i = 0; i < i_max+1; i++) {
        dst[offset2d(i, row_num, i_max+1)] = src[i];
    }
}

void copy_col_to_mat(double *dst, double *src, int col_num)
{
    int j;
    
    for (j = 0; j < j_max+1; j++) {
        dst[offset2d(col_num, j, i_max+1)] = src[j];
    }
}

/* returns the valus of L_x according to equation 10
argumetn list:
x_vals_mat - 1D array of the x valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus 
i, j = the point coordinate */
double L_x(double *x_vals_mat, double *alpha_vals_mat,
           double *phi_vals_mat, double *beta_vals_mat,
           double *gama_vals_mat, double *psi_vals_mat, int i, int j)
{
    double first_element, second_element, third_element, x_i_plus1_j, x_i_j,
    x_i_minus1_j, x_i_plus1_j_plus1, x_i_plus1_j_minus1, x_i_minus1_j_plus1,
    x_i_minus1_j_minus1, x_i_j_plus1, x_i_j_minus1;

    if (i == 0 || i == i_max || j == 0 || j == j_max) {
        return 0;
    }

    x_i_plus1_j = x_vals_mat[offset2d(i+1, j, i_max+1)];
    x_i_j = x_vals_mat[offset2d(i, j, i_max+1)];
    x_i_minus1_j = x_vals_mat[offset2d(i-1, j, i_max+1)];
    x_i_plus1_j_plus1 = x_vals_mat[offset2d(i+1,j+1, i_max+1)];
    x_i_plus1_j_minus1 = x_vals_mat[offset2d(i+1, j-1, i_max+1)];
    x_i_minus1_j_plus1 = x_vals_mat[offset2d(i-1, j+1, i_max+1)];
    x_i_minus1_j_minus1 = x_vals_mat[offset2d(i-1, j-1, i_max+1)];
    x_i_j_plus1 = x_vals_mat[offset2d(i, j+1, i_max+1)];
    x_i_j_minus1 = x_vals_mat[offset2d(i, j-1, i_max+1)];

    first_element = alpha_vals_mat[offset2d(i, j, i_max+1)] *
                    ((x_i_plus1_j - 2 * x_i_j + x_i_minus1_j) +
                    0.5 * phi_vals_mat[offset2d(i, j, i_max+1)] *
                    (x_i_plus1_j - x_i_minus1_j));
    second_element = 0.5 * beta_vals_mat[offset2d(i, j, i_max+1)] *
                     (x_i_plus1_j_plus1 - x_i_plus1_j_minus1 -
                     x_i_minus1_j_plus1 + x_i_minus1_j_minus1);
    third_element = gama_vals_mat[offset2d(i, j, i_max+1)] *
                    ((x_i_j_plus1 - 2 * x_i_j + x_i_j_minus1) +
                    0.5 * psi_vals_mat[offset2d(i, j, i_max+1)] *
                    (x_i_j_plus1 - x_i_j_minus1));

    return first_element - second_element + third_element;
}

/* returns the valus of L_y according to equation 11
argumetn list:
y_vals_mat - 1D array of the y valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
phi_vals_mat - 1D array of the phi valus
psi_vals_mat - 1D array of the psi valus 
i, j = the point coordinate */
double L_y(double *y_vals_mat, double *alpha_vals_mat,
           double *phi_vals_mat, double *beta_vals_mat,
           double *gama_vals_mat, double *psi_vals_mat, int i, int j)
{
    double first_element, second_element, third_element, y_i_plus1_j, y_i_j,
    y_i_minus1_j, y_i_plus1_j_plus1, y_i_plus1_j_minus1, y_i_minus1_j_plus1,
    y_i_minus1_j_minus1, y_i_j_plus1, y_i_j_minus1;

    if (i == 0 || i == i_max || j == 0 || j == j_max) {
        return 0;
    }

    y_i_plus1_j = y_vals_mat[offset2d(i+1, j, i_max+1)];
    y_i_j = y_vals_mat[offset2d(i, j, i_max+1)];
    y_i_minus1_j = y_vals_mat[offset2d(i-1, j, i_max+1)];
    y_i_plus1_j_plus1 = y_vals_mat[offset2d(i+1,j+1, i_max+1)];
    y_i_plus1_j_minus1 = y_vals_mat[offset2d(i+1, j-1, i_max+1)];
    y_i_minus1_j_plus1 = y_vals_mat[offset2d(i-1, j+1, i_max+1)];
    y_i_minus1_j_minus1 = y_vals_mat[offset2d(i-1, j-1, i_max+1)];
    y_i_j_plus1 = y_vals_mat[offset2d(i, j+1, i_max+1)];
    y_i_j_minus1 = y_vals_mat[offset2d(i, j-1, i_max+1)];

    first_element = alpha_vals_mat[offset2d(i, j, i_max+1)] *
                    ((y_i_plus1_j - 2 * y_i_j + y_i_minus1_j) +
                    0.5 * phi_vals_mat[offset2d(i, j, i_max+1)] *
                    (y_i_plus1_j - y_i_minus1_j));
    second_element = 0.5 * beta_vals_mat[offset2d(i, j, i_max+1)] *
                     (y_i_plus1_j_plus1 - y_i_plus1_j_minus1 -
                     y_i_minus1_j_plus1 + y_i_minus1_j_minus1);
    third_element = gama_vals_mat[offset2d(i, j, i_max+1)] *
                    ((y_i_j_plus1 - 2 * y_i_j + y_i_j_minus1) +
                    0.5 * psi_vals_mat[offset2d(i, j, i_max+1)] *
                    (y_i_j_plus1 - y_i_j_minus1));

    return first_element - second_element + third_element;
}

/* doing the first sweep according to equation 12;
returns 0 on success
argument list:
fx_vals_mat - 1D array for the fx valus 
fy_vals_mat - 1D array for the fy valus
x_vals_mat - 1D array of the x valus
y_vals_mat - 1D array of the y valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus */
int sweep1(double *fx_vals_mat, double *fy_vals_mat, double *x_vals_mat,
           double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat,
           double *beta_vals_mat, double *gama_vals_mat,
           double *psi_vals_mat, double *A, double *B,
           double *C,  double *D, double *temp_row)
{
    int j_index, success = 0;

    /* solving for each j */
    for (j_index = 0; j_index < j_max+1; j_index++) {
        LHS_sweep1(A, B, C, alpha_vals_mat, j_index);
        RHS_sweep1_x(D, x_vals_mat, alpha_vals_mat, phi_vals_mat,
                     beta_vals_mat, gama_vals_mat, psi_vals_mat, j_index);
        BC_sweep1(A, B, C, D);
        success = tridiag(A, B, C, D, temp_row, 0, i_max);
        if (success == 1) {
            break;
        }
        copy_row_to_mat(fx_vals_mat, temp_row, j_index);

        LHS_sweep1(A, B, C, alpha_vals_mat, j_index);
        RHS_sweep1_y(D, y_vals_mat, alpha_vals_mat, phi_vals_mat,
                     beta_vals_mat, gama_vals_mat, psi_vals_mat, j_index);
        BC_sweep1(A, B, C, D);
        success = tridiag(A, B, C, D, temp_row, 0, i_max);
        if (success == 1) {
            break;
        }
        copy_row_to_mat(fy_vals_mat, temp_row, j_index);
    }
    if (success == 1) {
        return 1;
    }

    return 0;
}

/* doing the first sweep according to equation 13;
returns 0 on success
argument list:
Cx_vals_mat - 1D array for the Cx valus 
Cy_vals_mat - 1D array for the Cy valus
fx_vals_mat - 1D array of the fx valus
fy_vals_mat - 1D array of the fy valus
gama_vals_mat - 1D array of the gama valus */
int sweep2(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat,
           double *fy_vals_mat, double *gama_vals_mat, double *A,
           double *B, double *C, double *D, double *temp_row)
{
    int i_index, success = 0;

    /* solving for each i */
    for (i_index = 0; i_index < i_max+1; i_index++) {
        LHS_sweep2(A, B, C, gama_vals_mat, i_index);
        RHS_sweep2_x(D, fx_vals_mat, i_index);
        BC_sweep2(A, B, C, D);
        /*test*/
        // for (int index = 0; index < j_max+1; index++) {
        //     printf("%g\n", B[index]);
        // }
        /*test*/
        success = tridiag(A, B, C, D, temp_row, 0, j_max);
        if (success == 1) {
            printf("1\n");
            break;
        }
        copy_col_to_mat(Cx_vals_mat, temp_row, i_index);

        LHS_sweep2(A, B, C, gama_vals_mat, i_index);
        RHS_sweep2_y(D, fy_vals_mat, i_index);
        BC_sweep2(A, B, C, D);
        /*test*/
        // for (int index = 0; index < j_max+1; index++) {
        //     printf("%g\n", Bx2[index]);
        // }
        /*test*/
        success = tridiag(A, B, C, D, temp_row, 0, j_max);
        if (success == 1) {
            printf("2\n");
            break;
        }
    copy_col_to_mat(Cy_vals_mat, temp_row, i_index);
    }
    if (success == 1) {
        return 1;
    }

    return 0;
}

/* populates the A, B, C vectors accoding to eq 15 
argument list:
A - 1D array for the A valus 
B - 1D array for the B valus 
C - 1D array for the C valus
alpha_vals_mat - 1D array of the alpha valus
j - the row number */
void LHS_sweep1(double *A, double *B, double *C, double *alpha_vals_mat, int j)
{
    int i;

    for (i = 0; i < i_max+1; i++) {
        A[i] = -alpha_vals_mat[offset2d(i, j, i_max+1)];
        B[i] = r + 2 * alpha_vals_mat[offset2d(i, j, i_max+1)];
        C[i] = -alpha_vals_mat[offset2d(i, j, i_max+1)];
    }
}

/* populates the A, B, C vectors accoding to eq 15 
argument list:
A - 1D array for the A valus 
B - 1D array for the B valus 
C - 1D array for the C valus
alpha_vals_mat - 1D array of the alpha valus
j - the row number */
void LHS_sweep2(double *A, double *B, double *C, double *gama_vals_mat, int i)
{
    int j;

    for (j = 0; j < j_max+1; j++) {
        A[j] = -gama_vals_mat[offset2d(i, j, i_max+1)];
        B[j] = r + 2 * gama_vals_mat[offset2d(i, j, i_max+1)];
        C[j] = -gama_vals_mat[offset2d(i, j, i_max+1)];
    }
}

/* populates the D vectors accoding to eq 15 
argument list:
D - 1D array for the D valus 
x_vals_mat - 1D array of the x valus
alpha_vals_mat - 1D array of the alpha valus
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus 
j - the row number */
void RHS_sweep1_x(double *D, double *x_vals_mat, double *alpha_vals_mat,
                double *phi_vals_mat, double *beta_vals_mat,
                double *gama_vals_mat, double *psi_vals_mat, int j)
{
    int i;
    
    for (i = 0; i < i_max+1; i++) {
        D[i] = r * omega * L_x(x_vals_mat, alpha_vals_mat, phi_vals_mat,
                               beta_vals_mat, gama_vals_mat,
                               psi_vals_mat, i, j);
    }
}

/* populates the D vectors accoding to eq 15 
argument list:
D - 1D array for the D valus 
fx_vals_mat - 1D array of the fx valus
i - the row number */
void RHS_sweep2_x(double *D,double *fx_vals_mat, int i)
{
    int j;
    
    for (j = 0; j < j_max+1; j++) {
        D[j] = fx_vals_mat[offset2d(i, j, i_max+1)];
    }
}

/* populates the D vectors accoding to eq 15 
argument list:
D - 1D array for the D valus 
y_vals_mat - 1D array of the y valus
alpha_vals_mat - 1D array of the alpha valus
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus 
j - the row number */
void RHS_sweep1_y(double *D, double *y_vals_mat, double *alpha_vals_mat,
                double *phi_vals_mat, double *beta_vals_mat,
                double *gama_vals_mat, double *psi_vals_mat, int j)
{
    int i;
    
    for (i = 0; i < i_max+1; i++) {
        D[i] = r * omega * L_y(y_vals_mat, alpha_vals_mat, phi_vals_mat,
                               beta_vals_mat, gama_vals_mat,
                               psi_vals_mat, i, j);
    }
}

/* populates the D vectors accoding to eq 15 
argument list:
D - 1D array for the D valus 
fy_vals_mat - 1D array of the fy valus
i - the row number */
void RHS_sweep2_y(double *D,double *fy_vals_mat, int i)
{
    int j;
    
    for (j = 0; j < j_max+1; j++) {
        D[j] = fy_vals_mat[offset2d(i, j, i_max+1)];
    }
}

/* applaying the boundary conditions on the A, B, C and D vectors according to eq 16
argument list:
A, B, C are the tridaig diaganols
D is the RHS vector */
void BC_sweep1(double *A, double *B, double *C, double *D)
{
    A[0] = 0;
    B[0] = 1;
    C[0] = 0;
    D[0] = 0;
    A[i_max] = 0;
    B[i_max] = 1;
    C[i_max] = 0;
    D[i_max] = 0;
}

/* applaying the boundary conditions on the A, B, C and D vectors according to eq 16
argument list:
A, B, C are the tridaig diaganols
D is the RHS vector */
void BC_sweep2(double *A, double *B, double *C, double *D)
{
    A[0] = 0;
    B[0] = 1;
    C[0] = 0;
    D[0] = 0;
    A[j_max] = 0;
    B[j_max] = 1;
    C[j_max] = 0;
    D[j_max] = 0;
}

/* a, b, c, are the vectors of the diagonal and the two
off-diagonals. The vector d is the RHS vector, the vector
u is the solution vector, "is" is the starting point, and
ie is the last point */
int tridiag(double *a, double *b, double *c, double *d, double *u, int is, int ie)
{
    int i;
    double beta;

    for (i = is + 1; i <= ie; i++)
    {
        if (b[i - 1] == 0.)
            return (1);
        beta = a[i] / b[i - 1];
        b[i] = b[i] - c[i - 1] * beta;
        d[i] = d[i] - d[i - 1] * beta;
    }

    u[ie] = d[ie] / b[ie];
    for (i = ie - 1; i >= is; i--)
    {
        u[i] = (d[i] - c[i] * u[i + 1]) / b[i];
    }
    return (0);
} 

/* doing the first sweep according to eq 9-17;
returns Vec2 with Lx_max, Ly_max (the residual)
argument list:
Cx_vals_mat - 1D array of the Cx valus 
Cy_vals_mat - 1D array of the Cy valus
fx_vals_mat - 1D array of the fx valus 
fy_vals_mat - 1D array of the fy valus
x_vals_mat_current - 1D array of the current x valus
x_vals_mat_next - 1D array of the next x valus
y_vals_mat_current - 1D array of the current y valus
y_vals_mat_next - 1D array of the next y valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus
A-D are temperary vectors for inverting the tri-diag matrices */
Vec2 step(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat,
         double *fy_vals_mat, double *x_vals_mat_current,
         double *x_vals_mat_next, double *y_vals_mat_current,
         double *y_vals_mat_next, double *alpha_vals_mat,
         double *phi_vals_mat, double *beta_vals_mat,
         double *gama_vals_mat, double *psi_vals_mat,
         double *A, double *B, double *C, double *D, double *temp_row)
{
    int success, i, j, index; 
    Vec2 ans = {.x = 0, .y = 0};

    alpha_beta_gama(alpha_vals_mat, beta_vals_mat, gama_vals_mat, x_vals_mat_current, y_vals_mat_current);

    success = sweep1(fx_vals_mat, fy_vals_mat, x_vals_mat_current,
           y_vals_mat_current, alpha_vals_mat, phi_vals_mat,
           beta_vals_mat, gama_vals_mat, psi_vals_mat, A, B,
           C, D, temp_row);
    if (success != 0) {
        ans.x = 1;
        return ans;
    }
    success = sweep2(Cx_vals_mat, Cy_vals_mat, fx_vals_mat,
                     fy_vals_mat, gama_vals_mat, A, B,
                     C, D, temp_row);
    if (success != 0) {
        ans.x = 2;
        return ans;
    }

    for (i = 0; i < i_max+1; i++) {
        for (j = 0; j < j_max+1; j++) {
            index = offset2d(i, j, i_max+1);
            x_vals_mat_next[index] = Cx_vals_mat[index] + x_vals_mat_current[index];
            y_vals_mat_next[index] = Cy_vals_mat[index] + y_vals_mat_current[index];
        }
    }

    ans.x = calculate_max_L_x(x_vals_mat_current, alpha_vals_mat,
                              phi_vals_mat, beta_vals_mat, gama_vals_mat,
                              psi_vals_mat);
    ans.y = calculate_max_L_y(y_vals_mat_current, alpha_vals_mat,
                              phi_vals_mat, beta_vals_mat, gama_vals_mat,
                              psi_vals_mat);

    return ans;
}

/* printing a 1D array to the commend line
argument list:
data - 1D array */
void mat_print(double *data)
{
    int j_index, i_index;
    for (j_index = 0; j_index < j_max+1; j_index++) {
        for (i_index = 0; i_index < i_max+1; i_index++) {
            printf("%g ", data[offset2d(i_index, j_index, i_max+1)]);
        }
        printf("\n");
    }
}

/* returning the absolut maximum valus at the Lx matrix
argument list:
x_vals_mat - 1D array of the x valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus */
double calculate_max_L_x(double *x_vals_mat, double *alpha_vals_mat,
           double *phi_vals_mat, double *beta_vals_mat,
           double *gama_vals_mat, double *psi_vals_mat)
{
    int j_index, i_index;
    double max_L_x = 0, current_L_x;

    for (j_index = 0; j_index < j_max+1; j_index++) {
        for (i_index = 0; i_index < i_max+1; i_index++) {
            current_L_x = L_x(x_vals_mat, alpha_vals_mat,
                              phi_vals_mat, beta_vals_mat,
                              gama_vals_mat, psi_vals_mat, i_index, j_index);
            if (fabs(current_L_x) > max_L_x) {
                max_L_x = fabs(current_L_x);
            }
        }
    }
    return max_L_x;
}

/* returning the absolut maximum valus at the Ly matrix
argument list:
y_vals_mat - 1D array of the y valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus */
double calculate_max_L_y(double *y_vals_mat, double *alpha_vals_mat,
           double *phi_vals_mat, double *beta_vals_mat,
           double *gama_vals_mat, double *psi_vals_mat)
{
    int j_index, i_index;
    double max_L_y = 0, current_L_y;

    for (j_index = 0; j_index < j_max+1; j_index++) {
        for (i_index = 0; i_index < i_max+1; i_index++) {
            current_L_y = L_y(y_vals_mat, alpha_vals_mat,
                              phi_vals_mat, beta_vals_mat,
                              gama_vals_mat, psi_vals_mat, i_index, j_index);
            if (fabs(current_L_y) > max_L_y) {
                max_L_y = fabs(current_L_y);
            }
        }
    }
    return max_L_y;
}
