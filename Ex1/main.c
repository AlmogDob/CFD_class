#include <stdio.h>
#define MATRIX_IMPLEMENTATION
#include "Matrix_Double.h"  /* I included this library for debuging */
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#define MAXDIR 100
#define MAXWORD 100
#define PI 3.14159265359
#define dprintSTRING(expr) printf(#expr " = %s\n", expr)    /* macro for easy debuging*/
#define dprintINT(expr) printf(#expr " = %d\n", expr)
#define dprintF(expr) printf(#expr " = %g\n", expr)
#define dprintD(expr) printf(#expr " = %g\n", expr)

void read_input(char *dir);
void output_solution_d(char *dir, double *data);
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
           double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat);
int sweep2(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat,
           double *fy_vals_mat, double *gama_vals_mat);
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
int step(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat,
         double *fy_vals_mat, double *x_vals_mat_current,
         double *x_vals_mat_next, double *y_vals_mat_current,
         double *y_vals_mat_next, double *alpha_vals_mat,
         double *phi_vals_mat, double *beta_vals_mat,
         double *gama_vals_mat, double *psi_vals_mat);

/* Input variables */
double t, delta_x, delta_y, XSF, YSF, x_int = 1, r, omega;
int i_max, j_max, i_TEL, i_LE, i_TEU, i_min = 0, j_min = 0;


int main(int argc, char const *argv[])
{
    /* decleraitons */
    char input_dir[MAXDIR], output_dir[MAXDIR];
    int i_index, j_index, success = 0;
    double *x_vals_mat_init, *y_vals_mat_init, *x_vals_mat_current,
           *y_vals_mat_current, *x_vals_mat_next, *y_vals_mat_next,
           *alpha_vals_mat, *beta_vals_mat, *gama_vals_mat, *psi_vals_mat,
           *phi_vals_mat, *fx_vals_mat, *fy_vals_mat, *Cx_vals_mat,
           *Cy_vals_mat;

    /*------------------------------------------------------------*/

    /* Geting the input and output directories */
    if (--argc != 2) {
        fprintf(stderr, "ERROR: not right usage\nUsage: main 'input dir' 'output dir'\n");
        return -1;
    }

    strncpy(input_dir, (*(++argv)), MAXDIR);
    strncpy(output_dir, (*(++argv)), MAXDIR);

    if (input_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "Error: input too long\n");
        return -1;
    }
    if (output_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "Error: output too long\n");
        return -1;
    }

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

    /*------------------------------------------------------------*/
    
    initialize(x_vals_mat_init, y_vals_mat_init, alpha_vals_mat, beta_vals_mat, gama_vals_mat,
               psi_vals_mat, phi_vals_mat);

    copy_mat(x_vals_mat_current, x_vals_mat_init);
    copy_mat(x_vals_mat_next, x_vals_mat_init);
    copy_mat(y_vals_mat_current, y_vals_mat_init);
    copy_mat(y_vals_mat_next, y_vals_mat_init);

    for (i_index = 0; i_index < 1; i_index++) {
        success = step(Cx_vals_mat, Cy_vals_mat, fx_vals_mat, fy_vals_mat, x_vals_mat_current,
            x_vals_mat_next, y_vals_mat_current, y_vals_mat_next, alpha_vals_mat,
            phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat);
        if (success == 1) {
            fprintf(stderr, "ERROR: Step - sweep 1\n");
            exit(1);
        }
        if (success == 2) {
            fprintf(stderr, "ERROR: Step - sweep 2\n");
            exit(2);
        }

        copy_mat(x_vals_mat_current, x_vals_mat_next);
        copy_mat(y_vals_mat_current, y_vals_mat_next);
    }
    /*------------------------------------------------------------*/

    Mat xmat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = x_vals_mat_init}; /* for debuging */
    Mat ymat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = y_vals_mat_init}; 
    Mat xCurrentMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = x_vals_mat_current}; 
    Mat yCurrentMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = y_vals_mat_current}; 
    Mat xnextMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = x_vals_mat_next}; 
    Mat ynextMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = y_vals_mat_next}; 
    Mat fxMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = fx_vals_mat}; 
    Mat fyMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = fy_vals_mat}; 
    Mat CxMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = Cx_vals_mat}; 
    Mat CyMat = {.cols = i_max+1, .rows = j_max+1, .stride = i_max+1, .elements = Cy_vals_mat}; 

    // for (i_index = 0; i_index < i_max+1; i_index++) {
    //     alpha_vals_mat[offset2d(i_index, j_max, i_max+1)] = 3;
    // }

    // MAT_PRINT(xmat);
    // MAT_PRINT(ymat);
    MAT_PRINT(CxMat);

    FILE *fp = fopen("x_mat_current.txt", "wt");
    mat_print_to_file(fp, xCurrentMat, "");
    fclose(fp);

    fp = fopen("y_mat_current.txt", "wt");
    mat_print_to_file(fp, yCurrentMat, "");
    fclose(fp);

    fp = fopen("x_mat_next.txt", "wt");
    mat_print_to_file(fp, xnextMat, "");
    fclose(fp);

    fp = fopen("y_mat_next.txt", "wt");
    mat_print_to_file(fp, ynextMat, "");
    fclose(fp);

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
        }
    }

    delta_x = 1.0/(i_LE - i_TEL);

    fclose(fp);
}

/* output data; double version
argument list:
dir - the directory of the output file.
data - the solution vector */
void output_solution_d(char *dir, double *data)
{
    FILE *fp = fopen(dir, "wt");
    
    data = (void *)data;

    fclose(fp);
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
argument list:
x_vals_mat - 1D array of the x valsuse 
y_vals_mat - 1D array of the y valsuse */
void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat)
{
    int i_index, i_min = 0, j_index, j_min = 0, index = 0, num_points_befor_circle,
    num_of_outer_segments, num_of_top_outer_segments;

    double x, x_i_minos_2, x_i_minos_1, y_j_minos_2, y_j_minos_1, y_imax_jmax, x_imax_jmax,
    delta_theta, R, theta = 0, length_top_outer, segment_length, current_x_vals = 0;

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
        x_vals_mat[offset2d(i_max, j_index, i_max+1)] = x_vals_mat[offset2d(i_max, i_min, i_max+1)];
        y_vals_mat[offset2d(i_min, j_index, i_max+1)] = -y_vals_mat[offset2d(i_max, j_index, i_max+1)];
        x_vals_mat[offset2d(i_min, j_index, i_max+1)] = x_vals_mat[offset2d(i_min, j_min, i_max+1)];
    }
    /* Outer boundary */
    y_imax_jmax = y_vals_mat[offset2d(i_max, j_max, i_max+1)];
    x_imax_jmax = x_vals_mat[offset2d(i_max, j_max, i_max+1)];
    R = y_imax_jmax;

    num_of_outer_segments = i_max;
    num_of_top_outer_segments = num_of_outer_segments/2;
    length_top_outer = x_imax_jmax + 0.5*PI*R; /* length of stright part and quarter of the circle */
    segment_length = length_top_outer/num_of_top_outer_segments;

    /* the stright line part */
    for (num_points_befor_circle = 0;
         num_points_befor_circle < num_of_top_outer_segments + 1;
         num_points_befor_circle++) {
            current_x_vals = x_imax_jmax - num_points_befor_circle*segment_length; 
            if (current_x_vals < 0) {
                break;
            }
            x_vals_mat[offset2d(i_max-num_points_befor_circle, j_max, i_max+1)] = current_x_vals;
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
argument list:
x - x positon 
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
        return (mat[offset2d(i, j+1, i_max+1)] -2*mat[offset2d(i, j, i_max)] + mat[offset2d(i, j-1, i_max+1)]) / (1); /* second order second derivitive */
    }
    if (diraction == 'i') {
        if (i == i_min || i == i_max) {
            return 0;
        }
        return (mat[offset2d(i+1, j, i_max+1)] -2*mat[offset2d(i, j, i_max)] + mat[offset2d(i-1, j, i_max+1)]) / (1); /* second order second derivitive */
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

    /* eq 4 */
    for (j = 0; j < j_max+1; j++) {
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
    }

    /* eq 5 */
    for (i = 0; i < i_max+1; i++) {
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
    
    for (j = 0; j < i_max+1; j++) {
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
           double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat)
{
    int i_index, j_index, success = 0;

    /* memory allocating */
    double *Ax = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        Ax[i_index] = 0;
    }
    double *Bx = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        Bx[i_index] = 0;
    }
    double *Cx = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        Cx[i_index] = 0;
    }
    double *Dx = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        Dx[i_index] = 0;
    }
    double *Ay = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        Ay[i_index] = 0;
    }
    double *By = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        By[i_index] = 0;
    }
    double *Cy = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        Cy[i_index] = 0;
    }
    double *Dy = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        Dy[i_index] = 0;
    }
    double *temp_row = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        temp_row[i_index] = 0;
    }

    /* solving for each j */
    for (j_index = 0; j_index < j_max+1; j_index++) {
        LHS_sweep1(Ax, Bx, Cx, alpha_vals_mat, j_index);
        RHS_sweep1_x(Dx, x_vals_mat, alpha_vals_mat, phi_vals_mat,
                     beta_vals_mat, gama_vals_mat, psi_vals_mat, j_index);
        BC_sweep1(Ax, Bx, Cx, Dx);
        success = tridiag(Ax, Bx, Cx, Dx, temp_row, 0, i_max);
        if (success == 1) {
            break;
        }
        copy_row_to_mat(fx_vals_mat, temp_row, j_index);

        LHS_sweep1(Ay, By, Cy, alpha_vals_mat, j_index);
        RHS_sweep1_y(Dy, y_vals_mat, alpha_vals_mat, phi_vals_mat,
                     beta_vals_mat, gama_vals_mat, psi_vals_mat, j_index);
        BC_sweep1(Ay, By, Cy, Dy);
        success = tridiag(Ay, By, Cy, Dy, temp_row, 0, i_max);
        if (success == 1) {
            break;
        }
        copy_row_to_mat(fy_vals_mat, temp_row, j_index);
    }
    if (success == 1) {
        return 1;
    }

    free(Ax);
    free(Bx);
    free(Cx);
    free(Dx);
    free(Ay);
    free(By);
    free(Cy);
    free(Dy);
    free(temp_row);

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
           double *fy_vals_mat, double *gama_vals_mat)
{
    int i_index, j_index, success = 0;

    /* memory allocating */
    double *Ax2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        Ax2[j_index] = 0;
    }
    double *Bx2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        Bx2[j_index] = 0;
    }
    double *Cx2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        Cx2[j_index] = 0;
    }
    double *Dx2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        Dx2[j_index] = 0;
    }
    double *Ay2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        Ay2[j_index] = 0;
    }
    double *By2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        By2[j_index] = 0;
    }
    double *Cy2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        Cy2[j_index] = 0;
    }
    double *Dy2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        Dy2[j_index] = 0;
    }
    double *temp_row2 = (double *)malloc(sizeof(double) * (j_max + 1));
    for (j_index = 0; j_index < j_max+1; j_index++) {   /* filling the array with zeros */
        temp_row2[j_index] = 0;
    }

    /* solving for each i */
    for (i_index = 0; i_index < i_max+1; i_index++) {
        LHS_sweep2(Ax2, Bx2, Cx2, gama_vals_mat, i_index);
        RHS_sweep2_x(Dx2, fx_vals_mat, i_index);
        BC_sweep2(Ax2, Bx2, Cx2, Dx2);
        /*test*/
        // for (int index = 0; index < j_max+1; index++) {
        //     printf("%g\n", Bx2[index]);
        // }
        /*test*/
        success = tridiag(Ax2, Bx2, Cx2, Dx2, temp_row2, 0, j_max);
        if (success == 1) {
            printf("1\n");
            break;
        }
        copy_col_to_mat(Cx_vals_mat, temp_row2, i_index);

        LHS_sweep2(Ay2, By2, Cy2, gama_vals_mat, i_index);
        RHS_sweep2_y(Dy2, fy_vals_mat, i_index);
        BC_sweep2(Ay2, By2, Cy2, Dy2);
        success = tridiag(Ay2, By2, Cy2, Dy2, temp_row2, 0, j_max);
        if (success == 1) {
            printf("2\n");
            break;
        }
        copy_col_to_mat(Cy_vals_mat, temp_row2, i_index);
    }
    if (success == 1) {
        return 1;
    }
    
    printf("helo\n");
    free(Ax2);
    free(Bx2);
    free(Cx2);
    free(Dx2);
    free(Ay2);
    free(By2);
    free(Cy2);
    free(Dy2);
    free(temp_row2);
        printf("here\n");

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
returns 0 on success
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
psi_vals_mat - 1D array of the psi valus */
int step(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat,
         double *fy_vals_mat, double *x_vals_mat_current,
         double *x_vals_mat_next, double *y_vals_mat_current,
         double *y_vals_mat_next, double *alpha_vals_mat,
         double *phi_vals_mat, double *beta_vals_mat,
         double *gama_vals_mat, double *psi_vals_mat)
{
    int success, i, j, index; 

    success = sweep1(fx_vals_mat, fy_vals_mat, x_vals_mat_current,
           y_vals_mat_current, alpha_vals_mat, phi_vals_mat,
           beta_vals_mat, gama_vals_mat, psi_vals_mat);
    if (success != 0) {
        return 1;
    }
    success = sweep2(Cx_vals_mat, Cy_vals_mat, fx_vals_mat, fy_vals_mat, gama_vals_mat);
    if (success != 0) {
        return 2;
    }

    for (i = 0; i < i_max+1; i++) {
        for (j = 0; j < j_max+1; j++) {
            index = offset2d(i, j, i_max+1);
            x_vals_mat_next[index] = Cx_vals_mat[index] + x_vals_mat_current[index];
            y_vals_mat_next[index] = Cy_vals_mat[index] + y_vals_mat_current[index];
        }
    }
    
    return 0;
}
