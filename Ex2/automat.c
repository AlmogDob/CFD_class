#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdarg.h>

#define dfprintINT(fp, expr) do{fprintf(fp, #expr "\n%d\n\n", expr);} while(0)     /* macro for easy debuging*/
#define dfprintD(fp, expr) do{fprintf(fp, #expr "\n%g\n\n", expr);} while(0)     /* macro for easy debuging*/

int create_empty_dir(char *parent_directory);
void print_command_to_file(FILE *fp, char *program, ...);
void create_input_file(char *file_name, int i_TEL, int i_LE, int i_TEU, 
                       int j_TEL, int j_LE, int j_TEU, double Mach_inf, 
                       double angle_of_attack_deg, double density, 
                       double environment_pressure, double delta_t,
                       double Gamma, double epse, double max_iteration);

int main()
{
    char parent_dir[] = "./auto_results";
    if (create_empty_dir("./auto_results") != 0) {
        return 1;
    }

    char temp_dir[BUFSIZ], temp1[BUFSIZ], temp_input[BUFSIZ], temp_num[BUFSIZ];

    strncpy(temp_dir, parent_dir, BUFSIZ);
    strncat(temp_dir, "/command_to_run.txt", BUFSIZ/2);
    FILE *fp = fopen(temp_dir, "wt");

    fprintf(fp, "make build_solver\n");

    for (int i = 1; i <= 20; i++) {
        strncpy(temp_dir, parent_dir, BUFSIZ);
        strncat(temp_dir, "/input", BUFSIZ/2);
        sprintf(temp1, "%d.txt", i);
        strncat(temp_dir, temp1, BUFSIZ/2);
        create_input_file(temp_dir, 11, 25, 39, 0, 0, 0, 1.5, 0, 1.225, 101325, (i*1)*1e-5, 1.4, 0.06, 1e6);


        strncpy(temp_dir, parent_dir, BUFSIZ);
        strncat(temp_dir, "/command_to_run.txt", BUFSIZ/2);

        strncpy(temp_input, parent_dir, BUFSIZ);
        strncat(temp_input, "/input", BUFSIZ/2);
        sprintf(temp1, "%d.txt", i);
        strncat(temp_input, temp1, BUFSIZ/2);

        sprintf(temp_num, "%d", i);

        print_command_to_file(fp,
                            "./solver",
                            temp_input,
                            "./mesh_output.txt",
                            "./auto_results",
                            temp_num,
                            NULL);

    }
    fprintf(fp, "make clean_solver\n");

    return 0;
}

/* if allready exisest, delet all the files inside 
returns 0 on success
this functin handls the errors so on fail just quit */
int create_empty_dir(char *parent_directory)
{
    char path_to_remove[BUFSIZ];

    if (mkdir(parent_directory, 0777) == -1) {
        if (errno == EEXIST) {
            DIR *dir = opendir(parent_directory);
            if (dir == NULL) {
                fprintf(stderr, "%s:%d: [Error] problem opening '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
                return 1;
            }
            struct dirent* entity;
            entity = readdir(dir);
            printf("\n");
            while (entity != NULL) {   /* directory is not empty */
                strncpy(path_to_remove, parent_directory, BUFSIZ);
                strncat(path_to_remove, "/", BUFSIZ/2);
                strncat(path_to_remove, entity->d_name, BUFSIZ);
                printf("%hhd: %s\n", entity->d_type, path_to_remove);
                if (entity->d_type == DT_REG) {
                    if (remove(path_to_remove) != 0) {
                        fprintf(stderr, "%s:%d: [Error] problem removing '%s': %s\n", __FILE__, __LINE__, path_to_remove, strerror(errno));
                        return 1;
                    }
                    printf("remove %s\n", path_to_remove);
                }
                entity = readdir(dir);
            }


            printf("\ndirectory already exist\n\n");

            closedir(dir);

            return 0;
        }

        fprintf(stderr, "%s:%d: [Error] problem making '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
        return 1;
    }
    return 0;
}

void print_command_to_file(FILE *fp, char *program, ...)
{
    fprintf(fp, "%s ", program);
    va_list args;
    va_start(args, program);
    char *temp_string;
    while ((temp_string = va_arg(args, char *)) != NULL) {
        fprintf(fp, "%s ", temp_string);
    }
    va_end(args);
    fprintf(fp, "\n");
}

void create_input_file(char *file_name, int i_TEL, int i_LE, int i_TEU, 
                       int j_TEL, int j_LE, int j_TEU, double Mach_inf, 
                       double angle_of_attack_deg, double density, 
                       double environment_pressure, double delta_t,
                       double Gamma, double epse, double max_iteration)
{
    FILE *input_fp = fopen(file_name, "wt");

    dfprintINT(input_fp, i_TEL);
    dfprintINT(input_fp, i_LE);
    dfprintINT(input_fp, i_TEU);
    dfprintINT(input_fp, j_TEL);
    dfprintINT(input_fp, j_LE);
    dfprintINT(input_fp, j_TEU);
    dfprintD(input_fp, Mach_inf);
    dfprintD(input_fp, angle_of_attack_deg);
    dfprintD(input_fp, density);
    dfprintD(input_fp, environment_pressure);
    dfprintD(input_fp, delta_t);
    dfprintD(input_fp, Gamma);
    dfprintD(input_fp, epse);
    dfprintD(input_fp, max_iteration);
}