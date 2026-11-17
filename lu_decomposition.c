#include<stdio.h>
#include<stdlib.h>
#include <windows.h>
#define Abs(a) ((a)>=0?(a):-(a))

static int total_variable(FILE *fp) {
    int n=0;
    char line[5000];
    while (fgets(line,5000,fp)!=NULL) {
        n++;
    }
    rewind(fp);
    return n;
}

static void lu_decomposition(FILE *fp, FILE *fo, double *soln) {
    int n= total_variable(fp);

    double **matrix= malloc(sizeof(double*)*n);
    if (!matrix) {
        return;
    }
    for (int i=0; i<n; i++) {
        matrix[i]=malloc(sizeof(double)*(n+1));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            fscanf(fp,"%lf\t",&matrix[i][j]);
        }
    }
    rewind(fp);

    for (int step=0; step<n-1; step++) {
        double col_max=Abs(matrix[step][step]);

        for (int i=step+1; i<n; i++) {
            if (Abs(matrix[i][step])>Abs(col_max)) {
                col_max=Abs(matrix[i][step]);
            }
        }

        if (col_max<1e-10) {
            printf("Singular Matrix. No unique Solution!\n");
            for (int ii=0; ii<n; ii++) {
                free(matrix[ii]);
            }
            free(matrix);
            return;
        }

        if (Abs(matrix[step][step])!=col_max) {
            for (int i=step+1; i<n; i++) {
                if (col_max==Abs(matrix[i][step])) {
                    double *temp= malloc(sizeof(double)*(n+1));
                    if (!temp) {
                        for (int ii=0; ii<n; ii++) {
                            free(matrix[ii]);
                        }
                        free(matrix);
                        return;
                    }

                    for (int it=0; it<n+1; it++) {
                        temp[it]=matrix[i][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[i][it]=matrix[step][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[step][it]=temp[it];
                    }
                    // pivoting done. now gaussian kora need
                    free(temp);
                    break;
                }
            }
        }

        // gaussian elimination here
        for (int i=step+1; i<n; i++) {
            double factor= matrix[i][step]/matrix[step][step];
            for(int j=step; j<n; j++) {
                if (j==step) {
                    matrix[i][j]= factor;
                }
                else {
                    matrix[i][j]+= -factor*matrix[step][j];
                }
            }
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            fprintf(fo,"%.6lf\t",matrix[i][j]);
        }
        if (i!=n-1)
            fprintf(fo,"\n");
    }

    double *y= malloc(sizeof(double)*n);
    if (!y) {
        for (int i=0; i<n; i++) {
            free(matrix[i]);
        }
        free(matrix);
        return;
    }
    for (int i=0; i<n; i++) {
        y[i]=0;
    }

    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<i; j++) {
            s+=matrix[i][j]*y[j];
        }
        y[i]=matrix[i][n]-s;
    }

    for (int i=0; i<n; i++) {
        soln[i]=0;
    }

    for (int i=n-1; i>=0; i--) {
        double s=0;
        for (int j=n-1; j>=i+1; j--) {
            s+=matrix[i][j]*soln[j];
        }
        soln[i]=(y[i]-s)/matrix[i][i];
    }

    // free matrix
    for (int i=0; i<n; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(y);
}

static void lu_decomposition_without_file(double **matrix, double *soln, const int n) {
    for (int step=0; step<n-1; step++) {
        double col_max=Abs(matrix[step][step]);

        for (int i=step+1; i<n; i++) {
            if (Abs(matrix[i][step])>Abs(col_max)) {
                col_max=Abs(matrix[i][step]);
            }
        }

        if (col_max<1e-10) {
            printf("Singular Matrix. No unique Solution!\n");
            return;
        }

        if (Abs(matrix[step][step])!=col_max) {
            for (int i=step+1; i<n; i++) {
                if (col_max==Abs(matrix[i][step])) {
                    double *temp= malloc(sizeof(double)*(n+1));
                    if (!temp) {
                        for (int ii=0; ii<n; ii++) {
                            free(matrix[ii]);
                        }
                        free(matrix);
                        return;
                    }

                    for (int it=0; it<n+1; it++) {
                        temp[it]=matrix[i][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[i][it]=matrix[step][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[step][it]=temp[it];
                    }
                    // pivoting done. now gaussian kora need
                    free(temp);
                    break;
                }
            }
        }

        // gaussian elimination here
        for (int i=step+1; i<n; i++) {
            double factor= matrix[i][step]/matrix[step][step];
            for(int j=step; j<n; j++) {
                if (j==step) {
                    matrix[i][j]= factor;
                }
                else {
                    matrix[i][j]+= -factor*matrix[step][j];
                }
            }
        }
    }

    if (Abs(matrix[n-1][n-1])<1e-10) {
        printf("Singular Matrix. No unique Solution!\n");
        return;
    }

    double *y= malloc(sizeof(double)*n);
    if (!y) {
        return;
    }
    for (int i=0; i<n; i++) {
        y[i]=0;
    }

    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<i; j++) {
            s+=matrix[i][j]*y[j];
        }
        y[i]=matrix[i][n]-s;
    }

    for (int i=0; i<n; i++) {
        soln[i]=0;
    }

    for (int i=n-1; i>=0; i--) {
        double s=0;
        for (int j=n-1; j>=i+1; j--) {
            s+=matrix[i][j]*soln[j];
        }
        soln[i]=(y[i]-s)/matrix[i][i];
    }

    free(y);
}

int main() {
    FILE *fp=fopen("5x6_augmented_matrix.txt","r");
    if (!fp) {
        return 1;
    }

    FILE *fo= fopen("output_matrix_5x6.txt","w");
    if (!fo) {
        fclose(fp);
        return 1;
    }

    // int n=2;
    // double temp[2][3]={
    //     {1,1,3},
    //     {1,-1,1}
    // };
    // double **matrix= malloc(sizeof(double*)*2);
    // for (int i=0; i<n+1; i++) {
    //     matrix[i]=malloc(sizeof(double));
    // }
    // for (int i=0; i<n; i++) {
    //     for (int j=0; j<n+1; j++) {
    //         matrix[i][j]=temp[i][j];
    //     }
    // }

    double *soln= malloc(sizeof(double)*total_variable(fp));
    if (!soln) {
        fclose(fp);
        return 1;
    }

    // double *soln= malloc(sizeof(double)*n);
    // if (!soln) {
    //     for (int i=0; i<n; i++) {
    //         free(matrix[i]);
    //     }
    //     free(matrix);
    // }

    lu_decomposition(fp,fo,soln);
    // lu_decomposition_without_file(matrix,soln,n);

    for (int i=0; i<total_variable(fp); i++) {
        printf("%.6lf\n",soln[i]);
    }
    // for (int i=0; i<n; i++) {
    //     printf("%.6lf\n",soln[i]);
    // }

    free(soln);
    fclose(fp);
    fclose(fo);
}