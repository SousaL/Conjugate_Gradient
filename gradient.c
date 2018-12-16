#include "gradient.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

MATRIX * init_matrix(int rows, int columns){
  /* Init matrix
    Args:
      rows:  pointer to matrix to store the result
      columns: pointer to matrix A
    Return:
      matrix
  */
  double *p_r = (double*)malloc(sizeof(double) * rows * columns);
  double **p_rows = (double**)malloc(sizeof(double*) * rows);
  int i = 0, j = 0;
  for(i = 0; i < rows; i++){
    p_rows[i] = p_r + columns * i;
    for(j = 0; j < columns; j++){
      p_rows[i][j] = 0.0;
    }
  }
  MATRIX * matrix = (MATRIX*)malloc(sizeof(MATRIX));
  matrix->m = p_rows;
  matrix->rows = rows;
  matrix->columns = columns;
  return matrix;
}

void print(MATRIX * matrix){
  /* Multiply the matrix
    Args:
      dst:  pointer to matrix to store the result
      matrix_a: pointer to matrix A
      matrix_b: pointer to matrix B
    Return:
  */
  if(matrix == NULL){
    printf("(MATRIX NULL)\n" );
    return;
  }
  int i = 0, j = 0;
  printf("rows: %d columns: %d\n", matrix->rows, matrix->columns);
  for(i = 0; i < matrix->rows; i++){
    for(j = 0; j < matrix->columns; j++){
      printf(" %f", matrix->m[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

MATRIX * zeros(int rows){
  /* Create a vector-column with zeros
    Args:
      rows:  number of rows
    Return:
      vector: pointer to vector with zeros
  */
  return init_matrix(rows, 1);
}

MATRIX * mult_matrix(MATRIX * a, MATRIX  * b){
  /* Multiply the matrix
    Args:
      matrix_a: pointer to matrix A
      matrix_b: pointer to matrix B
    Return:
      dst:  pointer to matrix to store the result

  */
  if(a->columns != b->rows){
    return NULL;
  }
  MATRIX * dst = init_matrix(a->rows,b->columns);

  int rows = a->rows;
  int columns = b->columns;
  int elements = a->columns;
  int i = 0, j = 0, k = 0;
  for(i = 0; i < rows; i++){
    for(j = 0; j < columns; j++){
      for(k = 0; k < elements; k++){
      //  printf("a[%d][%d] * b[%d][%d] + ",i,k,k,j);
        dst->m[i][j] += a->m[i][k] * b->m[k][j];
      }
    }
  }
  return dst;
}

MATRIX * diff_matrix(MATRIX * a, MATRIX  * b){
  /* Multiply the matrix
    Args:
      matrix_a: pointer to matrix A
      matrix_b: pointer to matrix B
    Return:
      dst:  pointer to matrix to store the result

  */
  if(a->columns != b->columns || a->rows != b->rows){
    return NULL;
  }
  MATRIX * dst = init_matrix(a->rows,a->columns);

  int rows = a->rows;
  int columns = a->columns;
  int i = 0, j = 0;
  for(i = 0; i < rows; i++){
    for(j = 0; j < columns; j++){
      //printf("a[%d][%d] - b[%d][%d] + \n",i,j,i,j);
        dst->m[i][j] = a->m[i][j] - b->m[i][j];
    }
  }
  return dst;
}

MATRIX * sum_matrix(MATRIX * a, MATRIX  * b){
  /* Multiply the matrix
    Args:
      matrix_a: pointer to matrix A
      matrix_b: pointer to matrix B
    Return:
      dst:  pointer to matrix to store the result

  */
  if(a->columns != b->columns || a->rows != b->rows){
    return NULL;
  }
  MATRIX * dst = init_matrix(a->rows,a->columns);

  int rows = a->rows;
  int columns = a->columns;
  int i = 0, j = 0;
  for(i = 0; i < rows; i++){
    for(j = 0; j < columns; j++){
      //printf("a[%d][%d] - b[%d][%d] + \n",i,j,i,j);
        dst->m[i][j] = a->m[i][j] + b->m[i][j];
    }
  }
  return dst;
}

MATRIX * transpose(MATRIX * a){
  /* transpose the matrix
    Args:
      matrix_a: pointer to matrix A
    Return:
      new:  pointer to matrix to store the result
    */
     MATRIX * new = init_matrix(a->columns, a->rows);
    int i = 0;
    int j = 0;
    for(i = 0; i < a->rows; i++){
      for(j = 0; j < a->columns; j++){
        new->m[j][i] = a->m[i][j];
      }
    }
    return new;
}

MATRIX * scalar(double scalar, MATRIX * a){
  /* scalar the matrix
    Args:
      scalar: number to Multiply each element at matrix
      matrix_a: pointer to matrix A
    Return:
      new:  pointer to matrix to store the result
    */
     MATRIX * new = init_matrix(a->rows, a->columns);
    int i = 0;
    int j = 0;
    for(i = 0; i < a->rows; i++){
      for(j = 0; j < a->columns; j++){
        new->m[i][j] = scalar * a->m[i][j];
      }
    }
    return new;
}

MATRIX * copy(MATRIX * a){
  /* transpose the matrix
    Args:
      matrix_a: pointer to matrix A
    Return:
      new:  pointer to matrix to store the result
    */
    if(a == NULL) return NULL;
    MATRIX * new = init_matrix(a->rows, a->columns);
    int i = 0;
    int j = 0;
    for(i = 0; i < a->rows; i++){
      for(j = 0; j < a->columns; j++){
        new->m[i][j] = a->m[i][j];
      }
    }
    return new;
}

double first_value(MATRIX * a){
  /* transpose the matrix
    Args:
      matrix_a: pointer to matrix A
    Return:
      new:  pointer to matrix to store the result
    */
    return a->m[0][0];
}

MATRIX * gradiente(MATRIX * A, MATRIX * b){
  int imax = 5000;
  double erro = 0.00001;
  int n = MAX(A->rows,A->columns);
  int i = 1;

  MATRIX * r, * tmp;
  MATRIX * x = zeros(n);
  tmp = mult_matrix(A, x);
  r = diff_matrix(b, tmp);
  MATRIX * d = copy(r);
  double sigma_novo = first_value(mult_matrix(transpose(r), r));
  double sigma0 = sigma_novo;
  double sigma_velho;
  double beta;

  while(i < imax && sigma_novo > (erro * erro * sigma0)){

    MATRIX * q = mult_matrix(A,d);
    double alpha = sigma_novo/(first_value(mult_matrix(transpose(d),q)));
    x = sum_matrix(x, scalar(alpha, d));
    if(i % 50 == 0){
      r = diff_matrix(b,mult_matrix(A,x));
    }
    else{
      r = diff_matrix(r,scalar(alpha,q));
    }

    sigma_velho = sigma_novo;


    sigma_novo = first_value(mult_matrix(transpose(r), r));
    beta = sigma_novo/sigma_velho;
    d = sum_matrix(r, scalar(beta,d));
    i++;
  }
  return x;
}
