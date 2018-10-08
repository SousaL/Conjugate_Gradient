#ifndef GRADIENT_H
#define GRADIENT_H



#define MAX(a,b) (((a)>(b))?(a):(b))
typedef struct matrix{
  double ** m;
  int rows;
  int columns;
} MATRIX;

MATRIX * init_matrix(int rows, int columns);

void print(MATRIX * matrix);

MATRIX * zeros(int rows);

MATRIX * mult_matrix(MATRIX * a, MATRIX  * b);

MATRIX * diff_matrix(MATRIX * a, MATRIX  * b);

MATRIX * sum_matrix(MATRIX * a, MATRIX  * b);

MATRIX * transpose(MATRIX * a);

MATRIX * scalar(double scalar, MATRIX * a);

MATRIX * copy(MATRIX * a);

double first_value(MATRIX * a);

MATRIX * gradiente(MATRIX * A, MATRIX * b);

#endif /* GRADIENT_H */
