#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "gradient.h"
#include "hb_io.h"
#include <time.h>
void le_arquivo ( char *input_file )
{
  int *colptr = NULL;
  int indcrd;
  char *indfmt = NULL;
  FILE *input;
  char *key = NULL;
  int khi;
  int klo;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero;
  int nrhs;
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  char *rhstyp = NULL;
  int *rowind = NULL;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;
  double *values = NULL;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  HB_VALUES_READ reads the values of an HB file.\n" );
  printf ( "\n" );
  printf ( "  Reading the file '%s'\n", input_file );

  input = fopen ( input_file, "rt" );

  if ( !input )
  {
    printf ( "\n" );
    printf ( "TEST05 - Warning!\n" );
    printf ( "  Error opening the file.\n" );
    return;
  }

  printf ( "  Reading the header.\n" );

  hb_header_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd,
    &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt,
    &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix );

  colptr = ( int * ) malloc ( ( ncol + 1 ) * sizeof ( int ) );

  if ( mxtype[2] == 'A' )
  {
    rowind = ( int * ) malloc ( nnzero * sizeof ( int ) );
  }
  else if ( mxtype[2] == 'E' )
  {
    rowind = ( int * ) malloc ( neltvl * sizeof ( int ) );
  }
  else
  {
    printf ( "\n" );
    printf ( "TEST05 - Warning!\n" );
    printf ( "  Illegal value of MXTYPE character 3.\n" );
    return;
  }

  printf ( "  Reading the structure.\n" );

  hb_structure_read ( input, ncol, mxtype, nnzero, neltvl,
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );

  if ( mxtype[2] == 'A' )
  {
    values = ( double * ) malloc ( nnzero * sizeof ( double ) );
  }
  else if ( mxtype[2] == 'E' )
  {
    values =  ( double * ) malloc ( neltvl * sizeof ( double ) );
  }
  else
  {
    printf ( "\n" );
    printf ( "TEST05 - Warning!\n" );
    printf ( "  Illegal value of MXTYPE character 3 = '%c'\n", mxtype[2] );
    return;
  }

  printf ( "  Reading the values.\n" );

  hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values );

  fclose ( input );
  // ncol = 5;
  // nrow = 5;
  hb_values_print ( ncol, colptr, mxtype, nnzero, neltvl, values );
  MATRIX * m = init_matrix(nrow, ncol);
  MATRIX * b = zeros(ncol);
  int l;
  srand(time(NULL));
  for (l = 0; l < ncol; l++) {
    b->m[l][0] = (double)400;
  }


  printf("%d %d \n", nrow, ncol);
  int i = 0;
  int j = 0;
  int k = 0;
  while(i < nnzero){
    if(i+1 == colptr[j]){
      j++;
      //puts("");

    }
    m->m[j-1][rowind[i]-1] = values[i];
    m->m[rowind[i]-1][j-1] = values[i];
    //printf("%d %d |  ", j-1, rowind[i]-1);
    //printf("%d ", rowind[i++]);
    i++;
  }
  //return;
  //print(m);
  //return;
  // m->m[0][0] = 5;
  // m->m[0][1] = 2;
  // m->m[0][2] = 1;
  // m->m[0][3] = 0;
  // m->m[0][4] = 0;
  //
  // m->m[1][0] = 2;
  // m->m[1][1] = 5;
  // m->m[1][2] = 5;
  // m->m[1][3] = 1;
  // m->m[1][4] = 0;
  //
  // m->m[2][0] = 1;
  // m->m[2][1] = 2;
  // m->m[2][2] = 5;
  // m->m[2][3] = 2;
  // m->m[2][4] = 1;
  //
  //
  // m->m[3][0] = 0;
  // m->m[3][1] = 1;
  // m->m[3][2] = 2;
  // m->m[3][3] = 5;
  // m->m[3][4] = 2;
  //
  // m->m[4][0] = 0;
  // m->m[4][1] = 0;
  // m->m[4][2] = 1;
  // m->m[4][3] = 2;
  // m->m[4][4] = 5;


  MATRIX * x = gradiente(m,b);
  print(x);
  MATRIX * b_p = mult_matrix(m, x);
  print(b_p);
  //printf("qweq***\n");
  free ( colptr );
  free ( rowind );
  free ( values );

  return;
}

int main(){
  le_arquivo("bcsstk01.rsa");
}
