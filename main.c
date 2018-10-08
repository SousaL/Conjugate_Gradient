#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "gradient.h"
#include "hb_io.h"

void test05 ( char *input_file )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests HB_VALUES_READ;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 January 2014

  Author:

    John Burkardt
*/
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

  hb_values_print ( ncol, colptr, mxtype, nnzero, neltvl, values );
  MATRIX * m = init_matrix(nrow, ncol);
  MATRIX * b = zeros(ncol);
  
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
    //printf("%d %d |  ", j, rowind[i]-1);
    //printf("%d ", rowind[i++]);
    i++;
  }
  //print(m);

  MATRIX * x = gradiente(m,b);
  print(x);
  free ( colptr );
  free ( rowind );
  free ( values );

  return;
}

int main(){
  test05("bcsstk01.rsa");
}
