#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "gradient.h"
#include "hb_io.h"
#include <time.h>
#include <mpi.h>

void le_arquivo ( char *input_file )
{
  int id, np;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

  MATRIX * m, * b;
  if(id == 0){
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

    input = fopen ( input_file, "rt" );

    if ( !input )
    {
      printf ( "\n" );
      printf ( "TEST05 - Warning!\n" );
      printf ( "  Error opening the file.\n" );
      return;
    }

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

        return;
      }

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
          return;
        }

        hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values );

        fclose ( input );
        // ncol = 5;
        // nrow = 5;
        //hb_values_print ( ncol, colptr, mxtype, nnzero, neltvl, values );
        m = init_matrix(nrow, ncol);
        b = zeros(ncol);
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


        free ( colptr );
        free ( rowind );
        free ( values );
        int dest;

        for(dest = 1; dest < np; dest++){
          MPI_Send(&m->rows, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
          MPI_Send(&m->columns, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
          MPI_Send(&(m->m[0][0]), m->rows * m->columns, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);

          MPI_Send(&ncol, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

      }
  }
  MPI_Status status;
  if(id > 0){
    int m_r, m_c, b_c;

    MPI_Recv(&m_r, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&m_c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    double *m_m = (double*)malloc(sizeof(double) * m_r * m_c);

    MPI_Recv(m_m, m_r * m_c, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    m = init_matrix(m_c, m_r);
    memcpy(&m->m[0][0], m_m, m_r * m_c * sizeof(double));

    MPI_Recv(&b_c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    b = zeros(b_c);
    int lo;
    for (lo = 0; lo < b_c; lo++) {
      b->m[lo][0] = (double)400;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MATRIX * x = gradiente(m,b);
  //print(x);
  MATRIX * b_p = mult_matrix(m, x);
  if(id == 0) print(b_p);
  //printf("qweq***\n");
  return;
}

int main(int argc, char **argv){
  MPI_Init(&argc, &argv);
  le_arquivo("bcsstk01.rsa");
  MPI_Finalize();
}
