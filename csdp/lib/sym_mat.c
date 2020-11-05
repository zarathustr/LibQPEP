/*
  Symmetrize a matrix in Fortran storage format.  
  */

#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

void sym_mat(A)
     struct blockmatrix A;
{
  int i;
  int j;
  int blk;
  double foo;
  double *ap;
  int n;

  /*
   * Loop through the blocks, symmetrizing one at a time.
   */

  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory)
	{
	case DIAG:
	  break;
	case MATRIX:
	  n=A.blocks[blk].blocksize;
	  ap=A.blocks[blk].data.mat;
#pragma omp parallel for schedule(dynamic,64) default(none) private(i,j,foo) shared(ap,n)
	  for (j=1; j<=n; j++)
	    for (i=1; i<j; i++)
	      {
		foo=(ap[ijtok(i,j,n)]+ap[ijtok(j,i,n)])/2.0;
		ap[ijtok(i,j,n)]=foo;
		ap[ijtok(j,i,n)]=foo;
	      };
	  break;
	case PACKEDMATRIX:
	default:
	  printf("sym_mat illegal block type \n");
	  exit(206);
	};
    };

}


