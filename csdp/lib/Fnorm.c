/*
  Compute the Frobenius norm of a symmetric matrix.  Note that this routine
  requires symmetry and saves memory accesses by only using the upper triangle
  of the matrix.
  */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "declarations.h"

double Fnorm(A)
     struct blockmatrix A;
{
  int blk;
  double nrm;
  double temp;
  int i,j;
  
  nrm=0;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory)
	{
	case DIAG:
	  temp=norm2(A.blocks[blk].blocksize,A.blocks[blk].data.vec+1);
	  nrm += temp*temp;
	  break;
	case MATRIX:
	  temp=0;
#pragma omp parallel for schedule(dynamic,64) default(none) private(i,j) shared(blk,A) reduction(+:temp)
          for (j=1; j<=A.blocks[blk].blocksize; j++)
#pragma omp simd
	    for (i=1; i<j; i++)
	      temp+=A.blocks[blk].data.mat[ijtok(i,j,A.blocks[blk].blocksize)]*A.blocks[blk].data.mat[ijtok(i,j,A.blocks[blk].blocksize)];
	  temp=2.0*temp;
#pragma omp simd
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    temp+=A.blocks[blk].data.mat[ijtok(i,i,A.blocks[blk].blocksize)]*A.blocks[blk].data.mat[ijtok(i,i,A.blocks[blk].blocksize)];
	  nrm += temp;
	  break;
	case PACKEDMATRIX:
	default:
	  printf("Fnorm illegal block type \n");
	  exit(206);
	};
    };

  nrm=sqrt(nrm);
  return(nrm);
}

/*
 * The Knorm is the sum of the Fnorms of the blocks.
 */

double Knorm(A)
     struct blockmatrix A;
{
  int blk;
  double nrm;
  double temp;

  nrm=0;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory)
	{
	case DIAG:
	  temp=norm2(A.blocks[blk].blocksize,A.blocks[blk].data.vec+1);
	  nrm += temp;
	  break;
	case MATRIX:
	  temp=norm2(A.blocks[blk].blocksize*A.blocks[blk].blocksize,
		     A.blocks[blk].data.mat);
	  nrm += temp;
	  break;
	case PACKEDMATRIX:
	default:
	  printf("Fnorm illegal block type \n");
	  exit(206);
	};
    };

  return(nrm);
}


double mat1norm(A)
     struct blockmatrix A;
{
  int blk;
  double nrm;
  double temp;

  nrm=0;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory)
	{
	case DIAG:
	  temp=norm1(A.blocks[blk].blocksize,A.blocks[blk].data.vec+1);
	  nrm += temp;
	  break;
	case MATRIX:
	  temp=norm1(1*A.blocks[blk].blocksize*A.blocks[blk].blocksize,
		     A.blocks[blk].data.mat);
	  nrm += temp;
	  break;
	case PACKEDMATRIX:
	default:
	  printf("Fnorm illegal block type \n");
	  exit(206);
	};
    };

  return(nrm);
}


double matinfnorm(A)
     struct blockmatrix A;
{
  int blk;
  int i;
  double nrm;

  nrm=0;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory)
	{
	case DIAG:
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    {
	      if (fabs(A.blocks[blk].data.vec[i]) > nrm)
		nrm=fabs(A.blocks[blk].data.vec[i]);
	    };
	  break;
	case MATRIX:
	  for (i=0; i<A.blocks[blk].blocksize*A.blocks[blk].blocksize; i++)
	    {
	      if (fabs(A.blocks[blk].data.vec[i]) > nrm)
		nrm=fabs(A.blocks[blk].data.vec[i]);
	    };
	  break;
	case PACKEDMATRIX:
	default:
	  printf("Fnorm illegal block type \n");
	  exit(206);
	};
    };

  return(nrm);
}

