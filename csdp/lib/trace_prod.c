/*
  Compute the trace of the product of two symmetric matrices.  Since this
  computation is very memory bound, we exploit the symmetry to cut the 
  memory accesses in half.  The basic formula is that 
 
  tr(AB)=sum(sum(A(i,j)*B(i,j),i=1..n),j=1..n)

  or more efficiently 

  tr(AB)=2*sum(sum(A(i,j)*B(i,j),i=1..j-1),j=1..n)+sum(A(i,i)*B(i,i),i=1..n)

  */

#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

double trace_prod(A,B)
     struct blockmatrix A;
     struct blockmatrix B;
{
  int blk;
  double sum;
  double temp;
  int i;
  int j;

  sum=0.0;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    sum += A.blocks[blk].data.vec[i]*B.blocks[blk].data.vec[i];
	  break;
	case MATRIX:
	  temp=0;
#pragma omp parallel for schedule(dynamic,64) default(none) private(i,j) shared(A,B,blk) reduction(+:temp)
	  for (j=1; j<=A.blocks[blk].blocksize; j++)
	    for (i=1; i<j; i++)
	      temp=temp+A.blocks[blk].data.mat[ijtok(i,j,A.blocks[blk].blocksize)]
		*B.blocks[blk].data.mat[ijtok(i,j,A.blocks[blk].blocksize)];
	  sum += 2.0*temp;
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	      sum += A.blocks[blk].data.mat[ijtok(i,i,A.blocks[blk].blocksize)]
                    *B.blocks[blk].data.mat[ijtok(i,i,A.blocks[blk].blocksize)];
	  break;
	case PACKEDMATRIX:
	default:
	  printf("trace_prod illegal block type \n");
	  exit(206);
	};
    };
  return(sum);
}

double trace(A)
     struct blockmatrix A;
{
  int blk;
  double sum;
  int i;


  sum=0.0;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    sum += A.blocks[blk].data.vec[i];
	  break;
	case MATRIX:
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    sum += A.blocks[blk].data.mat[ijtok(i,i,A.blocks[blk].blocksize)];
	  break;
	case PACKEDMATRIX:
	default:
	  printf("trace_prod illegal block type \n");
	  exit(206);
	};
    };
  return(sum);
}

