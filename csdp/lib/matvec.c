/*
 *  Compute y=A*x
 *
 *  We use 
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

void matvec(A,x,y)
     struct blockmatrix A;
     double *x;
     double *y;
{
  int blk,i,n;
  int p;
  double *ap;
  double scale1;
  double scale2;
  int inc;
  
  /*
   * Work through the blocks one at a time.
   */

  p=1;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    {
	      y[p]=A.blocks[blk].data.vec[i]*x[p];
	      p++;
	    };
	  break;
	case MATRIX:
	  /*
	   * Call dgemm to do the matrix multiplication.
	   */

	  n=A.blocks[blk].blocksize;
	  ap=A.blocks[blk].data.mat;
          inc=1;

	  scale1=1.0;
	  scale2=0.0;

#ifdef HIDDENSTRLEN
	  dgemv_("N",&n,&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc,1);
#else
#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
	  DGEMV("N",&n,&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc);
#else
	  dgemv("N",&n,&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc);
#endif
#else
#ifdef CAPSBLAS
	  DGEMV_("N",&n,&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc);
#else
	  dgemv_("N",&n,&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc);
#endif
#endif
#endif
	  
	  p=p+n;

	  break;
	case PACKEDMATRIX:
	default:
	  printf("matvec illegal block type \n");
	  exit(206);
	};

    };
}


/*
 * This version of the routine computes y=A*x assuming that A is symmetric.
 * Since DGEMV is memory bound, DSYMV is about twice as fast.
 */

void matvecsym(A,x,y)
     struct blockmatrix A;
     double *x;
     double *y;
{
  int blk,i,n;
  int p;
  double *ap;
  double scale1;
  double scale2;
  int inc;
  
  /*
   * Work through the blocks one at a time.
   */

  p=1;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    {
	      y[p]=A.blocks[blk].data.vec[i]*x[p];
	      p++;
	    };
	  break;
	case MATRIX:
	  /*
	   * Call dgemm to do the matrix multiplication.
	   */

	  n=A.blocks[blk].blocksize;
	  ap=A.blocks[blk].data.mat;
          inc=1;

	  scale1=1.0;
	  scale2=0.0;

#ifdef HIDDENSTRLEN
	  dsymv_("U",&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc,1);
#else
#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
	  DSYMV("U",&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc);	  
#else
  	  dsymv("U",&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc);	  
#endif
#else
#ifdef CAPSBLAS
  	  DSYMV_("U",&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc);
#else
  	  dsymv_("U",&n,&scale1,ap,&n,x+p,&inc,&scale2,y+p,&inc);
#endif
#endif
#endif
	  
	  p=p+n;

	  break;
	case PACKEDMATRIX:
	default:
	  printf("matvec illegal block type \n");
	  exit(206);
	};

    };
}


/*
 * This version of the routine computes y=A*x assuming that A is upper 
 * triangular.  Since DGEMV is memory bound, DTRMV is about twice as fast.
 */

void matvecR(A,x,y)
     struct blockmatrix A;
     double *x;
     double *y;
{
  int blk,i,n;
  int p;
  double *ap;
  int inc;

  /*
   * Copy x into y.  Figure out the overall length of x and y first.
   */

  n=0;
  for (blk=1; blk<=A.nblocks; blk++)
    n=n+A.blocks[blk].blocksize;
  
  for (p=1; p<=n; p++)
    y[p]=x[p];
  
  /*
   * Work through the blocks one at a time.
   */

  p=1;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    {
	      y[p]=A.blocks[blk].data.vec[i]*x[p];
	      p++;
	    };
	  break;
	case MATRIX:
	  /*
	   * Call dgemm to do the matrix multiplication.
	   */

	  n=A.blocks[blk].blocksize;
	  ap=A.blocks[blk].data.mat;
          inc=1;

#ifdef HIDDENSTRLEN
          dtrmv_("U","N","N",&n,ap,&n,y+p,&inc,1,1,1);
#else
#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
          DTRMV("U","N","N",&n,ap,&n,y+p,&inc);
#else
          dtrmv("U","N","N",&n,ap,&n,y+p,&inc);
#endif
#else
#ifdef CAPSBLAS
          DTRMV_("U","N","N",&n,ap,&n,y+p,&inc);
#else
  	  dtrmv_("U","N","N",&n,ap,&n,y+p,&inc);
#endif
#endif
#endif
	  
	  p=p+n;

	  break;
	case PACKEDMATRIX:
	default:
	  printf("matvec illegal block type \n");
	  exit(206);
	};

    };
}

/*
 * This version of the routine computes y=A'*x assuming that A is upper 
 * triangular.  Since DGEMV is memory bound, DTRMV is about twice as fast.
 */

void matvecRT(A,x,y)
     struct blockmatrix A;
     double *x;
     double *y;
{
  int blk,i,n;
  int p;
  double *ap;
  int inc;

  /*
   * Copy x into y.  Figure out the overall length of x and y first.
   */

  n=0;
  for (blk=1; blk<=A.nblocks; blk++)
    n=n+A.blocks[blk].blocksize;
  
  for (p=1; p<=n; p++)
    y[p]=x[p];
  
  /*
   * Work through the blocks one at a time.
   */

  p=1;
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
	  for (i=1; i<=A.blocks[blk].blocksize; i++)
	    {
	      y[p]=A.blocks[blk].data.vec[i]*x[p];
	      p++;
	    };
	  break;
	case MATRIX:
	  /*
	   * Call dgemm to do the matrix multiplication.
	   */

	  n=A.blocks[blk].blocksize;
	  ap=A.blocks[blk].data.mat;
          inc=1;

#ifdef HIDDENSTRLEN
          dtrmv_("U","T","N",&n,ap,&n,y+p,&inc,1,1,1);
#else
#ifdef NOUNDERBLAS
#ifdef CAPSBLAS
          DTRMV("U","T","N",&n,ap,&n,y+p,&inc);
#else
          dtrmv("U","T","N",&n,ap,&n,y+p,&inc);
#endif
#else
#ifdef CAPSBLAS
          DTRMV_("U","T","N",&n,ap,&n,y+p,&inc);
#else
  	  dtrmv_("U","T","N",&n,ap,&n,y+p,&inc);
#endif
#endif
#endif
	  
	  p=p+n;

	  break;
	case PACKEDMATRIX:
	default:
	  printf("matvec illegal block type \n");
	  exit(206);
	};

    };
}


