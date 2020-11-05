/*
 * Build an initial solution for an SDP problem.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "declarations.h"

void initsoln(n,k,C,a,constraints,pX0,py0,pZ0)
     int n;
     int k;
     struct blockmatrix C;
     double *a;
     struct constraintmatrix *constraints;
     struct blockmatrix *pX0;
     double **py0;
     struct blockmatrix *pZ0;
{
  int i,j,blk;
  double xi,eta;
  double minxi,mineta;
  double *normsofA;
  double normofC;
  double rescale;
  struct sparseblock *ptr;

  /*
   *  First, allocate storage for X0, y0, and Z0.
   */

  alloc_mat(C,pX0);
  alloc_mat(C,pZ0);
  *py0=(double *)malloc(sizeof(double)*(k+1));

  if (py0 == NULL)
    {
      printf("Storage allocation failed!\n");
      exit(205);
    };

  /*
   * Go ahead and set y to 0.
   */

  for (i=1; i<=k; i++)
    (*py0)[i]=0.0;

  /*
   * Allocate space to store the norms of the constraint matrices. 
   */

  normsofA=(double *)malloc(sizeof(double)*(k+1));
  if (normsofA == NULL)
    {
      printf("Storage allocation failed!\n");
      exit(205);
    };

  /*
   * Keep track of the minimum values of xi and eta.
   */

  minxi=1.0e300;
  mineta=1.0e300;
  
  /*
   * Next, work through the blocks one by one.  
   */

  for (blk=1; blk<=C.nblocks; blk++)
    {
      /*
       * Work out the norms of the objective and constraint blocks.
       */
      
      switch (C.blocks[blk].blockcategory)
	{
	case DIAG:
	  /*
 	   * Get the norm of this block of C.
	   */
	  
	  normofC=0;
	  for (i=1; i<=C.blocks[blk].blocksize; i++)
	    normofC=normofC+C.blocks[blk].data.vec[i]*C.blocks[blk].data.vec[i];
	  normofC=sqrt(normofC);
	  
	  /*
 	   * Get the norms of this block in constraints 1, 2, ...
	   */
	  
	  for (i=1; i<=k; i++)
	    {
	      normsofA[i]=0.0;
	      ptr=constraints[i].blocks;
	      while (ptr != NULL)
		{
		  if (ptr->blocknum == blk)
		    {
		      for (j=1; j<=ptr->numentries; j++)
			{
			  normsofA[i]=normsofA[i]+(ptr->entries[j])*(ptr->entries[j]);
			};
		      normsofA[i]=sqrt(normsofA[i]);
		      break;
		    };
		  ptr=ptr->next;
		};
	    };
	  break;
	case MATRIX:
	  /*
 	   * Get the norm of this block of C.
	   */
	  
	  normofC=0;
	  for (j=1; j<=C.blocks[blk].blocksize; j++)
	    for (i=1; i<=C.blocks[blk].blocksize; i++)
	      normofC=normofC+C.blocks[blk].data.mat[ijtok(i,j,C.blocks[blk].blocksize)]*C.blocks[blk].data.mat[ijtok(i,j,C.blocks[blk].blocksize)];
	  normofC=sqrt(normofC);
	  
	  /*
 	   * Get the norms of this block in constraints 1, 2, ...
	   */
	  
	  for (i=1; i<=k; i++)
	    {
	      normsofA[i]=0.0;
	      ptr=constraints[i].blocks;
	      while (ptr != NULL)
		{
		  if (ptr->blocknum == blk)
		    {
		      for (j=1; j<=ptr->numentries; j++)
			{
			  if (ptr->iindices[j] == ptr->jindices[j])
			    {
			      normsofA[i]=normsofA[i]+(ptr->entries[j])*(ptr->entries[j]);
			    }
			  else
			    {
			      normsofA[i]=normsofA[i]+2*(ptr->entries[j])*(ptr->entries[j]);
			    };
			};
		      normsofA[i]=sqrt(normsofA[i]);
		      break;
		    };
		  ptr=ptr->next;
		};
	    };
	  break;
	case PACKEDMATRIX:
	default:
	  printf("initsoln illegal block type \n");
	  exit(206);
	};
      
      /*
       * Now, setup X0 and Z0 for this block.
       */
  
      switch (C.blocks[blk].blockcategory)
	{
	case DIAG:
	  /*
	   * First, deal with the X0 block.
	   */
	  
	  xi=10;
	  
	  if (sqrt(C.blocks[blk].blocksize) > xi)
	    xi=sqrt(C.blocks[blk].blocksize);

	  for (i=1; i<=k; i++)
	    {
	      if (sqrt(C.blocks[blk].blocksize)*(1.0+fabs(a[i]))/(1+normsofA[i]) > xi)
		xi=sqrt(C.blocks[blk].blocksize)*(1.0+fabs(a[i]))/(1+normsofA[i]);
	    };

	  if (0==1)
	    printf("block=%d, xi=%e\n",blk,xi);

	  for (i=1; i<=C.blocks[blk].blocksize; i++)
	    (*pX0).blocks[blk].data.vec[i]=xi;

	  if (xi < minxi)
	    minxi=xi;
	  
	  /*
	   * Next, deal with the Z0 block.
	   */

	  eta=10;

	  if (sqrt(C.blocks[blk].blocksize) > eta)
	    eta=sqrt(C.blocks[blk].blocksize);

	  if (normofC > eta)
	    eta=normofC;
	  
	  for (i=1; i<=k; i++)
	    {
	      if (normsofA[i] > eta)
		eta=normsofA[i];
	    };

	  if (0==1)
	    printf("block=%d, eta=%e\n",blk,eta);

	  for (i=1; i<=C.blocks[blk].blocksize; i++)
	    (*pZ0).blocks[blk].data.vec[i]=eta;

	  if (eta < mineta)
	    mineta=eta;
	  
	  /*
	   * Done with current block.
	   */
	  break;
	case MATRIX:
	  /*
	   * First, deal with the X0 block.
	   */
	  
	  xi=10;
	  
	  if (sqrt(C.blocks[blk].blocksize) > xi)
	    xi=sqrt(C.blocks[blk].blocksize);

	  for (i=1; i<=k; i++)
	    {
	      if (C.blocks[blk].blocksize*(1.0+fabs(a[i]))/(1+normsofA[i]) > xi)
		xi=C.blocks[blk].blocksize*(1.0+fabs(a[i]))/(1+normsofA[i]);
	    };

	  if (0==1)
	    printf("block=%d, xi=%e\n",blk,xi);

	  for (j=1; j<=C.blocks[blk].blocksize; j++)
	    for (i=1; i<=C.blocks[blk].blocksize; i++)
	      (*pX0).blocks[blk].data.mat[ijtok(i,j,C.blocks[blk].blocksize)]=0.0;
	  
	  for (i=1; i<=C.blocks[blk].blocksize; i++)
	    (*pX0).blocks[blk].data.mat[ijtok(i,i,C.blocks[blk].blocksize)]=xi;

	  if (xi < minxi)
	    minxi=xi;
	  
	  /*
	   * Next, deal with the Z0 block.
	   */

	  eta=10;

	  if (sqrt(C.blocks[blk].blocksize) > eta)
	    eta=sqrt(C.blocks[blk].blocksize);

	  if (normofC > eta)
	    eta=normofC;
	  
	  for (i=1; i<=k; i++)
	    {
	      if (normsofA[i] > eta)
		eta=normsofA[i];
	    };

	  if (0==1)
	    printf("block=%d, eta=%e\n",blk,eta);

	  for (j=1; j<=C.blocks[blk].blocksize; j++)
	    for (i=1; i<=C.blocks[blk].blocksize; i++)
	      (*pZ0).blocks[blk].data.mat[ijtok(i,j,C.blocks[blk].blocksize)]=0.0;
	  
	  for (i=1; i<=C.blocks[blk].blocksize; i++)
	    (*pZ0).blocks[blk].data.mat[ijtok(i,i,C.blocks[blk].blocksize)]=eta;

	  if (eta < mineta)
	    mineta=eta;

	  /*
	   * Done with current block.
	   */
	  break;
	case PACKEDMATRIX:
	default:
	  printf("initsoln illegal block type \n");
	  exit(206);
	};
    };

  /*
   * Look for situations where X0 or Z0 is far too large and scale down 
   * as needed.
   */

  /*
   * First, look at A(X0) versus the right hand side.  We'll reuse 
   * normsofA as a work vector. If norm(A(X0) is 10 times bigger than 
   * norm(a), then scale it back down, but keep all the diagonal elements
   * at least 10.  
   */

  op_a(k,constraints,*pX0,normsofA);

  if (norm2(k,normsofA+1) > 10.0*norm2(k,a+1))
    {
      rescale=10*norm2(k,a+1)/norm2(k,normsofA+1);
      if (rescale*minxi < 10)
	rescale=10/minxi;
      if (0==1)
	printf("Rescaling X0 by factor of %e\n",rescale);	
      scalemat(rescale,*pX0,*pX0);
    };

  /*
   * Next, look at norm(C) versus norm(Z0) and rescale Z0 as with X0.
   */

  if (Fnorm(*pZ0) > 10*Fnorm(C))
    {
      rescale=10*Fnorm(C)/Fnorm(*pZ0);
      if (rescale*mineta < 10)
	rescale=10/mineta;
      if (0==1)
	printf("rescaling Z0 by %e\n",rescale);
      scalemat(rescale,*pZ0,*pZ0);
    };
  
  
  /*
   * Free up working storage.
   */

  free(normsofA);
  
}

