/*
 *  This is an easy to call version of the sdp routine.  It takes as
 *  input a problem (n,k,C,a,constraints,constant_offset), and an 
 *  initial solution (X,y,Z), allocates working storage, and calls sdp() 
 *  to solve the problem.  The solution is returned in X,y,Z,pobj,dobj, and 
 *  the return code from sdp is returned as the return value from easy_sdp.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "declarations.h"

int easy_sdp(n,k,C,a,constraints,constant_offset,pX,py,pZ,ppobj,pdobj)
     int n;
     int k;
     struct blockmatrix C;
     double *a;
     struct constraintmatrix *constraints;
     double constant_offset;
     struct blockmatrix *pX;
     double **py;
     struct blockmatrix *pZ;
     double *ppobj;
     double *pdobj;
{
  int ret;
  struct constraintmatrix fill;
  struct paramstruc params;
  struct blockmatrix work1;
  struct blockmatrix work2;
  struct blockmatrix work3;
  struct blockmatrix bestx;
  struct blockmatrix bestz;
  struct blockmatrix Zi;
  struct blockmatrix dZ;
  struct blockmatrix dX;
  struct blockmatrix cholxinv;
  struct blockmatrix cholzinv;
  double *workvec1;
  double *workvec2;
  double *workvec3;
  double *workvec4;
  double *workvec5;
  double *workvec6;
  double *workvec7;
  double *workvec8;
  double *diagO; 
  double *Fp;
  double *O;
  double *dy;
  double *dy1;
  double *rhs;
  double *besty;
  int printlevel;
  int ldam;
  struct sparseblock **byblocks;
  struct sparseblock *ptr;
  struct sparseblock *oldptr;
  int i;
  int j;
  struct sparseblock *p;
  struct sparseblock *q;
  struct sparseblock *prev=NULL;
  double gap;
  int nnz;
  int denseblocks;
  int numblocks;

   /*
    *  Initialize the parameters.
    */

   initparams(&params,&printlevel);


  /*
   *  Allocate working storage
   */

  alloc_mat(C,&work1);
  alloc_mat(C,&work2);
  alloc_mat(C,&work3);
  alloc_mat_packed(C,&bestx);
  alloc_mat_packed(C,&bestz);
  alloc_mat_packed(C,&cholxinv);
  alloc_mat_packed(C,&cholzinv);

  besty=(double *)malloc(sizeof(double)*(k+1));
   if (besty == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   if (n > k)
     {
       workvec1=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec1=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec1 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   if (n > k)
     {
       workvec2=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec2=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec2 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   if (n > k)
     {
       workvec3=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec3=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec3 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   if (n > k)
     {
       workvec4=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec4=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec4 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   if (n > k)
     {
       workvec5=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec5=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec5 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   if (n > k)
     {
       workvec6=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec6=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec6 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   if (n > k)
     {
       workvec7=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec7=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec7 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   if (n > k)
     {
       workvec8=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       workvec8=(double *)malloc(sizeof(double)*(k+1));
     };
   if (workvec8 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };


   if (n > k)
     {
       diagO=(double *)malloc(sizeof(double)*(n+1));
     }
   else
     {
       diagO=(double *)malloc(sizeof(double)*(k+1));
     };
   if (diagO == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };



   rhs=malloc(sizeof(double)*(k+1));
   if (rhs == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   dy=malloc(sizeof(double)*(k+1));
   if (dy == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   dy1=malloc(sizeof(double)*(k+1));
   if (dy1 == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   Fp=malloc(sizeof(double)*(k+1));
   if (Fp == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   /*
    * Check to make sure that there is at least one constraint.
    */

   if (k<= 0)
     {
       printf("Problem must have at least one constraint.\n");
       exit(206);
     };
   
   /*
    *  Work out the leading dimension for the array.  Note that we may not
    *  want to use k itself, for cache issues.
    */
   if ((k % 2) == 0)
     ldam=k+1;
   else
     ldam=k;

   O=malloc(sizeof(double)*ldam*ldam);
   if (O == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   alloc_mat(C,&Zi);
   alloc_mat(C,&dZ);
   alloc_mat(C,&dX);

   /*
    *  Fill in lots of details in the constraints data structure that haven't
    *  necessarily been done before now.
    */

   /*
    * Next, setup issparse and NULL out all nextbyblock pointers.
    */

   for (i=1; i<=k; i++)
     {
       p=constraints[i].blocks;
       while (p != NULL)
	 {
	   /*
	    * First, set issparse.  The 0.125 parameter was hand tuned.
	    */

	   if (((1.0*k*k*p->numentries*p->numentries) > 0.125*(1.0*k*p->blocksize*p->blocksize*p->blocksize)) && ((p->numentries) > 5))
	     {
	       p->issparse=0;
	     }
	   else
	     {
	       p->issparse=1;
	     };

	   if (C.blocks[p->blocknum].blockcategory == DIAG)
	     p->issparse=1;

	   /*
	    * Setup the cross links.
	    */
	   
	   p->nextbyblock=NULL;
	   p=p->next;
	 };
     };
   
   /*
    * Now, cross link.
    */
   for (i=1; i<=k; i++)
     {
       p=constraints[i].blocks;
       while (p != NULL)
	 {
	   if (p->nextbyblock == NULL)
	     {
	       /*
		* link in the remaining blocks.
		*/
	       for (j=i+1; j<=k; j++)
		 {
		   q=constraints[j].blocks;
		   
		   while (q != NULL)
		     {
		       if (q->blocknum == p->blocknum)
			 {
			   if (p->nextbyblock == NULL)
			     {
			       p->nextbyblock=q;
			       q->nextbyblock=NULL;
			       prev=q;
			     }
			   else
			     {
			       prev->nextbyblock=q;
			       q->nextbyblock=NULL;
			       prev=q;
			     };
			   break;
			 };
		       q=q->next;
		     };
		 };
	     };
	   p=p->next;
	 };
     };

   /*
    * If necessary, print out information on sparsity of blocks.
    */
   
   if (printlevel >= 1)
     {
       numblocks=0;
       denseblocks=0;
       for (i=1; i<=k; i++)
	 {
	   p=constraints[i].blocks;
	   while (p != NULL)
	     {
               if (printlevel >= 4)
                 printf("Sparsity of constraint blocks: %d,%d,%d,%d \n",i,p->blocknum,p->issparse,p->numentries);
               
               if (p->issparse == 0)
                 denseblocks=denseblocks+1;

               numblocks=numblocks+1;
	       p=p->next;
	     };
	 };

       if (printlevel >= 2)
         printf("Percentage of dense constraint blocks is %f\n",
                (100.0*denseblocks)/(1.0*numblocks));
     };
   
   /*
    * Allocate space for byblocks pointers.
    */

   byblocks=(struct sparseblock **)malloc((C.nblocks+1)*sizeof(struct sparseblock *));
   if (byblocks == NULL)
     {
       printf("Storage Allocation Failed!\n");
       exit(205);
     };

   for (i=1; i<=C.nblocks; i++)
     byblocks[i]=NULL;

   /*
    * Fill in byblocks pointers.
    */
   
   for (i=1; i<=k; i++)
     {
       ptr=constraints[i].blocks;
       while (ptr != NULL)
	 {
	   if (byblocks[ptr->blocknum]==NULL)
	     {
	       byblocks[ptr->blocknum]=ptr;
	     };
	   ptr=ptr->next;
	 };
     };

   /*
    *  Compute "fill".  This data structure tells us which elements in the
    *  block diagonal matrix have corresponding elements in one of the
    *  constraints, and which constraint this element first appears in.
    *
    */

   makefill(k,C,constraints,&fill,work1,printlevel);

   /*
    * Compute the nonzero structure of O. This is only used for debugging
    * at this point.
    */

   if (printlevel >= 2)
     {
       nnz=structnnz(n,k,C,constraints);
       printf("Structural density of O %d, %e \n",nnz,nnz*1.0/(k*k*1.0));
     };


   /*
    * Sort entries in the constraints.
    */

   sort_entries(k,C,constraints);

   /*
    * Check the symmetry of C.
    */

   checkc(n,C,printlevel);

   /*
    * Check constraints.
    */

   checkconstraints(n,k,C,constraints,printlevel);

   /*
    *  Now, call sdp().
    */

   ret=sdp(n,k,C,a,constant_offset,constraints,byblocks,fill,*pX,*py,*pZ,
	   cholxinv,cholzinv,ppobj,pdobj,work1,work2,work3,workvec1,
	   workvec2,workvec3,workvec4,workvec5,workvec6,workvec7,workvec8,
	   diagO,bestx,besty,bestz,Zi,O,rhs,dZ,dX,dy,dy1,Fp,
	   printlevel,params);

   if (printlevel >= 1)
     {
       if (ret==0)
	 printf("Success: SDP solved\n");
       if (ret==1)
	 printf("Success: SDP is primal infeasible\n");
       if (ret==2)
	 printf("Success: SDP is dual infeasible\n");
       if (ret==3)
	 printf("Partial Success: SDP solved with reduced accuracy\n");
       if (ret >= 4)
	 printf("Failure: return code is %d \n",ret);

       if (ret==1)
	 {
	   op_at(k,*py,constraints,work1);
	   addscaledmat(work1,-1.0,*pZ,work1);
	   printf("Certificate of primal infeasibility: a'*y=%.5e, ||A'(y)-Z||=%.5e\n",-1.0,Fnorm(work1));
	 };

       if (ret==2)
	 {
	   op_a(k,constraints,*pX,workvec1);
	   printf("Certificate of dual infeasibility: tr(CX)=%.5e, ||A(X)||=%.5e\n",trace_prod(C,*pX),norm2(k,workvec1+1));
	 };

       if ((ret==0) || (ret>=3))
	 {
	   if (printlevel >= 3)
	     {
	       printf("XZ Gap: %.7e \n",trace_prod(*pZ,*pX));
	       gap=*pdobj-*ppobj;
	       printf("Real Gap: %.7e \n",gap);
	     };

	   if (printlevel >= 1)
	     {
               gap=*pdobj-*ppobj;

	       printf("Primal objective value: %.7e \n",*ppobj);
	       printf("Dual objective value: %.7e \n",*pdobj);
	       printf("Relative primal infeasibility: %.2e \n",
		      pinfeas(k,constraints,*pX,a,workvec1));
	       printf("Relative dual infeasibility: %.2e \n",
		      dinfeas(k,C,constraints,*py,*pZ,work1,Fnorm(C)));
	       printf("Real Relative Gap: %.2e \n",gap/(1+fabs(*pdobj)+fabs(*ppobj)));
	       printf("XZ Relative Gap: %.2e \n",trace_prod(*pZ,*pX)/(1+fabs(*pdobj)+fabs(*ppobj)));

	       printf("DIMACS error measures: %.2e %.2e %.2e %.2e %.2e %.2e\n",
		      pinfeas(k,constraints,*pX,a,workvec1)*(1+norm2(k,a+1))/
		      (1+norminf(k,a+1)),
		      0.0,
		      dimacserr3(k,C,constraints,*py,*pZ,work1),
		      0.0,
		      gap/(1+fabs(*pdobj)+fabs(*ppobj)),
		      trace_prod(*pZ,*pX)/(1+fabs(*pdobj)+fabs(*ppobj)));
	     };
	 };

     };

   /*
    *  Now, free up all of the storage.
    */

   free_mat(work1);
   free_mat(work2);
   free_mat(work3);
   free_mat_packed(bestx);
   free_mat_packed(bestz);
   free_mat_packed(cholxinv);
   free_mat_packed(cholzinv);

   free_mat(Zi);
   free_mat(dZ);
   free_mat(dX);

   free(besty);
   free(workvec1);
   free(workvec2);
   free(workvec3);
   free(workvec4);
   free(workvec5);
   free(workvec6);
   free(workvec7);
   free(workvec8);
   free(rhs);
   free(dy);
   free(dy1);
   free(Fp);
   free(O);
   free(diagO);
   free(byblocks);
 
   /*
    * Free up the fill data structure.
    */

   ptr=fill.blocks;
   while (ptr != NULL)
     {
       free(ptr->entries);
       free(ptr->iindices);
       free(ptr->jindices);
       oldptr=ptr;
       ptr=ptr->next;
       free(oldptr);
     };


  /*
   * Finally, free the constraints array.
   */

   return(ret);

}

int structnnz(n,k,C,constraints)
     int n;
     int k;
     struct blockmatrix C;
     struct constraintmatrix *constraints;
{
  int i,j;
  int ii,jj;
  int nnz;
  struct sparseblock *ptri;
  struct sparseblock *ptrj;
  

  nnz=0;
  for (i=1; i<=k; i++)
    for (j=1; j<=k; j++)
      {
	ptri=constraints[i].blocks;

	while (ptri != NULL)
	  {

	    ptrj=constraints[j].blocks;
	    while (ptrj != NULL)
	      {

		if (ptri->blocknum == ptrj->blocknum)
		  {
		    if (C.blocks[ptri->blocknum].blockcategory==MATRIX)
		      {
			nnz++;
			goto NEXTJ;
		      }
		    else
		      { /* DIAG block */
			for (ii=1; ii<=ptri->numentries; ii++)
			  for (jj=1; jj<=ptrj->numentries; jj++)
			    {
			      if (ptri->iindices[ii]==ptrj->iindices[jj])
				{
				  nnz++;
				  goto NEXTJ;
				};
			    };
		      };
		  };
		ptrj=ptrj->next;
	      };

	    ptri=ptri->next;
	  }; /* end while */

      NEXTJ:;
      }; /* end nested fors */
 
  return(nnz);
}

int actnnz(n,lda,A)
     int n;
     int lda;
     double A[];
{
  int i,j;
  int nnz;
 
  nnz=0;
  for (i=1; i<=n; i++)
    {
      if (A[ijtok(i,i,lda)] != 0.0)
	nnz++;
      for (j=i+1; j<=n; j++)
	{
	  if (A[ijtok(i,j,lda)] != 0.0)
	    {
	      nnz++;
	      nnz++;
	    };
	};
    };

  return(nnz);
}

int bandwidth(n,lda,A)
     int n;
     int lda;
     double A[];
{
  int i;
  int j;
  int bw;

  bw=0;
  for (j=2; j<=n; j++)
    {
      for (i=1; i<=j-1; i++)
	{
	  if (A[ijtok(i,j,lda)] != 0.0)
	    {
	      if ((j-i) > bw)
		bw=j-i;
	      break;
	    };
	};
    };

  return(bw);
}

/*
 * Sanity checks for the C matrix data structure.  If this fails, we'll just
 * exit(206) to indicate that there was improper input.  Otherwise, we'll 
 * return 0;
 */

int checkc(int n,struct blockmatrix C,int printlevel)
{
  int i,j,k;
  int totalsize;

  totalsize=0;
  for (k=1; k<=C.nblocks; k++)
    {
      if (C.blocks[k].blockcategory==DIAG)
	{
	  if (printlevel > 5)
	    printf("blockcategory=diag\n");
	};
      if (C.blocks[k].blockcategory==MATRIX)
	{
	  if (printlevel > 5)
	    printf("blockcategory=matrix\n");
	};
      
      totalsize=totalsize+C.blocks[k].blocksize;

      if (C.blocks[k].blockcategory==MATRIX)
	{
	  for (i=1; i<=C.blocks[k].blocksize; i++)
	    {
	      for (j=1; j<=C.blocks[k].blocksize; j++)
		{
		  if (C.blocks[k].data.mat[ijtok(i,j,C.blocks[k].blocksize)] !=
		      C.blocks[k].data.mat[ijtok(j,i,C.blocks[k].blocksize)])
		    {
                      if (printlevel >= 1)
                        printf("C is not symmetric, %d, %d, %d\n",k,i,j);
		      exit(206);
		    }
		};
	    };
	};
    };

  if (totalsize != n)
    {
      if (printlevel >= 1)
        printf("Sum of block sizes does not equal n!\n");
      exit(206);
    };

  return(0);
}

/*
 * Sanity tests on the constraints data structure.
 */

int checkconstraints(n,k,C,constraints,printlevel)
     int n;
     int k;
     struct blockmatrix C;
     struct constraintmatrix *constraints;
     int printlevel;
{
  int i,j;
  struct sparseblock *p;

  for (i=1; i<=k; i++)
    {
      p=constraints[i].blocks;
      if (p==NULL)
	{
          if (printlevel >= 1)
            printf("Constraint %d is empty!\n",i);
	  exit(206);
	};
      while (p != NULL)
	{
	  if (p->constraintnum != i)
	    {
              if (printlevel >= 1)
                printf("p->constraintnum != i, i=%d \n",i);
	      exit(206);
	    };
	  if (p->blocksize != C.blocks[p->blocknum].blocksize)
	    {
              if (printlevel >= 1)
                printf("p->blocksize is wrong, constraint %d \n",i);
	      exit(206);
	    };
	  if (printlevel > 5)
	    printf("Constraint %d, block %d, entries %d\n",i,p->blocknum,p->numentries);
	  for (j=1; j<=p->numentries; j++)
	    {

              /*
               * For debugging, print out the constraint entry.
               */
              
	      if (printlevel >6)
		printf(" (%d, %d)=%lf\n",p->iindices[j],p->jindices[j],p->entries[j]);

              /*
               * Make sure that all indices are in range.
               */
              
	      if (p->iindices[j] > C.blocks[p->blocknum].blocksize)
		{
                  if (printlevel >= 1)
                    printf("i index is larger than blocksize!\n");
		  exit(206);
		};

	      if (p->jindices[j] > C.blocks[p->blocknum].blocksize)
		{
                  if (printlevel >= 1)
                    printf("j index is larger than blocksize!\n");
		  exit(206);
		};

	      if (p->iindices[j] < 1)
		{
                  if (printlevel >= 1)
                    printf("i index is less than 1!\n");
		  exit(206);
		};
              
	      if (p->jindices[j] < 1)
		{
                  if (printlevel >= 1)
                    printf("j index is less than 1!\n");
		  exit(206);
		};

              /*
               * Make sure that entries are sorted properly.
               */

              if (p->iindices[j] > p->jindices[j])
                {
                  if (printlevel >= 1)
                    {
                      printf("i index is greater than j index!\n");
                      printf("constraint=%d\n",i);
                      printf("iindex=%d\n",p->iindices[j]);
                      printf("jindex=%d\n",p->jindices[j]);
                    };

		  exit(206);
                };


              /*
               * Check for any duplicate entries that snuck in.
               */
              
              if (j < p->numentries)
                if ((p->iindices[j]==p->iindices[j+1]) &
                    (p->jindices[j]==p->jindices[j+1]))
                  {
                    if (printlevel >= 1)
                      {
                        printf("Duplicate entry!\n");
                        printf("constraint=%d\n",i);
                        printf("iindex=%d\n",p->iindices[j]);
                        printf("jindex=%d\n",p->jindices[j]);
                      };

                    exit(206);
                  };
              
	      if (C.blocks[p->blocknum].blockcategory==DIAG)
		{
		  if (p->iindices[j] != p->jindices[j])
		    {
                      if (printlevel >= 1)
                        printf("Off diagonal entry in diagonal block!\n");
		      exit(206);
		    };
		};
	    };
	  p=p->next;
	};
    };
  return(0);
} 
