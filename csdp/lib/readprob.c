/*
  Read in out a problem in SDPA sparse format.  Return 0 if ok, 1 if 
  failure.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "declarations.h"

void skip_to_end_of_line();
int get_line();
void countentry();
int addentry();

int read_prob(fname,pn,pk,pC,pa,pconstraints,printlevel)
     char *fname;
     int *pn;
     int *pk;
     struct blockmatrix *pC;
     double **pa;
     struct constraintmatrix **pconstraints;
     int printlevel;
     
{
  struct constraintmatrix *myconstraints;
  FILE *fid;
  int i,j;
  int buflen;
  char *buf;
  int c;
  int nblocks;
  int blksz;
  int blk;
  char *ptr1;
  char *ptr2;
  int matno;
  int blkno;
  int indexi;
  int indexj;
  double ent;
  int ret;
  struct sparseblock *p;
  int *isdiag;
  double *tempdiag;

  /*
   * Open the file for reading, and determine the length of the longest
   * line.
   */

  fid=fopen(fname,"r");
 
  if (fid == (FILE *) NULL)
    {
      if (printlevel >= 1)
        printf("Couldn't open problem file for reading! \n");
      exit(201);
    };

  /*
   * This constant allows for up to 500000*20=10,000,000 byte long lines.
   * You might see a line this long in a problem with 500000 constraints.
   */
  
  buflen=8000000;

  buf=(char *)malloc(buflen*sizeof(char));
  if (buf == NULL)
    {
      if (printlevel >= 1)
        printf("Storage allocation failed!\n");
      exit(205);
    };

  /*
   * In the first pass through the file, determine the size parameters,
   * and allocate space for everything.
   */

  fid=fopen(fname,"r");
 
  if (fid == (FILE *) NULL)
    {
      if (printlevel >= 1)
        printf("Couldn't open problem file for reading! \n");
      exit(201);
    };

  /*
   * First, read through the comment lines.
   */
 
  c=getc(fid);
  while ((c == '"') || (c == '*'))
    {
      skip_to_end_of_line(fid);
      c=getc(fid);
    };

  ungetc(c,fid);

  /*
   * Get the number of constraints (primal variables in SDPA terminology)
   */

  ret=get_line(fid,buf,buflen);
  if (ret == 0)
    {
      ret=sscanf(buf,"%d",pk);
      if ((ret!=1) || (*pk<=0))
	{
          if (printlevel >= 1)
            printf("Incorrect SDPA file.  Couldn't read mDIM\n");
	  fclose(fid);
	  return(1);
	};
    }
  else
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Couldn't read mDIM \n");
      fclose(fid);
      return(1);
    };

#ifndef BIT64
  /*
   * If operating in 32 bit mode, make sure that the dimension mDIM isn't
   * too big for 32 bits.  If we don't do this check, then integer overflow
   * won't be detected, and we'll allocate a bogus amount of storage for
   * O.
   */

  if (*pk > 23169)
    {
      if (printlevel >= 1)
        printf("This problem is too large to be solved in 32 bit mode!\n");
      exit(206);
    };
#endif

  /*
   * Read in the number of blocks.
   */
  ret=get_line(fid,buf,buflen);
  if (ret == 0)
    {
      ret=sscanf(buf,"%d",&nblocks);
      if ((ret != 1) || (nblocks <= 0))
	{
          if (printlevel >= 1)
            printf("Incorect SDPA file. Couldn't read nBLOCKS. \n");
	  fclose(fid);
	  return(1);
	};
    }
  else
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Couldn't read nBLOCKS. \n");
      fclose(fid);
      return(1);
    };

  /*
   * Keep track of which blocks have off diagonal entries. 
   */

  isdiag=(int *)malloc((nblocks+1)*sizeof(int));
  for (i=1; i<=nblocks; i++)
    isdiag[i]=1;

  /*
   * Allocate space for the C matrix.
   */

  pC->nblocks=nblocks;
  pC->blocks=(struct blockrec *)malloc((nblocks+1)*sizeof(struct blockrec));
  if (pC->blocks == NULL)
    {
      if (printlevel >= 1)
        printf("Storage allocation failed!\n");
      exit(205);
    }

  /*
   * Allocate space for the constraints.
   */

  myconstraints=(struct constraintmatrix *)malloc((*pk+1)*sizeof(struct constraintmatrix));

  if (myconstraints == NULL)
    {
      if (printlevel >= 1)
        printf("Storage allocation failed!\n");
      exit(205);
    };
  
  /*
   * Null out all pointers in constraints.
   */
  for (i=1; i<=*pk; i++)
    {
      myconstraints[i].blocks=NULL;
    };

  *pa=(double *)malloc((*pk+1)*sizeof(double));

  if (*pa == NULL)
    {
      if (printlevel >= 1)
        printf("Storage allocation failed!\n");
      exit(205);
    };

  /*
   * And read the block structure.
   */

  *pn=0;

  ret=get_line(fid,buf,buflen);
  if (ret == 0)
    {
      /*
       * Decode nblocks numbers out of the buffer.  Put the results in 
       * block_structure.
       */
      ptr1=buf;
      for (blk=1; blk<=nblocks; blk++)
	{
	  blksz=strtol(ptr1,&ptr2,10);
	  ptr1=ptr2;

	  /*
	   * negative numbers are used to indicate diagonal blocks.  First,
	   * update n.
	   */

	  *pn=*pn+abs(blksz);

	  /*
	   * Now, handle diagonal blocks and matrix blocks separately.
	   */
	  if (blksz < 0)
	    {
	      /*
	       * It's a diag block.
	       */
	      pC->blocks[blk].blocksize=abs(blksz);
	      pC->blocks[blk].blockcategory=DIAG;
	      pC->blocks[blk].data.vec=(double *)malloc((1+abs(blksz))*sizeof(double));
	      if (pC->blocks[blk].data.vec == NULL)
		{
                  if (printlevel >= 1)
                    printf("Storage allocation failed!\n");
		  exit(205);
		};
	      for (i=1; i<=abs(blksz); i++)
		pC->blocks[blk].data.vec[i]=0.0;

	    }
	  else
	    {
	      /*
	       * It's a matrix block.  In the I32LP64 model, we can't do array 
               * indexing for huge blocks using 32 bit integers.  46340 is 
               * the largest sized block that we can handle.  Trying to 
               * malloc a larger block could cause integer overflow and
               * result in a negative malloc size.
b	       */
              if (blksz > 46340)
                {
                  printf("This problem is too large to be solved in I32LP64 mode.\n");
                  exit(206);
                  
                };
              /*
               * Setup the block in C.
               */
	      pC->blocks[blk].blocksize=abs(blksz);
	      pC->blocks[blk].blockcategory=MATRIX;
	      pC->blocks[blk].data.mat=(double *)malloc((blksz*blksz)*sizeof(double));
	      if (pC->blocks[blk].data.mat == NULL)
		{
                  if (printlevel >= 1)
                    printf("Storage allocation failed!\n");
		  exit(205);
		};

	      for (j=1; j<=blksz; j++)
		for (i=1; i<=blksz; i++)
		  pC->blocks[blk].data.mat[ijtok(i,j,blksz)]=0.0;

	    };
		
	};

    }
  else
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Couldn't read block sizes.\n");
      fclose(fid);
      free(isdiag);
      return(1);
    };

  /*
   *  Read in the right hand side values.
   */

  ret=get_line(fid,buf,buflen);
  if (ret == 0)
    {
      /*
       * Decode k numbers out of the buffer.  Put the results in 
       * a.
       */
      ptr1=buf;
      for (i=1; i<=*pk; i++)
	{
	  (*pa)[i]=strtod(ptr1,&ptr2);
	  /*
	   * Check for a case where ptr2 didn't advance.  This indicates
	   * a strtod failure.
	   */
	  if (ptr1==ptr2)
	    {
              if (printlevel >= 1)
                printf("Incorect SDPA file. Can't read RHS values.\n");
	      fclose(fid);
	      free(isdiag);
	      return(1);
	    };
	  ptr1=ptr2;
	};
    }
  else
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Can't read values.\n");
      fclose(fid);
      free(isdiag);
      return(1);
    };

  /*
   *  Now, loop through the entries, 
   *  counting entries in the constraint matrices block by block.
   */

  ret=fscanf(fid,"%d %d %d %d %le ",&matno,&blkno,&indexi,&indexj,&ent);

  if (ret != 5)
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Return code from fscanf is %d, should be 5\n",ret);
      fclose(fid);
      free(isdiag);
      return(1);
    };

  do {

    /*
     * Check the validity of these values.
     */

    if ((matno < 0) || (matno > *pk) ||
	(blkno<1) || (blkno>nblocks) ||
	(indexi < 1) || (indexi > pC->blocks[blkno].blocksize) ||
	(indexj < 1) || (indexj > pC->blocks[blkno].blocksize))
      {
        if (printlevel >= 1)
          printf("Incorect SDPA file. Bad values in line: %d %d %d %d %e \n",
                 matno,blkno,indexi,indexj,ent);
	fclose(fid);
	free(isdiag);
	return(1);
      };

    if (matno != 0)
      {
	if (ent != 0.0)
	  countentry(myconstraints,matno,blkno,pC->blocks[blkno].blocksize);
      }
    else
      {
	/*
	 * An entry in C. ignore it for now.
	 */
      };
    ret=fscanf(fid,"%d %d %d %d %le",&matno,&blkno,&indexi,&indexj,&ent);
  } while (ret == 5);

  if ((ret != EOF) && (ret != 0))
    {
      if (printlevel >= 1)
        printf("Incorrect SDPA file, while reading entries.  ret=%d \n",ret);
      fclose(fid);
      free(isdiag);
      return(1);
    };

  fclose(fid);

  /*
   * Now, go through each of the blks in each of the constraint matrices,
   * and allocate space for the entries and indices.
   */
  for (i=1; i<=*pk; i++)
    {
      p=myconstraints[i].blocks;

      while (p != NULL)
	{
	  /*
	   * allocate storage for the entries in this block of this constraint.
	   */
	  p->entries=(double *)malloc((p->numentries+1)*sizeof(double));
          if (p->entries == NULL)
	    {
              if (printlevel >= 1)
                printf("Storage allocation failed!\n");
	      exit(205);
	    };

	  p->iindices=(int *)malloc((p->numentries+1)*sizeof(int));

          if (p->iindices == NULL)
	    {
              if (printlevel >= 1)
                printf("Storage allocation failed!\n");
	      exit(205);
	    };

	  p->jindices=(int *)malloc((p->numentries+1)*sizeof(int));

          if (p->jindices == NULL)
	    {
              if (printlevel >= 1)
                printf("Storage allocation failed!\n");
	      exit(205);
	    };

	  p->numentries=0;
	  p=p->next;
	};
    };


  /*
   *  In the final pass through the file, fill in the actual data.
   */

  zero_mat(*pC);
  
  /*
   * Open the file for reading, and then read in all of the actual 
   * matrix entries.
   * line.
   */
  fid=fopen(fname,"r");
 
  if (fid == (FILE *) NULL)
    {
      if (printlevel >= 1)
        printf("Couldn't open problem file for reading! \n");
      exit(201);
    };

  /*
   * First, read through the comment lines.
   */
 
  c=getc(fid);
  while ((c == '"') || (c == '*'))
    {
      skip_to_end_of_line(fid);
      c=getc(fid);
    };

  ungetc(c,fid);

  /*
   * Get the number of constraints (primal variables in SDPA terminology)
   */

  ret=get_line(fid,buf,buflen);
  if (ret == 0)
    {
      sscanf(buf,"%d",pk);
    }
  else
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Couldn't read mDIM \n");
      fclose(fid);
      free(isdiag);
      return(1);
    };

  /*
   * Read in the number of blocks.
   */
  ret=get_line(fid,buf,buflen);
  if (ret == 0)
    {
      sscanf(buf,"%d",&nblocks);
    }
  else
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Couldn't read nBLOCKS. \n");
      fclose(fid);
      free(isdiag);
      return(1);
    };

  /*
   * And read the block structure.
   */

  ret=get_line(fid,buf,buflen);
  if (ret != 0)
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Couldn't read block sizes.\n");
      fclose(fid);
      free(isdiag);
      return(1);
    };

  /*
   *  Read in the right hand side values.
   */

  ret=get_line(fid,buf,buflen);
  if (ret == 0)
    {
      /*
       * Decode k numbers out of the buffer.  Put the results in 
       * a.
       */
      ptr1=buf;
      for (i=1; i<=*pk; i++)
	{
	  (*pa)[i]=strtod(ptr1,&ptr2);
	  /*
	   * Check for a case where ptr2 didn't advance.  This indicates
	   * a strtod failure.
	   */
	  if (ptr1==ptr2)
	    {
              if (printlevel >= 1)
                printf("Incorect SDPA file. Can't read RHS values.\n");
	      fclose(fid);
	      free(isdiag);
	      return(1);
	    };
	  ptr1=ptr2;
	};
    }
  else
    {
      if (printlevel >= 1)
        printf("Incorect SDPA file. Can't read a values.\n");
      fclose(fid);
      free(isdiag);
      return(1);
    };

  /*
   * Now, read the actual entries.
   */
  ret=fscanf(fid,"%d %d %d %d %le ",&matno,&blkno,&indexi,&indexj,&ent);
  do {

    /*
     * No need for sanity checking the second time around.
     */

    /*
     * Mark this block as not diagonal if indexi!=indexj.
     */
    if ((indexi != indexj)  && (ent != 0.0))
      isdiag[blkno]=0;

    if (matno != 0)
      {
	if (ent != 0.0)
          {
	  ret=addentry(myconstraints,matno,blkno,indexi,indexj,ent);
        
          if (ret != 0)
            {
              if (printlevel >= 1)
                {
                  printf("Incorrect SDPA file. Duplicate entry.\n");
                  printf("matno=%d\n",matno);
                  printf("blkno=%d\n",blkno);
                  printf("indexi=%d\n",indexi);
                  printf("indexj=%d\n",indexj);
                };
              
              fclose(fid);
              free(isdiag);
              return(1);
            };
          };
      }
    else
      {
	/*
	 * An entry in C. 
	 */
	if (ent != 0.0)
	  {
	    blksz=pC->blocks[blkno].blocksize;
	    if (pC->blocks[blkno].blockcategory == DIAG)
	      {
                if (pC->blocks[blkno].data.vec[indexi] != 0.0)
                  {
                    /*
                     * We've got a duplicate entry in C!
                     */
                    if (printlevel >= 1)
                      {
                        printf("Incorrect SDPA file. Duplicate entry.\n");
                        printf("matno=%d\n",matno);
                        printf("blkno=%d\n",blkno);
                        printf("indexi=%d\n",indexi);
                        printf("indexj=%d\n",indexj);
                      };

                    fclose(fid);
                    free(isdiag);
                    return(1);
                  }
                else
                  pC->blocks[blkno].data.vec[indexi]=ent;
	      }
	    else
	      {
                if (pC->blocks[blkno].data.mat[ijtok(indexi,indexj,blksz)] != 0.0)
                  {
                    /*
                     * We've got a duplicate entry in C!
                     */
                    if (printlevel >= 1)
                      {
                        printf("Incorrect SDPA file. Duplicate entry.\n");
                        printf("matno=%d\n",matno);
                        printf("blkno=%d\n",blkno);
                        printf("indexi=%d\n",indexi);
                        printf("indexj=%d\n",indexj);
                      };

                    fclose(fid);
                    free(isdiag);
                    return(1);                    
                  }
                else
                  {
                    pC->blocks[blkno].data.mat[ijtok(indexi,indexj,blksz)]=ent;
                    pC->blocks[blkno].data.mat[ijtok(indexj,indexi,blksz)]=ent;
                  };
	      };
	  };
      };
    ret=fscanf(fid,"%d %d %d %d %le ",&matno,&blkno,&indexi,&indexj,&ent);
  } while (ret == 5);

  if ((ret != EOF) && (ret != 0))
    {
      if (printlevel >= 1)
        printf("Incorrect SDPA file. \n");
      fclose(fid);
      free(isdiag);
      return(1);
    };

  /*
   * At this point, we'll stop to recognize whether any of the blocks
   * are "hidden LP blocks"  and correct the block type if needed.
   */

  for (i=1; i<=nblocks; i++)
    {
      if ((pC->blocks[i].blockcategory != DIAG) && 
	  (isdiag[i]==1) && (pC->blocks[i].blocksize > 1))
	{
	  /*
	   * We have a hidden diagonal block!
	   */
	  if (printlevel >= 2)
	    {
	      printf("Block %d is actually diagonal.\n",i);
	    };
	  blksz=pC->blocks[i].blocksize;
	  tempdiag=(double *)malloc((blksz+1)*sizeof(double));
	  for (j=1; j<=blksz; j++)
	    tempdiag[j]=pC->blocks[i].data.mat[ijtok(j,j,blksz)];
	  free(pC->blocks[i].data.mat);
	  pC->blocks[i].data.vec=tempdiag;
	  pC->blocks[i].blockcategory=DIAG;
	};
    };

  /*
   * If the printlevel is high, print out info on constraints and block
   * matrix structure.
   */
  if (printlevel >= 3)
    {
      printf("Block matrix structure.\n");
      for (blk=1; blk<=pC->nblocks; blk++)
	{
	  if (pC->blocks[blk].blockcategory == DIAG)
	    printf("Block %d, DIAG, %d \n",blk,pC->blocks[blk].blocksize);
	  if (pC->blocks[blk].blockcategory == MATRIX)
	    printf("Block %d, MATRIX, %d \n",blk,pC->blocks[blk].blocksize);
	};
    };

  /*
   * Next, setup issparse and NULL out all nextbyblock pointers.
   */

  for (i=1; i<=*pk; i++)
    {
      p=myconstraints[i].blocks;
      while (p != NULL)
	{
	  p->nextbyblock=NULL;
	  p=p->next;
	};
    };
  
  /*
   * Free unneeded memory.
   */

  free(buf);
  free(isdiag);

  /*
   *  Put back all the returned values.
   */

  *pconstraints=myconstraints;
  
  fclose(fid);
  return(0);
}

/*
 *  This routine skips to the end of the current line of input from the 
 *  file fid.
 */

void skip_to_end_of_line(fid)
     FILE *fid;
{
  char c;
 
  c=getc(fid);
  while (c != '\n')
    c=getc(fid);
}

/*
 *  This routine reads a line of input into a buffer, and translates all
 *  occurences of "," "{" "}" "(" ")" to spaces.
 *
 */

int get_line(fid,buffer,bufsiz)
     FILE *fid;
     char *buffer;
     int bufsiz;
{
  int k;
  char c;

  if (fgets(buffer,bufsiz-1,fid) != NULL)
    {
      c=buffer[0];
      k=0;
      while ((c != '\n')  && (c!=0) && (k<bufsiz))
        {
          if (c==',')
            buffer[k]=' ';
          if (c=='(')
            buffer[k]=' ';
          if (c==')')
            buffer[k]=' ';
          if (c=='{')
            buffer[k]=' ';
          if (c=='}')
            buffer[k]=' ';
          k=k+1;
          c=buffer[k];
        };

      /*
       * return an error if the line is longer than the buffer!
       */
      
      if (k<bufsiz-5)
        return(0);
      else
        return(1);
    }
  else
    return(2);
}
  

void countentry(constraints,matno,blkno,blocksize)
     struct constraintmatrix *constraints;
     int matno;
     int blkno;
     int blocksize;
{
  struct sparseblock *p;
  struct sparseblock *q;
  
  p=constraints[matno].blocks;

  if (p == NULL)
    {
      /*
       * We haven't yet allocated any blocks.
       */
      p=(struct sparseblock *)malloc(sizeof(struct sparseblock));
      if (p==NULL)
	{
          printf("Storage allocation failed!\n");
	  exit(205);
	};
      p->constraintnum=matno;
      p->blocknum=blkno;
      p->numentries=1;
      p->next=NULL;
      p->entries=NULL;
      p->iindices=NULL;
      p->jindices=NULL;
      p->blocksize=blocksize;
      constraints[matno].blocks=p;
    }
  else
    {
      /*
       * We have some existing blocks.  See whether this block is already
       * in the chain.
       */

      while ((p->next) != NULL)
	{
	  if (p->blocknum == blkno)
	    {
	      /*
	       * Found the right block.
	       */
	      p->numentries=p->numentries+1;
	      return;
	    };
	  p=p->next;
	};
      /*
       * If we get here, we still have to check the last block in the
       * chain.
       */
      if (p->blocknum == blkno)
	{
	  /*
	   * Found the right block.
	   */
	  p->numentries=p->numentries+1;
	  return;
	};
      /*
       * If we get here, then the block doesn't exist yet.
       */
      q=(struct sparseblock *)malloc(sizeof(struct sparseblock));
      if (q==NULL)
	{
          printf("Storage allocation failed!\n");
	  exit(205);
	};
      /*
       * Fill in information for this block.
       */
      q->blocknum=blkno;
      q->constraintnum=matno;
      q->numentries=1;
      q->next=NULL;
      q->entries=NULL;
      p->iindices=NULL;
      p->jindices=NULL;
      q->blocksize=blocksize;
      /*
       * Now link it into the list.
       */
      p->next=q;
    };

}

int addentry(constraints,matno,blkno,indexi,indexj,ent)
     struct constraintmatrix *constraints;
     int matno;
     int blkno;
     int indexi;
     int indexj;
     double ent;
{
  struct sparseblock *p;
  int itemp;
    
  /*
   * Arrange things so that indexi <= indexj.
   */

  if (indexi > indexj)
    {
      itemp=indexi;
      indexi=indexj;
      indexj=itemp;
    };

  /* 
   * Find the appropriate block.
   */
  
  p=constraints[matno].blocks;
  
  if (p == NULL)
    {
      printf("Internal Error in readprob.c !\n");
      exit(206);
    }
  else
    {
      /*
       * We have some existing blocks.  See whether this block is already
       * in the chain.
       */

      while (p != NULL)
	{
	  if (p->blocknum == blkno)
	    {
	      /*
	       * Found the right block. 
	       */

	      p->numentries=(p->numentries)+1;
	      p->entries[(p->numentries)]=ent;
	      p->iindices[(p->numentries)]=indexi;
	      p->jindices[(p->numentries)]=indexj;

              /*
               * We've successfully added the entry.
               */
              
	      return(0);
	    };
	  p=p->next;
	};
      /*
       * If we get here, we have an internal error.
       */
      printf("Internal Error in CSDP readprob.c !\n");
      exit(206);
    };

  /*
   * Everything is good, so return 0.
   */

  return(0);

}


