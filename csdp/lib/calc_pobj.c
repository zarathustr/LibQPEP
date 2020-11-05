/*
  Compute the primal objective function value pobj=Trace(C*X)  Since we
  only need the trace, it makes more sense to just compute the diagonal
  of the product and sum the entries..
  */

#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

double calc_pobj(C,X,constant_offset)
     struct blockmatrix C;
     struct blockmatrix X;
     double constant_offset;
{
  double pobj;

  pobj=constant_offset+trace_prod(C,X);

  return(pobj);
	 
}
