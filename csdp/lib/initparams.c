/*
 *  Setup default values of the parameters.
 */

#include <stdio.h>
#include <strings.h>
#include "declarations.h"

void initparams(params,pprintlevel)
     struct paramstruc *params;
     int *pprintlevel;
{
  FILE *paramfile;
  int ret;
  double value;
  char parametername[30];
  char junk[2];

  /*
   * First, set default values for all parameters.
   */

  params->axtol=1.0e-8;
  params->atytol=1.0e-8;
  params->objtol=1.0e-8;
  params->pinftol=1.0e8;
  params->dinftol=1.0e8;
  params->maxiter=100;
  params->minstepfrac=0.90;
  params->maxstepfrac=0.97;
  params->minstepp=1.0e-8;
  params->minstepd=1.0e-8;
  params->usexzgap=1;
  params->tweakgap=0;
  params->affine=0;
  params->perturbobj=1;
  params->fastmode=0;
  *pprintlevel=0;

  /*
   * Attempt to open param.csdp.  If it doesn't open, then just use 
   * the default values.
   */

  paramfile=fopen("param.csdp","r");

  if (paramfile != NULL)
    {
      while (1==1)
        {
          /*
           * Skip over leading white space if there is any.  Ignore the return
           * code, since it's OK if there are no leading spaces.
           */

          ret=fscanf(paramfile,"%*[ \t\n]");

          /*
           * handle end of file.
           */

          if (ret == EOF)
            break;
          
          /*
           * Get the parameter name. 
           */
          
          ret=fscanf(paramfile,"%29[A-Za-z0-9]",parametername);

          if (ret != 1)
            {
              /*
               * We don't have a parameter name.  Go on to the next line.
               */
              ret=fscanf(paramfile,"%*[^\n]\n");
              continue;
            };
          
          /*
           * Skip over any trailing white space, the = sign, and any 
           * leading white space after =.  
           */

          ret=fscanf(paramfile,"%*[ \t]");
          ret=fscanf(paramfile,"%1[=]",junk);
          if (ret != 1)
            {
              printf("param.csdp line missing =.  Skipping to next line.\n");
              ret=fscanf(paramfile,"%*[^\n]\n");
              continue;
            };

          /*
           * Skip over any white space after equal sign.
           */
          
          ret=fscanf(paramfile,"%*[ \t]");
          
          /*
           * Get the value.
           */
          
          ret=fscanf(paramfile,"%le",&value);
          
          /*
           * Skip to the end of line.
           */
          ret=fscanf(paramfile,"%*[^\n]\n");
          
          /*
           * For debugging, print out the parametername and value.
           */

          /*

          printf("parameter name: %s\n",parametername);
          printf("value: %le\n",value);

          */
          
          /*
           * Now, adjust the parameter as needed.
           */
          
          if (strcasecmp(parametername,"axtol")==0)
            {
              params->axtol=value;
              continue;
            };

          if (strcasecmp(parametername,"atytol")==0)
            {
              params->atytol=value;
              continue;
            };

          if (strcasecmp(parametername,"objtol")==0)
            {
              params->objtol=value;
              continue;
            };

          if (strcasecmp(parametername,"pinftol")==0)
            {
              params->pinftol=value;
              continue;
            };

          if (strcasecmp(parametername,"dinftol")==0)
            {
              params->dinftol=value;
              continue;
            };

          if (strcasecmp(parametername,"maxiter")==0)
            {
              params->maxiter=value;
              continue;
            };

          if (strcasecmp(parametername,"minstepfrac")==0)
            {
              params->minstepfrac=value;
              continue;
            };

          if (strcasecmp(parametername,"maxstepfrac")==0)
            {
              params->maxstepfrac=value;
              continue;
            };

          if (strcasecmp(parametername,"minstepp")==0)
            {
              params->minstepp=value;
              continue;
            };

          if (strcasecmp(parametername,"minstepd")==0)
            {
              params->minstepd=value;
              continue;
            };

          if (strcasecmp(parametername,"usexzgap")==0)
            {
              params->usexzgap=value;
              continue;
            };

          if (strcasecmp(parametername,"tweakgap")==0)
            {
              params->tweakgap=value;
              continue;
            };

          if (strcasecmp(parametername,"affine")==0)
            {
              params->affine=value;
              continue;
            };

          if (strcasecmp(parametername,"printlevel")==0)
            {
              *pprintlevel=value;
              continue;
            };

          if (strcasecmp(parametername,"perturbobj")==0)
            {
              params->perturbobj=value;
              continue;
            };

          if (strcasecmp(parametername,"fastmode")==0)
            {
              params->fastmode=value;
              continue;
            };

          printf("param.csdp: unrecognized parameter, %s\n",
                 parametername);
        };

      /*
       * Close the parameter file.   
       */

      fclose(paramfile);

    };

  /*
   * Print out the parameters that will be used.  
   */

  if (*pprintlevel >= 2)
    {
      printf("params->axtol is %e \n",params->axtol);
      printf("params->atytol is %e \n",params->atytol);
      printf("params->objtol is %e \n",params->objtol);
      printf("params->pinftol is %e \n",params->pinftol);
      printf("params->dinftol is %e \n",params->dinftol);
      printf("params->maxiter is %d \n",params->maxiter);
      printf("params->minstepfrac is %e \n",params->minstepfrac);
      printf("params->maxstepfrac is %e \n",params->maxstepfrac);
      printf("params->minstepp is %e \n",params->minstepp);
      printf("params->minstepd is %e \n",params->minstepd);
      printf("params->usexzgap is %d \n",params->usexzgap);
      printf("params->tweakgap is %d \n",params->tweakgap);
      printf("params->affine is %d \n",params->affine);
      printf("params->printlevel is %d \n",*pprintlevel);
      printf("params->perturbobj is %e \n",params->perturbobj);
      printf("params->fastmode is %d \n",params->fastmode);
    };


}
