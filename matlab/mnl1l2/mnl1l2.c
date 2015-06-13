#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lbfgs.h"


#include "arithmetic_ansi.h"

#include <stdio.h>
#include <math.h>
#include "mex.h"

typedef struct {
  int ns, nc, ncls, n, display;
  lbfgsfloatval_t lambda;
  lbfgsfloatval_t *x;
  int *y;
  lbfgsfloatval_t *out;
  lbfgsfloatval_t *g;
  lbfgsfloatval_t *s;
  int *killsign;
  double *pniter;
} mnl1l2problem_t;

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *ww,
    lbfgsfloatval_t *gg,
    const int dd,
    const lbfgsfloatval_t step
    )
{
    int i,j;
    lbfgsfloatval_t fval = 0.0;
    lbfgsfloatval_t outmax, sumexp;

    mnl1l2problem_t* prob = (mnl1l2problem_t*)instance;

    /* Set gradient vector to zero */
    vecset(gg, 0, dd);

    for (i = 0;i < prob->n; i++) { 
      int ixtr = i*prob->ncls;

      /* Calculate the output values */
      for (j = 0;j < prob->ncls; j++) { 
	int ixepo = ixtr+j;
	vecdot(&prob->out[ixepo], &prob->x[ixepo*dd], ww, dd);

      }
      /* Maximum output value */
      vecmax(&outmax, &prob->out[ixtr], prob->ncls);
      

      /* Sum of exponentiated output values */
      sumexp = 0;
      for (j = 0; j < prob->ncls; j++) {
	int ixepo = ixtr+j;
	sumexp += exp(prob->out[ixepo] - outmax);
	
      }

      /* Objective (likelihood term) */
      fval += -prob->out[ixtr + prob->y[i]] + outmax + log(sumexp);

      /* Gradient (likelihood term) */
      for (j = 0; j < prob->ncls; j++) {
	int ixepo = ixtr+j;
	lbfgsfloatval_t coeff = -(j==prob->y[i]) + exp(prob->out[ixepo]-outmax-log(sumexp));
	vecadd(gg, &prob->x[ixepo*dd], coeff, dd);
      }
    }

    /* Regularization term */
    for (i=0; i<prob->nc; i++) {
      lbfgsfloatval_t nm, nmg;

      /* norm of the weight vector */
      vecnorm(&nm, &ww[prob->ns*i], prob->ns);

      /* norm of the gradient vector */
      vecnorm(&nmg, &gg[prob->ns*i], prob->ns);

      /* Objective value */
      fval += prob->lambda*nm;

      prob->killsign[i] = 0;
      if (nm==0) {
	/* At the origin */
	if (nmg <= prob->lambda) {
	  /* Trapped by the regularization */
	  vecset(&gg[prob->ns*i], 0, prob->ns);
	}
	else {
	  /* Not trapped  */
	  vecscale(&gg[prob->ns*i], (1-prob->lambda/nmg), prob->ns);
	}
      }
      else {
	if (nm<1e-6) {
	  prob->killsign[i] = 1;
	}
	/* Regular point */
	vecadd(&gg[prob->ns*i], &ww[prob->ns*i], prob->lambda/nm, prob->ns);

      }
      /* mexPrintf("killsign[%d]=%d\n",i, prob->killsign[i]);*/

    }



    /* Quadratic regularization
    fval += 0.5*prob->lambda*w2;
    vecadd(gg, ww, prob->lambda, n);
    */

    return fval;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t *s,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int dd,
    int k,
    int ls
    )
{
  mnl1l2problem_t *prob = (mnl1l2problem_t*)instance;
  lbfgsfloatval_t snorm;

  if (prob->display>1) {
    vecnorm(&snorm, s, dd);
    mexPrintf("Iteration %d: w=[%f %f] fval=%g gnorm=%g snorm=%g step=%g\n", k, x[0], x[1], fx, gnorm, snorm, step);
    fflush(0);
    /*mexEvalString("drawnow;");*/
  }
  else if (prob->display>0) {
    mexPrintf(".");
    fflush(0); 
    /*mexEvalString("drawnow;");*/
  }
    
    /* Copy the descent direction */
    veccpy(prob->s, s, dd);

    /* Copy the number of iterations */
    *prob->pniter = (double)k;

    return 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i, dd, ret, nTrials;
  double *ytmp, niter;
  lbfgsfloatval_t fval, *ww;
  mnl1l2problem_t problem;
  lbfgs_parameter_t param;
  const int *size;
  const char *code;
  const int noptfields    = 6;
  const char *optfields[] = {"W0","epsg","epsf","maxiter","lsmethod", "display"};
  mxArray *opt, *tmp;

  mexPrintf("mnl1l2.c (rev 0.2)\n");

  if (nrhs <3 || nlhs <1) {
    mexErrMsgTxt("Syntax: [W, z, status] = mnl1l2(X, Y, lambda, <opt>)");
  }

  size           = mxGetDimensions(prhs[0]);
  
  nTrials        = mxGetNumberOfElements(prhs[1]);

  if (size[0]==0 || size[1]==0 || size[2]==0) {
    mexErrMsgTxt("X must have the size of [ns, nc, ncls*n].");
  }

  if (nTrials==0) {
    mexErrMsgTxt("Y must be specified as {0,...,ncls-1}^n.");
  }

  if (((int)size[2]/nTrials)*nTrials != size[2]) {
    mexErrMsgTxt("size(X,3) must be divisible by length(Y).");
  }

  if (mxGetNumberOfElements(prhs[2])==0) {
    mexErrMsgTxt("Lambda must be specified.");
  }
  
  problem.ns     = size[0];
  problem.nc     = size[1];
  problem.ncls   = size[2]/nTrials;
  problem.n      = nTrials;
  problem.x      = mxGetPr(prhs[0]);

  dd             = problem.ns*problem.nc;
  
  ytmp           = mxGetPr(prhs[1]);
  problem.y      = malloc(sizeof(int)*problem.n);
  for (i=0; i< nTrials; i++) {
    problem.y[i] = (int)ytmp[i];
  }

  problem.lambda = *(mxGetPr(prhs[2]));
  problem.killsign = malloc(sizeof(int)*problem.nc);


  /* Load options */
  opt = mxCreateStructMatrix(1, 1, noptfields, optfields);
  if (nrhs > 3) {
    for (i=0; i<noptfields; i++) {
      if ((tmp=mxGetField(prhs[3], 0, optfields[i]))!=NULL) {
	mxSetField(opt, 0, optfields[i], tmp);
      }
    }
  } 
  
  
  /* 1st output arrgument: weight vectors */
  plhs[0]        = mxCreateDoubleMatrix(problem.ns, problem.nc, mxREAL);
  ww             = mxGetPr(plhs[0]);
  
  /* Load the default value (if it exists) */
  tmp = mxGetField(opt, 0, "W0");
  if (tmp!=NULL && mxGetM(tmp)==problem.ns && mxGetN(tmp)==problem.nc) {
    veccpy(ww, mxGetPr(tmp), dd);
  }
  else {
    vecset(ww, 0, dd);
  }


  /* 2nd output argument: output values */
  if (nlhs>1) {
    plhs[1]      = mxCreateDoubleMatrix(1, problem.ncls*problem.n, mxREAL);
    problem.out  = mxGetPr(plhs[1]);
  }
  else {
    problem.out  = lbfgs_malloc(problem.ncls*problem.n);
  }

  /* 3rd output argument: status */
  if (nlhs>2) {
    const char *fnames[] = {"ret","niter","fval", "g", "s"};
    mxArray *stmp        = mxCreateDoubleMatrix(problem.ns, problem.nc, mxREAL);
    mxArray *nitertmp    = mxCreateDoubleMatrix(1, 1, mxREAL);

    plhs[2]              = mxCreateStructMatrix(1, 1, 5, fnames);
    mxSetField(plhs[2], 0, "ret", mxCreateDoubleMatrix(1, 1, mxREAL));
    mxSetField(plhs[2], 0, "niter", nitertmp);
    mxSetField(plhs[2], 0, "fval", mxCreateDoubleMatrix(1, 1, mxREAL));
    mxSetField(plhs[2], 0, "g", mxCreateDoubleMatrix(problem.ns, problem.nc, mxREAL));
    mxSetField(plhs[2], 0, "s", stmp);

    problem.s            = mxGetPr(stmp);
    problem.pniter       = mxGetPr(nitertmp);
  }
  else {
    problem.s      = lbfgs_malloc(problem.ns*problem.nc);
    problem.pniter = &niter;
  }

  if (ww == NULL || problem.y == NULL || problem.out == NULL || problem.s == NULL) {
    mexErrMsgTxt("Failed to allocate a memory block for variables.");
  }

  mexPrintf("ns=%d nc=%d ncls=%d n=%d lambda=%f\n", \
	 problem.ns,\
	 problem.nc,\
	 problem.ncls,\
	 problem.n,\
	 problem.lambda);



  /**
   *   L-BFGS parameters
  */
  lbfgs_parameter_init(&param);


  /* Use default scaling */
  param.scaling = 1;

  /* Tolerance of inf norm of gradient */
  tmp = mxGetField(opt, 0, "epsg");
  if (tmp!=NULL && mxGetM(tmp)*mxGetN(tmp)>0) {
    param.epsg = *mxGetPr(tmp);
  }
  
  /* Tolerance of relative change in function value */
  tmp = mxGetField(opt, 0, "epsf");
  if (tmp!=NULL && mxGetM(tmp)*mxGetN(tmp)>0) {
    param.epsf = *mxGetPr(tmp);
  }

  /* Maximum number of iterations */
  tmp = mxGetField(opt, 0, "maxiter");
  if (tmp!=NULL && mxGetM(tmp)*mxGetN(tmp)>0) {
    param.max_iterations = (int)(*mxGetPr(tmp));
  }
  else {
    /* Never stops */
    param.max_iterations = 0;
  }

  /* Line search option */
  tmp = mxGetField(opt, 0, "lsmethod");
  if (tmp!=NULL && mxGetM(tmp)*mxGetN(tmp)>0) {
    param.linesearch = (int)(*mxGetPr(tmp));
  }
  else {
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
  }
  param.soc_nblocks = problem.nc;
  param.soc_kill    = problem.killsign; 


  /* Display option */
  tmp = mxGetField(opt, 0, "display");
  if (tmp!=NULL && mxGetM(tmp)*mxGetN(tmp)>0) {
    problem.display = (int)(*mxGetPr(tmp));
  }
  else {
    problem.display = 1;
  }

  /*
    Start the L-BFGS optimization; this will invoke the callback functions
    evaluate() and progress() when necessary.
  */
  ret = lbfgs(dd, ww, &fval, evaluate, progress, (void*)&problem, &param);

  /* Report the result. */
  switch(ret) {
  case LBFGS_SUCCESS:
    code = "LBFGS_SUCCESS";
    break;
  case LBFGS_PROGRESS_TOO_SLOW:
    code = "LBFGS_PROGRESS_TOO_SLOW";
    break;
  case LBFGSERR_MAXIMUMITERATION:
    code = "LBFGSERR_MAXIMUMITERATION";
    break;
  case LBFGS_ALREADY_MINIMIZED:
    code = "LBFGS_ALREADY_MINIMIZED";
    break;
  case LBFGSERR_ROUNDING_ERROR:
    code = "LBFGSERR_ROUNDING_ERROR";
    break;
  case LBFGSERR_MAXIMUMLINESEARCH:
    code = "LBFGSERR_MAXIMUMLINESEARCH";
    break;
  case LBFGSERR_INCREASEGRADIENT:
    code = "LBFGSERR_INCREASEGRADIENT";
    break;
  default:
    code = "Unknown";
  }

  if (problem.display>0) {
    mexPrintf("\n");
    mexPrintf("L-BFGS optimization terminated with status code = %d [%s]\n", ret, code);
    mexPrintf("  fval = %f, w = [%f %f %f %f]\n", fval, ww[0], ww[1], ww[2], ww[3]);
  }


  if (nlhs>2) {
    double *pret, *pfval, *gg;

    pret  = mxGetPr(mxGetField(plhs[2], 0, "ret"));
    pfval = mxGetPr(mxGetField(plhs[2], 0, "fval"));
    gg    = mxGetPr(mxGetField(plhs[2], 0, "g"));

    *pret  = ret;
    *pfval = evaluate((void*)&problem, ww, gg, dd, 0);
  }
  else {
    lbfgs_free(problem.s);
  }
    

  /* Free variables */
  if (nlhs<2) {
    lbfgs_free(problem.out);
  }

  /* Put NULLs so that mxDestroyArray does not destroy
     the structure entries */
  for (i=0;i<noptfields;i++) {
    mxSetField(opt, 0, optfields[i], NULL);
  }
  mxDestroyArray(opt);

  free(problem.killsign);
  free(problem.y);
  
  return;
}
        

