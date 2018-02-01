#include "math.h"
#include "stdlib.h"
#include "R.h"
#include "Rmath.h"
#include "stdio.h"
#include "surv.h"

// modify
void SRCT(double *x, double *t, double *delta, double *Sc, int *n, double *RCTV)
{
	/*  computes RCT(x,y)  */
  int    i, j;
  double jp=0, p1 = 0, p2 = 0;
  *RCTV = 0;

  /*  I need to change  */
  for(i=0;i<(*n);i++){
		for(j=0;j<(*n);j++){
            if((x[j]>x[i]) && (t[j]>t[i]))
			    jp += 1;
            if(x[j]>x[i])
                p1 += 1;
            if(t[j]>t[i])
                p2 += 1;
		}
		if (delta[i] == 1) {
			*RCTV += pow(jp / (*n) - p1*p2 / ((*n)*(*n)), 2) / pow(Sc[i], 3);
		}
		jp = 0;
		p1 = 0;
		p2 = 0;
	}
  *RCTV = *RCTV/(1.0*(*n));
  
  return;
}



