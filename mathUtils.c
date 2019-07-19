
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

/*******************************************************
*               Prototypes
********************************************************/

   double logbase (double , double );


/*******************************************************
*               Coding section
********************************************************/

/* Returns the
log base b of y. */

    double logbase (double y, double b)
    {
      double lg;
      lg = (double) log10(y)/log10(b);
      return(lg);
    }