/*
 * seq_art.h
 *
 *  Created on: 11.11.2013
 *      Author: wolfgang
 */

#ifndef SEQ_ART_H_
#define SEQ_ART_H_

#include <R.h>
#include <Rmath.h>

void rsim_unif(int n, double min, double max, double *x)
{
   int i;
   GetRNGstate();
   for (i = 0; i < n; ++i)
        x[i] = runif(min, max);
   PutRNGstate();
}


#endif /* SEQ_ART_H_ */
