/**************************************************************************/
/*   treba - probabilistic FSM and HMM training and decoding              */
/*   Copyright © 2013 Mans Hulden                                         */

/*   This file is part of treba.                                          */

/*   Treba is free software: you can redistribute it and/or modify        */
/*   it under the terms of the GNU General Public License version 2 as    */
/*   published by the Free Software Foundation.                           */

/*   Treba is distributed in the hope that it will be useful,             */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/*   GNU General Public License for more details.                         */

/*   You should have received a copy of the GNU General Public License    */
/*   along with treba.  If not, see <http://www.gnu.org/licenses/>.       */
/**************************************************************************/

/* Deterministic finite frequency-automaton based algorithms, including
   - state merging (ALERGIA, MDI, likelihood ratio based merges)
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "treba.h"

extern double chiinv(double p, int df);
extern double monte_carlo_multinomial_test(unsigned int *c0, double *p, int size, int iterations);
extern double approx_binomial_test(double p, int successes, int failures);
extern double flex_binomial_test(double p, int successes, int failures);

#define DFFA_TRANS(DFFA,STATE,SYMBOL) *((DFFA)->transitions + (STATE) * (DFFA)->alphabet_size + (SYMBOL))
#define DFFA_FREQ(DFFA,STATE,SYMBOL) *((DFFA)->transition_freqs + (STATE) * (DFFA)->alphabet_size + (SYMBOL))

#define DFFA_FINAL_FREQ(DFFA,STATE) *(((DFFA)->final_freqs + (STATE)))
#define DFFA_TOTAL_FREQ(DFFA,STATE) *(((DFFA)->total_freqs + (STATE)))

struct dffa_state {
    struct dffa_state *target;
    struct dffa_state *next_arc;
    int symbol;
    int freq;
    int number;
};

extern int g_t0;
extern PROB g_merge_prior;
int dffa_maxnumber;

void dffa_destroy(struct dffa *dffa) {
    if (dffa != NULL) {
        free(dffa->transitions);
        free(dffa->transition_freqs);
        free(dffa->final_freqs);
        free(dffa->total_freqs);
        free(dffa);
    }
}

struct dffa *dffa_copy(struct dffa *dffa) {
    struct dffa *dffa_new;
    dffa_new = calloc(1, sizeof(struct dffa));
    dffa_new->num_states = dffa->num_states;
    dffa_new->alphabet_size = dffa->alphabet_size;
    dffa_new->transitions = malloc(dffa->alphabet_size * dffa->num_states * sizeof(int));
    dffa_new->transition_freqs = malloc(dffa->alphabet_size * dffa->num_states * sizeof(int));
    dffa_new->final_freqs = malloc(dffa->num_states * sizeof(int));
    dffa_new->total_freqs = malloc(dffa->num_states * sizeof(int));
    memcpy(dffa_new->transitions, dffa->transitions, dffa->alphabet_size * dffa->num_states * sizeof(int));
    memcpy(dffa_new->transition_freqs, dffa->transition_freqs, dffa->alphabet_size * dffa->num_states * sizeof(int));
    memcpy(dffa_new->final_freqs, dffa->final_freqs, dffa->num_states * sizeof(int));
    memcpy(dffa_new->total_freqs, dffa->total_freqs, dffa->num_states * sizeof(int));
    return(dffa_new);
}

int *redarray, *redqueue, *bluearray, *bluequeue;
int lastred;

int dffa_find_incoming(struct dffa *dffa, int q, int *qf, int *a) {
    /* Returns the first (state,symbol)->q that it finds */
    /* Modifies both qf (source state) and a (symbol)    */
    /* qf must be red, q must be blue                    */

    int i, j;
    
    for (i = 0; i <= lastred; i++) {
	for (j = 0; j < dffa->alphabet_size; j++) {
	    if (DFFA_TRANS(dffa, redqueue[i], j) == q) {
		*qf = redqueue[i];
		*a = j;
		return 1;
	    }
	}
    }
    printf("No found incoming\n"); exit(1);
}

void dffa_stochastic_fold(struct dffa *dffa, int p, int q) {
    int a;
    if (p == q) {
	fprintf(stderr,"FATAL p = q in fold %i %i\n",p,q);
	exit(1);
	return;
    }
    if (DFFA_FINAL_FREQ(dffa, q) < 0) {
	//	fprintf(stderr,"<0 at %i %i\n",p,q); exit(1);
    }
    DFFA_FINAL_FREQ(dffa, p) += DFFA_FINAL_FREQ(dffa, q);
     for (a = 0; a < dffa->alphabet_size; a++) {
	if (DFFA_FREQ(dffa,q,a) > 0) {
	    if (DFFA_FREQ(dffa,p,a) > 0) {
		DFFA_FREQ(dffa,p,a) += DFFA_FREQ(dffa,q,a);
		dffa_stochastic_fold(dffa, DFFA_TRANS(dffa, p, a), DFFA_TRANS(dffa, q, a));
	    } else {
		DFFA_TRANS(dffa,p,a) = DFFA_TRANS(dffa,q,a);
		DFFA_FREQ(dffa,p,a) = DFFA_FREQ(dffa,q,a);
	    }
	}
    }
}

void dffa_update_freq(struct dffa *dffa, int p) {
    int i;
    DFFA_TOTAL_FREQ(dffa, p) = 0;
    for (i = 0; i < dffa->alphabet_size; i++) {
	DFFA_TOTAL_FREQ(dffa, p) += DFFA_FREQ(dffa,p,i);
    }
    DFFA_TOTAL_FREQ(dffa, p) += DFFA_FINAL_FREQ(dffa, p);
}

void dffa_stochastic_merge(struct dffa *dffa, int p, int q) {
    int qf, a;
    /* p = red state, q = blue state */
    /* find arc (qf, a) such that (qf,a) -> q */
    if (dffa_find_incoming(dffa, q, &qf, &a) == 0) {
	exit(0);
    }
    DFFA_TRANS(dffa, qf, a) = p;
    dffa_stochastic_fold(dffa, p, q);
}

void dffa_mark_accessible(struct dffa *dffa, int start, int *accessible, int *numaccessible) {
    /* Recursively mark accessible states. */
    int j;
    accessible[start] = 1;
    (*numaccessible)++;
    for (j = 0; j < dffa->alphabet_size; j++) {
	if (DFFA_FREQ(dffa, start, j) > 0 && !accessible[DFFA_TRANS(dffa,start,j)]) {
	    dffa_mark_accessible(dffa, DFFA_TRANS(dffa, start, j), accessible, numaccessible);
	}
    }    
}

/* Likelihood ratio test for merging qu and qv */
int dffa_chi_test_lr(struct dffa *dffa, int qu, int qv, double alpha) {
    int i;
    double ll1, ll2, ll, R, S, RS, Ri, Si, prior, criticalvalue;
    prior = 0.1;
    R = (double)(DFFA_TOTAL_FREQ(dffa, qu)) + ((double)dffa->alphabet_size * prior + prior);
    S = (double)(DFFA_TOTAL_FREQ(dffa, qv)) + ((double)dffa->alphabet_size * prior + prior);
    RS = R + S;
    for (i = 0, ll1 = 0, ll2 = 0; i < dffa->alphabet_size; i++) {
	Ri = (double)DFFA_FREQ(dffa, qu, i) + prior;
	Si = (double)DFFA_FREQ(dffa, qv, i) + prior;
	ll1 += Ri * log (((Ri+Si)/RS) / (Ri/R));
	ll2 += Si * log (((Ri+Si)/RS) / (Si/S));
    }
    Ri = (double)DFFA_FINAL_FREQ(dffa, qu) + prior;
    Si = (double)DFFA_FINAL_FREQ(dffa, qv) + prior;
    ll1 += Ri * log(((Ri+Si)/RS) / (Ri/R));
    ll2 += Si * log(((Ri+Si)/RS) / (Si/S));
    ll = -2 * (ll1 + ll2);
    //   fprintf(stderr, "LLe for %i and %i is %.17g [SS1: %i][SS2: %i]\n", qu, qv, ll, (int)R,(int)S);
    criticalvalue = chiinv(alpha, dffa->alphabet_size+1);    /* Set degrees of freedom */
    if (ll > criticalvalue) {
	return 0;
    }
    return 1;
}

int dffa_exact_multinomial_test(struct dffa *dffa, int qu, int qv, double alpha) {
    unsigned int i, *c0, *c1;
    double *p, p0, p1, prior;
    prior = 0.01;
    if (alpha + 0.2 < 1 && dffa_chi2_test(dffa, qu, qv, alpha + 0.2) == 1) {
	return 1;
    }
    if (alpha - 0.2 > 0 && dffa_chi2_test(dffa, qu, qv, alpha - 0.2) == 0) {
	return 0;
    }
    c0 = malloc((dffa->alphabet_size + 1) * sizeof(int));
    c1 = malloc((dffa->alphabet_size + 1) * sizeof(int));
    p = malloc((dffa->alphabet_size + 1) * sizeof(double));

    for (i = 0; i < dffa->alphabet_size; i++) {
	c0[i] = DFFA_FREQ(dffa, qu, i);
	c1[i] = DFFA_FREQ(dffa, qv, i);
	p[i] = (double)(DFFA_FREQ(dffa, qu, i) +  DFFA_FREQ(dffa, qv, i) + prior) / (double)(DFFA_TOTAL_FREQ(dffa, qu) + DFFA_TOTAL_FREQ(dffa, qv) + prior * dffa->alphabet_size + 1);
    }
    c0[i] = DFFA_FINAL_FREQ(dffa, qu);
    c1[i] = DFFA_FINAL_FREQ(dffa, qv);
    p[i] = (double)(DFFA_FINAL_FREQ(dffa, qu) +  DFFA_FINAL_FREQ(dffa, qv) + prior) / (double)(DFFA_TOTAL_FREQ(dffa, qu) + DFFA_TOTAL_FREQ(dffa, qv) + prior * dffa->alphabet_size + 1);
    p1 = monte_carlo_multinomial_test(c1, p, dffa->alphabet_size + 1, 10000);
    if (p1 < alpha) {
	free(c0); free(c1); free(p);
	return 0;
    }
    p0 = monte_carlo_multinomial_test(c0, p, dffa->alphabet_size + 1, 10000);
    if (p0 < alpha) {
	free(c0); free(c1); free(p);
	return 0;
    }
    free(c0); free(c1); free(p);
    //    fprintf(stderr, "p0: %.17g  p1: %.17g\n",p0,p1);
    return 1;
}

int dffa_binomial_test(struct dffa *dffa, int qu, int qv, double alpha) {
    double mle_p, prior;
    int i, successes_u, successes_v, failures_u, failures_v;
    prior = 0.01;
    mle_p = (double)(DFFA_FINAL_FREQ(dffa, qu) +  DFFA_FINAL_FREQ(dffa, qv) + prior) / (double)(DFFA_TOTAL_FREQ(dffa, qu) + DFFA_TOTAL_FREQ(dffa, qv) + prior * dffa->alphabet_size + 1);
    successes_u = DFFA_FINAL_FREQ(dffa,qu);
    successes_v = DFFA_FINAL_FREQ(dffa,qv);
    failures_u = DFFA_TOTAL_FREQ(dffa,qu) - successes_u;
    failures_v = DFFA_TOTAL_FREQ(dffa,qv) - successes_v;
    if (approx_binomial_test(mle_p,successes_u,failures_u) < alpha || approx_binomial_test(mle_p,successes_v,failures_v) < alpha) {
	return 0;
    }
    for (i = 0; i < dffa->alphabet_size; i++) {
	mle_p = (double)(DFFA_FREQ(dffa, qu, i) +  DFFA_FREQ(dffa, qv, i) + prior) / (double)(DFFA_TOTAL_FREQ(dffa, qu) + DFFA_TOTAL_FREQ(dffa, qv) + prior * dffa->alphabet_size + 1);
	successes_u = DFFA_FREQ(dffa,qu,i);
	successes_v = DFFA_FREQ(dffa,qv,i);
	failures_u = DFFA_TOTAL_FREQ(dffa,qu) - successes_u;
	failures_v = DFFA_TOTAL_FREQ(dffa,qv) - successes_v;
	 if (approx_binomial_test(mle_p,successes_u,failures_u) < alpha || approx_binomial_test(mle_p,successes_v,failures_v) < alpha) {
	     return 0;
	 }	 
    }
    return 1;
}

/* Exact binomial test for each transition + halting                               */
/* That is, we create the "merged" state as the MLE of qu + qv / total transitions */
/* and do a binomial test for each arc based on the MLE p where                    */      
/* success = a was followed, failure = some other arc was followed (or halt)       */

int dffa_exact_binomial_test(struct dffa *dffa, int qu, int qv, double alpha) {
    double mle_p, prior;
    int i, successes_u, successes_v, failures_u, failures_v, cutoff;
    cutoff = 15;
    prior = 0.01;

    mle_p = (double)(DFFA_FINAL_FREQ(dffa, qu) +  DFFA_FINAL_FREQ(dffa, qv) + prior) / (double)(DFFA_TOTAL_FREQ(dffa, qu) + DFFA_TOTAL_FREQ(dffa, qv) + prior * dffa->alphabet_size + 1);
    successes_u = DFFA_FINAL_FREQ(dffa,qu);
    successes_v = DFFA_FINAL_FREQ(dffa,qv);
    failures_u = DFFA_TOTAL_FREQ(dffa,qu) - successes_u;
    failures_v = DFFA_TOTAL_FREQ(dffa,qv) - successes_v;
    if (flex_binomial_test(mle_p,successes_u,failures_u) < alpha || approx_binomial_test(mle_p,successes_v,failures_v) < alpha) {
	return 0;
    }

    for (i = 0; i < dffa->alphabet_size; i++) {
	mle_p = (double)(DFFA_FREQ(dffa, qu, i) +  DFFA_FREQ(dffa, qv, i) + prior) / (double)(DFFA_TOTAL_FREQ(dffa, qu) + DFFA_TOTAL_FREQ(dffa, qv) + prior * dffa->alphabet_size + 1);
	successes_u = DFFA_FREQ(dffa,qu,i);
	successes_v = DFFA_FREQ(dffa,qv,i);
	failures_u = DFFA_TOTAL_FREQ(dffa,qu) - successes_u;
	failures_v = DFFA_TOTAL_FREQ(dffa,qv) - successes_v;
	if (flex_binomial_test(mle_p,successes_u,failures_u) < alpha || approx_binomial_test(mle_p,successes_v,failures_v) < alpha) {
	    return 0;
	}	

    }
    return 1;
}

/* Pearson's chi-squared test for merging qu and qv */
int dffa_chi2_test(struct dffa *dffa, int qu, int qv, double alpha) {
    int i;
    double chi1, chi2, R, S, Ri, Si, prior, criticalvalue, E;
    prior = 0.1;
    R = (double)(DFFA_TOTAL_FREQ(dffa, qu)) + ((double)dffa->alphabet_size * prior + prior);
    S = (double)(DFFA_TOTAL_FREQ(dffa, qv)) + ((double)dffa->alphabet_size * prior + prior);
    for (i = 0, chi1 = chi2 = 0; i < dffa->alphabet_size; i++) {
	Ri = (double)DFFA_FREQ(dffa, qu, i) + prior;
	Si = (double)DFFA_FREQ(dffa, qv, i) + prior;
	E = (Ri+Si)/(R+S) * R;
	chi1 += ((Ri-E)*(Ri-E))/E;
	E = (Ri+Si)/(R+S) * S;
	chi2 += ((Si-E)*(Si-E))/E;
	//printf("Seen transitions with %i: %.17g\n",i,Ri);
	//printf("Expected transitions with %i: %.17g\n",i,E);
    }
    Ri = (double)DFFA_FINAL_FREQ(dffa, qu) + prior;
    Si = (double)DFFA_FINAL_FREQ(dffa, qv) + prior;
    E = (Ri+Si)/(R+S) * R; 
    chi1 += ((Ri-E)*(Ri-E))/E;
    E = (Ri+Si)/(R+S) * S; 
    chi2 += ((Si-E)*(Si-E))/E;
    //  fprintf(stderr, "Chisquare for %i and %i is %.17g %.17g \n", qu, qv, chi1, chi2);
    criticalvalue = chiinv(alpha, dffa->alphabet_size);        /* Set degrees of freedom (slots - 1) */
    //fprintf(stderr, "Critical value is %.17g\n", criticalvalue);
    //exit(1);
    if (chi1 > criticalvalue) {
	return 0;
    }
    if (chi2 > criticalvalue) {
	return 0;
    }
    return 1;
}


/* Pearson's chi-squared test for merging qu and qv */
int dffa_chi_test(struct dffa *dffa, int qu, int qv, double alpha) {
    int i;
    double chisquare, ratio1, ratio2, R, S, t, Ri, Si, prior, criticalvalue;
    prior = 0.1;
    R = (double)(DFFA_TOTAL_FREQ(dffa, qu)) + ((double)dffa->alphabet_size * prior + prior);
    S = (double)(DFFA_TOTAL_FREQ(dffa, qv)) + ((double)dffa->alphabet_size * prior + prior);
    ratio1 = sqrt(S/R);
    ratio2 = sqrt(R/S);
    for (i = 0, chisquare = 0; i < dffa->alphabet_size; i++) {
	Ri = (double)DFFA_FREQ(dffa, qu, i) + prior;
	Si = (double)DFFA_FREQ(dffa, qv, i) + prior;
	t = ratio1 * Ri - ratio2 * Si;
	chisquare += (t * t) / (Ri + Si);
    }
    Ri = (double)DFFA_FINAL_FREQ(dffa, qu) + prior;
    Si = (double)DFFA_FINAL_FREQ(dffa, qv) + prior;
    t = ratio1 * Ri - ratio2 * Si;
    chisquare += (t * t) / (Ri + Si);
    //fprintf(stderr, "Chisquare for %i and %i is %.17g [SS1: %i][SS2: %i]\n", qu, qv, chisquare, (int)R,(int)S);
    criticalvalue = chiinv(alpha, dffa->alphabet_size+1);    /* Set degrees of freedom */
    
    if (chisquare > criticalvalue) {
	return 0;
    }
    return 1;
}

/* Hoeffding bound test for merging states (ALERGIA) */
int dffa_alergia_test(struct dffa *dffa, int qu, int qv, double alpha) {
    int a;
    int f1, n1, f2, n2;
    double gamma, bound;
   
    f1 = DFFA_FINAL_FREQ(dffa,qu);
    n1 = DFFA_TOTAL_FREQ(dffa,qu);
    f2 = DFFA_FINAL_FREQ(dffa,qv);
    n2 = DFFA_TOTAL_FREQ(dffa,qv);
    gamma = fabs(((double)f1)/((double)n1) - ((double)f2)/((double)n2));
    bound = ((sqrt(1.0/(double)n1) + sqrt(1.0/(double)n2)) * sqrt(log(2.0/alpha)))/1.41421356237309504880;
    if (gamma > bound)
	return 0;
    for (a = 0; a < dffa->alphabet_size; a++) {
	f1 = DFFA_FREQ(dffa,qu,a);
	n1 = DFFA_TOTAL_FREQ(dffa,qu);
	f2 = DFFA_FREQ(dffa,qv,a);
	n2 = DFFA_TOTAL_FREQ(dffa,qv);
	gamma = fabs(((double)f1)/((double)n1) - ((double)f2)/((double)n2));
	bound = ((sqrt(1.0/(double)n1) + sqrt(1.0/(double)n2)) * sqrt(log(2.0/alpha)))/1.41421356237309504880;
	if (gamma > bound)
	    return 0;
    }
    return 1;
}

int dffa_merge_compatible(struct dffa *dffa, int qu, int qv, double alpha, int (*mergetest)(), int recursive) {
    int i;
    if (mergetest(dffa,qu,qv,alpha) == 0) {
	return 0;
    }
    if (!recursive)
	return 1;
    for (i = 0; i < dffa->alphabet_size; i++) {
	if (DFFA_FREQ(dffa, qu, i) > 0 && DFFA_FREQ(dffa, qv, i) > 0) {
	    if (mergetest(dffa, DFFA_TRANS(dffa, qu, i), DFFA_TRANS(dffa, qv, i), alpha) == 0) {
	    	return 0;
	    }
	}
    }
    return 1;
}

double dffa_get_sequence_probability(struct dffa *dffa, int *sequence, int length, int initialstate) {
    int i, state;
    PROB totalprob, prob;
    for (i = 0, state = initialstate, totalprob = 0; i < length; i++) {
	if (DFFA_FREQ(dffa, state, sequence[i]) == 0) {
	    return 0;
	}
	prob = log(((double)DFFA_FREQ(dffa, state, sequence[i])) / ((double)DFFA_TOTAL_FREQ(dffa, state)));
	totalprob = totalprob + prob;
	state = DFFA_TRANS(dffa, state, sequence[i]);
    }
    if (DFFA_FINAL_FREQ(dffa, state) == 0) {
	return 0;
    }
    prob = log(((double)DFFA_FINAL_FREQ(dffa, state)) / ((double) DFFA_TOTAL_FREQ(dffa, state)));
    totalprob = totalprob + prob;
    return(totalprob);
}

double dffa_mdi_score(struct dffa *dffaold, struct dffa *dffanew, struct observations *o, int numobs) {
    int i, A, A1, A2, numaccessible, *accessible;
    double score, q, r;
    for (i = 0; i < dffaold->num_states; i++) {
	dffa_update_freq(dffaold, i);
	dffa_update_freq(dffanew, i);
    }
    for (q = 0, r = 0; o != NULL; o = o->next) {
	q += dffa_get_sequence_probability(dffaold, o->data, o->size, 0) * ((double)o->occurrences/(double)numobs) ;
	r += dffa_get_sequence_probability(dffanew, o->data, o->size, 0) * ((double)o->occurrences/(double)numobs) ;
    }
    score = q - r;
    numaccessible = 0;
    accessible = calloc(dffaold->num_states, sizeof(int));
    dffa_mark_accessible(dffaold, 0, accessible, &numaccessible);
    A1 = numaccessible;
    dffa_mark_accessible(dffanew, 0, accessible, &numaccessible);
    A2 = numaccessible;
    A = abs(A1 - A2);
    score = score / ((double)A);
    free(accessible);
    return(score);
}

double mdi_currentsum;
int mdi_statediff;
void dffa_mdi_parallel(struct dffa *dffa, int red, int blue, double rbprob, double bprob, int numobs) {
    /* Traverse and collect, starting from blue state bprob(x) and rbprob(x)  */
    /* bprob(x) = probability of string traversing only blue states           */
    /* rbprob(x) = probability of string traversing (and folding) red + blue  */
    int i;
    double px, fbprob, frbprob;
    //    fprintf(stderr,"Have red: %i [%.17g] blue: %i [%.17g]\n", red,exp(rbprob),blue, exp(bprob));

    if (red != -1) { mdi_statediff++; } /* If we arrive in red + blue, the new machine has 1 state fewer */
    if (DFFA_FINAL_FREQ(dffa, blue) > 0) {
	if (red != -1 && DFFA_FINAL_FREQ(dffa, red) > 0) {
	    frbprob = rbprob + log((double)(DFFA_FINAL_FREQ(dffa, blue) + DFFA_FINAL_FREQ(dffa, red)) / (double)(DFFA_TOTAL_FREQ(dffa, blue) + DFFA_TOTAL_FREQ(dffa, red)));
	} else {
	    rbprob = 0;
	}
	fbprob = bprob + log((double)DFFA_FINAL_FREQ(dffa, blue)/(double)DFFA_TOTAL_FREQ(dffa, blue));
	px = (double)DFFA_FINAL_FREQ(dffa, blue)/(double)numobs;
	mdi_currentsum += px * fabs(fbprob - frbprob);
    }
    for (i = 0; i < dffa->alphabet_size; i++) {
	if (DFFA_FREQ(dffa, blue, i) > 0) {
	    bprob += log((double)DFFA_FREQ(dffa, blue, i)/(double)DFFA_TOTAL_FREQ(dffa, blue));
	    if (red != -1 && DFFA_FREQ(dffa, red, i) > 0) {
		rbprob += log((double)(DFFA_FREQ(dffa, blue, i) + DFFA_FREQ(dffa, red, i)) / (double)(DFFA_TOTAL_FREQ(dffa, blue) + DFFA_TOTAL_FREQ(dffa, red)));
		dffa_mdi_parallel(dffa, DFFA_TRANS(dffa, red, i), DFFA_TRANS(dffa, blue, i), rbprob, bprob, numobs);
	    } else {
		dffa_mdi_parallel(dffa, -1, DFFA_TRANS(dffa, blue, i), rbprob, bprob, numobs);
	    }
	}
    }    
}

double dffa_mdi_score2(struct dffa *dffa, int red, int blue, int numobs) {
    double score;
    mdi_currentsum = 0.0 ;
    mdi_statediff = 0;
    dffa_mdi_parallel(dffa, red, blue, 0.0, 0.0, numobs);
    score = mdi_currentsum / (double) mdi_statediff;
    return(score);
}

int dffa_mdi_compatible(struct dffa *dffa, int qu, int qv, double alpha, struct observations *o, int numobs) {
    double score;
    score = dffa_mdi_score2(dffa, qu, qv, numobs);
    if (score < alpha)
	return 1;
    return 0;
}

struct dffa_state *find_add_trans(struct dffa_state *dffa, int symbol, int visits) {
    struct dffa_state *r;
    for (r = dffa; r != NULL; r = r->next_arc) {
	if (r->symbol == symbol) {
	    r->freq += visits;
	    return(r->target);
	}
    }
    r = calloc(1, sizeof(struct dffa_state));
    r->number = dffa->number;
    r->next_arc = dffa->next_arc;
    r->symbol = symbol;
    r->freq = visits;
    r->target = calloc(1, sizeof(struct dffa_state));
    r->target->number = dffa_maxnumber++;
    r->target->symbol = -1;
    dffa->next_arc = r;
    return(r->target);
}

struct dffa *dffa_accessible(struct dffa *dffa) {
    int i, j, *accessible, *accessiblemap, numaccessible;
    struct dffa *dffa_new;
    numaccessible = 0;
    accessible = calloc(dffa->num_states, sizeof(int));
    accessiblemap = calloc(dffa->num_states, sizeof(int));
    dffa_mark_accessible(dffa, 0, accessible, &numaccessible);
    fprintf(stderr, "Found %i accessible states\n", numaccessible);
    for (i = 0, j = 0; i < dffa->num_states; i++) {
	if (accessible[i]) {
	    accessiblemap[i] = j++;
	}
    }
    dffa_new = dffa_init(numaccessible, dffa->alphabet_size);
    for (i = 0; i < dffa->num_states; i++) {
	if (accessible[i]) {
	    DFFA_FINAL_FREQ(dffa_new, accessiblemap[i]) = DFFA_FINAL_FREQ(dffa, i);
	    for (j = 0; j < dffa->alphabet_size; j++) {
		if (DFFA_FREQ(dffa, i, j) > 0) {
		    DFFA_TRANS(dffa_new, accessiblemap[i], j) = accessiblemap[DFFA_TRANS(dffa, i, j)];
		    DFFA_FREQ(dffa_new, accessiblemap[i], j) = DFFA_FREQ(dffa, i, j);
		}
	    }
	    dffa_update_freq(dffa_new, accessiblemap[i]);
	}
    }
    free(accessible);
    free(accessiblemap);    
    return(dffa_new);
}

/* Generic state-merging algorithm.                                           */    
/* Starts with a trie of all strings, and repeatedly merges states.           */
/* States are divided into red, blue, and white states.                       */
/* Initially, only the root is red.                                           */
/* Blue states are always the descendants of all red states.                  */
/* The algorithm tries to merge any red state with a blue state               */
/* based on a compatibility test. If the test passes, the states are merged,  */
/* and a new set of blue states is calculated.                                */

/* Currently supported tests:                                                 */
/* ALERGIA:    calculates a Hoeffding bound for each transition frequency.    */
/* CHISQUARED: calculates the CHI^2-statistic from the (MLE of the) merged    */
/*             state and the (1) red state, and (2) blue state. If any of the */
/*             two states get a pvalue less than alpha the merge is rejected. */
/* LR:         Like chisquare, except a likelihood ratio test is used.        */
/* EXACTM:     Calculates the "exact" p-value between the MLE and the red and */
/*             the blue states.  This is done by Monte Carlo simulation of    */
/*             an exact multinomial test.                                     */
/* EXACTB:     Calculates the "exact" p-value between the MLE and the red and */
/*             the blue states (as a binomial separately for each pair of     */
/*             symbols).                                                      */
/* BINOMIAL:   Calculates the approximate binomial p-value.                   */


struct dffa *dffa_state_merge(struct observations *o, PROB alpha, int test, int recursive) {
    struct dffa *dffa;
    int qb, qr, i, j, lastblue, round;
    static int (*merge_test)() = &dffa_alergia_test;
    switch (test) {
    case MERGE_TEST_ALERGIA:
	merge_test = &dffa_alergia_test;
	break;
    case MERGE_TEST_CHISQUARED:
	merge_test = &dffa_chi2_test;
	g_t0 = 0;
	break;
    case MERGE_TEST_LR:
	merge_test = &dffa_chi_test_lr;
	g_t0 = 0;
	break;
    case MERGE_TEST_EXACT_M:
	g_t0 = 0;
	merge_test = &dffa_exact_multinomial_test;
	break;
    case MERGE_TEST_EXACT_B:
	g_t0 = 0;
	merge_test = &dffa_exact_binomial_test;
	break;
    case MERGE_TEST_BINOMIAL:
	g_t0 = 0;
	merge_test = &dffa_binomial_test;
    }
	    
    dffa = observations_to_dffa(o);
    redarray = calloc(dffa->num_states, sizeof(int));
    redqueue = calloc(dffa->num_states, sizeof(int));
    bluequeue = calloc(dffa->num_states, sizeof(int));
    bluearray = calloc(dffa->num_states, sizeof(int));
    redarray[0] = 1;
    redqueue[0] = 0;
    lastred = 0;
    lastblue = -1;
    /* Mark all descendants of 0 blue */
    for (j = 0; j < dffa->alphabet_size; j++) {
	if (DFFA_FREQ(dffa,0,j) > 0 && !redarray[DFFA_TRANS(dffa,0,j)]) {
	    //fprintf(stderr,"ADDING %i to blue\n", DFFA_TRANS(dffa,0,j));
	    bluequeue[++lastblue] = DFFA_TRANS(dffa,0,j);
	}
    }
    for (round = 1 ; ; round++) {
	//fprintf(stderr, "BLUE: ");
	for (i = 0 ; i <= lastblue; i++) {
	    //	    fprintf(stderr," [%i]", bluequeue[i]);
	}
	//fprintf(stderr,"\n");
	//fprintf(stderr, "RED: ");
	for (i = 0 ; i <= lastred; i++) {
	    //	    fprintf(stderr," [%i]", redqueue[i]);
	}

	//		fprintf(stderr,"\n****************\n");
	//fisher_yates_shuffle_array(bluequeue, lastblue);
	for (i = 0, qb = -1; i <= lastblue; i++) { /* Choose next blue state */
	    dffa_update_freq(dffa, bluequeue[i]);
	    if (DFFA_TOTAL_FREQ(dffa, bluequeue[i]) >= g_t0) {
		qb = bluequeue[i];
		//	fprintf(stderr, "Chose %i from blue\n", bluequeue[i]);
		break;
	    }
	}
	if (qb < 0) {
	    break;
	}
	for (i = 0, qr = -1; i <= lastred; i++) {
	    dffa_update_freq(dffa,redqueue[i]);
	    if (dffa_merge_compatible(dffa, redqueue[i], qb, alpha, merge_test, recursive)) {
		qr = redqueue[i];
		//fprintf(stderr, "Chose %i from red\n", redqueue[i]);
		break;
	    }
	}
	if (qb == qr) {
	    fprintf(stderr, "FATAL: %i is %i\n", qb, qr);
	}
	if (qr >= 0) {
	    dffa_stochastic_merge(dffa,qr,qb);
	} else {
	    redarray[qb] = 1;
	    redqueue[++lastred] = qb;
	}
	/* Mark all descendants of all red states blue (that aren't already red) */
	for (i = 0, lastblue = -1; i <= lastred; i++) {
	    for (j = 0; j < dffa->alphabet_size; j++) {
		if (DFFA_FREQ(dffa,redqueue[i],j) > 0 && !redarray[DFFA_TRANS(dffa,redqueue[i],j)] && bluearray[DFFA_TRANS(dffa,redqueue[i],j)] != round) {
		    bluequeue[++lastblue] = DFFA_TRANS(dffa,redqueue[i],j);
		    bluearray[DFFA_TRANS(dffa,redqueue[i],j)] = round;
		}
	    }
	}
    }
    dffa = dffa_accessible(dffa);
    free(bluequeue);
    free(redarray);
    free(redqueue);
    return(dffa);
}

struct dffa *dffa_mdi(struct observations *o, PROB alpha) {
    struct dffa *dffa;
    struct observations *obs;
    int qb, qr, i, j, lastblue, round, numobs;
    dffa = observations_to_dffa(o);
    for (numobs = 0, obs = o; obs != NULL; obs = obs->next) {
	numobs += obs->occurrences;
    }
    redarray = calloc(dffa->num_states, sizeof(int));
    redqueue = calloc(dffa->num_states, sizeof(int));
    bluequeue = calloc(dffa->num_states, sizeof(int));
    bluearray = calloc(dffa->num_states, sizeof(int));
    redarray[0] = 1;
    redqueue[0] = 0;
    lastred = 0;
    lastblue = -1;
    /* Mark all descendants of 0 blue */
    for (j = 0; j < dffa->alphabet_size; j++) {
	if (DFFA_FREQ(dffa,0,j) > 0 && !redarray[DFFA_TRANS(dffa,0,j)]) {
	    //fprintf(stderr,"ADDING %i to blue\n", DFFA_TRANS(dffa,0,j));
	    bluequeue[++lastblue] = DFFA_TRANS(dffa,0,j);
	}
    }
    for (round = 1 ; ; round++) {
	//	fprintf(stderr, "BLUE: ");
	//for (i = 0 ; i <= lastblue; i++) {
	//    fprintf(stderr," [%i]", bluequeue[i]);
	//}
	//fprintf(stderr,"\n");
	//fprintf(stderr, "RED: ");
	//for (i = 0 ; i <= lastred; i++) {
	//   fprintf(stderr," [%i]", redqueue[i]);
	//}

	//fprintf(stderr,"\n****************\n");
	//fisher_yates_shuffle_array(bluequeue, lastblue);
	for (i = 0, qb = -1; i <= lastblue; i++) { /* Choose next blue state */
	    dffa_update_freq(dffa, bluequeue[i]);
	    if (DFFA_TOTAL_FREQ(dffa, bluequeue[i]) >= g_t0) {
		qb = bluequeue[i];
		//fprintf(stderr, "Chose %i from blue\n", bluequeue[i]);
		break;
	    }
	}
	if (qb < 0) {
	    break;
	}

	for (i = 0, qr = -1; i <= lastred; i++) {
	    dffa_update_freq(dffa,redqueue[i]);
	    if (dffa_mdi_compatible(dffa, redqueue[i], qb, alpha, o, numobs)) {
		//fprintf(stderr,"Found compatible %i and %i\n", redqueue[i]  ,qb);
		qr = redqueue[i];
		//fprintf(stderr, "Chose %i from red\n", redqueue[i]);
		break;
	    }
	}
	if (qb == qr) {
	    fprintf(stderr, "FATAL: %i is %i\n", qb, qr);
	}
	if (qr >= 0) {
	    dffa_stochastic_merge(dffa,qr,qb);
	} else {
	    redarray[qb] = 1;
	    redqueue[++lastred] = qb;
	}
	/* Mark all descendants of all red states blue (that aren't already red) */
	for (i = 0, lastblue = -1; i <= lastred; i++) {
	    for (j = 0; j < dffa->alphabet_size; j++) {
		if (DFFA_FREQ(dffa,redqueue[i],j) > 0 && !redarray[DFFA_TRANS(dffa,redqueue[i],j)] && bluearray[DFFA_TRANS(dffa,redqueue[i],j)] != round) {
		    bluequeue[++lastblue] = DFFA_TRANS(dffa,redqueue[i],j);
		    bluearray[DFFA_TRANS(dffa,redqueue[i],j)] = round;
		}
	    }
	}
    }
    dffa = dffa_accessible(dffa);
    free(bluequeue);
    free(redarray);
    free(redqueue);
    return(dffa);
}

struct wfsa *dffa_to_wfsa(struct dffa *dffa) {
    struct wfsa *fsm;
    double prior;
    int i, j, k;
    prior = g_merge_prior;
    fsm = wfsa_init(dffa->num_states, dffa->alphabet_size);
    fprintf(stderr, "Converting to PFSA: %i states %i symbols\n", dffa->num_states, dffa->alphabet_size);
    for (i = 0; i < dffa->num_states; i++) {
	for (j = 0; j < dffa->num_states; j++) {
	    for (k = 0; k < dffa->alphabet_size; k++) {	
		*TRANSITION(fsm,i,k,j) = SMRZERO_LOG;
	    }
	}
    }
    for (i = 0; i < dffa->num_states; i++) {
	if (DFFA_FINAL_FREQ(dffa,i) > 0) {
	    *FINALPROB(fsm,i) = LOG( ((((double) DFFA_FINAL_FREQ(dffa,i)) + prior)) / ( (double) DFFA_TOTAL_FREQ(dffa,i) + ((dffa->alphabet_size + 1) * prior )));
	} else {
	    *FINALPROB(fsm,i) = LOG( prior / ( (double) DFFA_TOTAL_FREQ(dffa,i) + ((dffa->alphabet_size + 1) * prior )));
	}
	for (j = 0; j < dffa->alphabet_size; j++) {
	    if (DFFA_FREQ(dffa, i, j) > 0) {
		*TRANSITION(fsm, i, j, DFFA_TRANS(dffa, i, j)) = LOG( (((double)DFFA_FREQ(dffa,i,j)) + prior ) / ((double) DFFA_TOTAL_FREQ(dffa, i) + ((dffa->alphabet_size + 1) * prior)));
	    } else  {
		*TRANSITION(fsm, i, j, i) = LOG( prior / ((double) DFFA_TOTAL_FREQ(dffa, i) + ((dffa->alphabet_size + 1) * prior)));
	    }
	}
    }
    return(fsm);
}

struct dffa *dffa_init(int num_states, int alphabet_size) {
    struct dffa *dffa;
    dffa = malloc(sizeof(struct dffa));
    dffa->num_states = num_states;
    dffa->alphabet_size = alphabet_size;
    dffa->transitions = calloc(num_states * alphabet_size, sizeof(int));
    dffa->transition_freqs = calloc(num_states * alphabet_size, sizeof(int));
    dffa->final_freqs = calloc(num_states, sizeof(int));
    dffa->total_freqs = calloc(num_states, sizeof(int));
    return(dffa);
}

void dffa_chain_to_dffa(struct dffa_state *ds, struct dffa *dffa) {
    if (ds->number >= 0) {
	for ( ; ds != NULL; ds = ds->next_arc) {
	    if (ds->symbol >= 0) {
		DFFA_TRANS(dffa, ds->number, ds->symbol) = ds->target->number;
		DFFA_FREQ(dffa, ds->number, ds->symbol) = ds->freq;
		if (ds->target != NULL) {
		    dffa_chain_to_dffa(ds->target, dffa);
		}
	    } else {
		DFFA_FINAL_FREQ(dffa, ds->number) = ds->freq;
	    }
	    DFFA_TOTAL_FREQ(dffa, ds->number) += ds->freq;
	    ds->number = -1;
	}
    }
}

struct dffa *observations_to_dffa(struct observations *o) {
    struct dffa_state *pos, *root;
    struct dffa *dffa;
    int i, *obsdata, maxsym;
    root = calloc(1,sizeof(struct dffa_state));
    root->symbol = -1;
    root->number = 0;
    dffa_maxnumber = 1;
    for (maxsym = 0 ; o != NULL; o = o->next) {
	for (pos = root, i = 0, obsdata = o->data; i < o->size; i++) {
	    maxsym = obsdata[i] > maxsym ? obsdata[i] : maxsym;
	    pos = find_add_trans(pos, obsdata[i], o->occurrences);
	}
	pos = find_add_trans(pos, -1, o->occurrences);
    }
    dffa = dffa_init(dffa_maxnumber, maxsym + 1);
    dffa_chain_to_dffa(root, dffa);
    return(dffa);
}
