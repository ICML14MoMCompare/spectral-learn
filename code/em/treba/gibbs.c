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

/* Gibbs sampling algorithm and structures for PFSA and HMMs */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdint.h>

#include "treba.h"

extern int g_alphabet_size;

struct gibbs_state_chain *gibbs_init_fsm(struct observations *o, int num_states, int alphabet_size, int *obslen) {
    int i,j,*data;
    struct gibbs_state_chain *go;
    struct observations *tempo;

    /* Initialize the sequence of observations for Gibbs sampling */
    /* It's a simple array of two ints, like so:                  */

    /*        0     ...  obslen-1                                 */
    /*     --------     --------                                  */
    /*     |state |     |state |                                  */
    /*     |sym   |     |sym   |                                  */
    /*     --------     --------                                  */
    
    /* We chain together multiple observations into one big       */
    /* observation by adding the symbol # at the end of each o_i  */
    /* which transitions to state 0. State 0 after # is never     */
    /* resampled. Also, we augment the alphabet by one to         */
    /* accommodate for #. After sampling, # is used to calculate  */
    /* the halting probability at that state, and is removed when */
    /* constructing an automaton based on the samples taken.      */

    for (tempo = o, (*obslen) = 0; tempo != NULL; tempo = tempo->next, (*obslen)++) {
	(*obslen) += tempo->size;
    }
    (*obslen)++;
    go = malloc(sizeof(struct gibbs_state_chain) * (*obslen));

    for (j = 0; o != NULL; o = o->next) {
	data = o->data;
	/* First state is always 0 */
	(go+j)->state = 0;
	for (i = 0; i < o->size; i++, j++) {
	    if (i > 0) {
		(go+j)->state = rand_int_range(0, num_states);
	    }
	    /* Add transition */
	    (go+j)->sym = *(data+i);
	}
	(go+j)->state = rand_int_range(0, num_states);
	(go+j)->sym = alphabet_size;
	j++;
    }
    (go+j)->state = 0;
    (go+j)->sym = -1; /* sentinel */
    return(go);
}

struct gibbs_state_chain *gibbs_init_hmm(struct observations *o, int num_states, int alphabet_size, int *obslen) {
    int i,j,*data;
    struct gibbs_state_chain *go;
    struct observations *tempo;

    /* Initialize the sequence of observations for Gibbs sampling */
    /* It's a simple array of two ints, like so:                  */

    /*        0     ...  obslen-1                                 */
    /*     --------     --------                                  */
    /*     |state |     |state |                                  */
    /*     |sym   |     |sym   |                                  */
    /*     --------     --------                                  */
    
    /* We chain together multiple observations into one big       */
    /* observation by forcing each o to begin with state 0 (INIT) */
    /* and ending in state n (END).                               */

    for (tempo = o, (*obslen) = 0; tempo != NULL; tempo = tempo->next) {
	(*obslen) += (tempo->size + 2);
    }

    go = malloc(sizeof(struct gibbs_state_chain) * (*obslen));

    for (j = 0; o != NULL; o = o->next) {
	data = o->data;
	/* First state is always 0 */
	(go+j)->state = 0;
	(go+j)->sym = -1;
	j++;
	for (i = 0; i < o->size; i++, j++) {
	    (go+j)->state = rand_int_range(1, num_states - 2);
	    (go+j)->sym = *(data+i);
	}
	(go+j)->state = num_states - 1;
	(go+j)->sym = -1;
	j++;
    }
    return(go);
}

/* Macros to access Gibbs HMM counts easily */

#define Ccurr(SOURCE_STATE, SYMBOL, TARGET_STATE) (*((gibbs_counts) + (num_states * alphabet_size * (SOURCE_STATE) + (SYMBOL) * num_states + (TARGET_STATE))))

#define Csampled(SOURCE_STATE, SYMBOL, TARGET_STATE) (*((gibbs_sampled_counts) + (num_states * alphabet_size * (SOURCE_STATE) + (SYMBOL) * num_states + (TARGET_STATE))))

#define CsampledHMMemit(STATE, SYMBOL) (*((gibbs_sampled_counts_emit) + (alphabet_size * (STATE) + (SYMBOL))))
#define CsampledHMMtrans(SOURCE_STATE, TARGET_STATE) (*((gibbs_sampled_counts_trans) + (num_states * (SOURCE_STATE) + (TARGET_STATE))))

struct hmm *gibbs_counts_to_hmm(struct hmm *hmm, unsigned int *gibbs_sampled_counts_trans, unsigned int *gibbs_sampled_counts_emit, unsigned int *gibbs_counts_sampled_states, int alphabet_size, int num_states, double beta_t, double beta_e) {
    int i, j;
    PROB newprob;
    /* Construct HMM from GS sampled counts */
    for (i = 0; i < num_states; i++) {
	gibbs_counts_sampled_states[i] = 0; /* Clear */
	for (j = 0; j < num_states; j++) {
	    gibbs_counts_sampled_states[i] += CsampledHMMtrans(i,j);
	}
    }
    for (i = 0 ; i < num_states; i++) {
	for (j = 0; j < num_states; j++) {
	    if ((i == 0 && j == 0) || i == num_states - 1 || j == 0) {
		newprob = SMRZERO_LOG;
	    }  else {
		newprob = log2(((double)CsampledHMMtrans(i,j) + beta_t) / ((double)gibbs_counts_sampled_states[i] + (num_states - 1) * beta_t));
	    }
	    *HMM_TRANSITION_PROB(hmm, i, j) = newprob;
	}
    }
    for (i = 0; i < num_states; i++) {
	for (j = 0; j < alphabet_size; j++) {
	    if (i == 0 || i == num_states - 1) {
		*HMM_EMISSION_PROB(hmm, i, j) = SMRZERO_LOG; /* states 0 and n don't emit */
	    } else {
		*HMM_EMISSION_PROB(hmm, i, j) = log2(((double)CsampledHMMemit(i,j) + beta_e) / ((double)gibbs_counts_sampled_states[i] + alphabet_size * beta_e));
	    }
	}
    }
    return(hmm);
}

struct wfsa *gibbs_counts_to_wfsa(struct wfsa *fsm, unsigned int *gibbs_sampled_counts, unsigned int *gibbs_counts_sampled_states, int alphabet_size, int num_states, double beta, double ANbeta) {
    int i,j,k;
    /* Construct WFSA */
    for (i = 0 ; i < num_states; i++) {
	gibbs_counts_sampled_states[i] = 0; /* Clear */
	*FINALPROB(fsm,i) = 0;
	for (j = 0; j < num_states; j++) {
	    for (k = 0; k < alphabet_size; k++) {
		gibbs_counts_sampled_states[i] += Csampled(i,k,j);
	    }
	}
    }
    for (i = 0 ; i < num_states; i++) {
	for (j = 0; j < num_states; j++) {
	    for (k = 0; k < alphabet_size; k++) {
		if (k == alphabet_size - 1) {
		    *FINALPROB(fsm,i) += ((double)Csampled(i,k,j) + beta) / ((double)gibbs_counts_sampled_states[i] + ANbeta);
		} else {
		    *(TRANSITION(fsm,i,k,j)) = log2(((double)Csampled(i,k,j) + beta) / ((double)gibbs_counts_sampled_states[i] + ANbeta));
		}
	    }
	}
    }
    for (i = 0; i < num_states; i++) {
	*FINALPROB(fsm,i) = log2(*FINALPROB(fsm,i));
    }
    return(fsm);
}

PROB gibbs_sampler_fsm(struct wfsa *fsm, struct observations *o, double beta, int num_states, int maxiter, int burnin, int lag) {

    int i, j, k, l, obslen, alphabet_size, z, zprev, znext, a, aprev, indicator, newstate, samplecount, high, low, mid;

    uint32_t steps, *gibbs_counts, *gibbs_counts_states;
    unsigned int *gibbs_sampled_counts, *gibbs_counts_sampled_states;

    PROB ANbeta, g_k, g_sum, *current_prob, cointoss;
    struct gibbs_state_chain *chain;
	
    alphabet_size = g_alphabet_size + 1; /* Use extra symbol for end-of-word (#) */
    /* Build initial array of states */
    chain = gibbs_init_fsm(o, num_states, g_alphabet_size, &obslen);

    gibbs_counts = calloc(num_states * num_states * alphabet_size, sizeof(unsigned int));
    gibbs_counts_states = calloc(num_states, sizeof(unsigned int));

    gibbs_counts_sampled_states = calloc(num_states, sizeof(unsigned int));
    gibbs_sampled_counts = calloc(num_states * num_states * alphabet_size, sizeof(unsigned int));

    current_prob = calloc(num_states, sizeof(PROB));

    /* Accumulate initial counts from initial random sequence */
    for (i = 0; i < obslen-1; i++) {
	Ccurr( (chain+i)->state , (chain+i)->sym, (chain+i+1)->state )++;
	gibbs_counts_states[(chain+i)->state]++;
    }

    ANbeta = alphabet_size * num_states * beta;

    for (i = 0, samplecount = 1, steps = 0; i < maxiter; i++) {
	for (j = 0; j < obslen; j++) {
	    if (j == 0 || ((chain+j-1)->sym == g_alphabet_size)) { /* Don't resample "initial" states, i.e. first */
		continue;                                          /* state in chain, or states preceded by #     */
	    }
	    
	    a = (chain+j)->sym;          /* Current symbol  */
	    aprev = (chain+j-1)->sym;    /* Previous symbol */
	    z = (chain+j)->state;        /* Current state   */
	    zprev = (chain+j-1)->state;  /* Previous state  */
	    znext = (chain+j+1)->state;  /* Next state      */

	    /* Adjust counts before evaluating new state */
	    Ccurr(zprev,aprev,z)--;
	    Ccurr(z,a,znext)--;

	    gibbs_counts_states[z]--;
	    for (k = 0, g_sum = 0; k < num_states; k++) {
		/* Sample */
		indicator = (k == zprev && aprev == a && znext == k) ? 1 : 0;
		g_k = (((double)Ccurr(k,a,znext)) + beta + indicator) * (((double)Ccurr(zprev,aprev,k)) + beta) / ((double)gibbs_counts_states[k] + ANbeta);
		assert(g_k >= 0);
		g_sum += g_k;
		current_prob[k] = g_sum;
	    }

	    /* Do #num_states-way weighted coin toss to select new state        */
	    /* This is done by a binary search through the array current_prob[] */
	    /* which has stored the sums of all probabilities.                  */
	    /* That is, we search for the first element in current_prob[]       */
	    /* greater than our random number (0 >= r < g_sum)                  */

	    cointoss = rand_double() * g_sum;
	    for (low = 0, high = num_states - 1; low != high; ) {
		mid = (low + high) / 2;
		if (current_prob[mid] <= cointoss) {
		    low = mid + 1;
		} else {
		    high = mid;
		}
	    }
	    newstate = high;

	    /* Adjust transition counts now that new state is chosen */
	    Ccurr(zprev,aprev,newstate)++;
	    Ccurr(newstate,a,znext)++;

	    gibbs_counts_states[newstate]++;
	    (chain+j)->state = newstate;
	    steps++;
	}
	if (i >= burnin && (i - burnin) % lag == 0) {
	    /* Update samplecount: sample count table entries will be updated */
	    /* as we go along based on the timestamp                          */
	    for (l = 0 ; l < obslen - 1; l++) {
		Csampled( (chain+l)->state, (chain+l)->sym, (chain+l+1)->state )++;
	    }
	    samplecount++;
	}
	if (i > 0 && (i == 10 || i == 100 || i % 1000 == 0)) {
	    if (i == 10 || i == 100)
		lag = 10;
	    burnin = i;
	    fsm = gibbs_counts_to_wfsa(fsm, gibbs_sampled_counts, gibbs_counts_sampled_states, alphabet_size, num_states, beta, ANbeta);
	    fprintf(stderr, "%i\t%.17g\n", i, loglikelihood_all_observations_fsm(fsm, o));
	    for (l = 0 ; l < num_states * num_states * alphabet_size; l++) {
		gibbs_sampled_counts[l] = 0;
	    }
	}

    }

    fsm = gibbs_counts_to_wfsa(fsm, gibbs_sampled_counts, gibbs_counts_sampled_states, alphabet_size, num_states, beta, ANbeta);
    free(chain);
    free(gibbs_counts);
    free(gibbs_sampled_counts);
    free(gibbs_counts_states);
    free(gibbs_counts_sampled_states);
    free(current_prob);
    return(loglikelihood_all_observations_fsm(fsm, o)); /* log likelihood */
}


PROB gibbs_sampler_hmm(struct hmm *hmm, struct observations *o, double beta_e, double beta_t, int num_states, int maxiter, int burnin, int lag) {

    #define CcurrHMMemit(STATE, SYMBOL) (*((gibbs_counts_emit) + (alphabet_size * (STATE) + (SYMBOL))))
    #define CcurrHMMtrans(SOURCE_STATE, TARGET_STATE) (*((gibbs_counts_trans) + (num_states * (SOURCE_STATE) + (TARGET_STATE))))

    int i, j, k, l, obslen, alphabet_size, z, zprev, znext, a, indicator, newstate, samplecount, high, low, mid;

    uint32_t steps, *gibbs_counts_trans, *gibbs_counts_emit, *gibbs_counts_states;
    unsigned int *gibbs_sampled_counts_trans, *gibbs_sampled_counts_emit, *gibbs_counts_sampled_states;

    PROB g_k, g_sum, *current_prob, cointoss;
    struct gibbs_state_chain *chain;

    alphabet_size = g_alphabet_size;

    /* Build initial array of states */
    chain = gibbs_init_hmm(o, num_states, g_alphabet_size, &obslen);

    gibbs_counts_trans = calloc(num_states * num_states, sizeof(unsigned int));
    gibbs_counts_emit = calloc(num_states * alphabet_size, sizeof(unsigned int));
    gibbs_counts_states = calloc(num_states, sizeof(unsigned int));

    gibbs_counts_sampled_states = calloc(num_states, sizeof(unsigned int));
    gibbs_sampled_counts_trans = calloc(num_states * num_states, sizeof(unsigned int));
    gibbs_sampled_counts_emit = calloc(num_states * alphabet_size, sizeof(unsigned int));

    current_prob = calloc(num_states, sizeof(PROB));

    /* Accumulate initial counts from initial random sequence */
    for (i = 0; i < obslen-1; i++) {
	CcurrHMMtrans( (chain+i)->state, (chain+i+1)->state )++;
	gibbs_counts_states[(chain+i)->state]++;
    }
    for (i = 0; i < obslen-1; i++) {
      if ((chain+i)->sym >= 0) {
	  CcurrHMMemit( (chain+i)->state, (chain+i)->sym)++;
      }
    }

    for (i = 0, samplecount = 1, steps = 0; i < maxiter; i++) {
	for (j = 0; j < obslen; j++) {
	    if ((chain+j)->sym < 0) {  /* Don't resample states INIT(0) or END(last state) */
		continue;
	    }
	    
	    a = (chain+j)->sym;          /* Current symbol  */
	    z = (chain+j)->state;        /* Current state   */
	    zprev = (chain+j-1)->state;  /* Previous state  */
	    znext = (chain+j+1)->state;  /* Next state      */

	    /* Adjust counts before evaluating new state */
	    CcurrHMMemit(z,a)--;
	    CcurrHMMtrans(z,znext)--;
	    CcurrHMMtrans(zprev,z)--;
	    gibbs_counts_states[z]--;

	    for (k = 1, g_sum = 0; k < num_states - 1; k++) {
		/* Sample */
		indicator = (k == zprev && znext == k) ? 1 : 0;
		g_k = ((((double)CcurrHMMemit(k,a)) + beta_e) / (((double)gibbs_counts_states[k]) + alphabet_size * beta_e)) *
		    (((((double)CcurrHMMtrans(zprev,k)) + beta_t) * (((double)CcurrHMMtrans(k,znext)) + indicator + beta_t)) /
                      (((double)gibbs_counts_states[k]) + num_states * beta_t));

		assert(g_k >= 0);
		g_sum += g_k;
		current_prob[k] = g_sum;
	    }

	    /* Do #num_states-way weighted coin toss to select new state        */
	    /* This is done by a binary search through the array current_prob[] */
	    /* which has stored the sums of all probabilities.                  */
	    /* That is, we search for the first element in current_prob[]       */
	    /* greater than our random number (0 >= r < g_sum)                  */

	    cointoss = rand_double() * g_sum;
	    for (low = 1, high = num_states - 2; low != high; ) {
		mid = (low + high) / 2;
		if (current_prob[mid] <= cointoss) {
		    low = mid + 1;
		} else {
		    high = mid;
		}
	    }
	    newstate = high;

	    /* Adjust transition counts now that new state is chosen */
	    CcurrHMMtrans(newstate,znext)++;
	    CcurrHMMtrans(zprev,newstate)++;
	    CcurrHMMemit(newstate,a)++;
	    gibbs_counts_states[newstate]++;
	    (chain+j)->state = newstate;
	    steps++;
	}
	if (i >= burnin && (i - burnin) % lag == 0) {
	    /* Update samplecount: sample count table entries will be updated */
	    /* as we go along based on the timestamp                          */
	    for (l = 0 ; l < obslen - 1; l++) {
		CsampledHMMtrans((chain+l)->state, (chain+l+1)->state)++;
		if ((chain+l)->sym >= 0) {
		    CsampledHMMemit((chain+l)->state, (chain+l)->sym)++;
		}
	    }
	    samplecount++;
	}
	/* For parallel computing tests */
	if (i > 0 && (i == 10 || i == 100 || i % 1000 == 0)) {
	    if (i == 10 || i == 100)
		lag = 10;
	    burnin = i;
	    hmm = gibbs_counts_to_hmm(hmm, gibbs_sampled_counts_trans, gibbs_sampled_counts_emit, gibbs_counts_sampled_states, alphabet_size, num_states, beta_t, beta_e);
	    fprintf(stderr, "%i\t%.17g\n", i, loglikelihood_all_observations_hmm(hmm, o));
	    for (l = 0 ; l < num_states * num_states; l++) {
		gibbs_sampled_counts_trans[l] = 0;
	    }
	    for (l = 0 ; l < num_states * alphabet_size; l++) {
		gibbs_sampled_counts_emit[l] = 0;
	    }
	}
    }

    hmm = gibbs_counts_to_hmm(hmm, gibbs_sampled_counts_trans, gibbs_sampled_counts_emit, gibbs_counts_sampled_states, alphabet_size, num_states, beta_t, beta_e);

    free(chain);
    free(gibbs_counts_trans);
    free(gibbs_counts_emit);
    free(gibbs_sampled_counts_trans);
    free(gibbs_sampled_counts_emit);
    free(gibbs_counts_states);
    free(gibbs_counts_sampled_states);
    free(current_prob);
    return(loglikelihood_all_observations_hmm(hmm, o)); /* log likelihood */
}
