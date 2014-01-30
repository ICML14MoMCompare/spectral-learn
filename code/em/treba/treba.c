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

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <pthread.h>
#include <signal.h>
#include <getopt.h>

#include "treba.h"
#include "globals.h"
#include "fastlogexp.h"

#ifdef _WIN32
 #include <windows.h>
 #define srandom srand
 #define random rand
#endif

static char *usagestring = 
"treba [options] observationfile\n";
static char *helpstring =
"\n"
"Train/induce FSMs/HMMs, calculate probabilities and best paths of sequences.\n\n"

"Main options:\n\n"
" -T , --train=ALG        Train FSM with state merging, Baum-Welch, B-W +\n"
"                         deterministic annealing, Gibbs sampling, Viterbi\n"
"                         training, Viterbi+B-W,...\n"
"                         ALG is one of merge,mdi,bw,dabw,gs,vb,vit,vitbw\n"
"                         merge = State-merging algorithms (such as ALERGIA)\n"
"                         mdi = MDI algorithm\n"
"                         bw = Baum-Welch, dabw = Baum-Welch w/ det. annealing.\n"
"                         vit = Viterbi-Baum-Welch (Hard EM)\n"
"                         vitbw = vit until convergence followed by B-W\n"
"                         vb = Variational Bayes\n"
"                         gs = Gibbs sampling\n"
" -C , --cuda             Use CUDA parallellization (for Gibbs sampling).\n"
"                         Requires compiled support to cuda.\n"
" -H , --hmm              Operate on HMM (instead of PFSA).\n"
" -D , --decode=METHOD    Decode (find best path through automaton) with\n"
"                         forward, backward, or Viterbi.\n"
"                         METHOD one of f[,p],b[,p],vit[,p]; adding ,p specifies\n"
"                         probabilities to be printed as well as the path.\n"
" -L , --likelihood=TYPE  Calculate probability of sequences; forward\n"
"                         probability or best path (Viterbi). TYPE one of f,vit,b\n"
"                         (forward, Viterbi, backward)\n"
" -G , --generate=NUM     Generate (randomly) NUM words from HMM of HMM/PFSA\n"
" -M , --merge=ALG        Set merge test for merge-based learning algorithms.\n"
"                         ALG one of alergia,chi2,lr,binomial,exactm,exact\n"
" -i , --input-format=F   Set input probability/weight format from FSMs/HMMs\n"
"                         F one of real,log10,ln,log2,nlog10,nln,nlog2.\n"
"                         n prefix to FMT is negative. Default is real.\n"
" -o , --output-format=F  Set  probability/weight format from FSMs/HMMs\n"
"                         F one of real,log10,ln,log2,nlog10,nln,nlog2.\n"
"                         n prefix to FMT is negative. Default is real.\n"
" -f , --file=FILE        Specify FSM/PFSA/HMM to read from file.\n"
" -g , --initialize=TYPE  Specify type of initial random FSM. TYPE in\n"
"                         [b,d]NUMSTATES[,#NUMSYMS] b=Bakis,d=deterministic,n=ergodic\n"
" -u , --uniform-probs    Set uniform probabilities on initially generated FSM.\n"
"\n"
"Training options:\n\n"
" -A , --alpha=NUM        Main parameter for state merging algorithms\n"
"                         Default: 0.05 (ALERGIA/MDI)\n"
" -y , --t0=NUM           t0-parameter for ALERGIA/MDI\n"
" -b , --burnin=NUM       Burnin iterations for Gibbs sampler\n"
" -l , --lag=NUM          Lag of Gibbs sampler.\n"
" -d , --max-delta=NUM    Maximum change in log likelihood between iterations\n"
"                         for convergence. Default is 0.1.\n"
" -x , --max-iter         Maximum number of iterations to run. Default is 100000\n"
" -p , --prior=NUM        Pseudocounts/parameter to use in various algorithms\n"
"                         Default is 0.02 (Gibbs/Var. Bayes), 1.0 (Hard EM), \n"
"                         For HMMs, use two comma separated values NUM1,NUM2\n"
"                         NUM1 = state-to-state prior, NUM2 = emission prior\n"
" -r , --restarts=OPT     Number of restarts and iterations per restart before\n"
"                         beginning final run of B-W. Restart OPT is given as\n"
"                         numrestarts,iterations-per-restart\n"
" -R , --recursive-merge  Do merge tests recursively (for merging algorithms).\n"
" -a , --annealopts=PAR   Parameters for deterministic annealing.\n"
"                         PAR specified as betamin,betamax,alpha.\n"
" -t , --threads=NUM      Number of threads to launch in parallel for Baum-Welch\n"
"                         Can be specified as fraction of available CPUs c/NUM\n";

PROB g_loglikelihood = 0;

#ifdef USE_CUDA
 static char *versionstring = "treba v1.01 (compiled with CUDA support)";
 extern double gibbs_sampler_cuda_fsm(struct wfsa *fsm, struct observations *o, double beta, int num_states, int maxiter, int burnin, int lag);
 extern double gibbs_sampler_cuda_hmm(struct hmm *hmm, struct observations *o, double beta_t, double beta_e, int num_states, int maxiter, int burnin, int lag);
#else
 static char *versionstring = "treba v1.01 (compiled without CUDA support)";
#endif /* USE_CUDA */

PROB rand_double() {
    return (drand48());
}

/* Returns int from <= r < to */ 
int rand_int_range(int from, int to) {
    double ffrom = (double)from;
    double fto = (double)to;
    return from + (int)(fto * random() / (RAND_MAX + ffrom));
}

void wfsa_randomize_deterministic(struct wfsa *fsm, int uniform) {
    int i, j, k, tablesize;
    PROB randnum, randsum, *tableptr;
    for (i=0; i < fsm->num_states; i++) {
        randsum = 0;
        tablesize = fsm->num_states * fsm->alphabet_size;
        for (j = 0; j < fsm->alphabet_size; j++) {
            k = rand_int_range(0,fsm->num_states);
            tableptr = TRANSITION(fsm,i,j,k);
            randnum = (uniform == 0) ? rand_double() : 1;
            *tableptr = randnum;
            randsum += randnum;
        }
        randnum = (uniform == 0) ? rand_double() : 1;
        *(fsm->final_table + i) = randnum;
        randsum += randnum;
        for (tableptr = TRANSITION(fsm,i,0,0), j = 0; j < tablesize; j++) {
            *(tableptr+j) /= randsum;
        }
        *(fsm->final_table + i) /= randsum;
    }
}

void hmm_normalize(struct hmm *hmm) {
    int source, target, symbol;
    double sum;
    for (source = 1; source < hmm->num_states - 1; source++) {
	sum = 0.0;
	for (symbol = 0 ; symbol < hmm->alphabet_size; symbol++) {
	    sum += *HMM_EMISSION_PROB(hmm, source, symbol);
	}
	for (symbol = 0 ; symbol < hmm->alphabet_size; symbol++) {
	    *HMM_EMISSION_PROB(hmm, source, symbol) /= sum;
	}
    }

    for (target = 1; target < hmm->num_states; target++) {
	for (source = 0, sum = 0.0; source < hmm->num_states - 1; source++) {
	    sum += *HMM_TRANSITION_PROB(hmm, source, target);
	}
	for (source = 0; source < hmm->num_states - 1; source++) {
	    *HMM_TRANSITION_PROB(hmm, source, target) /= sum;
	}
    }
}

void hmm_randomize(struct hmm *hmm, int bakis, int uniform) {
    int source, target, symbol;
    for (source = 0; source < hmm->num_states; source++) {
	for (target = 0; target < hmm->num_states; target++) {
	    if (target == 0 || source == hmm->num_states - 1) {
		*HMM_TRANSITION_PROB(hmm, source, target) = 0.0;
	    } else if (bakis && target < source) {
		*HMM_TRANSITION_PROB(hmm, source, target) = 0.0;		    		
	    } else {
		*HMM_TRANSITION_PROB(hmm, source, target) = rand_double();
	    }
	}
    }
    for (source = 0; source < hmm->num_states; source++) {
	for (symbol = 0; symbol < hmm->alphabet_size; symbol++) {
	    if (source == 0 || source == hmm->num_states - 1) {
		*HMM_EMISSION_PROB(hmm, source, symbol) = 0.0;
	    } else {
		*HMM_EMISSION_PROB(hmm, source, symbol) = rand_double();
	    }
	}
    }
    hmm_normalize(hmm);
}

void wfsa_randomize_nondeterministic(struct wfsa *fsm, int bakis, int uniform) {
    int i, j, k, tablesize;
    PROB randnum, randsum, *tableptr;
    for (i = 0; i < fsm->num_states; i++) {
	randsum = 0;
	tablesize = fsm->num_states * fsm->alphabet_size;
	for (j = (bakis == 1 ? i : 0); j < fsm->num_states; j++) {
	    for (k = 0 ; k < fsm->alphabet_size; k++) {
		randnum = (uniform == 0) ? rand_double() : 1;
		*TRANSITION(fsm,i,k,j) = randnum;
		randsum += randnum;
	    }
	}
	randnum = (uniform == 0) ? rand_double() : 1;
	*(fsm->final_table + i) = randnum;
	randsum += randnum;
	for (tableptr = TRANSITION(fsm,i,0,0), j = 0; j < tablesize; j++) {
	    *(tableptr+j) /= randsum;
	}
	*(fsm->final_table + i) /= randsum;
    }
}

struct hmm *hmm_init(int num_states, int alphabet_size) {
    struct hmm *hmm;
    hmm = malloc(sizeof(struct hmm));
    hmm->num_states = num_states;
    hmm->alphabet_size = alphabet_size;
    hmm->transition_table = calloc(num_states * num_states, sizeof(PROB));
    hmm->emission_table = calloc(num_states * alphabet_size, sizeof(PROB));
    return(hmm);
}

struct wfsa *wfsa_init(int num_states, int alphabet_size) {
    struct wfsa *fsm;
    fsm = malloc(sizeof(struct wfsa));
    if (fsm == NULL) {
	fprintf(stderr, "Out of memory. Fatal.\n"); exit(1);
    }
    fsm->num_states = num_states;
    fsm->alphabet_size = alphabet_size;
    fsm->state_table = calloc(num_states * num_states * alphabet_size, sizeof(PROB));
    if (fsm->state_table == NULL) {
	fprintf(stderr, "Out of memory. Fatal.\n"); exit(1);
    }
    fsm->final_table = calloc(num_states, sizeof(PROB));
    return(fsm);
}

struct hmm *hmm_copy(struct hmm *hmm) {
    struct hmm *newhmm;
    newhmm = malloc(sizeof(struct hmm));
    newhmm->num_states = hmm->num_states;
    newhmm->alphabet_size = hmm->alphabet_size;
    newhmm->transition_table = malloc(hmm->num_states * hmm->num_states * sizeof(PROB));
    newhmm->emission_table = malloc(hmm->num_states * hmm->alphabet_size * sizeof(PROB));
    memcpy(newhmm->transition_table, hmm->transition_table, hmm->num_states * hmm->num_states * sizeof(PROB));
    memcpy(newhmm->emission_table, hmm->emission_table, hmm->num_states * hmm->alphabet_size * sizeof(PROB));
    return(newhmm);
}

struct wfsa *wfsa_copy(struct wfsa *fsm) {
    struct wfsa *newfsm;
    newfsm = malloc(sizeof(struct wfsa));
    newfsm->num_states = fsm->num_states;
    newfsm->alphabet_size = fsm->alphabet_size;
    newfsm->state_table = malloc(fsm->num_states * fsm->num_states * fsm->alphabet_size * sizeof(PROB));
    newfsm->final_table = malloc(fsm->num_states * sizeof(PROB));
    memcpy(newfsm->state_table, fsm->state_table, fsm->num_states * fsm->num_states * fsm->alphabet_size * sizeof(PROB));
    memcpy(newfsm->final_table, fsm->final_table, fsm->num_states * sizeof(PROB));
    return(newfsm);
}

void wfsa_destroy(struct wfsa *fsm) {
    free(fsm->state_table);
    free(fsm->final_table);
    free(fsm);
}

void hmm_destroy(struct hmm *hmm) {
    free(hmm->transition_table);
    free(hmm->emission_table);
    free(hmm);
}

void wfsa_to_log2(struct wfsa *fsm) {
    int i,j,k;
    for (i = 0; i < fsm->num_states; i++) {
	*FINALPROB(fsm,i) = input_convert(*FINALPROB(fsm,i));
    }
    for (i = 0; i < fsm->num_states; i++) {
	for (j = 0; j < fsm->alphabet_size; j++) {
	    for (k = 0; k < fsm->num_states; k++) {
		*TRANSITION(fsm,i,j,k) = input_convert(*TRANSITION(fsm,i,j,k));
	    }
	}
    }
}

void hmm_to_log2(struct hmm *hmm) {
    int i,j;
    for (i = 0; i < hmm->num_states; i++) {
	for (j = 0; j < hmm->num_states; j++) {
	    *HMM_TRANSITION_PROB(hmm,i,j) = input_convert(*HMM_TRANSITION_PROB(hmm,i,j));
	}
    }
    for (i = 0; i < hmm->num_states; i++) {
	for (j = 0; j < hmm->alphabet_size; j++) {
	    *HMM_EMISSION_PROB(hmm,i,j) = input_convert(*HMM_EMISSION_PROB(hmm,i,j));
	}
    }
}

struct wfsa *wfsa_read_file(char *filename) {
    char *wfsa_char_data, *w, *lastline;
    int elements, source, target, symbol, finalstate, maxstate, maxsymbol;
    PROB prob;
    struct wfsa *fsm;
    if ((wfsa_char_data = file_to_mem(filename)) == NULL) {
	exit(1);
    }
    /* Figure out alphabet size and number of states */
    for (w = wfsa_char_data, maxstate = 0, maxsymbol = 0; ; ) {
	lastline = w;
	elements = line_count_elements(&w);
	if (elements == -1) {
	    break;
	}
	if (elements == 0) {
	    continue; /* Comment line */
	}
	if (strstr(lastline, ">") != NULL) {
	    fprintf(stderr, "ERROR: Expecting FSA file: for HMMs use the --hmm flag.\n");
	    exit(EXIT_FAILURE);
	}
	switch (elements) {
	case 1:
	    sscanf(lastline, "%i", &finalstate);
	    maxstate = maxstate > finalstate ? maxstate : finalstate;
	    break;
	case 2:
	    sscanf(lastline, "%i %lg", &finalstate, &prob);
	    maxstate = maxstate > finalstate ? maxstate : finalstate;
	    break;
	case 3:
	    sscanf(lastline, "%i %i %i", &source, &target, &symbol);
	    maxstate = maxstate > source ? maxstate : source;
	    maxstate = maxstate > target ? maxstate : target;
	    maxsymbol = maxsymbol > symbol ? maxsymbol : symbol;
	    break;
	case 4:
	    sscanf(lastline, "%i %i %i %lg", &source, &target, &symbol, &prob);	    
	    maxstate = maxstate > source ? maxstate : source;
	    maxstate = maxstate > target ? maxstate : target;
	    maxsymbol = maxsymbol > symbol ? maxsymbol : symbol;
	    break;
	default:
	    perror("WFSA file format error");
	    free(wfsa_char_data);
	    exit(1);
	}
    }
    fsm = wfsa_init(maxstate+1, maxsymbol+1);
    for (w = wfsa_char_data, maxstate = 0, maxsymbol = 0; ; ) {
	lastline = w;
	elements = line_count_elements(&w);
	if (elements == 0) {
	    continue; /* Comment line */
	}
	if (elements == -1) {
	    break;
	}
	switch (elements) {
	case 1:
	    sscanf(lastline, "%i", &finalstate);
	    *FINALPROB(fsm, finalstate) = SMRONE_REAL;
	    break;
	case 2:
	    sscanf(lastline, "%i %lg", &finalstate, &prob);
	    *FINALPROB(fsm, finalstate) = prob;
	    break;
	case 3:
	    sscanf(lastline, "%i %i %i", &source, &target, &symbol);
	    *TRANSITION(fsm, source, symbol, target) = SMRONE_REAL;
	    break;
	case 4:
	    sscanf(lastline, "%i %i %i %lg", &source, &target, &symbol, &prob);
	    *TRANSITION(fsm, source, symbol, target) = prob;
	    break;
	default:
	    perror("WFSA file format error");
	    free(wfsa_char_data);
	    exit(1);
	}
    }
    free(wfsa_char_data);
    return(fsm);
}

struct hmm *hmm_read_file(char *filename) {
    char *hmm_char_data, *w, *lastline;
    int elements, source, target, symbol, maxstate, maxsymbol;
    PROB prob;
    struct hmm *hmm;
    if ((hmm_char_data = file_to_mem(filename)) == NULL) {
	exit(1);
    }
    /* Figure out alphabet size and number of states */
    for (w = hmm_char_data, maxstate = 0, maxsymbol = 0; ; ) {
	lastline = w;
	elements = line_count_elements(&w);
	if (elements == 0) {
	    continue; /* Comment line */
	}
	if (elements == -1) {
	    break;
	}
	switch (elements) {
	case 3:
	    sscanf(lastline, "%i %i %lg", &source, &symbol, &prob); /* Transition probability */
	    maxsymbol = maxsymbol > symbol ? maxsymbol : symbol;
	    maxstate = maxstate > source ? maxstate : source;
	    maxstate = maxstate > target ? maxstate : target;
	    maxsymbol = maxsymbol > symbol ? maxsymbol : symbol;
	    break;
	case 4:
	    sscanf(lastline, "%i > %i %lg", &source, &target, &prob); /* Emission probability */
	    maxstate = maxstate > source ? maxstate : source;
	    maxstate = maxstate > target ? maxstate : target;
	    break;
	default:
	    perror("HMM file format error");
	    free(hmm_char_data);
	    exit(1);
	}
    }

    hmm = hmm_init(maxstate+1, maxsymbol+1);
    for (w = hmm_char_data, maxstate = 0, maxsymbol = 0; ; ) {
	lastline = w;
	elements = line_count_elements(&w);
	if (elements == 0) {
	    continue; /* Comment line */
	}
	if (elements == -1) {
	    break;
	}
	switch (elements) {
	case 3:
	    sscanf(lastline, "%i %i %lg", &source, &symbol, &prob); /* Transition probability */
	    *HMM_EMISSION_PROB(hmm, source, symbol) = prob;
	    break;
	case 4:
	    sscanf(lastline, "%i > %i %lg", &source, &target, &prob); /* Emission probability */
	    *HMM_TRANSITION_PROB(hmm, source, target) = prob;
	    break;
	default:
	    perror("HMM file format error");
	    free(hmm_char_data);
	    exit(1);
	}
    }
    free(hmm_char_data);
    return(hmm);
}

struct wfsa *g_lastwfsa = NULL; /* Stores global pointer to last WFSA to spit out in case of SIGINT */
struct hmm *g_lasthmm = NULL; /* Stores global pointer to last WFSA to spit out in case of SIGINT */

void interrupt_sigproc() {
    fprintf(stderr, "Received SIGINT. Exiting.\n");
    if (g_lastwfsa != NULL) {
	wfsa_print(g_lastwfsa);
    }
    exit(1);
}

void interrupt_sigproc_hmm() {
    fprintf(stderr, "Received SIGINT. Exiting.\n");
    if (g_lasthmm != NULL) {
	hmm_print(g_lasthmm);
    }
    exit(1);
}

inline PROB log_add(PROB x, PROB y) {
    PROB temp, negdiff;
    PROB result;
    if (x == LOGZERO) return (y);
    if (y == LOGZERO) return (x);
    if (y > x) {
	temp = x;
	x = y;
	y = temp;
    }
    negdiff = y - x;
    if (negdiff <= -61) { return x; }
#ifdef LOG_LUT
    result = log1plus_table_interp(negdiff);
#elif LOG_LIB
    result = log2(exp2(negdiff)+1);
#else
    result = log1plus_minimax(negdiff);
#endif /* LOG_LUT */
    return(result+x);
}

PROB trellis_backward(struct trellis *trellis, int *obs, int length, struct wfsa *fsm) {
    int i, sourcestate, targetstate, symbol;
    PROB target_prob;
    for (i = 0; i <= length + 1; i++)
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++)
	    TRELLIS_CELL(sourcestate,i)->bp = LOGZERO;
    
    /* Fill last and penultimate column */
    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
	TRELLIS_CELL(targetstate,length+1)->bp = 0;
	TRELLIS_CELL(targetstate,length)->bp = 0 + *FINALPROB(fsm, targetstate);
    }
    /* Fill rest */
    for (i = length-1; i >= 0 ; i--) {
	symbol = obs[i];
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++) {
	    TRELLIS_CELL(sourcestate,i)->bp = LOGZERO;
	    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
		target_prob = *TRANSITION(fsm, sourcestate, symbol, targetstate);
		if (target_prob <= SMRZERO_LOG) { continue; }
		TRELLIS_CELL(sourcestate,i)->bp = log_add(TRELLIS_CELL(sourcestate,i)->bp, TRELLIS_CELL(targetstate,i+1)->bp + target_prob);
	    }
	}
    }
    return(TRELLIS_CELL(0,0)->bp);
}

PROB trellis_backward_hmm(struct trellis *trellis, int *obs, int length, struct hmm *hmm) {
    int i, sourcestate, targetstate, symbol;
    PROB target_prob;
    for (i = 0; i <= length + 1; i++) /* Clear all */
	for (sourcestate = 0; sourcestate < hmm->num_states; sourcestate++)
	    TRELLIS_CELL_HMM(sourcestate, i)->bp = LOGZERO;
    
    /* Fill last and penultimate column */
    TRELLIS_CELL_HMM(hmm->num_states - 1, length + 1)->bp = 0;  /* Set prob 1.0 for final */
    for (sourcestate = 0; sourcestate < hmm->num_states - 1; sourcestate++) {
	TRELLIS_CELL_HMM(sourcestate, length)->bp = *HMM_TRANSITION_PROB(hmm, sourcestate, hmm->num_states - 1);
    }
    /* Fill rest */
    for (i = length-1; i >= 0 ; i--) {
	symbol = obs[i];
	for (sourcestate = 0; sourcestate < hmm->num_states - 1; sourcestate++) {
	    if (sourcestate == 0 && i != 0)
		continue;
	    TRELLIS_CELL_HMM(sourcestate, i)->bp = LOGZERO;
	    for (targetstate = 0; targetstate < hmm->num_states - 1; targetstate++) {
		target_prob = *HMM_TRANSITION_PROB(hmm, sourcestate, targetstate) + *HMM_EMISSION_PROB(hmm, targetstate, symbol);
		if (target_prob <= SMRZERO_LOG) { continue; }
		TRELLIS_CELL_HMM(sourcestate, i)->bp = log_add(TRELLIS_CELL_HMM(sourcestate, i)->bp, TRELLIS_CELL_HMM(targetstate, i+1)->bp + target_prob);
	    }
	}
    }
    return(TRELLIS_CELL_HMM(0,0)->bp);
}

PROB trellis_viterbi(struct trellis *trellis, int *obs, int length, struct wfsa *fsm) {
    int i, sourcestate, targetstate, symbol, final_state;
    PROB target_prob, final_prob;
    
    for (i = 0; i <= length + 1; i++)
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++)
	    TRELLIS_CELL(sourcestate,i)->fp = LOGZERO;
    
    /* Calculate first transition */
    TRELLIS_CELL(0,0)->fp = 0;
    for (i = 0; i < 1 && i < length; i++) {
	symbol = obs[i];
	for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
	    target_prob = *TRANSITION(fsm, 0, symbol, targetstate);
	    if (target_prob > SMRZERO_LOG) {
		TRELLIS_CELL(targetstate,1)->fp = target_prob;
		TRELLIS_CELL(targetstate,1)->backstate = 0;
	    }
	}
    }
    /* Calculate remaining transitions */
    for (i = 1; i < length; i++) {
	symbol = obs[i];
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++) {
	    if (TRELLIS_CELL(sourcestate,i)->fp == LOGZERO) { continue; }
	    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
		target_prob = *TRANSITION(fsm, sourcestate, symbol, targetstate);	       
		if (target_prob <= SMRZERO_LOG) { continue; }
		if (TRELLIS_CELL(targetstate,(i+1))->fp == LOGZERO) {
		    TRELLIS_CELL(targetstate,(i+1))->fp = TRELLIS_CELL(sourcestate,i)->fp + target_prob;
		    TRELLIS_CELL(targetstate,(i+1))->backstate = sourcestate;
		} else {		    
		    if (TRELLIS_CELL(targetstate,(i+1))->fp < TRELLIS_CELL(sourcestate,i)->fp + target_prob) {
			TRELLIS_CELL(targetstate,(i+1))->fp = TRELLIS_CELL(sourcestate,i)->fp + target_prob;
			TRELLIS_CELL(targetstate,(i+1))->backstate = sourcestate;
		    }
		}
	    }
	}
    }
    
    /* Calculate final state probabilities */
    i = length;
    final_state = -1;
    for (targetstate = 0, final_prob = SMRZERO_LOG; targetstate < fsm->num_states; targetstate++) {
	if (TRELLIS_CELL(targetstate,i)->fp == LOGZERO) {
	    TRELLIS_CELL(targetstate,(i+1))->backstate = -1;
	    continue;
	}
	if (*FINALPROB(fsm, targetstate) > SMRZERO_LOG && TRELLIS_CELL(targetstate,i)->fp > SMRZERO_LOG) { 
	    TRELLIS_CELL(targetstate,(i+1))->fp = TRELLIS_CELL(targetstate,i)->fp + *FINALPROB(fsm, targetstate);
	} else {
	    continue;
	}
	if (TRELLIS_CELL(targetstate,(i+1))->fp > final_prob) {
	    final_prob = TRELLIS_CELL(targetstate,(i+1))->fp;
	    final_state = targetstate;
	}
    }
    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
	if (targetstate != final_state) {
	    TRELLIS_CELL(targetstate,(i+1))->backstate = -1;
	} else {
	    TRELLIS_CELL(targetstate,(i+1))->backstate = targetstate;
	}
    }
    return(final_prob);
}

PROB trellis_forward_hmm(struct trellis *trellis, int *obs, int length, struct hmm *hmm) {
    int i, sourcestate, targetstate, end_state;
    PROB target_prob, final_prob;
    
    end_state = hmm->num_states - 1;
    for (i = 0; i <= length + 1; i++)
	for (sourcestate = 0; sourcestate < hmm->num_states; sourcestate++)
	    TRELLIS_CELL_HMM(sourcestate,i)->fp = LOGZERO;
    
    /* Calculate first transition */
    TRELLIS_CELL_HMM(0,0)->fp = 0;
    for (targetstate = 1; targetstate < hmm->num_states ; targetstate++) {
	if ((targetstate == end_state && length != 0) || (length == 0 && targetstate != end_state)) {
	    continue;
	}
	if (targetstate == end_state)
	    target_prob = *HMM_TRANSITION_PROB(hmm, 0, targetstate);
	else
	    target_prob = *HMM_TRANSITION_PROB(hmm, 0, targetstate) + *HMM_EMISSION_PROB(hmm, targetstate, obs[0]);
	if (target_prob > SMRZERO_LOG) {
	    TRELLIS_CELL_HMM(targetstate, 1)->fp = target_prob;
	}
    }

    /* Calculate remaining transitions */
    for (i = 1; i <= length; i++) {
	for (sourcestate = 1; sourcestate < end_state; sourcestate++) {
	    if (TRELLIS_CELL_HMM(sourcestate, i)->fp == LOGZERO) { continue; }
	    for (targetstate = 0; targetstate < hmm->num_states; targetstate++) {
		if (i != length && targetstate == end_state) {
		    continue;
		}
		if (targetstate == hmm->num_states - 1 || i == length)
		    target_prob = *HMM_TRANSITION_PROB(hmm, sourcestate, targetstate);
		else
		    target_prob = *HMM_TRANSITION_PROB(hmm, sourcestate, targetstate) + *HMM_EMISSION_PROB(hmm, targetstate, obs[i]);

		if (target_prob <= SMRZERO_LOG) { continue; }
		TRELLIS_CELL_HMM(targetstate,(i+1))->fp = log_add(TRELLIS_CELL_HMM(sourcestate, i)->fp + target_prob, TRELLIS_CELL_HMM(targetstate,(i+1))->fp);
	    }
	}
    }
    final_prob = TRELLIS_CELL_HMM(end_state, i)->fp;
    return(final_prob);
}

PROB trellis_viterbi_hmm(struct trellis *trellis, int *obs, int length, struct hmm *hmm) {
    int i, sourcestate, targetstate, symbol, end_state;
    PROB target_prob, final_prob;
    
    end_state = hmm->num_states - 1;
    for (i = 0; i <= length + 1; i++)
	for (sourcestate = 0; sourcestate < hmm->num_states; sourcestate++)
	    TRELLIS_CELL_HMM(sourcestate,i)->fp = LOGZERO;
    
    /* Calculate first transition */
    TRELLIS_CELL_HMM(0,0)->fp = 0;
    for (targetstate = 1; targetstate < hmm->num_states ; targetstate++) {
	if (targetstate == end_state && length != 0) {
	    continue;
	}
	target_prob = *HMM_TRANSITION_PROB(hmm, 0, targetstate);
	if (target_prob > SMRZERO_LOG) {
	    TRELLIS_CELL_HMM(targetstate,1)->fp = target_prob;
	    TRELLIS_CELL_HMM(targetstate,1)->backstate = 0;
	}
    }

    /* Calculate remaining transitions */
    for (i = 1; i <= length; i++) {
	symbol = obs[i-1];
	for (sourcestate = 1; sourcestate < end_state; sourcestate++) {
	    if (TRELLIS_CELL_HMM(sourcestate, i)->fp == LOGZERO) { continue; }
	    for (targetstate = 0; targetstate < hmm->num_states; targetstate++) {
		if (i != length && targetstate == end_state) {
		    continue;
		}
		target_prob = *HMM_EMISSION_PROB(hmm, sourcestate, symbol) + *HMM_TRANSITION_PROB(hmm, sourcestate, targetstate);
		if (target_prob <= SMRZERO_LOG) { continue; }
		if (TRELLIS_CELL_HMM(targetstate, i+1)->fp == LOGZERO) {
		    TRELLIS_CELL_HMM(targetstate, i+1)->fp = TRELLIS_CELL_HMM(sourcestate, i)->fp + target_prob;
		    TRELLIS_CELL_HMM(targetstate, i+1)->backstate = sourcestate;
		} else {
		    if (TRELLIS_CELL_HMM(targetstate, i+1)->fp < TRELLIS_CELL_HMM(sourcestate,i)->fp + target_prob) {
			TRELLIS_CELL_HMM(targetstate, i+1)->fp = TRELLIS_CELL_HMM(sourcestate,i)->fp + target_prob;
			TRELLIS_CELL_HMM(targetstate, i+1)->backstate = sourcestate;
		    }
		}
	    }
	}
    }
    final_prob = TRELLIS_CELL_HMM(end_state, i)->fp;
    return(final_prob);
}

PROB trellis_forward_fsm(struct trellis *trellis, int *obs, int length, struct wfsa *fsm) {
    int i, sourcestate, targetstate, symbol;
    PROB target_prob, final_prob;

    for (i = 0; i <= length + 1; i++)
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++)
	    TRELLIS_CELL(sourcestate,i)->fp = LOGZERO;
    
    /* Calculate first transition */
    TRELLIS_CELL(0,0)->fp = 0;
    for (i = 0; i < 1 && i < length; i++) {
	symbol = obs[i];
	for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
	    target_prob = *TRANSITION(fsm, 0, symbol, targetstate);
	    if (target_prob > SMRZERO_LOG) {
		TRELLIS_CELL(targetstate,1)->fp = target_prob;
		TRELLIS_CELL(targetstate,1)->backstate = 0;
	    }
	}
    }
    /* Calculate remaining transitions */
    for (i = 1; i < length; i++) {
	symbol = obs[i];
	for (sourcestate = 0; sourcestate < fsm->num_states; sourcestate++) {
	    if (TRELLIS_CELL(sourcestate,i)->fp == LOGZERO) { continue; }
	    for (targetstate = 0; targetstate < fsm->num_states; targetstate++) {
		target_prob = *TRANSITION(fsm, sourcestate, symbol, targetstate);
		if (target_prob <= SMRZERO_LOG) { continue; }
		TRELLIS_CELL(targetstate,(i+1))->fp = log_add(TRELLIS_CELL(sourcestate,i)->fp + target_prob, TRELLIS_CELL(targetstate,(i+1))->fp);
	    }
	}
    }
    
    /* Calculate final state probabilities */
    i = length;
    for (targetstate = 0, final_prob = SMRZERO_LOG; targetstate < fsm->num_states; targetstate++) {
	if (TRELLIS_CELL(targetstate,i)->fp == LOGZERO) { continue; }
	if (*FINALPROB(fsm, targetstate) > SMRZERO_LOG && TRELLIS_CELL(targetstate,i)->fp > SMRZERO_LOG) {
	    TRELLIS_CELL(targetstate,(i+1))->fp = TRELLIS_CELL(targetstate,i)->fp + *FINALPROB(fsm, targetstate);
	} else {
	    continue;
	}
	final_prob = log_add(final_prob, TRELLIS_CELL(targetstate,(i+1))->fp);
    }
    return(final_prob);
}

struct trellis *trellis_init(struct observations *o, int num_states) {
    int olenmax;
    struct trellis *trellis;
    for (olenmax = 0 ; o != NULL; o = o->next) {
	olenmax = olenmax < o->size ? o->size : olenmax;
    }
    trellis = calloc((olenmax + 2) * num_states, sizeof(struct trellis));
    return(trellis);
}

void trellis_print(struct trellis *trellis, struct wfsa *fsm, int obs_len) {
    int i, j;
    for (j = fsm->num_states-1; j >= 0 ; j--) {
	for (i = 0; i <= obs_len + 1 ; i++) {
	    printf("FP:%4.4f BP:%4.4f B:%i\t",  EXP(TRELLIS_CELL(j,i)->fp), EXP(TRELLIS_CELL(j,i)->bp), TRELLIS_CELL(j,i)->backstate);
	}
	printf("\n");
    }
}

void forward_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len) {
    int i, j, beststate;
    PROB bestprob;
    for (i = 0; i <= obs_len+1; i++) {
	if (i == obs_len)
	    continue;
	bestprob = SMRZERO_LOG;
	beststate = -1;
	for (j = 0; j < fsm->num_states; j++) {
	    if (TRELLIS_CELL(j,i)->fp == LOGZERO) {
		continue;
	    }
	    if (TRELLIS_CELL(j,i)->fp > bestprob) {
		bestprob = TRELLIS_CELL(j,i)->fp;
		beststate = j;
	    }
	}
	printf("%i", beststate);
	if (i < obs_len) printf(" ");
    }
    printf("\n");
}

void forward_print_path_hmm(struct trellis *trellis, struct hmm *hmm, int obs_len) {
    int i, j, beststate;
    PROB bestprob;
    for (i = 0; i <= obs_len + 1; i++) {
	if (i == obs_len + 1)
	    continue;
	bestprob = SMRZERO_LOG;
	beststate = -1;
	for (j = 0; j < hmm->num_states - 1; j++) {
	    if (TRELLIS_CELL_HMM(j, i)->fp == LOGZERO) {
		continue;
	    }
	    if (TRELLIS_CELL_HMM(j, i)->fp > bestprob) {
		bestprob = TRELLIS_CELL_HMM(j, i)->fp;
		beststate = j;
	    }
	}
	printf("%i", beststate);
	if (i < obs_len + 1) printf(" ");
    }
    printf("%i", hmm->num_states - 1);
    printf("\n");
}

void backward_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len) {
    int i, j, beststate;
    PROB bestprob;
    for (i = 0; i <= obs_len; i++) {
	bestprob = SMRZERO_LOG;
	beststate = -1;
	for (j = 0; j < fsm->num_states; j++) {
	    if (TRELLIS_CELL(j,i)->bp == LOGZERO) {
		continue;
	    }
	    if (TRELLIS_CELL(j,i)->bp > bestprob) {
		bestprob = TRELLIS_CELL(j,i)->fp;
		beststate = j;
	    }
	}
	printf("%i", beststate);
	if (i < obs_len) printf(" ");
    }
    printf("\n");
}

void backward_print_path_hmm(struct trellis *trellis, struct hmm *hmm, int obs_len) {
    int i, j, beststate;
    PROB bestprob;
    for (i = 0; i <= obs_len + 1; i++) {
	bestprob = SMRZERO_LOG;
	beststate = -1;
	for (j = 0; j < hmm->num_states; j++) {
	    if (TRELLIS_CELL_HMM(j,i)->bp == LOGZERO) {
		continue;
	    }
	    if (TRELLIS_CELL_HMM(j,i)->bp > bestprob) {
		bestprob = TRELLIS_CELL_HMM(j,i)->bp;
		beststate = j;
	    }
	}
	printf("%i", beststate);
	if (i <= obs_len) printf(" ");
    }
    printf("\n");
}

void viterbi_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len) {
    int i, laststate, *path;
    path = malloc(sizeof(int) * (obs_len+1));
    for (i = 0; i < fsm->num_states; i++) {
	if (TRELLIS_CELL(i, obs_len+1)->backstate != -1) {
	    laststate = i;
	}
    }
    *(path+obs_len) = laststate;
    for (i = obs_len; i > 0; i--) {
	*(path+i-1) = TRELLIS_CELL(laststate, i)->backstate;
	laststate = TRELLIS_CELL(laststate, i)->backstate;
    }
    for (i = 0 ; i <= obs_len; i++) {
	printf("%i", path[i]);
	if (i < obs_len) {
	    printf(" ");
	}
    }
    printf("\n");
    free(path);
}

void viterbi_print_path_hmm(struct trellis *trellis, struct hmm *hmm, int obs_len) {
    int i, laststate, *path;
    path = malloc(sizeof(int) * (obs_len + 2));
    laststate = hmm->num_states - 1;
    *(path+obs_len+1) = laststate;
    for (i = obs_len + 1; i > 0; i--) {
	*(path+i-1) = TRELLIS_CELL_HMM(laststate, i)->backstate;
	laststate = TRELLIS_CELL_HMM(laststate, i)->backstate;
    }
    for (i = 0 ; i <= obs_len + 1; i++) {
	printf("%i", path[i]);
	if (i < obs_len + 1) {
	    printf(" ");
	}
    }
    printf("\n");
    free(path);
}

void viterbi(struct wfsa *fsm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB viterbi_prob;
    trellis = trellis_init(o, fsm->num_states);
    for (obs = o; obs != NULL; obs = obs->next) {
	viterbi_prob = trellis_viterbi(trellis, obs->data, obs->size, fsm);
	if (algorithm == DECODE_VITERBI_PROB)
	    printf("%.17g\t", output_convert(viterbi_prob));
	if (algorithm == LIKELIHOOD_VITERBI)
	    printf("%.17g\n", output_convert(viterbi_prob));
	if (algorithm == DECODE_VITERBI_PROB || algorithm == DECODE_VITERBI) {
	    if (viterbi_prob > SMRZERO_LOG) {
		viterbi_print_path(trellis, fsm, obs->size);
	    } else {
		printf("\n");
	    }
	}
    }
    free(trellis);
}

void viterbi_hmm(struct hmm *hmm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB viterbi_prob;
    trellis = trellis_init(o, hmm->num_states);
    for (obs = o; obs != NULL; obs = obs->next) {
	viterbi_prob = trellis_viterbi_hmm(trellis, obs->data, obs->size, hmm);
	if (algorithm == DECODE_VITERBI_PROB)
	    printf("%.17g\t", output_convert(viterbi_prob));
	if (algorithm == LIKELIHOOD_VITERBI)
	    printf("%.17g\n", output_convert(viterbi_prob));
	if (algorithm == DECODE_VITERBI_PROB || algorithm == DECODE_VITERBI) {
	    if (viterbi_prob > SMRZERO_LOG) {
		viterbi_print_path_hmm(trellis, hmm, obs->size);
	    } else {
		printf("\n");
	    }
	}
    }
    free(trellis);
}

PROB hmm_sum_transition_prob(struct hmm *hmm, int state) {
    /* Get sum of probabilities for transition in a state (in reals) */
    PROB sum;
    int i;
    sum = 0;
    for (i = 0; i < hmm->num_states; i++) {
	if (EXP(*HMM_TRANSITION_PROB(hmm, state, i) >= SMRZERO_LOG))
	    sum += EXP(*HMM_TRANSITION_PROB(hmm, state, i));	
    }
    return(sum);
}

PROB hmm_sum_emission_prob(struct hmm *hmm, int state) {
    /* Get sum of probabilities for emission in a state (in reals) */
    PROB sum;
    int i;
    sum = 0;
    for (i = 0; i < hmm->alphabet_size; i++) {
	if (EXP(*HMM_EMISSION_PROB(hmm, state, i) >= SMRZERO_LOG))
	    sum += EXP(*HMM_EMISSION_PROB(hmm, state, i));
    }
    return(sum);
}

PROB wfsa_sum_prob(struct wfsa *fsm, int state) {
    PROB sum;
    int i, s;
    /* Get sum of probabilities for a state (in reals) */
    sum = *FINALPROB(fsm,state) >= SMRZERO_LOG ? EXP(*FINALPROB(fsm,state)) : 0;
    for (i = 0; i < fsm->num_states; i++) {
	for (s = 0 ; s < fsm->alphabet_size; s++) {
	    if (*TRANSITION(fsm, state, s, i) >= SMRZERO_LOG) {
		sum += EXP(*TRANSITION(fsm, state, s, i));
	    }
	}
    }
    return(sum);
}

int wfsa_random_transition(struct wfsa *fsm, int state, int *symbol, PROB *prob) {
    /* Choose a random arc from state "state"         */
    /* Return target state, and put symbol in *symbol */
    /* If stop: symbol = -1                           */
    PROB thissum, r;
    int i, s;
    r = (PROB) random() / RAND_MAX;
    r = r * wfsa_sum_prob(fsm, state);
    thissum = 0;
    
    if (*FINALPROB(fsm,state) >= SMRZERO_LOG) {
	thissum += EXP(*FINALPROB(fsm, state));
	if (thissum >= r) {
	    *symbol = -1;
	    *prob = *FINALPROB(fsm,state);
	    return state;
	}
    }
    for (i = 0; i < fsm->num_states; i++) {
	for (s = 0 ; s < fsm->alphabet_size; s++) {
	    if ((*prob = *TRANSITION(fsm, state, s, i)) >= SMRZERO_LOG) {		
		thissum += EXP(*prob);
		if (thissum >= r) {
		    *symbol = s;
		    return(i);
		}
	    }
	}
    }
    perror("Inconsistent probabilities in FSM");
    exit(1);
}

int hmm_random_transition(struct hmm *hmm, int state, PROB *prob) {
    /* Choose a random arc from state "state"           */
    /* Return target state, and put probability in prob */
    PROB thissum, r;
    int i;
    r = (PROB) random() / RAND_MAX;
    r = r * hmm_sum_transition_prob(hmm, state);
    thissum = 0;
    
    for (i = 0; i < hmm->num_states; i++) {
	if ((*prob = *HMM_TRANSITION_PROB(hmm, state, i)) >= SMRZERO_LOG) {		
	    thissum += EXP(*prob);
	    if (thissum >= r) {
		return(i);
	    }
	}
    }
    perror("Inconsistent probabilities in HMM");
    exit(1);
}

int hmm_random_emission(struct hmm *hmm, int state, PROB *prob) {
    /* Choose a random symbol to emit from state  */
    /* Return symbol, and put probability in prob */
    PROB thissum, r;
    int i;
    r = (PROB) random() / RAND_MAX;
    r = r * hmm_sum_emission_prob(hmm, state);
    thissum = 0;
    
    for (i = 0; i < hmm->alphabet_size; i++) {
	if ((*prob = *HMM_EMISSION_PROB(hmm, state, i)) >= SMRZERO_LOG) {		
	    thissum += EXP(*prob);
	    if (thissum >= r) {
		return(i);
	    }
	}
    }
    perror("Inconsistent probabilities in HMM");
    exit(1);
}

void generate_words(struct wfsa *fsm, int numwords) {
    int i, j, k, state, symbol, *output = NULL, *stateseq = NULL;
    PROB prob, totalprob;

    output = malloc(sizeof(int) * g_gen_max_length);
    stateseq = malloc(sizeof(int) * g_gen_max_length);
    
    for (i = 0; i < numwords; i++) {
	stateseq[0] = 0;
	totalprob = LOGZERO;
	for (state = 0, j = 0; j < g_gen_max_length ; j++) {
	    state = wfsa_random_transition(fsm, state, &symbol, &prob);
	    totalprob = (totalprob == LOGZERO) ? prob : totalprob + prob;
	    stateseq[j+1] = state;
	    if (symbol >= 0) {
		output[j] = symbol;
	    } else {
		break;
	    }
        }
	printf("%.17g\t", output_convert(totalprob));
	if (j < g_gen_max_length) {
	    for (k = 0; k < j; k++) {
		printf("%i", output[k]);
		if (k < j-1) {
		    printf(" ");
		}
	    }
	    printf("\t");
	    for (k = 0; k <= j; k++) {
		printf("%i", stateseq[k]);
		if (k < j) {
		    printf(" ");
		}
	    }
	    printf("\n");

	} else {
	    i--;
	}
    }
    free(output);
    free(stateseq);
}

void generate_words_hmm(struct hmm *hmm, int numwords) {
    int i, j, k, state, symbol, *output = NULL, *stateseq = NULL;
    PROB tprob, eprob, totalprob;

    output = malloc(sizeof(int) * g_gen_max_length);
    stateseq = malloc(sizeof(int) * g_gen_max_length);
    

    for (i = 0; i < numwords; i++) {	
	stateseq[0] = 0;
	totalprob = LOGZERO;
	for (state = 0, j = 0; j < g_gen_max_length ; j++) {
	    state = hmm_random_transition(hmm, state, &tprob);
	    stateseq[j+1] = state;
	    if (state != hmm->num_states - 1) {
		symbol = hmm_random_emission(hmm, state, &eprob);
		output[j] = symbol;
	    } else {
		totalprob = (totalprob == LOGZERO) ? tprob : totalprob + tprob;
		break;
	    }
	    totalprob = (totalprob == LOGZERO) ? tprob + eprob : totalprob + tprob + eprob;
        }
	printf("%.17g\t", output_convert(totalprob));
	if (j < g_gen_max_length) {
	    for (k = 0; k < j; k++) {
		printf("%i", output[k]);
		if (k < j-1) {
		    printf(" ");
		}
	    }
	    printf("\t");
	    for (k = 0; k <= j + 1; k++) {
		printf("%i", stateseq[k]);
		if (k <= j) {
		    printf(" ");
		}
	    }
	    printf("\n");
	} else { 
	    i--; /* Try again, the string is too long... */
	}
    }
    free(output);
    free(stateseq);
}

PROB loglikelihood_all_observations_fsm(struct wfsa *fsm, struct observations *o) {
    struct observations *obs;
    struct trellis *trellis;
    PROB forward_prob;
    PROB ll;
    trellis = trellis_init(o, fsm->num_states);
    for (obs = o, ll = LOGZERO; obs != NULL; obs = obs->next) {
	forward_prob = obs->occurrences * trellis_forward_fsm(trellis, obs->data, obs->size, fsm);
	ll = ll == LOGZERO ? forward_prob : ll + forward_prob;
    }
    free(trellis);
    return(ll);
}

PROB loglikelihood_all_observations_hmm(struct hmm *hmm, struct observations *o) {
    struct observations *obs;
    struct trellis *trellis;
    PROB forward_prob;
    PROB ll;
    trellis = trellis_init(o, hmm->num_states);
    for (obs = o, ll = LOGZERO; obs != NULL; obs = obs->next) {
	forward_prob = obs->occurrences * trellis_forward_hmm(trellis, obs->data, obs->size, hmm);
	ll = ll == LOGZERO ? forward_prob : ll + forward_prob;
    }
    free(trellis);
    return(ll);
}

void forward_hmm(struct hmm *hmm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB forward_prob;
    trellis = trellis_init(o, hmm->num_states);
    for (obs = o; obs != NULL; obs = obs->next) {
	forward_prob = trellis_forward_hmm(trellis, obs->data, obs->size, hmm);
	if (algorithm == DECODE_FORWARD_PROB)
	    printf("%.17g\t", output_convert(forward_prob));
	if (algorithm == LIKELIHOOD_FORWARD)
	    printf("%.17g\n", output_convert(forward_prob));
	if (algorithm == DECODE_FORWARD_PROB || algorithm == DECODE_FORWARD) {
	    if (forward_prob > SMRZERO_LOG) {
		forward_print_path_hmm(trellis, hmm, obs->size);
	    } else {
		printf("\n");
	    }
	}
    }
    free(trellis);
}

void forward_fsm(struct wfsa *fsm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB forward_prob;
    trellis = trellis_init(o, fsm->num_states);
    for (obs = o; obs != NULL; obs = obs->next) {
	forward_prob = trellis_forward_fsm(trellis, obs->data, obs->size, fsm);
	if (algorithm == DECODE_FORWARD_PROB)
	    printf("%.17g\t", output_convert(forward_prob));
	if (algorithm == LIKELIHOOD_FORWARD)
	    printf("%.17g\n", output_convert(forward_prob));
	if (algorithm == DECODE_FORWARD_PROB || algorithm == DECODE_FORWARD) {
	    if (forward_prob > SMRZERO_LOG) {
		forward_print_path(trellis, fsm, obs->size);
	    } else {
		printf("\n");
	    }
	}
    }
    free(trellis);
}

void backward_fsm(struct wfsa *fsm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB backward_prob;
    trellis = trellis_init(o, fsm->num_states);
    for (obs = o; obs != NULL; obs = obs->next) {
	backward_prob = trellis_backward(trellis, obs->data, obs->size, fsm);
	if (algorithm == DECODE_BACKWARD_PROB)
	    printf("%.17g\t", output_convert(backward_prob));
	if (algorithm == LIKELIHOOD_BACKWARD)
	    printf("%.17g\n", output_convert(backward_prob));
	if (algorithm == DECODE_BACKWARD_PROB || algorithm == DECODE_BACKWARD) {
	    if (backward_prob > SMRZERO_LOG) {
		backward_print_path(trellis, fsm, obs->size);
	    } else {
		printf("\n");
	    }
	}
    }
    free(trellis);
}

void backward_hmm(struct hmm *hmm, struct observations *o, int algorithm) {
    struct observations *obs;
    struct trellis *trellis;
    PROB backward_prob;
    trellis = trellis_init(o, hmm->num_states);
    for (obs = o; obs != NULL; obs = obs->next) {
	backward_prob = trellis_backward_hmm(trellis, obs->data, obs->size, hmm);
	if (algorithm == DECODE_BACKWARD_PROB)
	    printf("%.17g\t", output_convert(backward_prob));
	if (algorithm == LIKELIHOOD_BACKWARD)
	    printf("%.17g\n", output_convert(backward_prob));
	if (algorithm == DECODE_BACKWARD_PROB || algorithm == DECODE_BACKWARD) {
	    if (backward_prob > SMRZERO_LOG) {
		backward_print_path_hmm(trellis, hmm, obs->size);
	    } else {
		printf("\n");
	    }
	}
    }
    free(trellis);
}

PROB train_viterbi(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta) {
    struct observations *obs;
    struct trellis *trellis;
    int i,j,k,iter, source, target, laststate, symbol, occurrences, *fsm_vit_counts, *fsm_vit_totalcounts, *fsm_vit_finalcounts;
    PROB viterbi_prob, loglikelihood, prevloglikelihood, newprob;
    trellis = trellis_init(o, fsm->num_states);
    fsm_vit_counts = malloc(sizeof(int) * fsm->num_states * fsm->num_states * fsm->alphabet_size);
    fsm_vit_totalcounts = malloc(sizeof(int) * fsm->num_states);
    fsm_vit_finalcounts = malloc(sizeof(int) * fsm->num_states);
    
    prevloglikelihood = 0;
    for (iter = 0 ; iter < maxiterations; iter++) {
        /* Clear counts */
        for (i = 0; i < fsm->num_states; i++) {
            fsm_vit_totalcounts[i] = fsm_vit_finalcounts[i] = 0;
        }
        for (i = 0; i < fsm->num_states * fsm->num_states * fsm->alphabet_size; i++) {
            fsm_vit_counts[i] = 0;
        }
        loglikelihood = 0;
        for (obs = o; obs != NULL; obs = obs->next) {
	    occurrences = obs->occurrences;
            viterbi_prob = trellis_viterbi(trellis, obs->data, obs->size, fsm);
            if (viterbi_prob <= SMRZERO_LOG) {
                continue;
            } else {
                loglikelihood += viterbi_prob * occurrences;
                /* Update final counts */
                for (i = 0, laststate = -1; i < fsm->num_states; i++) {
                    if (TRELLIS_CELL(i,(obs->size+1))->backstate != -1) {
                        laststate = i;
                        break;
                    }
                }
                if (laststate == -1) {
                    printf("Could not find last state\n");
                    continue;
                }
                fsm_vit_finalcounts[laststate] += occurrences;
                fsm_vit_totalcounts[laststate] += occurrences;
                /* Update arc counts */
                for (i = obs->size; i > 0; i--) {
                    target = laststate;
                    laststate = TRELLIS_CELL(laststate,i)->backstate;
                    source = laststate;
                    symbol = *((obs->data)+i-1);
                    *FSM_COUNTS(fsm_vit_counts, source, symbol, target) += occurrences;
                    fsm_vit_totalcounts[source]+= occurrences;
                }
            }
        }
	fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g\n", iter+1, loglikelihood, ABS(prevloglikelihood - loglikelihood));
	if (ABS(prevloglikelihood - loglikelihood) < maxdelta)  {
	    break;
	}
        /* Update WFSA */
        for (i=0; i < fsm->num_states; i++) {
            for (j = 0; j < fsm->alphabet_size; j++) {
                for (k = 0; k < fsm->num_states; k++) {
                    if (fsm_vit_totalcounts[i] > 0) {
                        newprob = LOG(*FSM_COUNTS(fsm_vit_counts,i,j,k) + g_viterbi_pseudocount) - LOG((PROB)fsm_vit_totalcounts[i] + g_viterbi_pseudocount * (PROB)fsm->alphabet_size * (PROB)fsm->num_states + g_viterbi_pseudocount);
                    } else {
                        newprob = SMRZERO_LOG;
                    }
                    *(TRANSITION(fsm,i,j,k)) = newprob;
                }
            }
        }
        for (i=0; i < fsm->num_states; i++) {
            if (fsm_vit_totalcounts[i] > 0) {
                newprob = LOG(fsm_vit_finalcounts[i] + g_viterbi_pseudocount) - LOG((PROB)fsm_vit_totalcounts[i] + g_viterbi_pseudocount * (PROB)fsm->alphabet_size * (PROB)fsm->num_states + g_viterbi_pseudocount);
            } else {
                newprob = SMRZERO_LOG;
            }
            *(fsm->final_table + i) = newprob;
        }
        prevloglikelihood = loglikelihood;
    }
    free(fsm_vit_counts);
    free(fsm_vit_totalcounts);
    free(fsm_vit_finalcounts);
    return(loglikelihood);
}

PROB train_viterbi_hmm(struct hmm *hmm, struct observations *o, int maxiterations, PROB maxdelta) {
    struct observations *obs;
    struct trellis *trellis;
    int i, j, iter, source, target, newsource, symbol, occurrences, *hmm_vit_counts_trans, *hmm_vit_counts_emit, *hmm_vit_totalcounts_trans, *hmm_vit_totalcounts_emit;
    PROB viterbi_prob, loglikelihood, prevloglikelihood, newprob;

    trellis = trellis_init(o, hmm->num_states);    

    hmm_vit_counts_trans = malloc(sizeof(int) * hmm->num_states * hmm->num_states);
    hmm_vit_counts_emit = malloc(sizeof(int) * hmm->num_states * hmm->alphabet_size);
    hmm_vit_totalcounts_trans = malloc(sizeof(int) * hmm->num_states);
    hmm_vit_totalcounts_emit = malloc(sizeof(int) * hmm->num_states);
    
    prevloglikelihood = 0;
    for (iter = 0 ; iter < maxiterations; iter++) {
        /* Clear counts */
        for (i = 0; i < hmm->num_states; i++) {
            hmm_vit_totalcounts_trans[i] = hmm_vit_totalcounts_emit[i] = 0;
        }
        for (i = 0; i < hmm->num_states * hmm->num_states; i++)
	    hmm_vit_counts_trans[i] = 0;
        for (i = 0; i < hmm->num_states * hmm->alphabet_size; i++)
	    hmm_vit_counts_emit[i] = 0;
       
        loglikelihood = 0;
        for (obs = o; obs != NULL; obs = obs->next) {
	    occurrences = obs->occurrences;
            viterbi_prob = trellis_viterbi_hmm(trellis, obs->data, obs->size, hmm);
            if (viterbi_prob <= SMRZERO_LOG) {
                continue;
            } else {
                loglikelihood += viterbi_prob * occurrences;
                /* Update trans count to final state */
		target = hmm->num_states - 1;
		source = TRELLIS_CELL_HMM(target, obs->size + 1)->backstate;		
		hmm_vit_totalcounts_trans[source] += occurrences;
		*HMM_TRANSITION_COUNTS(hmm_vit_counts_trans, source, target) += occurrences;
                /* Update counts following backpointers until zero */
                for (i = obs->size; i > 0; i--) {
		    symbol = *((obs->data)+i-1);
		    hmm_vit_totalcounts_emit[source] += occurrences;
		    *HMM_EMISSION_COUNTS(hmm_vit_counts_emit, source, symbol) += occurrences;
		    newsource = TRELLIS_CELL_HMM(source, i)->backstate;
		    target = source;
		    source = newsource;
		    hmm_vit_totalcounts_trans[source] += occurrences;
		    *HMM_TRANSITION_COUNTS(hmm_vit_counts_trans, source, target) += occurrences;
                }
            }
        }
	fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g\n", iter+1, loglikelihood, ABS(prevloglikelihood - loglikelihood));
	if (ABS(prevloglikelihood - loglikelihood) < maxdelta)  {
	    break;
	}
        /* Update HMM */
	for (i = 0; i < hmm->num_states; i++) {
	    for (j = 0; j < hmm->num_states; j++) {
		/* Transitions */
		if (i == hmm->num_states - 1) {
		    *HMM_TRANSITION_PROB(hmm, i, j) = SMRZERO_LOG;
		    continue;
		}
		if ((PROB) hmm_vit_totalcounts_trans[i] + g_viterbi_pseudocount > 0) {
		    newprob = LOG((PROB) *HMM_TRANSITION_COUNTS(hmm_vit_counts_trans, i, j) + (PROB)g_viterbi_pseudocount) - LOG((PROB) hmm_vit_totalcounts_trans[i] + (PROB)hmm->num_states * (PROB)g_viterbi_pseudocount);
		} else {
		    newprob = SMRZERO_LOG;
		}
		*HMM_TRANSITION_PROB(hmm, i, j) = newprob;
	    }
	    for (j = 0; j < hmm->alphabet_size; j++) {
		/* Emissions */
		if (i == hmm->num_states - 1 || i == 0) {
		    *HMM_EMISSION_PROB(hmm, i, j) = SMRZERO_LOG;
		    continue;
		}
		if ((PROB) hmm_vit_totalcounts_emit[i] + g_viterbi_pseudocount_emit > 0.0) {
		    newprob = LOG((PROB)*HMM_EMISSION_COUNTS(hmm_vit_counts_emit, i, j) + (PROB) g_viterbi_pseudocount_emit) - LOG( (PROB)hmm_vit_totalcounts_emit[i] + (PROB)hmm->alphabet_size * (PROB)g_viterbi_pseudocount_emit);
		} else {
		    newprob = SMRZERO_LOG;
		}
		*HMM_EMISSION_PROB(hmm, i, j) = newprob;
	    }
	}
        prevloglikelihood = loglikelihood;
    }
    free(hmm_vit_counts_trans);
    free(hmm_vit_counts_emit);
    free(hmm_vit_totalcounts_trans);
    free(hmm_vit_totalcounts_emit);
    return(loglikelihood);
}

inline void spinlock_lock(_Bool *ptr) {
    if (g_num_threads > 1)
    	while (__sync_lock_test_and_set(ptr, 1)) { }
}

inline void spinlock_unlock(_Bool *ptr) {
    if (g_num_threads > 1)
     	__sync_lock_release(ptr);
}

void *trellis_fill_bw(void *threadargs) {
    pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
    struct trellis *trellis;
    struct observations **obsarray, *obs;
    struct wfsa *fsm;
    PROB backward_prob, forward_prob, thisxi, beta;
    int i, t, symbol, source, target, minobs, maxobs, occurrences;

    trellis = ((struct thread_args *)threadargs)->trellis;
    obsarray = ((struct thread_args *)threadargs)->obsarray;
    minobs = ((struct thread_args *)threadargs)->minobs;
    maxobs = ((struct thread_args *)threadargs)->maxobs;
    fsm = (struct wfsa *)((struct thread_args *)threadargs)->fsmhmm;
    beta = ((struct thread_args *)threadargs)->beta;
    
    for (i = minobs; i <= maxobs; i++) {
	obs = *(obsarray+i);
	occurrences = obs->occurrences;
	/* E-step */
	backward_prob = trellis_backward(trellis, obs->data, obs->size, fsm);
	forward_prob = trellis_forward_fsm(trellis, obs->data, obs->size, fsm);
	pthread_mutex_lock(&mutex1);
	g_loglikelihood += backward_prob * occurrences;
	pthread_mutex_unlock(&mutex1);
	/* Traverse trellis and add */
	for (t = 0; t < obs->size; t++) {
	    symbol = obs->data[t];
	    for (source = 0; source < fsm->num_states; source++) {
		if (TRELLIS_CELL(source,t)->fp == LOGZERO) { continue; }
		for (target = 0; target < fsm->num_states; target++) {
		    if (TRELLIS_CELL(target,t+1)->bp == LOGZERO) { continue; }
		    if (*TRANSITION(fsm,source,symbol,target) <= SMRZERO_LOG) { continue; }
		    thisxi = TRELLIS_CELL(source,t)->fp + *TRANSITION(fsm,source,symbol,target) + TRELLIS_CELL(target,t+1)->bp;
		    thisxi = thisxi - backward_prob;
		    thisxi = g_train_da_bw == 0 ? thisxi : thisxi * beta;
		    thisxi += LOG(occurrences);
		    spinlock_lock(&fsm_counts_spin[source]);
		    *FSM_COUNTS(fsm_counts,source,symbol,target) = log_add(*FSM_COUNTS(fsm_counts,source,symbol,target), thisxi);
		    spinlock_unlock(&fsm_counts_spin[source]);
		}
	    }
	}
	/* Final states */
	for (source = 0; source < fsm->num_states; source++) {
	    target = source;
	    if (TRELLIS_CELL(source,t)->fp == LOGZERO)   { continue; }
	    if (TRELLIS_CELL(target,t+1)->bp == LOGZERO) { continue; }
	    thisxi = TRELLIS_CELL(source,t)->fp + *FINALPROB(fsm, source);
	    thisxi = thisxi - backward_prob ;
	    thisxi = g_train_da_bw == 0 ? thisxi : thisxi * beta;
	    thisxi += LOG(occurrences);

	    spinlock_lock(&fsm_counts_spin[source]);
	    fsm_finalcounts[source] = log_add(fsm_finalcounts[source], thisxi);
	    spinlock_unlock(&fsm_counts_spin[source]);
	}
    }
    return(NULL);
}

void *trellis_fill_bw_hmm(void *threadargs) {
    pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
    struct trellis *trellis;
    struct observations **obsarray, *obs;
    struct hmm *hmm;
    PROB backward_prob, forward_prob, thisxi, beta;
    int i, t, symbol, source, target, minobs, maxobs, occurrences;

    trellis = ((struct thread_args *)threadargs)->trellis;
    obsarray = ((struct thread_args *)threadargs)->obsarray;
    minobs = ((struct thread_args *)threadargs)->minobs;
    maxobs = ((struct thread_args *)threadargs)->maxobs;
    hmm = (struct hmm *)((struct thread_args *)threadargs)->fsmhmm;
    beta = ((struct thread_args *)threadargs)->beta;
    
    for (i = minobs; i <= maxobs; i++) {
	obs = *(obsarray+i);
	occurrences = obs->occurrences;
	/* E-step */
	backward_prob = trellis_backward_hmm(trellis, obs->data, obs->size, hmm);
	forward_prob = trellis_forward_hmm(trellis, obs->data, obs->size, hmm);
	pthread_mutex_lock(&mutex1);
	g_loglikelihood += backward_prob * occurrences;
	pthread_mutex_unlock(&mutex1);
	/* Traverse trellis and add */
	for (t = 0; t <= obs->size; t++) {
	    for (source = 0; source < hmm->num_states - 1; source++) {
		if (TRELLIS_CELL_HMM(source,t)->fp == LOGZERO) { continue; }
		/* Emission */
		if (source > 0 && t > 0) {
		    symbol = obs->data[t-1];
		    thisxi = TRELLIS_CELL_HMM(source, t)->fp + TRELLIS_CELL_HMM(source, t)->bp;
		    thisxi -= backward_prob;
		    thisxi += LOG(occurrences);
		    spinlock_lock(&hmm_counts_spin[source]);
		    *HMM_EMISSION_COUNTS(hmm_counts_emit, source, symbol) = log_add(*HMM_EMISSION_COUNTS(hmm_counts_emit, source, symbol), thisxi);
		    spinlock_unlock(&hmm_counts_spin[source]);
		}
		for (target = 1; target < hmm->num_states; target++) {
		    if (TRELLIS_CELL_HMM(target, t+1)->bp == LOGZERO) { continue; }
		    if (*HMM_TRANSITION_PROB(hmm, source, target) <= SMRZERO_LOG) { continue; }
		    if (t == obs->size) {
			thisxi = TRELLIS_CELL_HMM(source, t)->fp + *HMM_TRANSITION_PROB(hmm, source, target) + TRELLIS_CELL_HMM(target, t+1)->bp;
		    } else {
			thisxi = TRELLIS_CELL_HMM(source, t)->fp + *HMM_TRANSITION_PROB(hmm, source, target) + *HMM_EMISSION_PROB(hmm, target, obs->data[t]) + TRELLIS_CELL_HMM(target, t+1)->bp;
		    }
		    thisxi -= backward_prob;
		    thisxi = g_train_da_bw == 0 ? thisxi : thisxi * beta;
		    thisxi += LOG(occurrences);
		    spinlock_lock(&hmm_counts_spin[source]);
		    *HMM_TRANSITION_COUNTS(hmm_counts_trans, source, target) = log_add(*HMM_TRANSITION_COUNTS(hmm_counts_trans, source, target), thisxi);
		    spinlock_unlock(&hmm_counts_spin[source]);
		}
	    }
	}
    }
    return(NULL);
}

PROB train_baum_welch_hmm(struct hmm *hmm, struct observations *o, int maxiterations, PROB maxdelta, int vb) {
    struct trellis *trellis, *trellisarray[32];
    struct thread_args *threadargs[32];
    struct observations **obsarray;
    int i, source, target, symbol, iter, numobs, obsperthread;
    PROB newprob, prevloglikelihood, da_beta = 1.0;
    pthread_t threadids[32];
    
    if (g_train_da_bw) { da_beta = g_betamin; }
    obsarray = observations_to_array(o, &numobs);
    
    /* Each thread gets its own trellis */
    for (i = 0; i < g_num_threads; i++) {
	trellis = trellis_init(o, hmm->num_states);
	trellisarray[i] = trellis;
	threadargs[i] = malloc(sizeof(struct thread_args));
    }

    /* These are accessed by all threads through a spinlock */
    hmm_counts_trans = malloc(hmm->num_states * hmm->num_states * sizeof(PROB));
    hmm_counts_emit = malloc(hmm->num_states * hmm->alphabet_size * sizeof(PROB));
    hmm_totalcounts_trans = malloc(hmm->num_states * sizeof(PROB));
    hmm_totalcounts_emit = malloc(hmm->num_states * sizeof(PROB));
    hmm_counts_spin = calloc(hmm->num_states, sizeof(_Bool));
    
    prevloglikelihood = 0;

    obsperthread = (int) ((double) numobs/(double) g_num_threads);
    /* Helper thread parameters (divide observations equally among all) */
    for (i = 1; i < g_num_threads; i++) {
	threadargs[i]->minobs = (i-1) * obsperthread;
	threadargs[i]->maxobs = i * obsperthread - 1;
	threadargs[i]->trellis = trellisarray[i];
	threadargs[i]->obsarray = obsarray;
	threadargs[i]->fsmhmm = hmm;
    }
    /* Main thread parameters */
    threadargs[0]->minobs = (i-1) * obsperthread;
    threadargs[0]->maxobs = numobs - 1;
    threadargs[0]->trellis = trellisarray[0];
    threadargs[0]->obsarray = obsarray;
    threadargs[0]->fsmhmm = hmm;
   
    for (iter = 0 ; iter < maxiterations ; iter++) {
	g_loglikelihood = 0;
	for (i = 0; i < hmm->num_states * hmm->num_states; i++) { hmm_counts_trans[i] = LOGZERO; }
	for (i = 0; i < hmm->num_states * hmm->alphabet_size ; i++) { hmm_counts_emit[i] = LOGZERO; }
	for (i = 1; i < g_num_threads; i++) {
	    /* Launch threads */
	    threadargs[i]->beta = da_beta;
	    pthread_create(&threadids[i], NULL, &trellis_fill_bw_hmm, threadargs[i]);
	}

	/* Run main thread Baum-Welch */
	threadargs[0]->beta = da_beta;
	trellis_fill_bw_hmm(threadargs[0]);
	/* Wait for all to finish */
	for (i = 1; i < g_num_threads; i++) {
	    pthread_join(threadids[i],NULL);
	}

	if (!g_train_da_bw)
	    fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g\n", iter+1, g_loglikelihood, ABS(prevloglikelihood - g_loglikelihood));
	else 
	    fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g beta: %.17g\n", iter+1, g_loglikelihood, ABS(prevloglikelihood - g_loglikelihood), da_beta);
	    
	if (ABS(prevloglikelihood - g_loglikelihood) < maxdelta)  {
	    if (g_train_da_bw == 1 && da_beta < g_betamax) {
		da_beta *= g_alpha;
		if (da_beta > g_betamax) {
		    da_beta = g_betamax;
		}
	    } else {
		break;
	    }
	}

	/* Modify HMM (M-step) */
	signal(SIGINT, SIG_IGN); /* Disable interrupts to prevent corrupted HMM in case of SIGINT while updating */

	/* Clear and sum counts */
	for (i = 0; i < hmm->num_states ; i++) { 
	    hmm_totalcounts_trans[i] =  hmm_totalcounts_emit[i] = LOGZERO;
	}
	/* Is totalcounts for emit always the same as totalcounts for transition (except at state 0)? */
	for (source = 0; source < hmm->num_states - 1; source++) {
	    for (symbol = 0; symbol < hmm->alphabet_size; symbol++) {		
		hmm_totalcounts_emit[source] = log_add(*HMM_EMISSION_COUNTS(hmm_counts_emit, source, symbol), hmm_totalcounts_emit[source]);
	    }
	    for (target = 1; target < hmm->num_states; target++) {
		hmm_totalcounts_trans[source] = log_add(*HMM_TRANSITION_COUNTS(hmm_counts_trans, source, target), hmm_totalcounts_trans[source]);
	    }
	}

	/* Re-estimate emissions */
	for (source = 0; source < hmm->num_states; source++) {
	    for (symbol = 0; symbol < hmm->alphabet_size; symbol++) {
		newprob = *HMM_EMISSION_COUNTS(hmm_counts_emit, source, symbol);
		if (newprob == LOGZERO) { newprob = SMRZERO_LOG; }
		if (vb)
		    *HMM_EMISSION_PROB(hmm, source, symbol) = (digamma(EXP(newprob) + g_gibbs_beta) - digamma(EXP(hmm_totalcounts_emit[source]) + hmm->alphabet_size * g_gibbs_beta_emission))/M_LN2;
		else
		    *HMM_EMISSION_PROB(hmm, source, symbol) = newprob - hmm_totalcounts_emit[source];
	    }
	}
	/* Transitions */
	for (source = 0; source < hmm->num_states; source++) {
	    for (target = 0; target < hmm->num_states; target++) {
		newprob = *HMM_TRANSITION_COUNTS(hmm_counts_trans, source, target);
		if (newprob == LOGZERO) { newprob = SMRZERO_LOG; }
		if (vb)
		    *HMM_TRANSITION_PROB(hmm, source, target) = (digamma(EXP(newprob) + g_gibbs_beta) - digamma(EXP(hmm_totalcounts_trans[source]) + (hmm->num_states) * g_gibbs_beta_emission))/M_LN2;
		else
		    *HMM_TRANSITION_PROB(hmm, source, target) = newprob - hmm_totalcounts_trans[source];
	    }
	}
	g_lasthmm = hmm;                               /* Put fsm into global var to recover in case of SIGINT */
	signal(SIGINT, (void *)interrupt_sigproc_hmm); /* Re-enable interrupt */
	prevloglikelihood = g_loglikelihood;
    }
    for (i = 0; i < g_num_threads; i++) {
	free(trellisarray[i]);
	free(threadargs[i]);
    }
    free(hmm_counts_trans);
    free(hmm_counts_emit);
    free(hmm_totalcounts_trans);
    free(hmm_totalcounts_emit);
    free(hmm_counts_spin);
    free(obsarray);
    return(g_loglikelihood);
}


PROB train_baum_welch(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta, int vb) {
    struct trellis *trellis, *trellisarray[32];
    struct thread_args *threadargs[32];
    struct observations **obsarray;
    int i, source, target, symbol, iter, numobs, obsperthread;
    PROB newprob, prevloglikelihood, da_beta = 1.0, numstatetrans;
    pthread_t threadids[32];
    
    if (g_train_da_bw) { da_beta = g_betamin; }
    obsarray = observations_to_array(o, &numobs);
    
    /* Each thread gets its own trellis */
    for (i = 0; i < g_num_threads; i++) {
	trellis = trellis_init(o, fsm->num_states);
	trellisarray[i] = trellis;
	threadargs[i] = malloc(sizeof(struct thread_args));
    }

    /* These are accessed by all threads through a spinlock */
    fsm_counts = malloc(fsm->num_states * fsm->num_states * fsm->alphabet_size * sizeof(PROB));
    fsm_totalcounts = malloc(fsm->num_states * sizeof(PROB));
    fsm_finalcounts = malloc(fsm->num_states * sizeof(PROB));
    fsm_counts_spin = calloc(fsm->num_states, sizeof(_Bool));
    
    prevloglikelihood = 0;

    obsperthread = (int) ((double) numobs/(double) g_num_threads);
    /* Helper thread parameters (divide observations equally among all) */
    for (i = 1; i < g_num_threads; i++) {
	threadargs[i]->minobs = (i-1) * obsperthread;
	threadargs[i]->maxobs = i * obsperthread - 1;
	threadargs[i]->trellis = trellisarray[i];
	threadargs[i]->obsarray = obsarray;
	threadargs[i]->fsmhmm = fsm;
    }
    /* Main thread parameters */
    threadargs[0]->minobs = (i-1) * obsperthread;
    threadargs[0]->maxobs = numobs - 1;
    threadargs[0]->trellis = trellisarray[0];
    threadargs[0]->obsarray = obsarray;
    threadargs[0]->fsmhmm = fsm;

   
    for (iter = 0 ; iter < maxiterations ; iter++) {
	g_loglikelihood = 0;
	for (i = 0; i < fsm->num_states * fsm->num_states * fsm->alphabet_size ; i++) { fsm_counts[i] = LOGZERO; }
	for (i = 0; i < fsm->num_states ; i++) { fsm_finalcounts[i] = LOGZERO; }
	for (i = 1; i < g_num_threads; i++) {
	    /* Launch threads */
	    threadargs[i]->beta = da_beta;
	    pthread_create(&threadids[i], NULL, &trellis_fill_bw, threadargs[i]);
	}

	/* Run main thread Baum-Welch */
	threadargs[0]->beta = da_beta;
	trellis_fill_bw(threadargs[0]);
	/* Wait for all to finish */
	for (i = 1; i < g_num_threads; i++) {
	    pthread_join(threadids[i],NULL);
	}

	if (!g_train_da_bw)
	    fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g\n", iter+1, g_loglikelihood, ABS(prevloglikelihood - g_loglikelihood));
	else 
	    fprintf(stderr, "iteration %i loglikelihood=%.17g delta: %.17g beta: %.17g\n", iter+1, g_loglikelihood, ABS(prevloglikelihood - g_loglikelihood), da_beta);
	    
	if (ABS(prevloglikelihood - g_loglikelihood) < maxdelta)  {
	    if (g_train_da_bw == 1 && da_beta < g_betamax) {
		da_beta *= g_alpha;
		if (da_beta > g_betamax) {
		    da_beta = g_betamax;
		}
	    } else {
		break;
	    }
	}

	/* Modify WFSA (M-step) */
	signal(SIGINT, SIG_IGN); /* Disable interrupts to prevent corrupted WFSA in case of SIGINT while updating */	    
        
	numstatetrans = fsm->num_states * fsm->alphabet_size + 1;
	/* Sum counts */
	for (i = 0; i < fsm->num_states ; i++) { fsm_totalcounts[i] = LOGZERO; }
	for (source = 0; source < fsm->num_states; source++) {
	    for (symbol = 0; symbol < fsm->alphabet_size; symbol++) {
		for (target = 0; target < fsm->num_states; target++) {
		    fsm_totalcounts[source] = log_add(*FSM_COUNTS(fsm_counts,source,symbol,target), fsm_totalcounts[source]);
		}
	    }
	}
	for (source = 0; source < fsm->num_states; source++) {
	    fsm_totalcounts[source] = log_add(fsm_finalcounts[source], fsm_totalcounts[source]);
	}

	for (source = 0; source < fsm->num_states; source++) {
	    for (symbol = 0; symbol < fsm->alphabet_size; symbol++) {
		for (target = 0; target < fsm->num_states; target++) {
		    newprob = *FSM_COUNTS(fsm_counts,source,symbol,target);
		    if (newprob == LOGZERO) { newprob = SMRZERO_LOG; }
		    /* Variational Bayes: apply digamma function to count */
		    if (vb) {
			*TRANSITION(fsm,source,symbol,target) = (digamma(EXP(newprob) + g_gibbs_beta) - digamma(EXP(fsm_totalcounts[source]) + numstatetrans * g_gibbs_beta))/M_LN2;
		    } else {
			*TRANSITION(fsm,source,symbol,target) = newprob - fsm_totalcounts[source];
		    }
		}
	    }
	}
	for (source = 0; source < fsm->num_states; source++) {
	    newprob = fsm_finalcounts[source];
	    if (newprob == LOGZERO) { newprob = SMRZERO_LOG; }
	    /* Variational Bayes */
	    if (vb) {
		*FINALPROB(fsm,source) = (digamma(EXP(newprob) + g_gibbs_beta) - digamma(EXP(fsm_totalcounts[source]) + numstatetrans * g_gibbs_beta))/M_LN2;
	    } else {
		*FINALPROB(fsm,source) = newprob - fsm_totalcounts[source];
	    }
	}
	g_lastwfsa = fsm;                          /* Put fsm into global var to recover in case of SIGINT */
	signal(SIGINT, (void *)interrupt_sigproc); /* Re-enable interrupt */
	prevloglikelihood = g_loglikelihood;
    }

    for (i = 0; i < g_num_threads; i++) {
	free(trellisarray[i]);
	free(threadargs[i]);
    }
    free(fsm_counts);
    free(fsm_totalcounts);
    free(fsm_finalcounts);
    free(fsm_counts_spin);
    free(obsarray);
    return(g_loglikelihood);
}

PROB train_bw_hmm(struct hmm *hmm, struct observations *o, int maxiterations, PROB maxdelta) {
    struct hmm *besthmm;
    int i;
    PROB thisll, minll;
    if (g_random_restarts == 0)
	return(train_baum_welch_hmm(hmm, o, maxiterations, maxdelta, g_bw_vb));

    /* Random restarts */
    minll = -DBL_MAX;
    besthmm = NULL;
    for (i = 0; i < g_random_restarts; i++) {
	fprintf(stderr, "===Running BW restart #%i===\n", i+1);
	thisll = train_baum_welch_hmm(hmm, o, g_random_restart_iterations, maxdelta, g_bw_vb);
	if (thisll > minll) {
	    minll = thisll;
	    if (besthmm != NULL)
		hmm_destroy(besthmm);
	    besthmm = hmm_copy(hmm);
	}
	if (i < g_random_restarts + 1) {
	    hmm = hmm_init(g_num_states, g_alphabet_size);
	    if (g_generate_type == GENERATE_BAKIS)
		hmm_randomize(hmm,1,0);
	    else
		hmm_randomize(hmm,0,0);
	    hmm_to_log2(hmm);
	}
    }
    fprintf(stderr, "===Running final BW===\n");
    return(train_baum_welch_hmm(besthmm, o, maxiterations, maxdelta, g_bw_vb));
}

PROB train_bw(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta) {
    struct wfsa *bestfsm;
    int i;
    PROB thisll, minll;
    if (g_random_restarts == 0)
	return(train_baum_welch(fsm,o,maxiterations,maxdelta,g_bw_vb));

    /* Random restarts */
    minll = -DBL_MAX;
    bestfsm = NULL;
    for (i = 0; i < g_random_restarts; i++) {
	fprintf(stderr, "===Running BW restart #%i===\n", i+1);
	thisll = train_baum_welch(fsm, o, g_random_restart_iterations, maxdelta,g_bw_vb);
	if (thisll > minll) {
	    minll = thisll;
	    if (bestfsm != NULL)
		wfsa_destroy(bestfsm);
	    bestfsm = wfsa_copy(fsm);
	}
	if (i < g_random_restarts + 1) {
	    fsm = wfsa_init(g_num_states, g_alphabet_size);
	    if (g_generate_type == GENERATE_NONDETERMINISTIC)
		wfsa_randomize_nondeterministic(fsm,0,0);
	    if (g_generate_type == GENERATE_DETERMINISTIC)
		wfsa_randomize_deterministic(fsm,0);
	    if (g_generate_type == GENERATE_BAKIS)
		wfsa_randomize_nondeterministic(fsm,1,0);
	    wfsa_to_log2(fsm);
	}
    }
    fprintf(stderr, "===Running final BW===\n");
    return(train_baum_welch(bestfsm, o, maxiterations, maxdelta, g_bw_vb));
}

PROB train_viterbi_bw(struct wfsa *fsm, struct observations *o) {
    PROB ll;
    train_viterbi(fsm,o,g_maxiterations,g_maxdelta);
    ll = train_baum_welch(fsm,o,g_maxiterations,g_maxdelta,g_bw_vb);
    return(ll);
}

PROB train_viterbi_bw_hmm(struct hmm *hmm, struct observations *o) {
    PROB ll;
    train_viterbi_hmm(hmm, o, g_maxiterations, g_maxdelta);
    ll = train_baum_welch_hmm(hmm, o, g_maxiterations, g_maxdelta, g_bw_vb);
    return(ll);
}

int main(int argc, char **argv) {
    int opt, option_index = 0, algorithm = 0, numelem, obs_alphabet_size, use_cuda = 0, use_hmm = 0, statemergetest = MERGE_TEST_ALERGIA, recursive_merge_test = 0;
    char *fsmfile = NULL, optionchar;
    PROB ll;
    struct wfsa *fsm = NULL;
    struct hmm *hmm = NULL;
    struct observations *o = NULL;
    srandom((unsigned int)time((time_t *)NULL));
    srand48((unsigned int)time((time_t *)NULL));
    statemergetest = MERGE_TEST_ALERGIA;

    log1plus_taylor_init();
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    g_numcpus = sysinfo.dwNumberOfProcessors;
#else
    g_numcpus = (int) sysconf(_SC_NPROCESSORS_ONLN);
#endif

static struct option long_options[] =
	{
	    {"annealing-params",required_argument, 0, 'a'},
	    {"burnin",          required_argument, 0, 'b'},
	    {"max-delta",       required_argument, 0, 'd'},
	    {"file",            required_argument, 0, 'f'},
	    {"initialize",      required_argument, 0, 'g'},
	    {"help",                  no_argument, 0, 'h'},
	    {"lag",             required_argument, 0, 'l'},
	    {"input-format",    required_argument, 0, 'i'},
	    {"output-format",   required_argument, 0, 'o'},
	    {"prior",           required_argument, 0, 'p'},
	    {"restarts",        required_argument, 0, 'r'},
	    {"threads",         required_argument, 0, 't'},
	    {"uniform-probs",         no_argument, 0, 'u'},
	    {"version",               no_argument, 0, 'v'},
	    {"max-iterations",  required_argument, 0, 'x'},
	    {"t0",              required_argument, 0, 'y'},
	    {"alpha",           required_argument, 0, 'A'},
	    {"cuda",                  no_argument, 0, 'C'},
	    {"decode",          required_argument, 0, 'D'},
	    {"generate",        required_argument, 0, 'G'},
	    {"hmm",                   no_argument, 0, 'H'},
	    {"likelihood",      required_argument, 0, 'L'},
	    {"merge-test",      required_argument, 0, 'M'},
	    {"recursive-merge",       no_argument, 0, 'R'},
	    {"train",           required_argument, 0, 'T'},
	    {0, 0, 0, 0}
	};

 while ((opt = getopt_long(argc, argv, "a:b:d:f:g:hl:i:o:p:r:t:uvx:y:A:CD:G:HL:M:RT:", long_options, &option_index)) != -1) {
	switch(opt) {
	case 'v':
	    printf("This is %s\n", versionstring);
	    exit(0);
	case 'h':
	    printf("%s", helpstring);
	    exit(0);
	case 'u':
	    g_initialize_uniform = 1;
	    break;
	case 'a':
	    numelem = sscanf(optarg,"%lg,%lg,%lg",&g_betamin,&g_betamax,&g_alpha);
	    if (numelem < 3) {
		fprintf(stderr, "-a option requires betamin,betamax,alpha\n"); 
		exit(1);
	    }
	    break;
	case 'r':
	    numelem = sscanf(optarg,"%i,%i", &g_random_restarts, &g_random_restart_iterations);
	    if (numelem < 1) {
		fprintf(stderr, "-r option requires #num restarts[,#num iterations per restart]\n"); 
		exit(1);
	    }
	    break;
	case 'g':	    
	    if (*optarg == 'b' || *optarg == 'd' || *optarg == 'n') {
		numelem = sscanf(optarg,"%c%i,%i", &optionchar, &g_num_states, &g_alphabet_size);
		switch(optionchar) {
		case 'b': g_generate_type = GENERATE_BAKIS; break;
		case 'd': g_generate_type = GENERATE_DETERMINISTIC; break;
		case 'n': g_generate_type = GENERATE_NONDETERMINISTIC; break;
		default: fprintf(stderr, "-g option requires b|d|n\n"); exit(EXIT_FAILURE);
		}
		if (numelem < 2) { printf("Error in option -g\n"); exit(1); }
		if (numelem < 3) { g_alphabet_size = -1; }
	    } else { /* Default is nondeterministic, needs no specification */
		g_generate_type = GENERATE_NONDETERMINISTIC;
		numelem = sscanf(optarg,"%i,%i", &g_num_states, &g_alphabet_size);
		if (numelem < 2)
		    g_alphabet_size = -1;
	    }
	    break;
	case 't':
	    if (strncmp(optarg,"c/",2) == 0) {
		g_num_threads = g_numcpus / atoi(optarg+2);
	    } else if (strncmp(optarg,"c",1) == 0) {
		g_num_threads = g_numcpus - atoi(optarg+1);
	    } else {
		g_num_threads = atoi(optarg);
	    }
	    if (g_num_threads <= 0) { g_num_threads = 1; }
	    break;
	case 'y':
	    g_t0 = atoi(optarg);
	case 'd':
	    g_maxdelta = strtod(optarg, NULL);
	    break;
	case 'R':
	    recursive_merge_test = 1;
	    break;
	case 'M':
	    if (strcmp(optarg,"alergia") == 0)  { statemergetest = MERGE_TEST_ALERGIA;    }
	    if (strcmp(optarg,"chi2") == 0)     { statemergetest = MERGE_TEST_CHISQUARED; }
	    if (strcmp(optarg,"lr") == 0)       { statemergetest = MERGE_TEST_LR;         }
	    if (strcmp(optarg,"binomial") == 0) { statemergetest = MERGE_TEST_BINOMIAL;   }
	    if (strcmp(optarg,"exactm") == 0)   { statemergetest = MERGE_TEST_EXACT_M;    }
	    if (strcmp(optarg,"exactb") == 0)   { statemergetest = MERGE_TEST_EXACT_B;    }
	    break;
	case 'i':
	    if (strcmp(optarg,"real") == 0)   {  g_input_format = FORMAT_REAL;   }
	    if (strcmp(optarg,"log10") == 0)  {  g_input_format = FORMAT_LOG10;  }
	    if (strcmp(optarg,"ln") == 0)     {  g_input_format = FORMAT_LN;     }
	    if (strcmp(optarg,"log2") == 0)   {  g_input_format = FORMAT_LOG2;   }
	    if (strcmp(optarg,"nlog10") == 0) {  g_input_format = FORMAT_NLOG10; }
	    if (strcmp(optarg,"nln") == 0)    {  g_input_format = FORMAT_NLN;    }
	    if (strcmp(optarg,"nlog2") == 0)  {  g_input_format = FORMAT_NLOG2;  }
	    break;
	case 'o':
	    if (strcmp(optarg,"real") == 0)   {  g_output_format = FORMAT_REAL;   }
	    if (strcmp(optarg,"log10") == 0)  {  g_output_format = FORMAT_LOG10;  }
	    if (strcmp(optarg,"ln") == 0)     {  g_output_format = FORMAT_LN;     }
	    if (strcmp(optarg,"log2") == 0)   {  g_output_format = FORMAT_LOG2;   }
	    if (strcmp(optarg,"nlog10") == 0) {  g_output_format = FORMAT_NLOG10; }
	    if (strcmp(optarg,"nln") == 0)    {  g_output_format = FORMAT_NLN;    }
	    if (strcmp(optarg,"nlog2") == 0)  {  g_output_format = FORMAT_NLOG2;  }
	    break;
	case 'f':
	    fsmfile = strdup(optarg);
	    break;
	case 'A':
	    g_merge_alpha = strtod(optarg, NULL);
	    break;
	case 'p':
	    //g_viterbi_pseudocount = atoi(optarg);
	    numelem = sscanf(optarg,"%lg,%lg",&g_gibbs_beta,&g_gibbs_beta_emission);
	    g_viterbi_pseudocount = g_gibbs_beta;
	    if (numelem > 1)
		g_viterbi_pseudocount_emit = g_gibbs_beta_emission;
	    g_merge_prior = g_gibbs_beta;
	    break;
	case 'b':
	    g_gibbs_burnin = atoi(optarg);
	    break;
	case 'l':
	    g_gibbs_lag = atoi(optarg);
	    break;
	case 'x':
	    g_maxiterations = atoi(optarg);
	    break;
	case 'C':
#ifdef USE_CUDA
	    use_cuda = 1;
#else
	    perror("CUDA support not compiled");
            exit(EXIT_FAILURE);
#endif /* USE_CUDA */
	    break;
	case 'H':
	    use_hmm = 1;
	    break;
	case 'G':
	    g_generate_words = atoi(optarg);
	    algorithm = GENERATE_WORDS;
	    break;
	case 'T':
	    if (strcmp(optarg,"merge") == 0) { algorithm = TRAIN_MERGE;          }
	    if (strcmp(optarg,"mdi") == 0)   { algorithm = TRAIN_MDI;            }
	    if (strcmp(optarg,"vit") == 0)   { algorithm = TRAIN_VITERBI;        }
	    if (strcmp(optarg,"bw") == 0)    { algorithm = TRAIN_BAUM_WELCH;     }
	    if (strcmp(optarg,"gs") == 0)    { algorithm = TRAIN_GIBBS_SAMPLING; }
	    if (strcmp(optarg,"vb") == 0) {
		algorithm = TRAIN_VARIATIONAL_BAYES;
		g_bw_vb = 1;
	    }
	    if (strcmp(optarg,"vitbw") == 0)  { algorithm = TRAIN_VITERBI_BW;    }
	    if (strcmp(optarg,"dabw") == 0) {
	        algorithm = TRAIN_DA_BAUM_WELCH;
		g_train_da_bw = 1;
	    }
	    break;
	case 'L':
	    if (strcmp(optarg,"vit") == 0) { algorithm = LIKELIHOOD_VITERBI;  }
	    if (strcmp(optarg,"f") == 0)   { algorithm = LIKELIHOOD_FORWARD;  }
	    if (strcmp(optarg,"b") == 0)   { algorithm = LIKELIHOOD_BACKWARD; }
	    break;
	case 'D':
	    if (strcmp(optarg,"vit") == 0)   { algorithm = DECODE_VITERBI;       }
	    if (strcmp(optarg,"vit,p") == 0) { algorithm = DECODE_VITERBI_PROB;  }
	    if (strcmp(optarg,"f") == 0)     { algorithm = DECODE_FORWARD;       }
	    if (strcmp(optarg,"f,p") == 0)   { algorithm = DECODE_FORWARD_PROB;  }
	    if (strcmp(optarg,"b") == 0)     { algorithm = DECODE_BACKWARD;      }
	    if (strcmp(optarg,"b,p") == 0)   { algorithm = DECODE_BACKWARD_PROB; }
	    break;
	}
    }
    argc -= optind;
    argv += optind;

    if (use_hmm && (algorithm == TRAIN_MERGE || algorithm == TRAIN_MDI)) {
	perror("HMMs are not supported for this training/inference algorithm");
	exit(EXIT_FAILURE);
    }
    if (argc < 1 && ((algorithm && algorithm != GENERATE_WORDS) || (g_alphabet_size < 0 && fsmfile == NULL))) {
	fprintf(stderr, "Missing observation filename\n");
	fprintf(stderr, "Usage: %s",usagestring);
	exit(EXIT_FAILURE);
    }
    if (argc > 0) {
	if ((o = observations_read(argv[0])) == NULL) {
	    perror("Error reading observations file");	    
	    exit(EXIT_FAILURE);
	}
	obs_alphabet_size = observations_alphabet_size(o);
    }
    if (fsmfile == NULL && g_generate_type == 0 && (algorithm != TRAIN_MERGE && algorithm != TRAIN_MDI)) {
	perror("You must either specify a FSM file with -f, or initialize a random FSM with -g");
	exit(EXIT_FAILURE);
    }
    if (g_alphabet_size < 0 && o != NULL) {
	g_alphabet_size = obs_alphabet_size;
    }
    if (g_generate_type > 0) {
	if (!use_hmm) {
	    fsm = wfsa_init(g_num_states, g_alphabet_size);
	}
	if (use_hmm) {
	    hmm = hmm_init(g_num_states, g_alphabet_size);
	}
	if (algorithm != TRAIN_GIBBS_SAMPLING) { /* For Gibbs, we don't initialize anything, just get size */
	    if (!use_hmm) {
		if (g_generate_type == GENERATE_NONDETERMINISTIC)
		    wfsa_randomize_nondeterministic(fsm, 0, g_initialize_uniform);
		if (g_generate_type == GENERATE_DETERMINISTIC)
		    wfsa_randomize_deterministic(fsm, g_initialize_uniform);
		if (g_generate_type == GENERATE_BAKIS)
		    wfsa_randomize_nondeterministic(fsm, 1, g_initialize_uniform);
		wfsa_to_log2(fsm);
	    } else {
		if (g_generate_type == GENERATE_NONDETERMINISTIC)
		    hmm_randomize(hmm, 0, 0);
		if (g_generate_type == GENERATE_BAKIS) {
		    hmm_randomize(hmm, 1, 0);		
		    printf("BAKIS\n");
		}
		hmm_to_log2(hmm);
	    }
	}
    }

    if (fsmfile != NULL) {
	if (!use_hmm)
	    fsm = wfsa_read_file(fsmfile);
	else
	    hmm = hmm_read_file(fsmfile);
	if (g_input_format != FORMAT_LOG2) {
	    if (!use_hmm)
		wfsa_to_log2(fsm);
	    else
		hmm_to_log2(hmm);
	}
    }
    if (o != NULL && ((fsm != NULL && fsm->alphabet_size < obs_alphabet_size) | (hmm != NULL && hmm->alphabet_size < obs_alphabet_size))) {
	fprintf(stderr, "Error: the observations file has symbols outside the FSA alphabet.\n");
	fprintf(stderr, "FSA alphabet size: %i  Observations alphabet size %i.\n", fsm->alphabet_size, obs_alphabet_size);
	exit(1);
    }

    log1plus_init();
	
    switch (algorithm) {
    case GENERATE_WORDS:
	if (!use_hmm)
	    generate_words(fsm, g_generate_words);
	else
	    generate_words_hmm(hmm, g_generate_words);
	break;
    case DECODE_VITERBI:
    case DECODE_VITERBI_PROB:
    case LIKELIHOOD_VITERBI:
	if (!use_hmm)
	    viterbi(fsm, o, algorithm);
	else
	    viterbi_hmm(hmm, o, algorithm);
	break;
    case DECODE_FORWARD:
    case DECODE_FORWARD_PROB:
    case LIKELIHOOD_FORWARD:
	if (!use_hmm)
	    forward_fsm(fsm, o, algorithm);
	else
	    forward_hmm(hmm, o, algorithm);
	break;
    case DECODE_BACKWARD:
    case DECODE_BACKWARD_PROB:
    case LIKELIHOOD_BACKWARD:
	if (!use_hmm)
	    backward_fsm(fsm, o, algorithm);
	else
	    backward_hmm(hmm, o, algorithm);
	break;
    case TRAIN_VARIATIONAL_BAYES:
    case TRAIN_BAUM_WELCH:
    case TRAIN_DA_BAUM_WELCH:
	o = observations_sort(o);
	o = observations_uniq(o);
	if (!use_hmm) {
	    train_bw(fsm, o, g_maxiterations, g_maxdelta);
	    wfsa_print(fsm);
	} else {
	    train_bw_hmm(hmm, o, g_maxiterations, g_maxdelta);
	    hmm_print(hmm);
	}
	break;
    case TRAIN_VITERBI:
	o = observations_sort(o);
	o = observations_uniq(o);
	if (!use_hmm) {
	    train_viterbi(fsm, o, g_maxiterations, g_maxdelta);
	    wfsa_print(fsm);
	} else {
	    train_viterbi_hmm(hmm, o, g_maxiterations, g_maxdelta);
	    hmm_print(hmm);
	}
	break;
    case TRAIN_VITERBI_BW:
	o = observations_sort(o);
	o = observations_uniq(o);
	if (!use_hmm) {
	    train_viterbi_bw(fsm,o);
	    wfsa_print(fsm);
	} else {
	    train_viterbi_bw_hmm(hmm, o);
	    hmm_print(hmm);
	}
	break;
    case TRAIN_MERGE:
	o = observations_sort(o);
	o = observations_uniq(o);
	fsm = dffa_to_wfsa(dffa_state_merge(o, g_merge_alpha, statemergetest, recursive_merge_test));
	ll = loglikelihood_all_observations_fsm(fsm,o);
	fprintf(stderr, "loglikelihood=%.17g\n", ll);
	wfsa_print(fsm);
	break;
    case TRAIN_MDI:
	o = observations_sort(o);
	o = observations_uniq(o);
	fsm = dffa_to_wfsa(dffa_mdi(o, g_merge_alpha));
	ll = loglikelihood_all_observations_fsm(fsm,o);
	fprintf(stderr, "loglikelihood=%.17g\n", ll);
	wfsa_print(fsm);
	break;
    case TRAIN_GIBBS_SAMPLING:
	o = observations_sort(o);
	if (!use_hmm) {
	    if (!use_cuda) {
		ll = gibbs_sampler_fsm(fsm, o, g_gibbs_beta, g_num_states, g_maxiterations, g_gibbs_burnin, g_gibbs_lag);
	    } else {
		#ifdef USE_CUDA
		fprintf(stderr, "Launching CUDA.\n");
		ll = gibbs_sampler_cuda_fsm(fsm, o, g_gibbs_beta, g_num_states, g_maxiterations, g_gibbs_burnin, g_gibbs_lag);
		#endif
	    }
	    fprintf(stderr, "loglikelihood=%.17g\n", ll);
	    wfsa_print(fsm);
	} else {
	    if (!use_cuda) {
		ll = gibbs_sampler_hmm(hmm, o, g_gibbs_beta_emission, g_gibbs_beta, g_num_states, g_maxiterations, g_gibbs_burnin, g_gibbs_lag);
	    } else {
		#ifdef USE_CUDA
		fprintf(stderr, "Launching CUDA.\n");
		ll = gibbs_sampler_cuda_hmm(hmm, o, g_gibbs_beta_emission, g_gibbs_beta, g_num_states, g_maxiterations, g_gibbs_burnin, g_gibbs_lag);
		#endif
	    }
	    fprintf(stderr, "loglikelihood=%.17g\n", ll);
	    hmm_print(hmm);
	}
	break;
    default:
	if (use_hmm)
	    hmm_print(hmm);
	else
	    wfsa_print(fsm);
    }
    log1plus_free();
    if (fsm != NULL)
	wfsa_destroy(fsm);
    if (hmm != NULL)
	hmm_destroy(hmm);
    if (o != NULL)
	observations_destroy(o);
    exit(0);
}
