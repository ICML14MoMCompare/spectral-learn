/**************************************************************************/
/*   treba - probabilistic FSM and HMM training and decoding              */
/*   Copyright © 2012 Mans Hulden                                         */

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

#define LIKELIHOOD_FORWARD      1
#define LIKELIHOOD_BACKWARD     2
#define LIKELIHOOD_VITERBI      3
#define DECODE_FORWARD          4
#define DECODE_BACKWARD         5
#define DECODE_VITERBI          6
#define DECODE_FORWARD_PROB     7
#define DECODE_BACKWARD_PROB    8
#define DECODE_VITERBI_PROB     9
#define TRAIN_BAUM_WELCH        10
#define TRAIN_DA_BAUM_WELCH     11
#define TRAIN_VITERBI           12
#define TRAIN_VITERBI_BW        13
#define TRAIN_VARIATIONAL_BAYES 14
#define TRAIN_GIBBS_SAMPLING    15
#define TRAIN_MERGE             16
#define TRAIN_MDI               17
#define GENERATE_WORDS          18

/* State merging tests */
#define MERGE_TEST_ALERGIA      1 /* Alergia (Hoeffding bound test)                           */
#define MERGE_TEST_CHISQUARED   2 /* Chi squared                                              */
#define MERGE_TEST_LR           3 /* Likelihood ratio                                         */
#define MERGE_TEST_BINOMIAL     4 /* Approx binomial test                                     */
#define MERGE_TEST_EXACT_M      5 /* "Exact" (monte-carlo) multinomial test of p-value (slow) */
#define MERGE_TEST_EXACT_B      6 /* Exact binomial test of p-value (slow)                    */

#define GENERATE_NONDETERMINISTIC 1
#define GENERATE_DETERMINISTIC    2
#define GENERATE_UNIFORM          3
#define GENERATE_BAKIS            4

#define SMRONE_REAL  1

#define FORMAT_REAL     0
#define FORMAT_LOG10    1
#define FORMAT_LOG2     2
#define FORMAT_LN       3
#define FORMAT_NLOG10   4
#define FORMAT_NLOG2    5
#define FORMAT_NLN      6

#define LOG(X)        (log2((X)))
#define EXP(X)        (exp2((X)))
#define SMRZERO_LOG  -DBL_MAX
#define LOGZERO       DBL_MAX
#define ABS(X)        (fabs((X)))
typedef double PROB;

/* Auxiliary macros to access trellis, FSMs, and FSM counts */
#define FINALPROB(FSM, STATE) ((FSM)->final_table + (STATE))
#define TRELLIS_CELL(STATE, TIME) ((trellis) + (fsm->num_states) * (TIME) + (STATE))
#define TRELLIS_CELL_HMM(STATE, TIME) ((trellis) + (hmm->num_states) * (TIME) + (STATE))
#define TRANSITION(FSM, SOURCE_STATE, SYMBOL, TARGET_STATE) ((FSM)->state_table + (FSM)->num_states * (FSM)->alphabet_size * (SOURCE_STATE) + (SYMBOL) * (FSM)->num_states + (TARGET_STATE))
#define FSM_COUNTS(FSMC, SOURCE_STATE, SYMBOL, TARGET_STATE) ((FSMC) + (fsm->num_states * fsm->alphabet_size * (SOURCE_STATE) + (SYMBOL) * fsm->num_states + (TARGET_STATE)))

#define HMM_TRANSITION_COUNTS(HMMC, SOURCE_STATE, TARGET_STATE) ((HMMC) + (hmm->num_states * (SOURCE_STATE) + (TARGET_STATE)))
#define HMM_EMISSION_COUNTS(HMMC, STATE, SYMBOL) ((HMMC) + (hmm->alphabet_size * (STATE) + (SYMBOL)))

#define HMM_TRANSITION_PROB(HMM, SOURCE_STATE, TARGET_STATE) ((HMM)->transition_table + (HMM)->num_states * (SOURCE_STATE) + (TARGET_STATE))
#define HMM_EMISSION_PROB(HMM, STATE, SYMBOL) ((HMM)->emission_table + (HMM)->alphabet_size * (STATE) + (SYMBOL))

//PROB smrzero = SMRZERO_LOG;
//PROB smrone  = 0;

PROB *fsm_counts, *fsm_totalcounts, *fsm_finalcounts;
_Bool *fsm_counts_spin;
PROB *hmm_counts_trans, *hmm_counts_emit, *hmm_totalcounts_trans, *hmm_totalcounts_emit;
_Bool *hmm_counts_spin;

//pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;


struct thread_args {
    struct trellis *trellis;
    struct observations **obsarray;
    int minobs;
    int maxobs;
    void *fsmhmm;
    PROB beta;
};

struct observations *g_obsarray;

struct wfsa {
    int num_states;
    int alphabet_size;
    PROB *state_table;
    PROB *final_table;
};

struct hmm {
    int num_states;
    int alphabet_size;
    PROB *transition_table;
    PROB *emission_table;
};

struct observations {
    int size;
    int *data;
    int occurrences;
    struct observations *next;
};

struct trellis {
    PROB fp;
    PROB bp;
    int backstate;
};


 /* In io.c */

/* Input and output conversion */
PROB output_convert(PROB x);
PROB input_convert(PROB x);

/* Functions to handle file and text input */
char *file_to_mem(char *name);
int char_in_array(char c, char *array);
int line_count_elements(char **ptr);
char *line_to_int_array(char *ptr, int **line, int *size);

void hmm_print(struct hmm *hmm);


/* Gibbs */

struct gibbs_state_chain {
    int state;
    int sym;
};

PROB gibbs_sampler_fsm(struct wfsa *fsm, struct observations *o, double beta, int num_states, int maxiter, int burnin, int lag);
PROB gibbs_sampler_hmm(struct hmm *hmm, struct observations *o, double beta_e, double beta_t, int num_states, int maxiter, int burnin, int lag);
struct hmm *gibbs_counts_to_hmm(struct hmm *hmm, unsigned int *gibbs_sampled_counts_trans, unsigned int *gibbs_sampled_counts_emit, unsigned int *gibbs_counts_sampled_states, int alphabet_size, int num_states, double beta_t, double beta_e);
struct wfsa *gibbs_counts_to_wfsa(struct wfsa *fsm, unsigned int *gibbs_sampled_counts, unsigned int *gibbs_counts_sampled_states, int alphabet_size, int num_states, double beta, double ANbeta);
struct gibbs_state_chain *gibbs_init_fsm(struct observations *o, int num_states, int alphabet_size, int *obslen);
struct gibbs_state_chain *gibbs_init_hmm(struct observations *o, int num_states, int alphabet_size, int *obslen);
PROB gibbs_sampler_fsm(struct wfsa *fsm, struct observations *o, double beta, int num_states, int maxiter, int burnin, int lag);
PROB gibbs_sampler_hmm(struct hmm *hmm, struct observations *o, double beta_e, double beta_t, int num_states, int maxiter, int burnin, int lag);

void interrupt_sigproc(void);

inline void spinlock_lock(_Bool *ptr);
inline void spinlock_unlock(_Bool *ptr);

PROB rand_double();
int rand_int_range(int from, int to);

/* WFSA functions */
struct wfsa *wfsa_read_file(char *filename);
void wfsa_print(struct wfsa *fsm);
void wfsa_randomize_deterministic(struct wfsa *fsm, int uniform);
void wfsa_randomize_nondeterministic(struct wfsa *fsm, int bakis, int uniform);
struct wfsa *wfsa_init(int num_states, int alphabet_size);
struct wfsa *wfsa_copy(struct wfsa *fsm);
void wfsa_destroy(struct wfsa *fsm);
void wfsa_to_log2(struct wfsa *fsm);
PROB wfsa_sum_prob(struct wfsa *fsm, int state);
int wfsa_random_transition(struct wfsa *fsm, int state, int *symbol, PROB *prob);

/* Generation functions */
void generate_words(struct wfsa *fsm, int numwords);

/* Observation file/array functions */
int obssortcmp(struct observations **a, struct observations **b);
int observations_alphabet_size(struct observations *ohead);
struct observations **observations_to_array(struct observations *ohead, int *numobs);
struct observations *observations_uniq(struct observations *ohead);
struct observations *observations_sort(struct observations *ohead);
void observations_destroy(struct observations *ohead);
struct observations *observations_read(char *filename);


PROB loglikelihood_all_observations_fsm(struct wfsa *fsm, struct observations *o);
PROB loglikelihood_all_observations_hmm(struct hmm *hmm, struct observations *o);

/* Trellis functions */
PROB trellis_backward(struct trellis *trellis, int *obs, int length, struct wfsa *fsm);
PROB trellis_viterbi(struct trellis *trellis, int *obs, int length, struct wfsa *fsm);
PROB trellis_forward_fsm(struct trellis *trellis, int *obs, int length, struct wfsa *fsm);
PROB trellis_forward_hmm(struct trellis *trellis, int *obs, int length, struct hmm *hmm);
struct trellis *trellis_init(struct observations *o, int num_states);
void trellis_print(struct trellis *trellis, struct wfsa *fsm, int obs_len);

/* Trellis path printing functions */
void forward_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len);
void backward_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len);
void viterbi_print_path(struct trellis *trellis, struct wfsa *fsm, int obs_len);

/* Main decoding and likelihood calculations */
void viterbi(struct wfsa *fsm, struct observations *o, int algorithm);
void forward_fsm(struct wfsa *fsm, struct observations *o, int algorithm);
void forward_hmm(struct hmm *hmm, struct observations *o, int algorithm);
void backward_fsm(struct wfsa *fsm, struct observations *o, int algorithm);
void backward_hmm(struct hmm *hmm, struct observations *o, int algorithm);

/* Main training functions */
PROB train_viterbi(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta);
PROB train_viterbi_hmm(struct hmm *hmm, struct observations *o, int maxiterations, PROB maxdelta);
PROB train_baum_welch(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta, int vb);
PROB train_bw(struct wfsa *fsm, struct observations *o, int maxiterations, PROB maxdelta);
PROB train_viterbi_bw(struct wfsa *fsm, struct observations *o);
void *trellis_fill_bw(void *threadargs);

int main(int argc, char **argv);

/* dffa.c */

struct dffa {
    int num_states;
    int alphabet_size;
    int *transitions;
    int *transition_freqs;
    int *final_freqs;
    int *total_freqs;
};

struct wfsa *dffa_to_wfsa(struct dffa *dffa);
struct dffa *dffa_state_merge(struct observations *o, PROB alpha, int test, int recursive);
struct dffa *dffa_mdi(struct observations *o, PROB alpha);
struct dffa *observations_to_dffa(struct observations *o);
struct dffa *dffa_init(int num_states, int alphabet_size);
int dffa_chi2_test(struct dffa *dffa, int qu, int qv, double alpha);
