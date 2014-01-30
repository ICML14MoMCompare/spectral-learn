/*******************/
/* Default options */
/*******************/

int g_numcpus = 1;
int g_maxiterations = 20000;
PROB g_maxdelta = 0.1;
int g_input_format = FORMAT_REAL;
int g_output_format = FORMAT_REAL;
PROB g_viterbi_pseudocount = 1;
PROB g_viterbi_pseudocount_emit = 1;
int g_random_restarts = 0;
int g_random_restart_iterations = 3;
int g_generate_type = 0;
int g_gen_max_length = 100000;
int g_num_states = 5;
int g_alphabet_size = -1;
int g_initialize_uniform = 0;
int g_train_da_bw = 0;
int g_generate_words = 0;
/* Thread variables */
int g_num_threads = 1;
/* Deterministic annealing default parameters */
PROB g_betamin = 0.02;
PROB g_betamax = 1;
PROB g_alpha = 1.01;
/* Collapsed Gibbs sampling concentration, lag, burnin parameters, also Variational bayes priors (Dirichlet parameter) */
PROB g_gibbs_beta = 0.02;          /* For PFSA transitions, HMM transitions, VB transition prior */
PROB g_gibbs_beta_emission = 0.01; /* For HMM emissions, VB emission prior */
int g_gibbs_lag = 10;
int g_gibbs_burnin = 100;
/* Flag whether to adjust counts in BW by VB strategy */
int g_bw_vb = 0;
int g_t0 = 3;              /* Min number of visits to a state for a state to be mergeable in state-merging */
PROB g_merge_alpha = 0.05; /* The alpha parameter for ALERGIA and MDI */
PROB g_merge_prior = 0.02; /* Prior for avoiding missing transitions/final states with 0 prob in state-merging */
