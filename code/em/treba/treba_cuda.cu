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
#include <curand.h>
#include <curand_kernel.h>

#define THREADS_PER_BLOCK 256

#define UDIV_UP(a, b) (((a) + (b) - 1) / (b))
#define ALIGN_UP(a, b) (UDIV_UP(a, b) * (b))

#define CudaCheckError() __cudaCheckError(__FILE__,__LINE__)
inline void __cudaCheckError(const char *file, const int line) {
    cudaError err = cudaGetLastError();
    if (cudaSuccess != err) {
        fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
        exit(-1);
    }
}

struct gibbs_state_chain {
    int state;
    int sym;
};

/* Macro for accessing counts on host */
/* FSM */
#define Ccurr(SOURCE_STATE, SYMBOL, TARGET_STATE) (*((gibbs_counts) + (num_states * alphabet_size * (SOURCE_STATE) + (SYMBOL) * num_states + (TARGET_STATE))))

/* HMM */
#define CcurrHMMemit(STATE, SYMBOL) (*((gibbs_counts_emit) + (alphabet_size * (STATE) + (SYMBOL))))
#define CcurrHMMtrans(SOURCE_STATE, TARGET_STATE) (*((gibbs_counts_trans) + (num_states * (SOURCE_STATE) + (TARGET_STATE))))

/* Macro for accessing current counts on CUDA device */
/* FSM */
#define Ccurrdevice(SOURCE_STATE, SYMBOL, TARGET_STATE) (((d_gibbs_counts) + (dc_num_states * dc_alphabet_size * (SOURCE_STATE) + (SYMBOL) * dc_num_states + (TARGET_STATE))))

/* HMM */
#define CcurrdeviceEMIT(STATE, SYMBOL) ((d_gibbs_counts_emit + (dc_alphabet_size * (STATE) + (SYMBOL))))
#define CcurrdeviceTRANS(SOURCE, TARGET) ((d_gibbs_counts_trans + (dc_num_states * (SOURCE) + (TARGET))))

__global__ void init_random_generator(curandState *state, unsigned long seed) {
    int id;
    id = blockIdx.x * blockDim.x + threadIdx.x;

    /* From the curand manual 4.2, p. 20: */
    /* State setup can be an expensive operation. One way to speed up the setup is to use different seeds */
    /* for each thread and a constant sequence number of 0. This can be especially helpful if many        */ 
    /* generators need to be created. While faster to set up, this method provides less guarantees about  */
    /* the mathematical properties of the generated sequences. If there happens to be a bad interaction   */
    /* between the hash function that initializes the generator state from the seed and the periodicity   */ 
    /* of the generators, there might be threads with highly correlated outputs for some seed values.     */
    /* We do not know of any problem values; if they do exist they are likely to be rare.                 */

    /* In other words, the slow method is to use the same seed but different sequence numbers (id):       */
    /* curand_init(seed, id, 0, &state[id]);                                                              */
    /* Here, we trust their word that problem values are rare and use the (much) faster setup:            */
    curand_init(seed+id, 0, 0, &state[id]);
}

/* Update sampled (global) counts, i.e. add current counts to sampled counts */
__global__ void gibbs_sampler_update_kernel_fsm(int numthreads, unsigned int *d_gibbs_counts, unsigned int *d_gibbs_sampled_counts) {
    int id;
    id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= numthreads) { return; }
    d_gibbs_sampled_counts[id] += d_gibbs_counts[id];
}

__global__ void gibbs_sampler_update_kernel_hmm(int numthreads, unsigned int *d_gibbs_counts_trans, unsigned int *d_gibbs_counts_emit, unsigned int *d_gibbs_sampled_counts_trans, unsigned int *d_gibbs_sampled_counts_emit, int transtablesize, int emittablesize) {
    int id;
    id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= numthreads) { return; }
    if (id < transtablesize)
	d_gibbs_sampled_counts_trans[id] += d_gibbs_counts_trans[id];
    if (id < emittablesize)
	d_gibbs_sampled_counts_emit[id] += d_gibbs_counts_emit[id];
}

/* Constants (to avoid passing these as variables to the kernels) */
__constant__ int dc_num_states, dc_alphabet_size;
__constant__ unsigned int *dc_resamplable;
__constant__ float dc_beta, dc_beta_t, dc_beta_e, dc_ANbeta, *dc_weight_list;

/* Main Gibbs sampler kernels. Every thread chooses one new state k in the chain of */
/* observations proportional to the probability of the new chain that uses k.       */
/* We launch two kernels in series: one for choosing odd-numbered states, and       */
/* another for choosing even-numbered ones (since adjacent states depend on each    */
/* other, we can't sample _all_ states in parallel, only non-dependent ones.)       */

/* The process is the same for FSMs and HMMs:                                         */
/* (1) Calculate probabilities for all changes z_i -> k (for all possible k) at o_i   */
/* (2) Do a weighted selection based on the probabilities and change the state z_i to */
/*     the selected new state.                                                        */ 
/* (3) Change the transition (and emission for HMM) counts to reflect the new choice. */

__global__ void gibbs_sampler_kernel_fsm(int numthreads, struct gibbs_state_chain *mychain, unsigned int *d_gibbs_counts_states, unsigned int *d_gibbs_counts, curandState *globalrandstate, short int evenodd) {
    unsigned int id, a, aprev, z, zprev, znext, k, low, high, mid, idorig;
    float g_sum, g_k, ind1, ind2, ind0, cointoss;
    float *weight_ptr;
    id = idorig = blockIdx.x * blockDim.x + threadIdx.x;
    weight_ptr = dc_weight_list + id;

    if (id >= numthreads) { return; }
    id = id * 2 + evenodd;
    mychain = mychain + dc_resamplable[id];

    aprev = (mychain-1)->sym;    /* Previous symbol */
    zprev = (mychain-1)->state;  /* Previous state  */
    a = mychain->sym;            /* Current symbol  */
    z = mychain->state;          /* Current state   */
    znext = (mychain+1)->state;  /* Next state      */
    
    /* Find probabilities of changing current state z to some state k for all k */
    for (k = 0, g_sum = 0; k < dc_num_states; k++) {
	ind0 = ((zprev == k && aprev == a && znext == k) ? 1.0 : 0.0);
	ind1 = (z == k ? 1.0 : 0.0);
	ind2 = ((zprev == k && aprev == a && z == znext) ? 1.0 : 0.0);

	g_k = (((float)*Ccurrdevice(k,a,znext)) - ind1 -ind2 + dc_beta) * (((float)*Ccurrdevice(zprev,aprev,k)) - ind1 - ind2 + dc_beta + ind0) / ((float)d_gibbs_counts_states[k] - ind1 + dc_ANbeta);
	g_sum += g_k;
        weight_ptr[numthreads * k] = g_sum;
    }

    /* Do a binary search for the first element in weight_list    */
    /* larger than cointoss. This diverges, but is still slightly */
    /* faster than a linear (more coalesced) search.              */

    cointoss = curand_uniform(&globalrandstate[idorig]) * g_sum;

    for (low = 0, high = dc_num_states - 1; low != high; ) {
    	mid = (low + high) / 2;
    	if (weight_ptr[numthreads * mid] <= cointoss) {
    	    low = mid + 1;
    	} else {
    	    high = mid;
    	}
    }
    k = high;

    /* Update counts */
    if (k != z) {
	atomicAdd(Ccurrdevice(zprev,aprev,k),1);
	atomicSub(Ccurrdevice(zprev,aprev,z),1);
	atomicAdd(Ccurrdevice(k,a,znext),1);
	atomicSub(Ccurrdevice(z,a,znext),1);
	atomicAdd(&d_gibbs_counts_states[k],1);
	atomicSub(&d_gibbs_counts_states[z],1);
	mychain->state = k;	
    }
}

__global__ void gibbs_sampler_kernel_hmm(int numthreads, struct gibbs_state_chain *mychain, unsigned int *d_gibbs_counts_states, unsigned int *d_gibbs_counts_trans, unsigned int *d_gibbs_counts_emit, curandState *globalrandstate, short int evenodd) {
    unsigned int id, a, z, zprev, znext, k, low, high, mid, idorig;
    float g_sum, g_k, ind0, ind1, ind2, ind3, cointoss;
    float *weight_ptr;
    id = idorig = blockIdx.x * blockDim.x + threadIdx.x;
    weight_ptr = dc_weight_list + id;

    if (id >= numthreads) { return; }
    id = id * 2 + evenodd;
    mychain = mychain + dc_resamplable[id];

    zprev = (mychain-1)->state;  /* Previous state  */
    a = mychain->sym;            /* Current symbol  */
    z = mychain->state;          /* Current state   */
    znext = (mychain+1)->state;  /* Next state      */
    
    /* Find probabilities of changing current state z_i to some state k for all k */
    /* (except state 0 (INIT) and state n - 1 (END)                               */
    for (k = 1, g_sum = 0; k < dc_num_states - 1; k++) {

	ind0 = z == k ? 1.0 : 0.0;
	ind1 = (z == zprev && znext == k) ? 1.0 : 0.0;
	ind2 = (zprev == k && z == znext) ? 1.0 : 0.0;
	ind3 = (zprev == k && k == znext) ? 1.0 : 0.0;

        g_k = ((((float)*CcurrdeviceEMIT(k,a)) + dc_beta_e - ind0) / (((float)d_gibbs_counts_states[k]) - ind0 + dc_alphabet_size * dc_beta_e)) *
	      (((((float)*CcurrdeviceTRANS(zprev,k)) + dc_beta_t - ind0 - ind1) * (((float)*CcurrdeviceTRANS(k,znext)) - ind0 - ind2 + ind3 + dc_beta_t)) /
              (((float)d_gibbs_counts_states[k]) - ind0 + dc_num_states * dc_beta_t));

	g_sum += g_k;
        weight_ptr[numthreads * k] = g_sum;
    }

    /* Do a binary search for the first element in weight_list    */
    /* larger than cointoss. This diverges, but is still slightly */
    /* faster than a linear (more coalesced) search.              */

    cointoss = curand_uniform(&globalrandstate[idorig]) * g_sum;

    for (low = 1, high = dc_num_states - 2; low != high; ) {
    	mid = (low + high) / 2;
    	if (weight_ptr[numthreads * mid] <= cointoss) {
    	    low = mid + 1;
    	} else {
    	    high = mid;
    	}
    }
    k = high;

    /* Update current counts */

    if (k != z) {
	atomicAdd(CcurrdeviceTRANS(zprev,k),1);
	atomicSub(CcurrdeviceTRANS(zprev,z),1);

	atomicAdd(CcurrdeviceTRANS(k,znext),1);
	atomicSub(CcurrdeviceTRANS(z,znext),1);

	atomicAdd(CcurrdeviceEMIT(k,a),1);
	atomicSub(CcurrdeviceEMIT(z,a),1);

	atomicAdd(&d_gibbs_counts_states[k],1);
	atomicSub(&d_gibbs_counts_states[z],1);

	mychain->state = k;
    }
}

#define PROB double

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

/* nvcc really is a C++ compiler (i.e. g++), so all interaction with C code needs to be declared extern "C" */

extern "C" {
    extern int g_alphabet_size, g_gibbs_burnin, g_maxiterations;
    extern struct gibbs_state_chain *gibbs_init_fsm(struct observations *o, int num_states, int alphabet_size, int *obslen);
    extern struct gibbs_state_chain *gibbs_init_hmm(struct observations *o, int num_states, int alphabet_size, int *obslen);
    double gibbs_sampler_cuda_fsm(struct wfsa *fsm, struct observations *o, double beta, int num_states, int maxiter, int burnin, int lag);
    double gibbs_sampler_cuda_hmm(struct hmm *hmm, struct observations *o, double beta_t, double beta_e, int num_states, int maxiter, int burnin, int lag);
    extern struct wfsa *gibbs_counts_to_wfsa(struct wfsa *fsm, unsigned int *gibbs_sampled_counts, unsigned int *gibbs_counts_sampled_states, int alphabet_size, int num_states, double beta, double ANbeta);
    extern struct hmm *gibbs_counts_to_hmm(struct hmm *hmm, unsigned int *gibbs_sampled_counts_trans, unsigned int *gibbs_sampled_counts_emit, unsigned int *gibbs_counts_sampled_states, int alphabet_size, int num_states, double beta_t, double beta_e);
    extern PROB loglikelihood_all_observations_fsm(struct wfsa *fsm, struct observations *o);
    extern PROB loglikelihood_all_observations_hmm(struct hmm *hmm, struct observations *o);
}

double gibbs_sampler_cuda_fsm(struct wfsa *fsm, struct observations *o, double beta, int num_states, int maxiter, int burnin, int lag) {
    int alphabet_size, obslen, i,j, bdimeven, tdimeven, bdimodd, tdimodd, budim, tudim, samplecount;
    struct gibbs_state_chain *chain;

    unsigned int *gibbs_counts, *gibbs_sampled_counts, *gibbs_counts_sampled_states, *gibbs_counts_states, *resamplable, numthreadseven, numthreadsodd, updatethreads, chainlength;
    unsigned int *d_gibbs_counts, *d_gibbs_sampled_counts, *d_gibbs_counts_states, *d_resamplable;
    float *d_weight_list, fANbeta, fbeta;
    double ANbeta;
    struct gibbs_state_chain *d_chain;

    curandState *d_random_state;

    cudaSetDevice(0);
    /* We don't use shared memory, might as well use cache */
    cudaThreadSetCacheConfig(cudaFuncCachePreferL1);

    alphabet_size = g_alphabet_size + 1; /* Use extra symbol for end-of-word (#) */

    /* Init chain and counts locally */
    chain = gibbs_init_fsm(o, num_states, g_alphabet_size, &obslen);
    gibbs_counts = (unsigned int *) calloc(num_states * num_states * alphabet_size, sizeof(unsigned int));
    gibbs_sampled_counts = (unsigned int *) calloc(num_states * num_states * alphabet_size, sizeof(unsigned int));
    gibbs_counts_states = (unsigned int *) calloc(num_states, sizeof(unsigned int));
    gibbs_counts_sampled_states = (unsigned int *) calloc(num_states, sizeof(unsigned int));

    for (i = 0; i < obslen-1; i++) {
	Ccurr( (chain+i)->state , (chain+i)->sym, (chain+i+1)->state )++;
	gibbs_counts_states[(chain+i)->state]++;
    }
    resamplable = (unsigned int *) malloc(obslen * sizeof(unsigned int));
    /* Create array that indexes only the resamplable states              */
    /* That is, initial states (0 and states with incoming # not included */
    for (i = 0, j = 0; j < obslen; j++) {
	if (j == 0 || (chain+j-1)->sym == g_alphabet_size) { /* Don't resample "initial" states, i.e. first */
	    continue;                                        /* state in chain, or states preceded by #     */
	} else {
	    resamplable[i] = j;
	    i++;
	}
    }
    chainlength = i;

    //    fprintf(stderr, "CUDA: planning to launch %i threads\n", chainlength);

    /* Init constants */
    fbeta = (float)beta;
    fANbeta = alphabet_size * num_states * beta;
    cudaMemcpyToSymbol(dc_alphabet_size, &alphabet_size, sizeof(int));
    cudaMemcpyToSymbol(dc_num_states, &num_states, sizeof(int));
    cudaMemcpyToSymbol(dc_ANbeta, &fANbeta, sizeof(float));
    cudaMemcpyToSymbol(dc_beta, &fbeta, sizeof(float));

    /* Move chain and counts to device */
    cudaMalloc(&d_gibbs_counts, num_states * num_states * alphabet_size * sizeof(unsigned int));
    cudaMalloc(&d_gibbs_sampled_counts, num_states * num_states * alphabet_size * sizeof(unsigned int));
    cudaMalloc(&d_gibbs_counts_states, num_states * sizeof(unsigned int));
    cudaMalloc(&d_chain, obslen * sizeof(struct gibbs_state_chain));
    cudaMalloc(&d_resamplable, chainlength * sizeof(unsigned int));
    cudaMalloc(&d_weight_list, num_states * (chainlength / 2 + (chainlength % 2)) * sizeof(float));

    cudaMemcpyToSymbol(dc_weight_list, &d_weight_list, sizeof(float *));
    cudaMemcpyToSymbol(dc_resamplable, &d_resamplable, sizeof(unsigned int *));

    cudaMemcpy(d_gibbs_counts, gibbs_counts, num_states * num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_gibbs_sampled_counts, gibbs_sampled_counts, num_states * num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_gibbs_counts_states, gibbs_counts_states, num_states * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_chain, chain, obslen * sizeof(struct gibbs_state_chain), cudaMemcpyHostToDevice);
    cudaMemcpy(d_resamplable, resamplable, chainlength * sizeof(unsigned int), cudaMemcpyHostToDevice);

    /* We need to sample odd and even states in the chain separately to avoid */
    /* messing up the counts. So for each iteration we launch two consecutive */
    /* kernels, one to resample even-numbered states, and one to resample     */
    /* odd-numbered ones.                                                     */

    /* Even thread number = |chain|/2 + (|chain| % 2)  */
    /* Odd thread number = |chain|/2                   */
    /* Even threads access chain location: id*2        */
    /* Odd threads access chain location: id*2+1       */

    tdimeven = THREADS_PER_BLOCK;
    tdimodd = THREADS_PER_BLOCK;
    numthreadseven = chainlength/2 + (chainlength % 2);
    numthreadsodd =  chainlength/2 ;
    bdimeven = (int) ceil((double) numthreadseven / (double) tdimeven);
    bdimodd = (int) ceil((double) numthreadsodd / (double) tdimeven);
        
    if (bdimeven < 2) { tdimeven = ALIGN_UP(numthreadseven, 32); }
    if (bdimodd < 2)  { tdimodd = ALIGN_UP(numthreadsodd, 32);   }

    //fprintf(stderr,"Dimensions: EVEN: [%i x %i] ODD: [%i x %i]\n",bdimeven,tdimeven, bdimodd, tdimodd);

    /* Number of threads to do total count updates */
    updatethreads = num_states * num_states * alphabet_size;

    tudim = THREADS_PER_BLOCK;
    budim = (int) ceil((double) updatethreads / (double) tudim);
    if (budim < 2) { tudim = ALIGN_UP(updatethreads, 32); }
    /* fprintf(stderr,"Update dimensions: [%i x %i]\n",budim,tudim); */
    /* fprintf(stderr,"Update dimensions EVEN: [%i x %i]\n",bdimeven,tdimeven); */
    /* fprintf(stderr,"Update dimensions ODD: [%i x %i]\n",bdimodd,tdimodd); */
    /* fprintf(stderr,"EVEN: %i ODD: %i\n",numthreadseven, numthreadsodd); */

    /* fprintf(stderr,"Initing random: %i x %i\n", bdimeven*2, tdimeven); */
    /* fflush(stdout); */
    cudaMalloc(&d_random_state, (bdimeven * 2 * tdimeven) * sizeof(curandState));
    CudaCheckError();

    init_random_generator<<<bdimeven*2,tdimeven>>>(d_random_state, time(NULL));
    CudaCheckError();
   
    for (i = 0, samplecount = 0; i < g_maxiterations; i++) {
	gibbs_sampler_kernel_fsm<<<bdimeven,tdimeven>>>(numthreadseven, d_chain, d_gibbs_counts_states, d_gibbs_counts, d_random_state, 0);
	CudaCheckError();
	cudaDeviceSynchronize();
	gibbs_sampler_kernel_fsm<<<bdimodd,tdimodd>>>(numthreadsodd, d_chain, d_gibbs_counts_states, d_gibbs_counts, d_random_state, 1);
	CudaCheckError();
	cudaDeviceSynchronize();
	if (i > burnin && (i - burnin) % lag == 0) {
	    gibbs_sampler_update_kernel_fsm<<<budim,tudim>>>(updatethreads, d_gibbs_counts, d_gibbs_sampled_counts);
	    CudaCheckError();
	    cudaDeviceSynchronize();
	    samplecount++;
	}
	/* For parallel computing tests */
	if (i > 0 && (i == 10 || i == 100 || i % 1000 == 0)) {
	    if (i == 10 || i == 100)
		lag = 10;
	    burnin = i;

	    //fprintf(stderr, "Iteration: %i  Samples collected: %i\n", i, samplecount);
	    cudaMemcpy(gibbs_sampled_counts, d_gibbs_sampled_counts, num_states * num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyDeviceToHost);	    
	    ANbeta = (double) (alphabet_size * num_states * beta);
	    fsm = gibbs_counts_to_wfsa(fsm, gibbs_sampled_counts, gibbs_counts_sampled_states, alphabet_size, num_states, beta, ANbeta);
	    fprintf(stderr, "%i\t%.17g\n", i, loglikelihood_all_observations_fsm(fsm, o));
	    for (j = 0; j <  num_states * num_states * alphabet_size; j++) {
		gibbs_sampled_counts[j] = 0;
	    }
	    cudaMemcpy(d_gibbs_sampled_counts, gibbs_sampled_counts, num_states * num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyHostToDevice);
	}
    }

    /* Move collected counts back to host mem */
    cudaMemcpy(gibbs_sampled_counts, d_gibbs_sampled_counts, num_states * num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyDeviceToHost);

    /* Build WFSA from collected counts */
    ANbeta = (double) (alphabet_size * num_states * beta);
    fsm = gibbs_counts_to_wfsa(fsm, gibbs_sampled_counts, gibbs_counts_sampled_states, alphabet_size, num_states, beta, ANbeta);
    
    cudaFree(d_gibbs_counts);
    cudaFree(d_gibbs_counts_states);
    cudaFree(d_chain);
    cudaFree(d_resamplable);
    cudaFree(d_weight_list);
    cudaFree(d_random_state);
    cudaDeviceReset();
    return(loglikelihood_all_observations_fsm(fsm, o));
}

double gibbs_sampler_cuda_hmm(struct hmm *hmm, struct observations *o, double beta_t, double beta_e, int num_states, int maxiter, int burnin, int lag) {
    int alphabet_size, obslen, i,j, bdimeven, tdimeven, bdimodd, tdimodd, budim, tudim, samplecount;
    struct gibbs_state_chain *chain;

    unsigned int *gibbs_counts_trans, *gibbs_counts_emit, *gibbs_sampled_counts_trans, *gibbs_sampled_counts_emit, *gibbs_counts_sampled_states, *gibbs_counts_states, *resamplable, numthreadseven, numthreadsodd, updatethreads, chainlength;
    unsigned int *d_gibbs_counts_trans, *d_gibbs_counts_emit, *d_gibbs_sampled_counts_trans, *d_gibbs_sampled_counts_emit, *d_gibbs_counts_states, *d_resamplable;
    float *d_weight_list, fbeta_t, fbeta_e;
    struct gibbs_state_chain *d_chain;

    curandState *d_random_state;

    alphabet_size = g_alphabet_size;
    cudaSetDevice(0);
    /* We don't use shared memory, might as well use cache */
    cudaThreadSetCacheConfig(cudaFuncCachePreferL1);

    /* Init chain and counts locally */
    chain = gibbs_init_hmm(o, num_states, g_alphabet_size, &obslen);
    gibbs_counts_trans = (unsigned int *) calloc(num_states * num_states, sizeof(unsigned int));
    gibbs_counts_emit = (unsigned int *) calloc(num_states * alphabet_size, sizeof(unsigned int));
    gibbs_sampled_counts_trans = (unsigned int *) calloc(num_states * num_states, sizeof(unsigned int));
    gibbs_sampled_counts_emit = (unsigned int *) calloc(num_states * alphabet_size, sizeof(unsigned int));
    gibbs_counts_states = (unsigned int *) calloc(num_states, sizeof(unsigned int));
    gibbs_counts_sampled_states = (unsigned int *) calloc(num_states, sizeof(unsigned int));

    for (i = 0; i < obslen-1; i++) {
	CcurrHMMtrans( (chain+i)->state, (chain+i+1)->state )++;
	gibbs_counts_states[(chain+i)->state]++;
    }
    for (i = 0; i < obslen-1; i++) {
      if ((chain+i)->sym >= 0) {
	  CcurrHMMemit( (chain+i)->state, (chain+i)->sym)++;
      }
    }
    resamplable = (unsigned int *) malloc(obslen * sizeof(unsigned int));
    /* Create array that indexes only the resamplable states              */
    /* That is, initial states (0 and states with incoming # not included */
    for (i = 0, j = 0; j < obslen; j++) {
	if ((chain+j)->sym < 0)
	    continue;            /* Don't resample INIT or END states */
	resamplable[i] = j;
	i++;
    }
    chainlength = i;

    //fprintf(stderr, "CUDA: planning to launch %i threads\n", chainlength);

    /* Init constants */
    fbeta_t = (float)beta_t;
    fbeta_e = (float)beta_e;

    cudaMemcpyToSymbol(dc_alphabet_size, &alphabet_size, sizeof(int));
    cudaMemcpyToSymbol(dc_num_states, &num_states, sizeof(int));
    cudaMemcpyToSymbol(dc_beta_t, &fbeta_t, sizeof(float));
    cudaMemcpyToSymbol(dc_beta_e, &fbeta_e, sizeof(float));

    /* Move chain and counts to device */
    cudaMalloc(&d_gibbs_counts_trans, num_states * num_states * sizeof(unsigned int));
    cudaMalloc(&d_gibbs_counts_emit, num_states * alphabet_size * sizeof(unsigned int));
    cudaMalloc(&d_gibbs_sampled_counts_trans, num_states * num_states * sizeof(unsigned int));
    cudaMalloc(&d_gibbs_sampled_counts_emit, num_states * alphabet_size * sizeof(unsigned int));
    cudaMalloc(&d_gibbs_counts_states, num_states * sizeof(unsigned int));
    cudaMalloc(&d_chain, obslen * sizeof(struct gibbs_state_chain));
    cudaMalloc(&d_resamplable, chainlength * sizeof(unsigned int));
    cudaMalloc(&d_weight_list, num_states * (chainlength / 2 + (chainlength % 2)) * sizeof(float));

    cudaMemcpyToSymbol(dc_weight_list, &d_weight_list, sizeof(float *));
    cudaMemcpyToSymbol(dc_resamplable, &d_resamplable, sizeof(unsigned int *));

    cudaMemcpy(d_gibbs_counts_trans, gibbs_counts_trans, num_states * num_states * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_gibbs_counts_emit, gibbs_counts_emit, num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_gibbs_sampled_counts_trans, gibbs_sampled_counts_trans, num_states * num_states * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_gibbs_sampled_counts_emit, gibbs_sampled_counts_emit, num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_gibbs_counts_states, gibbs_counts_states, num_states * sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_chain, chain, obslen * sizeof(struct gibbs_state_chain), cudaMemcpyHostToDevice);
    cudaMemcpy(d_resamplable, resamplable, chainlength * sizeof(unsigned int), cudaMemcpyHostToDevice);

    /* We need to sample odd and even states in the chain separately to avoid */
    /* messing up the counts. So for each iteration we launch two consecutive */
    /* kernels, one to resample even-numbered states, and one to resample     */
    /* odd-numbered ones.                                                     */

    /* Even thread number = |chain|/2 + (|chain| % 2)  */
    /* Odd thread number = |chain|/2                   */
    /* Even threads access chain location: id*2        */
    /* Odd threads access chain location: id*2+1       */

    tdimeven = THREADS_PER_BLOCK;
    tdimodd = THREADS_PER_BLOCK;
    numthreadseven = chainlength/2 + (chainlength % 2);
    numthreadsodd =  chainlength/2 ;
    bdimeven = (int) ceil((double) numthreadseven / (double) tdimeven);
    bdimodd = (int) ceil((double) numthreadsodd / (double) tdimeven);
        
    if (bdimeven < 2) { tdimeven = ALIGN_UP(numthreadseven, 32); }
    if (bdimodd < 2)  { tdimodd = ALIGN_UP(numthreadsodd, 32);   }

    //fprintf(stderr,"Dimensions: EVEN: [%i x %i] ODD: [%i x %i]\n",bdimeven,tdimeven, bdimodd, tdimodd);

    /* Number of threads to do total count updates */
    updatethreads = num_states * (num_states > alphabet_size ? num_states : alphabet_size);

    tudim = THREADS_PER_BLOCK;
    budim = (int) ceil((double) updatethreads / (double) tudim);
    if (budim < 2) { tudim = ALIGN_UP(updatethreads, 32); }
    //fprintf(stderr,"Update dimensions: [%i x %i]\n",budim,tudim);
    cudaMalloc(&d_random_state, (bdimeven * 2 * tdimeven) * sizeof(curandState));
    CudaCheckError();
    init_random_generator<<<bdimeven*2,tdimeven>>>(d_random_state, time(NULL));
   
    for (i = 0, samplecount = 0; i < g_maxiterations; i++) {
	gibbs_sampler_kernel_hmm<<<bdimeven,tdimeven>>>(numthreadseven, d_chain, d_gibbs_counts_states, d_gibbs_counts_trans, d_gibbs_counts_emit, d_random_state, 0);
	cudaDeviceSynchronize();
	gibbs_sampler_kernel_hmm<<<bdimodd,tdimodd>>>(numthreadsodd, d_chain, d_gibbs_counts_states, d_gibbs_counts_trans, d_gibbs_counts_emit, d_random_state, 1);
	cudaDeviceSynchronize();
	if (i > burnin && (i - burnin) % lag == 0) {
	    gibbs_sampler_update_kernel_hmm<<<budim,tudim>>>(updatethreads, d_gibbs_counts_trans, d_gibbs_counts_emit, d_gibbs_sampled_counts_trans, d_gibbs_sampled_counts_emit, num_states * num_states, num_states * alphabet_size);
	    cudaDeviceSynchronize();
	    samplecount++;
	}
	if (i > 0 && i % 100 == 0) {
	    //fprintf(stderr, "Iteration: %i  Samples collected: %i\n", i, samplecount);
	}
	/* For parallel computing tests */
	if (i > 0 && (i == 10 || i == 100 || i % 1000 == 0)) {
	    if (i == 10 || i == 100)
		lag = 10;
	    burnin = i;

	    //fprintf(stderr, "Iteration: %i  Samples collected: %i\n", i, samplecount);
	    cudaMemcpy(gibbs_sampled_counts_trans, d_gibbs_sampled_counts_trans, num_states * num_states * sizeof(unsigned int), cudaMemcpyDeviceToHost);
	    cudaMemcpy(gibbs_sampled_counts_emit, d_gibbs_sampled_counts_emit, num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyDeviceToHost);

	    hmm = gibbs_counts_to_hmm(hmm, gibbs_sampled_counts_trans, gibbs_sampled_counts_emit, gibbs_counts_sampled_states, alphabet_size, num_states, beta_t, beta_e);
	    fprintf(stderr, "%i\t%.17g\n", i, loglikelihood_all_observations_hmm(hmm, o));
	    for (j = 0; j <  num_states * num_states; j++) {
		gibbs_sampled_counts_trans[j] = 0;
	    }
	    for (j = 0; j <  num_states * alphabet_size; j++) {
		gibbs_sampled_counts_emit[j] = 0;
	    }
	    cudaMemcpy(d_gibbs_sampled_counts_trans, gibbs_sampled_counts_trans, num_states * num_states * sizeof(unsigned int), cudaMemcpyHostToDevice);
	    cudaMemcpy(d_gibbs_sampled_counts_emit, gibbs_sampled_counts_emit, num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyHostToDevice);
	}
    }

    /* Move collected counts back to host mem */
    cudaMemcpy(gibbs_sampled_counts_trans, d_gibbs_sampled_counts_trans, num_states * num_states * sizeof(unsigned int), cudaMemcpyDeviceToHost);
    cudaMemcpy(gibbs_sampled_counts_emit, d_gibbs_sampled_counts_emit, num_states * alphabet_size * sizeof(unsigned int), cudaMemcpyDeviceToHost);

    /* Build HMM from collected counts */
    hmm = gibbs_counts_to_hmm(hmm, gibbs_sampled_counts_trans, gibbs_sampled_counts_emit, gibbs_counts_sampled_states, alphabet_size, num_states, beta_t, beta_e);
    
    cudaFree(d_gibbs_counts_trans);
    cudaFree(d_gibbs_counts_emit);
    cudaFree(d_gibbs_counts_states);
    cudaFree(d_chain);
    cudaFree(d_resamplable);
    cudaFree(d_weight_list);
    cudaFree(d_random_state);
    cudaDeviceReset();
    return(loglikelihood_all_observations_hmm(hmm, o));
}
