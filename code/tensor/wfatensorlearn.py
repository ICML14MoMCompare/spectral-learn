import ctypes as C
import numpy as np
from math import log,e
import hankelmatrixcreator
import sys
import time
import iohelpers
import math
import modelconversion
from scipy.sparse.linalg import lsqr
from scipy import sparse
import copy

DEBUG = False
VERBOSE = True

FAILURE_CONST = -100000.0

class TensorWFA:

    def __init__(self, n_symbols):
        self.n_symbols = n_symbols
        self.hankels_learnt = False

    def _estimate_hankels(self, data, prefixdict, suffixdict):
        self.h_pands,self.symbol_hankels,self.hp_pandsigma,self.hbar_pands,self.hbar_pandsigma,self.hbar_sigmaands = hankelmatrixcreator.construct_tensor_hankels(data, prefixdict, suffixdict, self.n_symbols, 100)
        
        if VERBOSE:
            print "Finished Hankel estimation"

    def compute_XY(self, hbar_pands, hbar_pandsigma, hbar_sigmaands, symbol_hankels, num_symbols, num_components):

        if VERBOSE:
            print "Constructing Q matrices.."
        self.qp,_,self.qs = sparse.linalg.svds(hbar_pands,num_components)
        self.qs_constructed = True

        qp = sparse.csr_matrix(np.mat((self.qp[:,0:num_components]).T))
        qs = sparse.csr_matrix(np.mat((self.qs[0:num_components, :])))
        qp.eliminate_zeros()
        qp.prune()
        qs.eliminate_zeros()
        qs.prune()

        if VERBOSE:
            print "Computing N.."
        hsquiggleps = qp*hbar_pands*qs.T
        hsquiggleps = sparse.csr_matrix(np.linalg.inv(hsquiggleps.A))
        N = qs.T*(hsquiggleps)*qp
        

        N.eliminate_zeros()
        N.prune()

        if VERBOSE:
            print "Computing X.."
        Xs = hbar_sigmaands*N*hbar_pandsigma
        X = np.empty((num_symbols, num_symbols), order='F', dtype=np.float64)
        X[:,:] = Xs.A
        self.rank = np.linalg.matrix_rank(X)
        if self.rank < num_components:
            print "Rank Deficient!!"
            return [],[], False

        if VERBOSE:
            print "Computing Y..."

        leftside = hbar_sigmaands*N
        leftside.eliminate_zeros()
        leftside.prune()
        rightside = N*hbar_pandsigma
        rightside.eliminate_zeros()
        rightside.prune()

        Y = np.empty((num_symbols, num_symbols, num_symbols))
        for sym in range(num_symbols):
            Y[:,sym,:] = (leftside*(symbol_hankels[sym])*rightside).A

        return X,Y, True

    def learn_tensor(self,  data, prefixdict, suffixdict, num_components):

        starttime = time.clock()

        if not self.hankels_learnt:
            begintime = time.clock()
            self._estimate_hankels(data, prefixdict, suffixdict)
            self.inittime = time.clock()-begintime
            self.hankels_learnt = True

        #adding aliases with "self" prefix for readability.
        h_pands = self.h_pands
        symbol_hankels = self.symbol_hankels
        hp_pandsigma = self.hp_pandsigma
        hbar_pands = self.hbar_pands
        hbar_pandsigma = self.hbar_pandsigma  
        hbar_sigmaands = self.hbar_sigmaands 
        num_symbols = self.n_symbols

        X,Y, success = self.compute_XY(hbar_pands, hbar_pandsigma, hbar_sigmaands, symbol_hankels, num_symbols, num_components)

        if not success:
            return success

        if VERBOSE:
            print "Performing tensor decomposition.."

        lib = C.CDLL('./cpp/libTensor.so')
        f_learn_tensor = lib.learn_tensor
        f_learn_tensor.argtypes = [np.ctypeslib.ndpointer(dtype = np.float64), np.ctypeslib.ndpointer(dtype = np.float64), C.c_int, np.ctypeslib.ndpointer(dtype = np.float64), np.ctypeslib.ndpointer(dtype = np.float64), C.c_int, C.c_int]
        f_learn_tensor.restype = C.c_double

        Otilde = np.empty((num_symbols, num_components), order='F', dtype=np.float64)
        gamma = np.empty((num_components), order='F', dtype=np.float64)

        res = f_learn_tensor(X,Y,num_components,Otilde,gamma,num_symbols, 1)

        if res == FAILURE_CONST:
            return False

        if VERBOSE:
            print "Building unnormalized model.."

        Otilde = np.mat(Otilde)

        # hbar_pandsigma = np.mat(self.hbar_pandsigma.toarray())
        # hbar_sigmaands = np.mat(self.hbar_sigmaands.toarray())
        # hbar_pands = np.mat(self.hbar_pands.toarray())

        Otildepinv = np.linalg.pinv(Otilde)
        spOtildepinv = sparse.csr_matrix(Otildepinv)
        spOtildepinv.eliminate_zeros()
        spOtildepinv.prune()
        Otildep = hbar_pandsigma*Otildepinv.T
        Otildes = Otildepinv*hbar_sigmaands

        alphatilde = np.mat(Otildep[0,:]).T
        stopvec = np.mat(Otildes[:,0])

        spOtildeppinv = sparse.csr_matrix(np.linalg.pinv(Otildep))
        spOtildeppinv.eliminate_zeros()
        spOtildeppinv.prune()

        Ttilde = self._compute_T(spOtildeppinv, num_components)

        Dgamma = np.mat(np.diag(gamma))
        Dgamma = Dgamma*Dgamma
        Ds = spOtildeppinv*hp_pandsigma*spOtildepinv.T
        Dspinv = np.mat(Ds.A)
        Ds = np.linalg.pinv(np.mat(Ds.A))
        Dsigma = np.eye(alphatilde.shape[0])-np.diag(stopvec)
        Dsigmapinv = np.linalg.pinv(Dsigma)
        Beta = Ds*(np.linalg.pinv(Ttilde))*Dgamma*stopvec

        for i in range(stopvec.shape[0]):
            stopvec[i] = Beta[i]/(1+Beta[i])

        alpha = np.empty((num_components), order='F', dtype=np.float64)
        alpha[:] = (alphatilde.T*Dsigmapinv*Dspinv).A

        O = np.empty((num_components, num_symbols), order='F', dtype=np.float64)
        O[:,:] = (Otilde*Dsigma).A.T
        T = np.empty((num_components, num_components), order='F', dtype=np.float64)
        T[:] = (Ds*Ttilde*Dspinv*Dsigmapinv).A


        if DEBUG:
            print "O before PNFA projection: ", O
            print "T before PNFA projection: ", T
            print "ainf before simplex projection: ", stopvec
            print "a before simplex projection: ", alpha

        O,T,alpha, stopvec = self._project_to_probabilities(O,T,alpha, stopvec.T.A[0], num_components, num_symbols)
        stopvec = np.mat(stopvec).T
        alpha = np.mat(alpha).T
        O = O.T

        if DEBUG:
            print "O after PNFA projection: ", O
            print "T after PNFA projection: ", T
            print "ainf after simplex projection: ", stopvec
            print "a after simplex projection: ", alpha


        self.initvec, self.ainf, self.wordstopvec, self.As = self.convert_hmm_to_wfa(O, T, alpha, stopvec)
        self.a = self.initvec.copy()
        self.buildtime = time.clock() - starttime

        return True

    def _compute_sparse_pinv(self, sparsemat, rank):
        u,s,v = sparse.linalg.svds(sparsemat, rank)
        u = np.mat(u)
        s = np.mat(np.diag(s))
        v = np.mat(v)
        pinv =  v.T*np.linalg.inv(s)*u.T
        return np.mat(pinv)

    def _compute_T(self, Otildepinv, num_components):

        A = (Otildepinv*self.h_pands).T
        A.eliminate_zeros()
        A.prune()
        B = (Otildepinv*self.hbar_pands).T
        B.eliminate_zeros()
        B.prune()
        T = np.empty((num_components, num_components))

        for i in range(num_components):
            T[i,:] = lsqr(A,B[:,i].A)[0]

        return T

    def _project_to_probabilities(self, O, T, alpha, ainf, num_components, num_symbols):
        lib = C.CDLL('../cpp/libTensor.so')
        f_simplex_proj = lib.do_simplex_projection
        f_simplex_proj.argtypes = [np.ctypeslib.ndpointer(dtype = np.float64), C.c_double, C.c_int]
        f_simplex_proj.restype = C.c_int

        f_probmat_proj = lib.do_probmat_projection
        f_probmat_proj.argtypes = [np.ctypeslib.ndpointer(dtype = np.float64), np.ctypeslib.ndpointer(dtype = np.float64), C.c_int, C.c_int]
        f_probmat_proj.restype = C.c_int

        f_simplex_proj(alpha, 1.0, num_components)
        zeros = np.zeros((num_components), order='F', dtype=np.float64)
        probO = np.zeros((num_components, num_symbols+1), order='F', dtype=np.float64)
        probO[:,0:num_symbols] = O
        probO[:,num_symbols:num_symbols+1] = np.mat(ainf).T
        f_probmat_proj(probO, zeros, num_components, num_symbols)
        O = probO[:,0:num_symbols]
        ainf = probO[:,num_symbols:num_symbols+1]
        ainf = ainf.T
        f_probmat_proj(T, zeros, num_components, num_components)

        return O, T, alpha, ainf


    def get_symbol_prediction(self):

        predictedsymbol = -1
        maxscore = np.finfo(float).eps
        for symbol in range(self.n_symbols):
            symbolscore = self.get_obs_prob(symbol)
            if symbolscore > maxscore:
                predictedsymbol = symbol
                maxscore = symbolscore

        stopscore = float(self.a.T*self.wordstopvec)

        if stopscore > maxscore:
            predictedsymbol = self.n_symbols

        return predictedsymbol


    def get_WER(self, testdata):
        errors = 0
        numpredictions = 0

        for seq in testdata:
            for obs in seq:
                numpredictions += 1
                predsymbol = self.get_symbol_prediction()
                self.update(obs)
                if predsymbol != obs:
                    errors += 1
            predsymbol = self.get_symbol_prediction()
            numpredictions += 1
            if predsymbol != self.n_symbols:
                errors += 1
            self.reset()

        return float(errors)/float(numpredictions)


    #provides average log-likelihood score
    def score(self, data): 
        loglike = 0

        for seq in data:
            seqloglike = 0
            self.reset()
            for obs in seq:
                seqloglike = seqloglike + log(self.get_obs_prob(obs))
                self.update(obs)
            loglike += seqloglike

        return loglike/(float(len(data)))

    #updates a/start/state vector after seeing symbol
    def update(self, obs):
        bomat = self.As[obs]
        numerator = self.a.T*bomat
        denom = numerator*self.ainf

        self.a = (numerator/denom).T

    #resets state vector
    def reset(self):
        self.a = self.initvec.copy()

    #returns the probablilty of a particular observation given current state
    def get_obs_prob(self, obs):
        prob = (self.a.T)*(self.As[obs])*self.ainf
        prob = min(prob,1)
        prob = max(prob,np.finfo(float).eps)
        return prob

    #returns the probability of an entire sequence, or "word"
    def get_word_prob(self,seq):
        seqprob = 0
        for obs in seq:
            prob = self.get_obs_prob(obs)
            if prob <= np.finfo(float).eps:
                return np.finfo(float).eps
            seqprob += log(prob)
            self.update(obs)

        endprob = float(self.a.T*self.wordstopvec)
        if endprob <= np.finfo(float).eps:
            return np.finfo(float).eps
       
        seqprob += log(endprob)
        self.reset()

        if not math.isnan(seqprob):
            return e**seqprob
        else:
            return np.finfo(float).eps


    def scorepautomac(self, testdata, truprobs):
        modelprobs = np.zeros((len(truprobs)))
        probsum = 0
        i = 0
        for seq in testdata:
            prob = self.get_word_prob(seq)
            modelprobs[i] = prob
            probsum += prob
            i += 1

        modelprobs /= float(probsum)

        i = 0
        scoresum = 0
        for truprob in truprobs:
            if modelprobs[i] < np.finfo(float).eps:
                modelprobs[i] =  np.finfo(float).eps
            scoresum += truprob*log(modelprobs[i],2)
            i += 1

        return 2.0**(-1.0*float(scoresum))

    def get_perplexity(self, testdata):
        modelprobs = np.zeros((len(testdata)))
        probsum = 0
        i = 0
        for seq in testdata:
            prob = self.get_word_prob(seq)
            modelprobs[i] = prob
            probsum += prob
            i += 1

        scoresum = 0
        for i in range(len(modelprobs)):
            if modelprobs[i] < np.finfo(float).eps:
                modelprobs[i] =  np.finfo(float).eps
            scoresum += log(modelprobs[i],2)

        scoresum /= float(len(testdata))

        return 2.0**(-1.0*float(scoresum))

    def convert_hmm_to_wfa(self, O, T, alpha, stopvec):

        As = {}
        for symbol in range(O.shape[0]):
            As[symbol] = np.mat(np.diag(O[symbol,:]))*T

        return alpha, np.mat(np.ones(T.shape[0])).T, stopvec, As


if __name__ == '__main__':

    PAUTOMACPATH = "/home/williamleif/Dropbox/icml2014-experiments/datasets/PAutomaC-competition_sets/"
    RESULTS_DIR = "/home/williamleif/Dropbox/icml2014-experiments/results/tensor/"
    MODEL_DIR = "/home/williamleif/Dropbox/icml2014-experiments/models/"

    metric = sys.argv[1]
    problem = sys.argv[2]
    n_symbols = sys.argv[3]
    n_symbols = int(n_symbols)

    if problem != "tree" and problem != "timeseries":

        traindata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.train")
        testdata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.test")

        if metric == "KL":
            groundtruth = iohelpers.parse_groundtruth_file(PAUTOMACPATH+problem+".pautomac_solution.txt")
        else:
            validdata = traindata[15000:20000]
            traindata = traindata[0:15000]

      
        maxbasissize = int(sys.argv[4])

        prefixdict, suffixdict = hankelmatrixcreator.top_k_string_bases(traindata,maxbasissize,n_symbols)
        basislength = len(prefixdict)

        bestsize = 0
        avruntime = 0
        nummodelsmade = 0
        bestscore = 0
        wfa = TensorWFA(n_symbols)
        begintime = time.clock()
        for i in range(3, 6):

            success = wfa.learn_tensor(traindata, prefixdict, suffixdict, i)

            if i == 3:
                inittime = wfa.inittime

            if not success:
                break

            if metric == "WER":
                score = wfa.get_WER(validdata)
            else:
                score = wfa.scorepautomac(testdata,groundtruth)

            if bestsize == 0:
                bestscore = score
                bestsize = i
                bestwfa = copy.deepcopy(wfa)
            elif score < bestscore and abs(score-1000) > 0.1:
                bestscore = score
                bestsize = i
                bestwfa = copy.deepcopy(wfa)

            print "Model size: ", i, " Score: ", score
            avruntime += wfa.buildtime
            nummodelsmade += 1

        # if metric == "WER":
        #     bestscore = bestwfa.get_WER(testdata)

        # iohelpers.write_results(RESULTS_DIR+"tensor-pautomac="+problem+"-"+metric+".txt", problem,"size= "+str(bestsize)+", basis size="+str(basislength), metric, bestscore, avruntime/float(nummodelsmade))
        # iohelpers.write_pnfa_model(MODEL_DIR+"tensor-"+str(bestsize)+"-pautomac="+problem+"-"+metric+".fsm", bestwfa)

        runtime = time.clock()-begintime
        fp = open("/home/williamleif/Dropbox/icml2014-experiments/results/runtimes/tensor", "w")
        fp.write("Init time: "+str(inittime)+" Runtime: "+str(runtime))

    else:
        
        RESULTS_DIR = "/home/williamleif/Dropbox/icml2014-experiments/results/real/"

        if problem == "tree":
            traindata = iohelpers.parse_file("/home/williamleif/Dropbox/icml2014-experiments/datasets/treebankdata.obs")
            validdata = traindata[0:5000]
            testdata = traindata[5000:10000]
            traindata = traindata[10000:len(traindata)]

        maxbasissize = int(sys.argv[4])

        prefixdict, suffixdict = hankelmatrixcreator.top_k_string_bases(traindata,maxbasissize,n_symbols)
        basislength = len(prefixdict)

        bestsize = 0
        avruntime = 0
        nummodelsmade = 0
        bestscore = 0
        wfa = TensorWFA(n_symbols)
        for i in range(3, n_symbols+1):

            success = wfa.learn_tensor(traindata, prefixdict, suffixdict, i)

            if not success:
                break

            if metric == "WER":
                score = wfa.get_WER(validdata)
            else:
                score = wfa.get_perplexity(validdata)

            if bestsize == 0:
                bestscore = score
                bestsize = i
                bestwfa = copy.deepcopy(wfa)
            elif score < bestscore and abs(score-1000) > 0.1:
                bestscore = score
                bestsize = i
                bestwfa = copy.deepcopy(wfa)

            print "Model size: ", i, " Score: ", score
            avruntime += wfa.buildtime
            nummodelsmade += 1

        if metric == "WER":
            bestscore = bestwfa.get_WER(testdata)
        else:
            bestscore = bestwfa.get_perplexity(testdata)

        iohelpers.write_results(RESULTS_DIR+"tensor-"+metric+".txt", problem,"size= "+str(bestsize)+", basis size="+str(basislength), metric, bestscore, 0)
        iohelpers.write_pnfa_model(MODEL_DIR+"tensor-"+str(bestsize)+"-pautomac="+problem+"-"+metric+".fsm", bestwfa)

