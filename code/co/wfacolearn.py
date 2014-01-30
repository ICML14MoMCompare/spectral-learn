import numpy as np
from math import log,sqrt,e
import ctypes as C
import hankelmatrixcreator
import iohelpers
import sys
import copy
import time
import sppy
import math

class CoWFA:

    def __init__(self, hankelmat, symbol_hankels,p):
        self.hankelmat = hankelmat
        self.symbol_hankels = symbol_hankels
        self.As = {}
        self.n_symbols = len(symbol_hankels.keys())
        self.p = p
        self.buildtime = 0


    def fit(self, tau, maxK, probabilistic=True):

        MAXDIM = 50
        starttime = time.clock()

        #load C library
        lib = C.CDLL('./cpp/libADMM.so')
        numprefs = self.hankelmat.shape[0]
        numsuffs = self.hankelmat.shape[1]

        #NOTE use of FORTRAN or 'F' order when passing things to C. The Eigen library uses this order.
        hankel_sigma = np.empty((self.n_symbols*numprefs,numsuffs),dtype=np.float64,order='F')
        i = 0
        #convert collection of H_symbol matrices into stacked H_sigma matrix
        for mat in self.symbol_hankels.values():
            hankel_sigma[numprefs*i:numprefs*(i+1),0:numsuffs] = mat
            i += 1

        #initialize large stacked A matrix and set ainf to the prefix estimate of h.
        A = np.zeros((self.n_symbols*numprefs,numsuffs),dtype=np.float64,order='F')
        ainf = np.zeros(numprefs,dtype=np.float64,order='F')
        ainf[:] = self.p

        if probabilistic:
            f_learn_pnfa = lib.admm_pnfa_learn
            f_learn_pnfa.argtypes = [np.ctypeslib.ndpointer(dtype = np.float64),np.ctypeslib.ndpointer(dtype = np.float64),np.ctypeslib.ndpointer(dtype = np.float64),np.ctypeslib.ndpointer(dtype = np.float64),C.c_int,C.c_int, C.c_double, C.c_int, C.c_int]
            f_learn_pnfa.restype = C.c_int
            #call C ADMM routine
            f_learn_pnfa(self.hankelmat,hankel_sigma,ainf,A,numprefs,self.n_symbols,tau,maxK, MAXDIM)
        else:
            f_learn_wa = lib.admm_wfa_learn
            f_learn_wa.argtypes = [np.ctypeslib.ndpointer(dtype = np.float64),np.ctypeslib.ndpointer(dtype = np.float64),np.ctypeslib.ndpointer(dtype = np.float64),C.c_int,C.c_int, C.c_double, C.c_int, C.c_int]
            f_learn_wa.restype = C.c_int
            f_learn_wa(self.hankelmat,hankel_sigma,A,numprefs,self.n_symbols,tau,maxK, MAXDIM)
        
        self.wordstopvec = np.mat(ainf).T
        self.stopvec = np.mat(np.ones((numprefs))).T
        
        self.initvec = np.mat(np.zeros((numprefs))).T
        self.initvec[0] = 1
        self.a = self.initvec.copy()
        i = 0
        insigindices = A < np.finfo(float).eps
        A[insigindices] = 0
        #extract single A_symbol operators from stacked solution
        for symbol in self.symbol_hankels.keys():
            self.As[symbol] = np.mat(A[numprefs*i:numprefs*(i+1),0:numsuffs])
            i += 1

        self.buildtime = time.clock()-starttime

    #returns average log-likelihood of sequence in the data
    def score(self, data):
        loglike = 0

        for seq in data:
            seqloglike = 0
            self.reset()
            for obs in seq:
                seqloglike = seqloglike + log(self.get_obs_prob(obs))
                self.update(obs)
            seqloglike = seqloglike + log(self.a.T*self.wordstopvec)
            loglike += seqloglike

        return loglike/(float(len(data)))

    #returns the probability assigned to a word.
    def get_word_prob(self,seq):

        seqprob = 0
        for obs in seq:
            prob = self.get_obs_prob(obs)
            if prob <= 0:
                return np.finfo(float).eps
            seqprob += log(prob)
            self.update(obs)

        endprob = float(self.a.T*self.wordstopvec)
        if endprob <= 0:
            return np.finfo(float).eps
       
        seqprob += log(endprob)
        self.reset()

        if not math.isnan(seqprob):
            return e**seqprob
        else:
            return np.finfo(float).eps

    def get_symbol_prediction(self):

        predictedsymbol = -1
        maxscore = np.finfo(float).eps
        for symbol in range(len(self.symbol_hankels)):
            symbolscore = self.get_obs_prob(symbol)
            if symbolscore > maxscore:
                predictedsymbol = symbol
                maxscore = symbolscore

        stopscore = float(self.a.T*self.wordstopvec)

        if stopscore > maxscore:
            predictedsymbol = len(self.symbol_hankels)


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
            if predsymbol != len(self.symbol_hankels):
                errors += 1
            self.reset()

        return float(errors)/float(numpredictions)

    #updates normalized internal state
    def update(self, symbol):
        Amat = np.mat(self.As[symbol])
        numerator = self.a.T*Amat
        denom = numerator*self.stopvec

        self.a = (numerator/denom).T

    #resets to inital state
    def reset(self):
        self.a = self.initvec.copy()

    #gets probability of single symbol/observation conditioned on current internal state
    def get_obs_prob(self, symbol):
        prob = (self.a.T)*(self.As[symbol])*self.stopvec
        prob = min(prob,1)
        prob = max(prob,np.finfo(float).eps)
        return prob

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
            scoresum += truprob*log(modelprobs[i],2)
            i += 1

        return 2.0**(-1.0*float(scoresum))

    def check_if_pnfa(self):
        negativeentries = False

        for Amat in self.As.values():
            if not np.all(Amat > 0):
                negativeentries = True
                break

        residuals = np.zeros((self.hankelmat.shape[0]))

        for i in range(self.hankelmat.shape[0]):
            rowsum = 0
            for Amat in self.As.values():
                rowsum += np.sum(Amat[i,:])
            rowsum += float(self.wordstopvec[i][0])

            residuals[i] = rowsum - 1

        return negativeentries, residuals

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


if __name__ == '__main__':

    PAUTOMACPATH = "/home/williamleif/Dropbox/icml2014-experiments/datasets/PAutomaC-competition_sets/"
    RESULTS_DIR = "/home/williamleif/Dropbox/icml2014-experiments/results/"
    MODEL_DIR = "/home/williamleif/Dropbox/icml2014-experiments/models/"

    wfatype = sys.argv[1]
    metric = sys.argv[2]
    problem = sys.argv[3]
    n_symbols = sys.argv[4]
    n_symbols = int(n_symbols)
    defmaxbasissize = int(sys.argv[5])

    maxbasissize = defmaxbasissize
    if len(sys.argv) == 6:
        maxbasissize = int(sys.argv[5])

    if wfatype == "PNFA":
        probabilistic = True
        RESULTS_DIR += "co-PNFA/"
        shift = 1
    elif wfatype == "WFA":
        probabilistic = False
        RESULTS_DIR += "co-WFA/"
        shift = -1
    else:
        sys.stderr.write("WFA type must be WFA or PNFA\n")
        sys.exit(0)

    if metric != "WER" and metric != "KL":
        sys.stderr.write("Metrix must be WER or KL\n")

    if problem != "tree" and problem != "timeseries":

        traindata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.train")
        testdata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.test")

      
        maxK = 500
        m = 1/sqrt(float(len(traindata)))
        begintime = time.clock()
        basisdict = hankelmatrixcreator.top_k_basis(traindata,maxbasissize,n_symbols, basis_length=4)
        hankelmat,symbolhankels,p = hankelmatrixcreator.construct_hankel_matrices_for_co(traindata,basisdict, n_symbols ,basis_length=4)
        _,s,_ = np.linalg.svd(hankelmat,full_matrices=0)
        hankelsing = s[0]
        maxtaupow = int(log(len(basisdict)*m*hankelsing, 10))
        inittime = time.clock()-begintime

        if metric == "KL":
            groundtruth = iohelpers.parse_groundtruth_file(PAUTOMACPATH+problem+".pautomac_solution.txt")
        else:
            validdata = traindata[15000:20000]
            traindata = traindata[0:15000]

        Cs = []
        for i in range(-4+shift,min(maxtaupow,3)+1+shift):
            Cs.append(5*(10**i))

        besttauorder = 0
        avruntime = 0
        nummodelsmade = 0
        for tauorder in Cs:
            wfa = CoWFA(hankelmat,symbolhankels,p)
            wfa.fit(tauorder , maxK, probabilistic)

            if metric == "WER":
                score = wfa.get_WER(validdata)
            else:
                score = wfa.scorepautomac(testdata,groundtruth)


            if besttauorder == 0:
                bestscore = score
                besttauorder = tauorder
            elif score < bestscore and score != 1000:
                bestscore = score
                besttauorder = tauorder

            print "Tau: ",tauorder, " Score: ", score

            avruntime += wfa.buildtime
            nummodelsmade += 1

        tauorder = besttauorder/5.0
        besttau = 0
        for i in range(1,10):
            tau = i*tauorder
            wfa = CoWFA(hankelmat,symbolhankels,p)
            wfa.fit(tau , maxK, probabilistic)

            if metric == "WER":
                score = wfa.get_WER(testdata)
            else:
                score = wfa.scorepautomac(testdata,groundtruth)

            if besttau == 0:
                bestscore = score
                besttau = tau
                bestwfa = copy.deepcopy(wfa)
            elif score < bestscore and abs(score-1000)<0.1:
                bestscore = score
                besttau = tau
                bestwfa = copy.deepcopy(wfa)

            print "Tau: ",tau, " Score: ", score

            avruntime += wfa.buildtime
            nummodelsmade += 1

        if metric == "WER":
            bestscore = bestwfa.get_WER(testdata)

        iohelpers.write_results(RESULTS_DIR+"co-"+wfatype+"-"+str(maxbasissize)+"-pautomac="+problem+"-"+metric+".txt", problem, wfatype+", "+"tau= "+str(besttau)+", basis size="+str(len(basisdict)), metric, bestscore, avruntime/float(nummodelsmade))

        if probabilistic:
            iohelpers.write_pnfa_model(MODEL_DIR+"co-"+wfatype+"-"+str(maxbasissize)+"-pautomac="+problem+"-"+metric+".fsm", bestwfa)


    else:
        

        RESULTS_DIR = "/home/williamleif/Dropbox/icml2014-experiments/results/real/"

        if problem == "tree":
            traindata = iohelpers.parse_file("/home/williamleif/Dropbox/icml2014-experiments/datasets/treebankdata.obs")
            validdata = traindata[0:5000]
            testdata = traindata[5000:10000]
            traindata = traindata[10000:len(traindata)]


        maxK = 500
        m = 1/sqrt(float(len(traindata)))
        basisdict = hankelmatrixcreator.top_k_basis(traindata,maxbasissize,n_symbols, basis_length=4)
        hankelmat,symbolhankels,p = hankelmatrixcreator.construct_hankel_matrices_for_co(traindata,basisdict, n_symbols ,basis_length=4)
        _,s,_ = np.linalg.svd(hankelmat,full_matrices=0)
        hankelsing = s[0]
        maxtaupow = int(log(len(basisdict)*m*hankelsing, 10))


        Cs = []
        for i in range(-4+shift,min(maxtaupow,3)+1+shift):
            Cs.append(5*(10**i))

        besttauorder = 0
        avruntime = 0
        nummodelsmade = 0
        for tauorder in Cs:
            wfa = CoWFA(hankelmat,symbolhankels,p)
            wfa.fit(tauorder , maxK, probabilistic)

            if metric == "WER":
                score = wfa.get_WER(validdata)
            else:
                score = wfa.get_perplexity(validdata)


            if besttauorder == 0:
                bestscore = score
                besttauorder = tauorder
            elif score < bestscore and score != 1000:
                bestscore = score
                besttauorder = tauorder

            print "Tau: ",tauorder, " Score: ", score

            avruntime += wfa.buildtime
            nummodelsmade += 1

        tauorder = besttauorder/5.0
        besttau = 0
        for i in range(1,10):#should be range to 10
            tau = i*tauorder
            wfa = CoWFA(hankelmat,symbolhankels,p)
            wfa.fit(tau , maxK, probabilistic)

            if metric == "WER":
                score = wfa.get_WER(testdata)
            else:
                score = wfa.get_perplexity(validdata)

            if besttau == 0:
                bestscore = score
                besttau = tau
                bestwfa = copy.deepcopy(wfa)
            elif score < bestscore and abs(score-1000)<0.1:
                bestscore = score
                besttau = tau
                bestwfa = copy.deepcopy(wfa)

            print "Tau: ",tau, " Score: ", score

            avruntime += wfa.buildtime
            nummodelsmade += 1

        if metric == "WER":
            bestscore = bestwfa.get_WER(testdata)
        else:
            bestscore = bestwfa.get_perplexity(testdata)
            
        iohelpers.write_results(RESULTS_DIR+"co-"+wfatype+"-"+str(maxbasissize)+"-"+metric+".txt", problem, wfatype+", "+"tau= "+str(besttau)+", basis size="+str(len(basisdict)), metric, bestscore, 0)

        if probabilistic:
            iohelpers.write_pnfa_model(MODEL_DIR+"co-"+wfatype+"-"+str(maxbasissize)+"-tree"+"-"+metric+".fsm", bestwfa)



