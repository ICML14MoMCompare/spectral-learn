import numpy as np
from math import log,e
import hankelmatrixcreator
import sppy
from sppy import linalg
import sys
import time
import iohelpers
import math
from scipy import sparse
import copy

VERBOSE = True

class SparseSpectralWFA:

    def __init__(self, n_symbols, prefixdict, suffixdict, basis_length, maxDim=80):

        self.n_symbols = n_symbols
        self.prefixdict = prefixdict
        self.suffixdict = suffixdict
        self.basis_length = basis_length
        self.buildtime = 0
        self.Asyms = []
        self.maxDim = maxDim


    def fit(self,data, n_components, substring=True):
        starttime = time.clock()
        self.n_components = n_components
        self.Asyms = []
        for i in range(self.n_symbols):
            self.Asyms.append(np.empty((n_components,n_components)))


        #getting hankel matrix estimates
        if VERBOSE:
            print "Getting hankel estimates.."
        if substring:
            hp,hs, H, self.symbol_hankels = hankelmatrixcreator.construct_substring_hankels(data,self.prefixdict, self.suffixdict, self.n_symbols)
        else:
            hp, hs,H, self.symbol_hankels = hankelmatrixcreator.construct_string_hankels(data,self.prefixdict, self.suffixdict, self.n_symbols)
        self.hp = hp
        self.hs = hs

        if VERBOSE:
            print "Performing SVD..."
        u,s,v = linalg.rsvd(H,self.maxDim)
        self.u = u.copy()
        self.v = v.copy()
        self.s = s.copy()

        self.inittime = time.clock()-starttime
    
        u = np.mat(u[:,0:n_components])
        v = np.mat(v[:,0:n_components]).T
        s = np.mat(np.diag(s[0:n_components]))

        leftopmat = sppy.csarray(((np.linalg.inv(s))*u.T).A)
        rightopmat = sppy.csarray(v.T.A)

        if VERBOSE:
            print "Computing operators..."
        for sym in range(self.n_symbols):
        	self.Asyms[sym] = np.mat((leftopmat.dot(self.symbol_hankels[sym]).dot(rightopmat)).toarray())

        #computing stopping and starting vectors
        spU = sppy.csarray(u.A)
        spV = sppy.csarray(v.A)
        self.ainf = sppy.csarray((np.linalg.inv(s)).A).dot(spU.T).dot(hp)
        self.a = spV.dot(hs)
        self.ainf = np.mat(self.ainf.toarray().T[0]).T
        self.a = np.mat(self.a.toarray().T[0]).T

        A = np.mat(np.zeros((self.n_components, self.n_components)))
        for sym in range(self.n_symbols):
            A += self.Asyms[sym]
        if substring:
            #transform starting vector since substring not prefix counts used in estimation
            self.a = (self.a.T*(np.mat(np.eye(self.n_components))-A)).T
            self.wordainf = (np.mat(np.eye(self.n_components))-A)*self.ainf
        else:
            self.wordainf = self.ainf.copy()
            self.ainf = np.linalg.inv(np.eye(self.n_components)-A)*self.wordainf

        #save a deep copy of the start vector for resets
        self.inita = self.a.copy()

        self.buildtime = time.clock()-starttime

    def resize(self,data,n_components, substring=True):
        starttime = time.clock()
        self.n_components = n_components
        hp = self.hp
        hs = self.hs

        for sym in range(self.n_symbols):
            self.Asyms[sym] = np.zeros((self.n_components, self.n_components))
    
    	u = np.mat(self.u[:,0:n_components])
        v = np.mat(self.v[:,0:n_components]).T
        s = np.mat(np.diag(self.s[0:n_components]))

        leftopmat = sppy.csarray(((np.linalg.inv(s))*u.T).A)
        rightopmat = sppy.csarray(v.T.A)

        for sym in range(self.n_symbols):
        	self.Asyms[sym] = np.mat((leftopmat.dot(self.symbol_hankels[sym]).dot(rightopmat)).toarray())

        #computing stopping and starting vectors
        spU = sppy.csarray(u.A)
        spV = sppy.csarray(v.A)
        self.ainf = sppy.csarray((np.linalg.inv(s)).A).dot(spU.T).dot(hp)
        self.a = spV.dot(hs)
        self.ainf = np.mat(self.ainf.toarray().T[0]).T
        self.a = np.mat(self.a.toarray().T[0]).T

        A = np.mat(np.zeros((self.n_components, self.n_components)))
        for sym in range(self.n_symbols):
            A += self.Asyms[sym]
        if substring:
            #transform starting vector since substring not prefix counts used in estimation
            self.a = (self.a.T*(np.mat(np.eye(self.n_components))-A)).T
            self.wordainf = (np.mat(np.eye(self.n_components))-A)*self.ainf
        else:
            self.wordainf = self.ainf.copy()
            self.ainf = np.linalg.inv(np.eye(self.n_components)-A)*self.wordainf

        #save a deep copy of the start vector for resets
        self.inita = self.a.copy()

        self.buildtime = time.clock()-starttime


    def get_symbol_prediction(self):

        predictedsymbol = -1
        maxscore = np.finfo(float).eps
        for symbol in range(self.n_symbols):
            symbolscore = self.get_obs_prob(symbol)
            if symbolscore > maxscore:
                predictedsymbol = symbol
                maxscore = symbolscore

        stopscore = float(self.a.T*self.wordainf)

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
        bomat = self.Asyms[obs]
        numerator = self.a.T*bomat
        denom = numerator*self.ainf

        self.a = (numerator/denom).T

    #resets state vector
    def reset(self):
        self.a = self.inita.copy()

    #returns the probablilty of a particular observation given current state
    def get_obs_prob(self, obs):
        prob = (self.a.T)*(self.Asyms[obs])*self.ainf
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

        endprob = float(self.a.T*self.wordainf)
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



if __name__ == '__main__':

    PAUTOMACPATH = "/home/williamleif/Dropbox/icml2014-experiments/datasets/PAutomaC-competition_sets/"
    RESULTS_DIR = "/home/williamleif/Dropbox/icml2014-experiments/results/"

    esttype = sys.argv[1]
    metric = sys.argv[2]
    problem = sys.argv[3]
    n_symbols = sys.argv[4]
    n_symbols = int(n_symbols)
    maxbasissize =int(sys.argv[5])  

    if esttype == "substring":
        substring = True
        RESULTS_DIR += 'spectral-substring/'
    elif esttype == 'string':
        substring = False
        RESULTS_DIR += 'spectral-string/'
    else:
        sys.stderr.write("Estimation type must be substring or string.")
        sys.exit(0)


    if problem != "tree" and problem != "timeseries":

        traindata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.train")
        testdata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.test")

        if metric == "KL":
            groundtruth = iohelpers.parse_groundtruth_file(PAUTOMACPATH+problem+".pautomac_solution.txt")
        else:
            validdata = traindata[15000:20000]
            traindata = traindata[0:15000]


        if substring:
            basisdict = hankelmatrixcreator.top_k_basis(traindata,maxbasissize,n_symbols, 4)
            basislength = len(basisdict)
        else:
            prefixdict, suffixdict = hankelmatrixcreator.top_k_string_bases(traindata,maxbasissize,n_symbols)
            basislength = len(prefixdict)

        if substring:
            wfa = SparseSpectralWFA(n_symbols, basisdict, basisdict, 4)
        else:
            wfa = SparseSpectralWFA(n_symbols, prefixdict,suffixdict, 100)
        bestsize = 0
        avruntime = 0
        nummodelsmade = 0
        sizes = []
        begintime = time.clock()
        sizes.extend(range(10,31,10))
        for i in sizes:
            #if i == 5:
            if i == 10:
                wfa.fit(traindata, i,substring)
                inittime = wfa.inittime
            else:
            	wfa.resize(traindata,i, substring)

            if metric == "WER":
                score = wfa.get_WER(validdata)
            else:
                score = wfa.scorepautomac(testdata,groundtruth)

            if bestsize == 0:
                bestscore = score
                bestsize = i
            elif score < bestscore:
                bestscore = score
                bestsize = i

            print "Model size: ", i, " Score: ", score
            avruntime += wfa.buildtime
            nummodelsmade += 1

        runtime = time.clock()-begintime

        if bestsize == 5:
            bestsize = 10;

        for i in range(int(bestsize)-9, int(bestsize)+10):

        		
            wfa.resize(traindata,i, substring)

            if metric == "WER":
                score = wfa.get_WER(testdata)
            else:
                score = wfa.scorepautomac(testdata,groundtruth)

            if bestsize == 0:
                bestwfa = copy.deepcopy(wfa)
                bestsize = i
            elif score < bestscore or math.isnan(bestscore):
                bestscore = score
                bestsize = i
            avruntime += wfa.buildtime
            nummodelsmade += 1

            print "Model size: ", i, " Score: ", score

        if metric == "WER":
            wfa.resize(traindata,i,substring)
            bestscore = wfa.get_WER(testdata)

        iohelpers.write_results(RESULTS_DIR+"spectral-"+esttype+"-pautomac="+problem+"-"+metric+".txt", problem, esttype+", "+"size= "+str(bestsize)+", basis size="+str(basislength), metric, bestscore, avruntime/float(nummodelsmade))

    else:

        RESULTS_DIR = "/home/williamleif/Dropbox/icml2014-experiments/results/real/"

        if problem == "tree":
            traindata = iohelpers.parse_file("/home/williamleif/Dropbox/icml2014-experiments/datasets/treebankdata.obs")
            validdata = traindata[0:5000]
            testdata = traindata[5000:10000]
            traindata = traindata[10000:len(traindata)]

        if substring:
            basisdict = hankelmatrixcreator.top_k_basis(traindata,maxbasissize,n_symbols, 4)
            basislength = len(basisdict)
        else:
            prefixdict, suffixdict = hankelmatrixcreator.top_k_string_bases(traindata,maxbasissize,n_symbols)
            basislength = len(prefixdict)

        if substring:
            wfa = SparseSpectralWFA(n_symbols, basisdict, basisdict, 4)
        else:
            wfa = SparseSpectralWFA(n_symbols, prefixdict,suffixdict, 100)
        bestsize = 0
        avruntime = 0
        nummodelsmade = 0
        sizes = [5]
        sizes.extend(range(10,71,10))
        for i in sizes:
            if i == 5:
                wfa.fit(traindata, i,substring)
            else:
                wfa.resize(traindata,i, substring)

            if metric == "WER":
                score = wfa.get_WER(validdata)
            else:
                score = wfa.get_perplexity(validdata)

            if bestsize == 0:
                bestscore = score
                bestsize = i
            elif score < bestscore:
                bestscore = score
                bestsize = i

            print "Model size: ", i, " Score: ", score
            avruntime += wfa.buildtime
            nummodelsmade += 1

        if bestsize == 5:
            bestsize = 10;

        for i in range(int(bestsize)-9, int(bestsize)+10):

                
            wfa.resize(traindata,i, substring)

            if metric == "WER":
                score = wfa.get_WER(testdata)
            else:
                score = wfa.get_perplexity(validdata)

            if bestsize == 0:
                bestscore = score
                bestsize = i
            elif score < bestscore or math.isnan(bestscore):
                bestscore = score
                bestsize = i
            avruntime += wfa.buildtime
            nummodelsmade += 1

            print "Model size: ", i, " Score: ", score

        wfa.resize(traindata,bestsize, substring)
        if metric == "WER":
            bestscore = wfa.get_WER(testdata)
        else:
            bestscore = wfa.get_perplexity(testdata)


        iohelpers.write_results(RESULTS_DIR+"spectral-"+esttype+"-"+metric+".txt", problem, esttype+", "+"size= "+str(bestsize)+", basis size="+str(basislength), metric, bestscore, 0)



