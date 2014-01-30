import numpy as np
from math import log,sqrt,e
import iohelpers
import sys
import copy
import time
import subprocess
import math

LIKEITSCONST = 100.0
MAXVALIDITERSCONST = 1000.0
NOTIMPROVCONST = 3
IMPROVEPS = 0.1
VERBOSE = True

class EmWFA:

	def __init__(self, n_components, n_symbols):
		self.n_components = n_components
		self.n_symbols = n_symbols

	def _parse_model(self, trebaoutput):

		trebaoutputlines = trebaoutput.split("\n")
		self.As = {}
		for symbol in range(self.n_symbols):
			self.As[symbol] = np.empty((self.n_components, self.n_components))

		self.wordstopvec = np.zeros((self.n_components))

		for line in trebaoutputlines:
			entries = line.split(' ')

			if len(entries) == 4:
				source_state = int(entries[0])
				target_state = int(entries[1])


				symbol = int(entries[2])
				prob = float(entries[3])
				Asym = self.As[symbol]
				Asym[source_state, target_state] = prob
			elif len(entries) == 2:
				source_state = int(entries[0])
				stopprob = float(entries[1])
				self.wordstopvec[source_state] = stopprob

		self.stopvec = np.mat(np.ones((self.n_components))).T
		self.wordstopvec = np.mat(self.wordstopvec).T
		self.initvec = np.mat(np.zeros((self.n_components))).T
		self.initvec[0] = 1
		self.a = self.initvec.copy()

	def realfit(self, obsFile, testdata, validdata):
		#must use wall time since subprocess being called..
		starttime = time.time()
		fp = open(".emtemp.fsm", "w")

		likeits = int(LIKEITSCONST/float(self.n_components))
		numvalidits = int(MAXVALIDITERSCONST/float(self.n_components))

		trebaoutput = subprocess.check_output(["treba","--train=bw","--initialize="+str(self.n_components), "--max-delta=0.5", "--restarts=5,"+str(likeits), "--max-iter=1", obsFile])
		fp.write(trebaoutput)
		fp.flush()

		lastwer = 0
		num_not_improv = 0
		for i in range(numvalidits):
			trebaoutput = subprocess.check_output(["treba","--train=bw","--file=.emtemp.fsm", "--max-iter="+str(likeits), obsFile])
			fp = open(".emtemp.fsm", "w")
			fp.write(trebaoutput)
			fp.flush()
			fp.close()
			self._parse_model(trebaoutput)

			kl = self.get_perplexity(validdata[0:1000])
			wer = self.get_WER(validdata[0:1000])

			if VERBOSE:
				print "WER: ", wer, " KL:", kl

			if lastwer == 0:
				lastwer = wer
				lastkl = kl
				bestwer = wer
				bestkl = kl
				self.bestwfa = copy.deepcopy(self)
			elif lastkl - kl < IMPROVEPS and lastwer - wer < IMPROVEPS:
				num_not_improv += 1
			else:
				num_not_improv = 0
				if lastkl > kl:
					lastkl = kl
				if lastwer > wer:
					lastwer = wer

				if kl < bestkl:
					bestkl = kl
				if wer < bestwer:
					bestwer = wer
					self.bestwfa = copy.deepcopy(self)

			if num_not_improv >= NOTIMPROVCONST:
				break

		self.kl = bestkl
		self.wer = bestwer

		self.buildtime = time.time()-starttime
		self._parse_model(trebaoutput)

	def fit(self, obsFile, testdata, validdata,groundtruth):
		#must use wall time since subprocess being called..
		starttime = time.time()
		fp = open(".emtemp.fsm", "w")

		likeits = int(LIKEITSCONST/float(self.n_components))
		numvalidits = int(MAXVALIDITERSCONST/float(self.n_components))

		trebaoutput = subprocess.check_output(["treba","--train=bw","--initialize="+str(self.n_components), "--max-delta=0.5", "--restarts=5,"+str(likeits), "--max-iter=1", obsFile])
		fp.write(trebaoutput)
		fp.flush()

		lastwer = 0
		num_not_improv = 0
		for i in range(numvalidits):
			trebaoutput = subprocess.check_output(["treba","--train=bw","--file=.emtemp.fsm", "--max-iter="+str(likeits), obsFile])
			fp = open(".emtemp.fsm", "w")
			fp.write(trebaoutput)
			fp.flush()
			fp.close()
			self._parse_model(trebaoutput)

			kl = self.scorepautomac(testdata,groundtruth)
			wer = self.get_WER(validdata)

			if VERBOSE:
				print "WER: ", wer, " KL:", kl

			if lastwer == 0:
				lastwer = wer
				lastkl = kl
				bestwer = wer
				bestkl = kl
				self.bestwfa = copy.deepcopy(self)
			elif lastkl - kl < IMPROVEPS and lastwer - wer < IMPROVEPS:
				num_not_improv += 1
			else:
				num_not_improv = 0
				if lastkl > kl:
					lastkl = kl
				if lastwer > wer:
					lastwer = wer

				if kl < bestkl:
					bestkl = kl
				if wer < bestwer:
					bestwer = wer
					self.bestwfa = copy.deepcopy(self)

			if num_not_improv >= NOTIMPROVCONST:
				break

		self.kl = bestkl
		self.wer = bestwer

		self.buildtime = time.time()-starttime
		self._parse_model(trebaoutput)


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
			if modelprobs[i] < np.finfo(float).eps:
				modelprobs[i] =  np.finfo(float).eps            
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

if __name__ == "__main__":

	PAUTOMACPATH = "/home/williamleif/Dropbox/icml2014-experiments/datasets/PAutomaC-competition_sets/"
	RESULTS_DIR = "/home/williamleif/Dropbox/icml2014-experiments/results/em/"

	problem = sys.argv[1]
	n_symbols = sys.argv[2]
	n_symbols = int(n_symbols)

	if problem != "tree" and problem != "timeseries":
		testdata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.test")
		validdata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.train")[15000:20000]
		iohelpers.clean_pautomacfile_for_em(PAUTOMACPATH+problem+".pautomac.train", PAUTOMACPATH+problem+".pautomac.em")

		groundtruth = iohelpers.parse_groundtruth_file(PAUTOMACPATH+problem+".pautomac_solution.txt")

		avruntime = 0
		nummodelsmade = 0
		klsize = 0
		sizes = [5]
		sizes.extend(range(10,41,10))
		for i in sizes:
			wfa = EmWFA(i, n_symbols)

			wfa.fit(PAUTOMACPATH+problem+".pautomac.em",testdata,validdata, groundtruth)

			kl = wfa.kl
			wer = wfa.bestwfa.get_WER(testdata)

			if i == 30:
				thirtykl = wfa.kl
				thirtywer = wfa.wer

			if klsize == 0:
				bestkl = kl
				bestwer = wer
				klsize = i
				wersize = i
			else:
				if kl < bestkl:
					bestkl = kl
					klsize = i
				if wer < bestwer:
					bestwer = wer
					wersize = i

			avruntime += wfa.buildtime
			nummodelsmade += 1

			print "Model size: ", i, " KL: ", kl, " WER: ", wer

		iohelpers.write_results(RESULTS_DIR+"em-pautomac="+problem+".txt", problem,  "KL size:"+str(klsize)+" WER size: "+str(wersize)+" 30 KL: "+str(thirtykl)+" 30 WER: "+str(thirtywer), "KL, WER", str(bestkl)+","+str(bestwer), avruntime/float(nummodelsmade))

	else:
		if problem == "tree":
			traindata = iohelpers.parse_file("/home/williamleif/Dropbox/icml2014-experiments/datasets/treebankdata.obs")
			validdata = traindata[0:5000]
			testdata = traindata[5000:10000]
			traindata = traindata[10000:len(traindata)]

		fp = open("treetemp.obs", "w")
		fp.write("\n".join([" ".join([str(j) for j in i]) for i in traindata]))
		fp.close()

		
		avruntime = 0
		nummodelsmade = 0
		klsize = 0
		sizes = [5]
		sizes.extend(range(10,41,10))
		for i in sizes:
			wfa = EmWFA(i, n_symbols)

			wfa.realfit("treetemp.obs",testdata,validdata)

			kl = wfa.kl
			wer = wfa.bestwfa.get_WER(testdata)

			if i == 30:
				thirtykl = wfa.kl
				thirtywer = wfa.wer

			if klsize == 0:
				bestkl = kl
				bestwer = wer
				klsize = i
				wersize = i
			else:
				if kl < bestkl:
					bestkl = kl
					klsize = i
				if wer < bestwer:
					bestwer = wer
					wersize = i

			avruntime += wfa.buildtime
			nummodelsmade += 1

			print "Model size: ", i, " KL: ", kl, " WER: ", wer

		iohelpers.write_results(RESULTS_DIR+"em-pautomac="+problem+".txt", problem,  "KL size:"+str(klsize)+" WER size: "+str(wersize)+" 30 KL: "+str(thirtykl), "KL, WER", str(bestkl)+","+str(bestwer), avruntime/float(nummodelsmade))