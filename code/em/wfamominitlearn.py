from wfaemlearn import EmWFA
import numpy as np 
import time
import sys
import iohelpers
import hankelmatrixcreator
import subprocess
import copy

NUMVALIDITERS = 100
NOTIMPROVCONST = 3
KLIMPROVEPS = 0.1
WERIMPROVEEPS = 0.001

VERBOSE = True

class MomentInitEmWFA(EmWFA):

	def __init__(self, n_components, n_symbols, initFile):
		self.n_components = n_components
		self.n_symbols = n_symbols
		self.initFile = initFile

	def fit(self, validata, obsFile):
		starttime = time.time()

		lastwer = 0
		num_not_improv = 0
		fp = open(".mominittemp.fsm", "w")
		for i in range(NUMVALIDITERS):

			if not fp.closed:
				fp.close()

			if i == 0:
				trebaoutput = subprocess.check_output(["treba","--train=bw","--file="+self.initFile, "--max-iter=5",obsFile])
			else:
				trebaoutput = subprocess.check_output(["treba","--train=bw","--file=.mominittemp.fsm", "--max-iter=5", obsFile])
			fp = open(".mominittemp.fsm", "w")
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
			elif lastkl - kl < KLIMPROVEPS and lastwer - wer < WERIMPROVEEPS:
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

			if num_not_improv >= NOTIMPROVCONST:
				break

		self.kl = bestkl
		self.wer = bestwer
		self.buildtime = time.time()-starttime
		self._parse_model(trebaoutput)

	def realfit(self, obsFile, validdata):
		starttime = time.time()

		lastwer = 0
		num_not_improv = 0
		fp = open(".mominittemp.fsm", "w")
		for i in range(NUMVALIDITERS):

			if not fp.closed:
				fp.close()

			if i == 0:
				trebaoutput = subprocess.check_output(["treba","--train=bw","--file="+self.initFile, "--max-iter=5",obsFile])
			else:
				trebaoutput = subprocess.check_output(["treba","--train=bw","--file=.mominittemp.fsm", "--max-iter=5", obsFile])
			fp = open(".mominittemp.fsm", "w")
			fp.write(trebaoutput)
			fp.flush()
			fp.close()
			self._parse_model(trebaoutput)

			kl = self.get_perplexity(validdata)
			wer = self.get_WER(validdata)

			if VERBOSE:
				print "WER: ", wer, " KL:", kl

			if lastwer == 0:
				lastwer = wer
				lastkl = kl
				bestwer = wer
				bestkl = kl
				self.bestwfa = copy.deepcopy(self)
			elif lastkl - kl < KLIMPROVEPS and lastwer - wer < WERIMPROVEEPS:
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


if __name__ == "__main__":

	PAUTOMACPATH = "/home/williamleif/Dropbox/icml2014-experiments/datasets/PAutomaC-competition_sets/"
	MODELPATH = "/home/williamleif/Dropbox/icml2014-experiments/models/"
	RESULTS_DIR = "/home/williamleif/Dropbox/icml2014-experiments/results/momentinit/"

	problem = sys.argv[1]
	n_symbols = sys.argv[2]
	n_symbols = int(n_symbols)
	modelfile = sys.argv[3]
	dimension = sys.argv[4]
	dimension = int(dimension)

	if problem != "tree":

		testdata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.test")
		validdata = iohelpers.parse_file(PAUTOMACPATH+problem+".pautomac.train")
		iohelpers.clean_pautomacfile_for_em(PAUTOMACPATH+problem+".pautomac.train", PAUTOMACPATH+problem+".pautomac.em")

		groundtruth = iohelpers.parse_groundtruth_file(PAUTOMACPATH+problem+".pautomac_solution.txt")

		wfa = MomentInitEmWFA(dimension, n_symbols, MODELPATH+modelfile)

		wfa.fit(PAUTOMACPATH+problem+".pautomac.em")

		wer = wfa.wer
		kl = wfa.kl

	else:

		traindata = iohelpers.parse_file("/home/williamleif/Dropbox/icml2014-experiments/datasets/treebankdata.obs")
		validdata = traindata[0:5000]
		testdata = traindata[5000:10000]
		traindata = traindata[10000:len(traindata)]

		fp = open("treetemp.obs", "w")
		fp.write("\n".join([" ".join([str(j) for j in i]) for i in traindata]))
		fp.close()

		wfa = MomentInitEmWFA(dimension, n_symbols, MODELPATH+modelfile)

		wfa.realfit("treetemp.obs",validdata)

		kl = wfa.kl
		wer = wfa.bestwfa.get_WER(testdata)



	print "KL: ", kl, " WER: ", wer

	iohelpers.write_results(RESULTS_DIR+"em"+modelfile, problem, "size= "+str(dimension), "KL,WER", str(kl)+", "+str(wer), wfa.buildtime)