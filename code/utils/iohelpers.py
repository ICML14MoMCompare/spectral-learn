def parse_file(file):
    fp = open(file,"r")

    fp.readline()
    data = [[int(i) for i in line.split()[1:]] for line in fp]

    return data

def parse_groundtruth_file(file):
    fp = open(file,"r")
    fp.readline()
    truprobs = [float(line) for line in fp]
    return truprobs

def write_results(file, problem_descript, algorithm_descript, metric, score, runtime):
	f = open(file, 'w')

	f.write("Problem: "+problem_descript+"\n")
	f.write("Algorithm: "+algorithm_descript+"\n")
	f.write(metric+": "+str(score)+"\n")
	f.write("Runtime: "+str(runtime)+"\n")
	f.close()

def write_pnfa_model(file, pnfa):

	pnfafile = open(file, 'w')
	numstates = pnfa.As[0].shape[0]
	for source_state in range(numstates):
		for target_state in range(numstates):
			for symbol in range(pnfa.n_symbols):
				Asym = pnfa.As[symbol]
				prob = float(Asym[source_state, target_state])
				pnfafile.write(str(source_state)+" "+str(target_state)+" "+str(symbol)+" "+str(prob)+"\n")

		pnfafile.write(str(source_state)+" "+str(float(pnfa.wordstopvec[source_state,0]))+"\n")

def clean_pautomacfile_for_em(filein, fileout, truncate=False):

	infp = open(filein, 'r')
	outfp = open(fileout, "w")

	infp.readline()

	if truncate:
		for line in infp.readlines()[0:15000]:
			entrys = line.split(' ')
			entrys = entrys[1:]
			line = ' '.join(entrys)
			outfp.write(line)
	else:
		for line in infp.readlines():
			entrys = line.split(' ')
			entrys = entrys[1:]
			line = ' '.join(entrys)
			outfp.write(line)
