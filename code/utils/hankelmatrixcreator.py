import numpy as np
from scipy import sparse
import sppy

DEBUG = True
EXPECTED_SPARSENESS = 0.0001

def top_k_basis(data, k, num_symbols, basis_length=4):

	count_dict = {}
	#iterating over all sequences in data
	for seq in data:
		for i in range(len(seq)):
			for j in range(i+1,len(seq)):
				if j - i > basis_length:
					break

				subseq = tuple(seq[i:j])

				if subseq in count_dict:
					count_dict[subseq] = count_dict[subseq] + 1
				else:
					count_dict[subseq] = 1

	topk = []
	topk.append(tuple([]))
	topk.extend(sorted(count_dict, key=count_dict.get, reverse=True)[0:k])

	basisdict = {}
	index = 0
	for item in topk:
		basisdict[item] = index
		index += 1

	print "Finished building basis"
	return basisdict

def top_k_string_bases(data, k, num_symbols):
	prefix_count_dict = {}
	suffix_count_dict = {}
	for seq in data:
		for i in range(1,len(seq)+1):
			prefix = tuple(seq[0:i])

			if prefix in prefix_count_dict:
				prefix_count_dict[prefix] = prefix_count_dict[prefix] + 1
			else:
				prefix_count_dict[prefix] = 1

			if i < len(seq):
				suffix = tuple(seq[i:len(seq)])
				if suffix in suffix_count_dict:
					suffix_count_dict[suffix] = suffix_count_dict[suffix] + 1
				else:
					suffix_count_dict[suffix] = 1

	topkprefix = []
	topkprefix.append(tuple([]))
	topkprefix.extend(sorted(prefix_count_dict, key=prefix_count_dict.get, reverse=True)[0:k])

	prefixdict = {}
	index = 0
	for item in topkprefix:
		prefixdict[item] = index
		index += 1

	topksuffix = []
	topksuffix.append(tuple([]))
	topksuffix.extend(sorted(prefix_count_dict, key=prefix_count_dict.get, reverse=True)[0:k])

	suffixdict = {}
	index = 0
	for item in topksuffix:
		suffixdict[item] = index
		index += 1

	print "Finished building basis"
	return prefixdict, suffixdict

def single_symbol_basis(num_symbols):

	basislist = []

	basislist.append(tuple([]))

	for i in range(num_symbols):
		dumlist = []
		dumlist.append(i)

		basislist.append(tuple(dumlist))

	basisdict = {}
	index = 0
	for item in basislist:
		basisdict[item] = index
		index += 1

	return basisdict


def construct_hankel_matrices_for_co(data,basisdict, num_symbols, basis_length=4):

	symbol_hankels = {}

	for i in range(num_symbols):
		symbol_hankels[i] = np.zeros((len(basisdict),len(basisdict)),dtype=np.float64,order='F')

	hankelmat = np.zeros((len(basisdict),len(basisdict)),dtype=np.float64,order='F')
	p = np.zeros((len(basisdict)), dtype=np.float64,order='F')
	prefix_counts = np.zeros(len(basisdict))
	suffix_counts = np.zeros(len(basisdict))

	for seq in data:
		#iterate over prefix start
		for i in range(len(seq)):
			#iterate over suffix start
			for j in range(i,len(seq)+1):
				#if prefix length greater than basis size then exit
				if j - i > basis_length:
					break
				#iterate over suffix end
				for k in range(j,len(seq)+1):
					#if suffix length greater than basis size then exit
					if k - j > basis_length:
						break

					prefix = tuple(seq[i:j])

					#if suffix empty string at end of word then special case
					if k == len(seq) and j == len(seq):
						suffix = tuple([])
					else:
						suffix = tuple(seq[j:k])

					if prefix in basisdict and suffix in basisdict:
						prefixind = basisdict[prefix]
						suffixind = basisdict[suffix]
						hankelmat[prefixind, suffixind] += 1
						prefix_counts[prefixind] += 1
						suffix_counts[suffixind] += 1

						if i == 0 and k == len(seq) and j == len(seq):
							p[prefixind] += 1

					#index must be one smaller when looking at prefix,symbol,suffix tuples.
					if j < len(seq) and k < len(seq):
						sym = seq[j]
						symprefix = tuple(seq[i:j])

						#if symsuffix empty string then special case
						if k + 1 == len(seq) and j + 1 == len(seq):
							symsuffix = tuple([])
						else:
							symsuffix = tuple(seq[j+1:k+1])

						if symprefix in basisdict and symsuffix in basisdict:
							symprefixind = basisdict[symprefix]
							symsuffixind = basisdict[symsuffix]
							symhankelmat = symbol_hankels[sym]
							symhankelmat[symprefixind, symsuffixind] += 1
	
	# hankelmat = np.mat(np.diag(prefix_counts))*np.mat(hankelmat)*np.mat(np.diag(suffix_counts))
	hankelmat /= float(len(data))
	p /= float(len(data))

	for sym in range(num_symbols):
		# symbol_hankels[sym] = np.mat(np.diag(prefix_counts))*np.mat(symbol_hankels[sym])*np.mat(np.diag(suffix_counts))
		symbol_hankels[sym] /= float(len(data))

	return hankelmat,symbol_hankels,hankelmat[0,:]


def construct_substring_hankels(data,prefixdict, suffixdict, n_symbols, basis_length=4):

		kappa = 5

		hankelmat = sppy.csarray((len(prefixdict),len(suffixdict)))
		hankelmat.reserve(int((len(prefixdict)*len(suffixdict))*EXPECTED_SPARSENESS))
		prefix_counts = np.zeros(len(prefixdict))
		suffix_counts = np.zeros(len(prefixdict))
		symbol_hankels = {}
		for sym in range(n_symbols):
		    symbol_hankels[sym] = sppy.csarray((len(prefixdict),len(suffixdict)))
		    symbol_hankels[sym].reserve(int((len(prefixdict)*len(suffixdict))*EXPECTED_SPARSENESS))


		#iterating over sequences
		for seq in data:
		    #iterating over prefix start positions
		    for i in range(len(seq)):
		        #iterating over prefix end positions
		        for j in range(i,len(seq)+1):
		            #break if prefix longer than anything in basis
		            if j - i > basis_length:
		                break
		            for k in range(j,len(seq)+1):
		                #break if suffix longer than anything in basis
		                if k - j > basis_length:
		                    break
		             
		                prefix = tuple(seq[i:j])

		                #if suffix empty string at end of word then special case
		                if k == len(seq) and j == len(seq):
		                    suffix = tuple([])
		                else:
		                    suffix = tuple(seq[j:k])

		                if prefix in prefixdict and suffix in suffixdict:
							prefixind = prefixdict[prefix]
							suffixind = suffixdict[suffix]
							hankelmat[prefixind, suffixind] += 1
							prefix_counts[prefixind] += 1
							suffix_counts[suffixind] += 1

		                if j < len(seq) and k < len(seq):

		                    sym = seq[j]

		                    symprefix = tuple(seq[i:j])
		                    #special case when suffix is just empty string
		                    if k + 1 == len(seq) and j + 1 == len(seq):
		                        symsuffix = tuple([])
		                    else:
		                        symsuffix = tuple(seq[j+1:k+1])

		                    if symprefix in prefixdict and symsuffix in suffixdict:
		                        symprefixind = prefixdict[symprefix]
		                        symsuffixind = suffixdict[symsuffix]
		                        symbol_hankel = symbol_hankels[sym]
		                        symbol_hankel[symprefixind, symsuffixind] += 1



		prefix_counts = np.sqrt(float(len(data))/(prefix_counts+kappa))
		suffix_counts = np.sqrt(float(len(data))/(suffix_counts+kappa))
		prefix_scale_mat = sparse.lil_matrix((len(prefixdict), len(prefixdict)))
		prefix_scale_mat.setdiag(prefix_counts)
		prefix_scale_mat = prefix_scale_mat.tocsr()
		prefix_scale_mat.eliminate_zeros()
		prefix_scale_mat.prune()
		prefix_scale_mat = sppy.csarray.fromScipySparse(prefix_scale_mat)
		suffix_scale_mat = sparse.lil_matrix((len(prefixdict), len(prefixdict)))
		suffix_scale_mat.setdiag(suffix_counts)
		suffix_scale_mat = suffix_scale_mat.tocsr()
		suffix_scale_mat.eliminate_zeros()
		suffix_scale_mat.prune()
		suffix_scale_mat = sppy.csarray.fromScipySparse(suffix_scale_mat)
		prefix_scale_mat = prefix_scale_mat
		suffix_scale_mat = suffix_scale_mat


		hankelmat.compress()
		hankelmat = prefix_scale_mat.dot(hankelmat).dot(suffix_scale_mat)

		for sym in range(n_symbols):
		    symbol_hankels[sym] = prefix_scale_mat.dot(symbol_hankels[sym]).dot(suffix_scale_mat)

		return hankelmat[0,:], hankelmat[:,0], hankelmat, symbol_hankels

def construct_string_hankels(data, prefixdict, suffixdict, n_symbols, basis_length=100):
		kappa = 5
		hankelmat = sppy.csarray((len(prefixdict),len(suffixdict)))
		hankelmat.reserve(int((len(prefixdict)*len(suffixdict))*EXPECTED_SPARSENESS))
		prefix_counts = np.zeros(len(prefixdict))
		suffix_counts = np.zeros(len(prefixdict))
		symbol_hankels = {}
		for sym in range(n_symbols):
		    symbol_hankels[sym] = sppy.csarray((len(prefixdict),len(suffixdict)))
		    symbol_hankels[sym].reserve(int((len(prefixdict)*len(suffixdict))*EXPECTED_SPARSENESS))


		#iterating over sequences
		for seq in data:
		    #iterating over prefix start positions
		    for i in range(len(seq)+1):
		        #break if prefix longer than anything in basis
		        if  i > basis_length:
		            break
		        #break if suffix longer than anything in basis
		        if len(seq) - i > basis_length:
		            break
		     
		        prefix = tuple(seq[0:i])

		        #if suffix empty string at end of word then special case
		        if i == len(seq):
		            suffix = tuple([])
		        else:
		            suffix = tuple(seq[i:len(seq)])

		        if prefix in prefixdict and suffix in suffixdict:
		            prefixind = prefixdict[prefix]
		            suffixind = suffixdict[suffix]
		            hankelmat[prefixind, suffixind] += 1


		        if i < len(seq):
		            sym = seq[i]

		            symprefix = tuple(seq[0:i])
		            #special case when suffix is just empty string
		            if i + 1 == len(seq):
		                symsuffix = tuple([])
		            else:
		                symsuffix = tuple(seq[i+1:len(seq)])

		            if symprefix in prefixdict and symsuffix in suffixdict:
		                symprefixind = prefixdict[symprefix]
		                symsuffixind = suffixdict[symsuffix]
		                symbol_hankel = symbol_hankels[sym]
		                symbol_hankel[symprefixind, symsuffixind] += 1


		prefix_counts = np.sqrt(float(len(data))/(prefix_counts+kappa))
		suffix_counts = np.sqrt(float(len(data))/(suffix_counts+kappa))
		prefix_scale_mat = sparse.lil_matrix((len(prefixdict), len(prefixdict)))
		prefix_scale_mat.setdiag(prefix_counts)
		prefix_scale_mat = prefix_scale_mat.tocsr()
		prefix_scale_mat.eliminate_zeros()
		prefix_scale_mat.prune()
		prefix_scale_mat = sppy.csarray.fromScipySparse(prefix_scale_mat)
		suffix_scale_mat = sparse.lil_matrix((len(prefixdict), len(prefixdict)))
		suffix_scale_mat.setdiag(suffix_counts)
		suffix_scale_mat = suffix_scale_mat.tocsr()
		suffix_scale_mat.eliminate_zeros()
		suffix_scale_mat.prune()
		suffix_scale_mat = sppy.csarray.fromScipySparse(suffix_scale_mat)
		prefix_scale_mat = prefix_scale_mat
		suffix_scale_mat = suffix_scale_mat

		hankelmat.compress()

		hankelmat = prefix_scale_mat.dot(hankelmat).dot(suffix_scale_mat)

		for sym in range(n_symbols):
		    #symbol_hankels[sym] = symbol_hankels[sym]*(1.0/(float(len(data)))) 
		    symbol_hankels[sym] = prefix_scale_mat.dot(symbol_hankels[sym]).dot(suffix_scale_mat)

		#NOTE: passing only one h valid since suffix basis and prefix basis identical
		# return hankelmat[1:hankelmat.shape[0],0],hankelmat[1:hankelmat.shape[0],1:hankelmat.shape[1]]
		return hankelmat[0,:], hankelmat[:,0], hankelmat, symbol_hankels

def construct_tensor_hankels(data, prefixdict, suffixdict, num_symbols, max_basis_length):
	basissize = len(prefixdict)
	kappa = 5

	h_pands = sparse.lil_matrix((basissize, basissize))
	hp_pandsigma = sparse.lil_matrix((basissize, num_symbols))
	symbol_hankels = {}
	# prefix_counts = np.zeros(basissize)
	# suffix_counts = np.zeros(basissize)

	h_pands = sparse.lil_matrix((basissize, basissize))
	hp_pandsigma = sparse.lil_matrix((basissize, num_symbols))
	symbol_hankels = {}

	for i in range(num_symbols):
		symbol_hankels[i] = sparse.lil_matrix((basissize, basissize))

	for seq in data:
		#iterate over prefix start
		for i in range(len(seq)):
			#iterate over suffix start
			for j in range(i,len(seq)+1):
				#if prefix length greater than basis size then exit
				if j - i > max_basis_length:
					break
				#iterate over suffix end
				for k in range(j,len(seq)+1):
					#if suffix length greater than basis size then exit
					if k - j > max_basis_length:
						break

					prefix = tuple(seq[i:j])

					#if suffix empty string at end of word then special case
					if k == len(seq) and j == len(seq):
						suffix = tuple([])
					else:
						suffix = tuple(seq[j:k])

					if prefix in prefixdict and suffix in suffixdict:
						prefixind = prefixdict[prefix]
						suffixind = suffixdict[suffix]
						# prefix_counts[prefixind] += 1
						# suffix_counts[suffixind] += 1
						h_pands[prefixind, suffixind] += 1

					#index must be one smaller when looking at prefix,symbol,suffix tuples.
					if j < len(seq) and k < len(seq):
						sym = seq[j]
						symprefix = tuple(seq[i:j])

						#if symsuffix empty string then special case
						if k + 1 == len(seq) and j + 1 == len(seq):
							symsuffix = tuple([])
						else:
							symsuffix = tuple(seq[j+1:k+1])

						if symprefix in prefixdict:
							symprefixind = prefixdict[symprefix]
							hp_pandsigma[symprefixind,sym] += 1


						if symprefix in prefixdict and symsuffix in suffixdict:
							symprefixind = prefixdict[symprefix]
							symsuffixind = suffixdict[symsuffix]
							symhankelmat = symbol_hankels[sym]
							symhankelmat[symprefixind, symsuffixind] += 1

	# prefix_counts = np.sqrt(float(len(data))/(prefix_counts+kappa))
	# suffix_counts = np.sqrt(float(len(data))/(suffix_counts+kappa))
	# prefix_scale_mat = sparse.lil_matrix((basissize, basissize))
	# prefix_scale_mat.setdiag(prefix_counts)
	# prefix_scale_mat = prefix_scale_mat.tocsr()
	# prefix_scale_mat.eliminate_zeros()
	# prefix_scale_mat.prune()
	# suffix_scale_mat = sparse.lil_matrix((basissize, basissize))
	# suffix_scale_mat.setdiag(suffix_counts)
	# suffix_scale_mat = suffix_scale_mat.tocsr()
	# suffix_scale_mat.eliminate_zeros()
	# suffix_scale_mat.prune()
	
	h_pands = sparse.csr_matrix(h_pands)
	h_pands.eliminate_zeros()
	h_pands.prune()
	h_pands  /= float(len(data))
	hp_pandsigma = sparse.csr_matrix(hp_pandsigma)
	hp_pandsigma /= float(len(data))
	hp_pandsigma.eliminate_zeros()
	hp_pandsigma.prune()

	for sym in range(num_symbols):
		symbol_hankels[sym] /= float(len(data))

	hbar_pands = sparse.csr_matrix((basissize, basissize))
	hbar_sigmaands = sparse.lil_matrix((num_symbols, basissize))
	hbar_pandsigma = sparse.lil_matrix((basissize, num_symbols))

	for sym in range(num_symbols):
		symhankelmat = symbol_hankels[sym]
		hbar_pands = hbar_pands+symhankelmat
		#hbar_sigmaands[sym,:] += symhankelmat.sum(0).A[0]
		#hbar_pandsigma[:,sym] += symhankelmat.sum(1).T.A[0]
		hbar_sigmaands[sym,:] = hbar_sigmaands[sym,:]+symhankelmat.sum(0)
		hbar_pandsigma[:,sym] = hbar_pandsigma[:,sym]+symhankelmat.sum(1)

	for sym in range(num_symbols):
		symbol_hankels[sym] = sparse.csr_matrix(symbol_hankels[sym])
		symbol_hankels[sym].eliminate_zeros()
		symbol_hankels[sym].prune()

	hbar_pands.eliminate_zeros();
	hbar_pands.prune()
	hbar_pandsigma = sparse.csr_matrix(hbar_pandsigma)
	hbar_pandsigma.eliminate_zeros()
	hbar_pandsigma.prune()
	hbar_sigmaands = sparse.csr_matrix(hbar_sigmaands)
	hbar_sigmaands.eliminate_zeros()
	hbar_sigmaands.prune()

	if DEBUG:
		print "hbar_pandsigma sum: ", hbar_pandsigma.sum()
		print "hbar_sigmaands sum: ", hbar_sigmaands.sum()
		print "hbar_pands sum: ", hbar_pands.sum()

		tensorsum = 0
		for sym in range(num_symbols):
			tensorsum += symbol_hankels[sym].sum()
		print "Tensor sum: ", tensorsum

	return h_pands,symbol_hankels,hp_pandsigma,hbar_pands,hbar_pandsigma, hbar_sigmaands













