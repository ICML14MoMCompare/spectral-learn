# MoMs-for-StochasticLanguages

Collection of three method of moments based algorithms for learning stochastic languages.
Accompaniment to the ICML 2014 paper "Methods of Moments for Learning Stochastic Languages: Unified Presentation and Empirical Comparison",
by Borja Balle, William L. Hamilton, and Joelle Pineau.
The learning algorithms produce weighted finite automata that can be used to make predictions over string.
For more details, the paper can be found here: http://www.cs.upc.edu/~bballe/papers/icml14-mom.pdf

### Includes: 
  Spectral learning algorithm (using both string and substring estimates).
  Convex optimization learning algorithm.
  Tensor decomposition based learning algorithm.
  Wrapper for the Treba EM library.
  12 synthetic datasets from the PAutomac competition (http://ai.cs.umbc.edu/icgi2012/challenge/Pautomac/download.php).
  One real-world NLP dataset from the parts-of-speech tags in the Penn-Treebank2 dataset (http://www.cis.upenn.edu/~treebank/).

### Prerequisites:

The code requires Python 2.7, SciPy (http://scipy.org/), the Python SpPy (http://pythonhosted.org/sppy/).
Necessary C++ libraries are provided.
The latest GCC compiler is recommended.

## Installation:

For the spectral methods:
  No installation required.

For the CO method: 
  Navigate to the co/cpp directory.
  Run "make" command.
  
For the tensor method:
  Navigate to the tensor/cpp directory.
  Run "make" command.

Using the code (with settings from the paper):

For all methods:
  Modify the PAUTOMACPATH, RESULTSPATH, etc. variables at the start of the main methods sections as necessary.
  
For the spectral method:
  Navigate to the spectral directory.
  Run "python wfaspectrallearn [esttype] [metric] [problem-id] [n-symbols] [basissize]
    estype is either "string" or "substring"
    metric is either "WER" or "KL" (perplexity)
    problem-id is the PAutomaC problem id number or "tree" for the Treebank2 dataset.
    n-symbols is the number of symbols in the target alphabet.
    basissize is the number of prefixes/suffixes to use in estimation.

For the convex optimization method:
  Navigate to the co directory.
  Run "python wfacolearn.py [modeltype] [metric] [problem-id] [n-symbols] [basissize]
    modeltype should be set to "WFA" (other settings not used in the paper).
    other args same as above.
    
For the tensor method:
  Navigate to the tensor directory.
  Run "python wfatensorlearn.py [metric] [problem-id] [n-symbols] [basisize]
  *this method automatically outputs a model description in the directory specified by MODELDIR in the main method.

For the EM wrapper code:
  Navigate to the em directory:
  Run "python wfaemlearn.py [problem-id] [n-symbols]

For the tensor-initialized EM:
  Navigate to the em directory:
  Run "python wfamominitlearn.py [metric] [problem-id] [numsymbols] [path-to-tensor-model-file]
