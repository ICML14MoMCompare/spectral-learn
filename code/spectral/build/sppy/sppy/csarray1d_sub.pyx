# cython: profile=False
from cython.operator cimport dereference as deref, preincrement as inc 
import numpy 
cimport numpy
import cython 
from libcpp.vector cimport vector
numpy.import_array()


      
cdef class csarray1d_signed_char:
    def __cinit__(self, shape):
        """
        Create a new dynamic array.
        """
        cdef int shapeVal         
        
        if type(shape) == tuple and len(shape) == 1: 
            shapeVal = shape[0]
        elif type(shape) == int: 
            shapeVal = shape
        else: 
            raise ValueError("Invalid input: " + str(shape))
            
        self.thisPtr = new SparseVectorExt[signed char](shapeVal) 
            
    def __dealloc__(self): 
        """
        Deallocate the SparseVectorExt object.  
        """
        del self.thisPtr
        
    def __getNDim(self): 
        """
        Return the number of dimensions of this array. 
        """
        return 1 
        
    def __getShape(self):
        """
        Return the shape of this array
        """
        return (self.thisPtr.size(), )
        
    def __getSize(self): 
        """
        Return the size of this array, that is the number of elements. 
        """
        return self.thisPtr.size()   
        
    def getnnz(self): 
        """
        Return the number of non-zero elements in the array. 
        """
        return self.thisPtr.nonZeros()

    def __setitem__(self, ind, val):
        """
        Set elements of the array. If i is integers then the corresponding 
        value in the array is set. 
        """

        if type(ind) == numpy.ndarray : 
            self.put(val, ind)
        elif type(ind) == int:
            ind = int(ind) 
            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid index " + str(ind))   
            self.thisPtr.insertVal(ind, val) 
        else:
            raise ValueError("Invalid index " + str(ind))  
            
            

    def put(self, val, numpy.ndarray[numpy.int_t, ndim=1] inds not None): 
        """
        Insert a value or array of values into the array. 
        """
        cdef unsigned int ix 
        self.reserve(len(inds))
        
        if type(val) == numpy.ndarray: 
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val[ix])
        else:
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val)
                
    def reserve(self, int n): 
        """
        Reserve n nonzero entries  
        """
        self.thisPtr.reserve(n)
        
    def toarray(self): 
        """
        Convert this sparse matrix into a numpy array. 
        """
        cdef numpy.ndarray[double, ndim=1, mode="c"] result = numpy.zeros(self.shape, numpy.float)
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result[inds[i]] += self.thisPtr.coeff(inds[i])   
            
        return result 

    def nonzero(self): 
        """
        Return a tuple of arrays corresponding to nonzero elements. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds = numpy.zeros(self.getnnz(), dtype=numpy.int64) 
        
        if self.getnnz() != 0:
            self.thisPtr.nonZeroInds(&inds[0])
        
        return (inds, )

    def __getitem__(self, ind):
        """
        Get a value or set of values from the array (denoted A). Currently 3 types of parameters 
        are supported. If i is an integer then the corresponding element of the array 
        is returned. If i is an arrays of ints then we return the corresponding 
        values of A[i[k]] (note: i,j must be sorted in ascending order). If i
         is a slice e.g. a[1:5] then we return the submatrix corresponding to 
        the slice. 
        """        
        if type(ind) == tuple: 
            if len(ind)==1: 
                ind = ind[0]
            else: 
                raise ValueError("Invalid input " + str(ind))
        
        
        if type(ind) == numpy.ndarray or type(ind) == slice:
            indList = []            
            if type(ind) == numpy.ndarray: 
                indList.append(ind) 
            elif type(ind) == slice: 
                if ind.start == None: 
                    start = 0
                else: 
                    start = ind.start
                if ind.stop == None: 
                    stop = self.shape[0]
                else:
                    stop = ind.stop  
                indArr = numpy.arange(start, stop)
                indList.append(indArr)
            
            return self.subArray(indList[0])
        else:
            #Deal with negative indices
            if ind < 0: 
                ind += self.thisPtr.rows()

            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid row index " + str(ind))       
            return self.thisPtr.coeff(ind)     

    def subArray(self, numpy.ndarray[numpy.int_t, ndim=1, mode="c"] inds): 
        """
        Explicitly perform an array slice to return a submatrix with the given
        indices. Only works with ascending ordered indices. This is similar 
        to using numpy.ix_. 
        """
        cdef numpy.ndarray[int, ndim=1, mode="c"] indsC 
        cdef csarray1d_signed_char result = csarray1d_signed_char(inds.shape[0])     
        
        indsC = numpy.ascontiguousarray(inds, dtype=numpy.int32) 
        
        if inds.shape[0] != 0: 
            self.thisPtr.slice(&indsC[0], indsC.shape[0], result.thisPtr) 
        return result 

    def sum(self): 
        """
        Sum all of the elements in this array. 
        """      
        return self.thisPtr.sumValues()

    def mean(self,): 
        """
        Find the mean value of this array. 
        """
        if self.thisPtr.size() != 0:
            return self.sum()/float(self.thisPtr.size())
        else: 
            return float("nan")

    def __abs__(self): 
        """
        Return a matrix whose elements are the absolute values of this array. 
        """
        cdef csarray1d_signed_char result = csarray1d_signed_char(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[signed char](self.thisPtr.abs())
        return result 

    def copy(self): 
        """
        Return a copied version of this array. 
        """
        cdef csarray1d_signed_char result = csarray1d_signed_char(self.shape)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[signed char](deref(self.thisPtr))
        return result 

         
    def __mul__(self, double x):
        """
        Return a new array multiplied by a scalar value x. 
        """
        cdef csarray1d_signed_char result = self.copy() 
        result.thisPtr.scalarMultiply(x)
        return result 

    def min(self): 
        """
        Find the minimum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef signed char minVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            minVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) < minVal: 
                minVal = self.thisPtr.coeff(inds[i])
            
        return minVal 

    def max(self): 
        """
        Find the maximum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef signed char maxVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            maxVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) > maxVal: 
                maxVal = self.thisPtr.coeff(inds[i])
            
        return maxVal 

        
    def var(self): 
        """
        Return the variance of the elements of this array. 
        """
        cdef double mean = self.mean() 
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef double result = 0
        
        if self.size == 0: 
            return float("nan")
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result += (self.thisPtr.coeff(inds[i]) - mean)**2
        
        result += (self.size - self.getnnz())*mean**2
        result /= float(self.size)
        
        return result 
    
    def std(self): 
        """
        Return the standard deviation of the array elements. 
        """
        return numpy.sqrt(self.var())

    def __neg__(self): 
        """
        Return the negation of this array. 
        """
        cdef csarray1d_signed_char result = csarray1d_signed_char(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[signed char](self.thisPtr.negate())
        return result 

    def __add__(csarray1d_signed_char self, csarray1d_signed_char a): 
        """
        Add two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot add matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_signed_char result = csarray1d_signed_char(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[signed char](self.thisPtr.add(deref(a.thisPtr)))
        return result    
        
    def __sub__(csarray1d_signed_char self, csarray1d_signed_char a): 
        """
        Subtract two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot subtract matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_signed_char result = csarray1d_signed_char(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[signed char](self.thisPtr.subtract(deref(a.thisPtr)))
        return result 

    def ones(self): 
        """
        Fill the array with ones. 
        """
        self.thisPtr.fill(1)

    def hadamard(self, csarray1d_signed_char a): 
        """
        Find the element-wise matrix (hadamard) product. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot elementwise multiply matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef csarray1d_signed_char result = csarray1d_signed_char(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[signed char](self.thisPtr.hadamard(deref(a.thisPtr)))
        return result 

    def dotCsarray1d(self, csarray1d_signed_char a): 
        if self.shape != a.shape: 
            raise ValueError("Cannot compute dot product of matrices of shapes " + str(self.shape) + " and " + str(a.shape))
            
        return self.thisPtr.dot(deref(a.thisPtr))



    shape = property(__getShape)
    size = property(__getSize)
    ndim = property(__getNDim)        


cdef class csarray1d_short:
    def __cinit__(self, shape):
        """
        Create a new dynamic array.
        """
        cdef int shapeVal         
        
        if type(shape) == tuple and len(shape) == 1: 
            shapeVal = shape[0]
        elif type(shape) == int: 
            shapeVal = shape
        else: 
            raise ValueError("Invalid input: " + str(shape))
            
        self.thisPtr = new SparseVectorExt[short](shapeVal) 
            
    def __dealloc__(self): 
        """
        Deallocate the SparseVectorExt object.  
        """
        del self.thisPtr
        
    def __getNDim(self): 
        """
        Return the number of dimensions of this array. 
        """
        return 1 
        
    def __getShape(self):
        """
        Return the shape of this array
        """
        return (self.thisPtr.size(), )
        
    def __getSize(self): 
        """
        Return the size of this array, that is the number of elements. 
        """
        return self.thisPtr.size()   
        
    def getnnz(self): 
        """
        Return the number of non-zero elements in the array. 
        """
        return self.thisPtr.nonZeros()

    def __setitem__(self, ind, val):
        """
        Set elements of the array. If i is integers then the corresponding 
        value in the array is set. 
        """

        if type(ind) == numpy.ndarray : 
            self.put(val, ind)
        elif type(ind) == int:
            ind = int(ind) 
            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid index " + str(ind))   
            self.thisPtr.insertVal(ind, val) 
        else:
            raise ValueError("Invalid index " + str(ind))  
            
            

    def put(self, val, numpy.ndarray[numpy.int_t, ndim=1] inds not None): 
        """
        Insert a value or array of values into the array. 
        """
        cdef unsigned int ix 
        self.reserve(len(inds))
        
        if type(val) == numpy.ndarray: 
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val[ix])
        else:
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val)
                
    def reserve(self, int n): 
        """
        Reserve n nonzero entries  
        """
        self.thisPtr.reserve(n)
        
    def toarray(self): 
        """
        Convert this sparse matrix into a numpy array. 
        """
        cdef numpy.ndarray[double, ndim=1, mode="c"] result = numpy.zeros(self.shape, numpy.float)
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result[inds[i]] += self.thisPtr.coeff(inds[i])   
            
        return result 

    def nonzero(self): 
        """
        Return a tuple of arrays corresponding to nonzero elements. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds = numpy.zeros(self.getnnz(), dtype=numpy.int64) 
        
        if self.getnnz() != 0:
            self.thisPtr.nonZeroInds(&inds[0])
        
        return (inds, )

    def __getitem__(self, ind):
        """
        Get a value or set of values from the array (denoted A). Currently 3 types of parameters 
        are supported. If i is an integer then the corresponding element of the array 
        is returned. If i is an arrays of ints then we return the corresponding 
        values of A[i[k]] (note: i,j must be sorted in ascending order). If i
         is a slice e.g. a[1:5] then we return the submatrix corresponding to 
        the slice. 
        """        
        if type(ind) == tuple: 
            if len(ind)==1: 
                ind = ind[0]
            else: 
                raise ValueError("Invalid input " + str(ind))
        
        
        if type(ind) == numpy.ndarray or type(ind) == slice:
            indList = []            
            if type(ind) == numpy.ndarray: 
                indList.append(ind) 
            elif type(ind) == slice: 
                if ind.start == None: 
                    start = 0
                else: 
                    start = ind.start
                if ind.stop == None: 
                    stop = self.shape[0]
                else:
                    stop = ind.stop  
                indArr = numpy.arange(start, stop)
                indList.append(indArr)
            
            return self.subArray(indList[0])
        else:
            #Deal with negative indices
            if ind < 0: 
                ind += self.thisPtr.rows()

            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid row index " + str(ind))       
            return self.thisPtr.coeff(ind)     

    def subArray(self, numpy.ndarray[numpy.int_t, ndim=1, mode="c"] inds): 
        """
        Explicitly perform an array slice to return a submatrix with the given
        indices. Only works with ascending ordered indices. This is similar 
        to using numpy.ix_. 
        """
        cdef numpy.ndarray[int, ndim=1, mode="c"] indsC 
        cdef csarray1d_short result = csarray1d_short(inds.shape[0])     
        
        indsC = numpy.ascontiguousarray(inds, dtype=numpy.int32) 
        
        if inds.shape[0] != 0: 
            self.thisPtr.slice(&indsC[0], indsC.shape[0], result.thisPtr) 
        return result 

    def sum(self): 
        """
        Sum all of the elements in this array. 
        """      
        return self.thisPtr.sumValues()

    def mean(self,): 
        """
        Find the mean value of this array. 
        """
        if self.thisPtr.size() != 0:
            return self.sum()/float(self.thisPtr.size())
        else: 
            return float("nan")

    def __abs__(self): 
        """
        Return a matrix whose elements are the absolute values of this array. 
        """
        cdef csarray1d_short result = csarray1d_short(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[short](self.thisPtr.abs())
        return result 

    def copy(self): 
        """
        Return a copied version of this array. 
        """
        cdef csarray1d_short result = csarray1d_short(self.shape)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[short](deref(self.thisPtr))
        return result 

         
    def __mul__(self, double x):
        """
        Return a new array multiplied by a scalar value x. 
        """
        cdef csarray1d_short result = self.copy() 
        result.thisPtr.scalarMultiply(x)
        return result 

    def min(self): 
        """
        Find the minimum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef short minVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            minVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) < minVal: 
                minVal = self.thisPtr.coeff(inds[i])
            
        return minVal 

    def max(self): 
        """
        Find the maximum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef short maxVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            maxVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) > maxVal: 
                maxVal = self.thisPtr.coeff(inds[i])
            
        return maxVal 

        
    def var(self): 
        """
        Return the variance of the elements of this array. 
        """
        cdef double mean = self.mean() 
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef double result = 0
        
        if self.size == 0: 
            return float("nan")
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result += (self.thisPtr.coeff(inds[i]) - mean)**2
        
        result += (self.size - self.getnnz())*mean**2
        result /= float(self.size)
        
        return result 
    
    def std(self): 
        """
        Return the standard deviation of the array elements. 
        """
        return numpy.sqrt(self.var())

    def __neg__(self): 
        """
        Return the negation of this array. 
        """
        cdef csarray1d_short result = csarray1d_short(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[short](self.thisPtr.negate())
        return result 

    def __add__(csarray1d_short self, csarray1d_short a): 
        """
        Add two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot add matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_short result = csarray1d_short(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[short](self.thisPtr.add(deref(a.thisPtr)))
        return result    
        
    def __sub__(csarray1d_short self, csarray1d_short a): 
        """
        Subtract two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot subtract matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_short result = csarray1d_short(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[short](self.thisPtr.subtract(deref(a.thisPtr)))
        return result 

    def ones(self): 
        """
        Fill the array with ones. 
        """
        self.thisPtr.fill(1)

    def hadamard(self, csarray1d_short a): 
        """
        Find the element-wise matrix (hadamard) product. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot elementwise multiply matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef csarray1d_short result = csarray1d_short(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[short](self.thisPtr.hadamard(deref(a.thisPtr)))
        return result 

    def dotCsarray1d(self, csarray1d_short a): 
        if self.shape != a.shape: 
            raise ValueError("Cannot compute dot product of matrices of shapes " + str(self.shape) + " and " + str(a.shape))
            
        return self.thisPtr.dot(deref(a.thisPtr))



    shape = property(__getShape)
    size = property(__getSize)
    ndim = property(__getNDim)        


cdef class csarray1d_int:
    def __cinit__(self, shape):
        """
        Create a new dynamic array.
        """
        cdef int shapeVal         
        
        if type(shape) == tuple and len(shape) == 1: 
            shapeVal = shape[0]
        elif type(shape) == int: 
            shapeVal = shape
        else: 
            raise ValueError("Invalid input: " + str(shape))
            
        self.thisPtr = new SparseVectorExt[int](shapeVal) 
            
    def __dealloc__(self): 
        """
        Deallocate the SparseVectorExt object.  
        """
        del self.thisPtr
        
    def __getNDim(self): 
        """
        Return the number of dimensions of this array. 
        """
        return 1 
        
    def __getShape(self):
        """
        Return the shape of this array
        """
        return (self.thisPtr.size(), )
        
    def __getSize(self): 
        """
        Return the size of this array, that is the number of elements. 
        """
        return self.thisPtr.size()   
        
    def getnnz(self): 
        """
        Return the number of non-zero elements in the array. 
        """
        return self.thisPtr.nonZeros()

    def __setitem__(self, ind, val):
        """
        Set elements of the array. If i is integers then the corresponding 
        value in the array is set. 
        """

        if type(ind) == numpy.ndarray : 
            self.put(val, ind)
        elif type(ind) == int:
            ind = int(ind) 
            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid index " + str(ind))   
            self.thisPtr.insertVal(ind, val) 
        else:
            raise ValueError("Invalid index " + str(ind))  
            
            

    def put(self, val, numpy.ndarray[numpy.int_t, ndim=1] inds not None): 
        """
        Insert a value or array of values into the array. 
        """
        cdef unsigned int ix 
        self.reserve(len(inds))
        
        if type(val) == numpy.ndarray: 
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val[ix])
        else:
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val)
                
    def reserve(self, int n): 
        """
        Reserve n nonzero entries  
        """
        self.thisPtr.reserve(n)
        
    def toarray(self): 
        """
        Convert this sparse matrix into a numpy array. 
        """
        cdef numpy.ndarray[double, ndim=1, mode="c"] result = numpy.zeros(self.shape, numpy.float)
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result[inds[i]] += self.thisPtr.coeff(inds[i])   
            
        return result 

    def nonzero(self): 
        """
        Return a tuple of arrays corresponding to nonzero elements. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds = numpy.zeros(self.getnnz(), dtype=numpy.int64) 
        
        if self.getnnz() != 0:
            self.thisPtr.nonZeroInds(&inds[0])
        
        return (inds, )

    def __getitem__(self, ind):
        """
        Get a value or set of values from the array (denoted A). Currently 3 types of parameters 
        are supported. If i is an integer then the corresponding element of the array 
        is returned. If i is an arrays of ints then we return the corresponding 
        values of A[i[k]] (note: i,j must be sorted in ascending order). If i
         is a slice e.g. a[1:5] then we return the submatrix corresponding to 
        the slice. 
        """        
        if type(ind) == tuple: 
            if len(ind)==1: 
                ind = ind[0]
            else: 
                raise ValueError("Invalid input " + str(ind))
        
        
        if type(ind) == numpy.ndarray or type(ind) == slice:
            indList = []            
            if type(ind) == numpy.ndarray: 
                indList.append(ind) 
            elif type(ind) == slice: 
                if ind.start == None: 
                    start = 0
                else: 
                    start = ind.start
                if ind.stop == None: 
                    stop = self.shape[0]
                else:
                    stop = ind.stop  
                indArr = numpy.arange(start, stop)
                indList.append(indArr)
            
            return self.subArray(indList[0])
        else:
            #Deal with negative indices
            if ind < 0: 
                ind += self.thisPtr.rows()

            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid row index " + str(ind))       
            return self.thisPtr.coeff(ind)     

    def subArray(self, numpy.ndarray[numpy.int_t, ndim=1, mode="c"] inds): 
        """
        Explicitly perform an array slice to return a submatrix with the given
        indices. Only works with ascending ordered indices. This is similar 
        to using numpy.ix_. 
        """
        cdef numpy.ndarray[int, ndim=1, mode="c"] indsC 
        cdef csarray1d_int result = csarray1d_int(inds.shape[0])     
        
        indsC = numpy.ascontiguousarray(inds, dtype=numpy.int32) 
        
        if inds.shape[0] != 0: 
            self.thisPtr.slice(&indsC[0], indsC.shape[0], result.thisPtr) 
        return result 

    def sum(self): 
        """
        Sum all of the elements in this array. 
        """      
        return self.thisPtr.sumValues()

    def mean(self,): 
        """
        Find the mean value of this array. 
        """
        if self.thisPtr.size() != 0:
            return self.sum()/float(self.thisPtr.size())
        else: 
            return float("nan")

    def __abs__(self): 
        """
        Return a matrix whose elements are the absolute values of this array. 
        """
        cdef csarray1d_int result = csarray1d_int(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[int](self.thisPtr.abs())
        return result 

    def copy(self): 
        """
        Return a copied version of this array. 
        """
        cdef csarray1d_int result = csarray1d_int(self.shape)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[int](deref(self.thisPtr))
        return result 

         
    def __mul__(self, double x):
        """
        Return a new array multiplied by a scalar value x. 
        """
        cdef csarray1d_int result = self.copy() 
        result.thisPtr.scalarMultiply(x)
        return result 

    def min(self): 
        """
        Find the minimum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef int minVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            minVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) < minVal: 
                minVal = self.thisPtr.coeff(inds[i])
            
        return minVal 

    def max(self): 
        """
        Find the maximum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef int maxVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            maxVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) > maxVal: 
                maxVal = self.thisPtr.coeff(inds[i])
            
        return maxVal 

        
    def var(self): 
        """
        Return the variance of the elements of this array. 
        """
        cdef double mean = self.mean() 
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef double result = 0
        
        if self.size == 0: 
            return float("nan")
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result += (self.thisPtr.coeff(inds[i]) - mean)**2
        
        result += (self.size - self.getnnz())*mean**2
        result /= float(self.size)
        
        return result 
    
    def std(self): 
        """
        Return the standard deviation of the array elements. 
        """
        return numpy.sqrt(self.var())

    def __neg__(self): 
        """
        Return the negation of this array. 
        """
        cdef csarray1d_int result = csarray1d_int(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[int](self.thisPtr.negate())
        return result 

    def __add__(csarray1d_int self, csarray1d_int a): 
        """
        Add two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot add matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_int result = csarray1d_int(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[int](self.thisPtr.add(deref(a.thisPtr)))
        return result    
        
    def __sub__(csarray1d_int self, csarray1d_int a): 
        """
        Subtract two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot subtract matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_int result = csarray1d_int(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[int](self.thisPtr.subtract(deref(a.thisPtr)))
        return result 

    def ones(self): 
        """
        Fill the array with ones. 
        """
        self.thisPtr.fill(1)

    def hadamard(self, csarray1d_int a): 
        """
        Find the element-wise matrix (hadamard) product. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot elementwise multiply matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef csarray1d_int result = csarray1d_int(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[int](self.thisPtr.hadamard(deref(a.thisPtr)))
        return result 

    def dotCsarray1d(self, csarray1d_int a): 
        if self.shape != a.shape: 
            raise ValueError("Cannot compute dot product of matrices of shapes " + str(self.shape) + " and " + str(a.shape))
            
        return self.thisPtr.dot(deref(a.thisPtr))



    shape = property(__getShape)
    size = property(__getSize)
    ndim = property(__getNDim)        


cdef class csarray1d_long:
    def __cinit__(self, shape):
        """
        Create a new dynamic array.
        """
        cdef int shapeVal         
        
        if type(shape) == tuple and len(shape) == 1: 
            shapeVal = shape[0]
        elif type(shape) == int: 
            shapeVal = shape
        else: 
            raise ValueError("Invalid input: " + str(shape))
            
        self.thisPtr = new SparseVectorExt[long](shapeVal) 
            
    def __dealloc__(self): 
        """
        Deallocate the SparseVectorExt object.  
        """
        del self.thisPtr
        
    def __getNDim(self): 
        """
        Return the number of dimensions of this array. 
        """
        return 1 
        
    def __getShape(self):
        """
        Return the shape of this array
        """
        return (self.thisPtr.size(), )
        
    def __getSize(self): 
        """
        Return the size of this array, that is the number of elements. 
        """
        return self.thisPtr.size()   
        
    def getnnz(self): 
        """
        Return the number of non-zero elements in the array. 
        """
        return self.thisPtr.nonZeros()

    def __setitem__(self, ind, val):
        """
        Set elements of the array. If i is integers then the corresponding 
        value in the array is set. 
        """

        if type(ind) == numpy.ndarray : 
            self.put(val, ind)
        elif type(ind) == int:
            ind = int(ind) 
            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid index " + str(ind))   
            self.thisPtr.insertVal(ind, val) 
        else:
            raise ValueError("Invalid index " + str(ind))  
            
            

    def put(self, val, numpy.ndarray[numpy.int_t, ndim=1] inds not None): 
        """
        Insert a value or array of values into the array. 
        """
        cdef unsigned int ix 
        self.reserve(len(inds))
        
        if type(val) == numpy.ndarray: 
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val[ix])
        else:
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val)
                
    def reserve(self, int n): 
        """
        Reserve n nonzero entries  
        """
        self.thisPtr.reserve(n)
        
    def toarray(self): 
        """
        Convert this sparse matrix into a numpy array. 
        """
        cdef numpy.ndarray[double, ndim=1, mode="c"] result = numpy.zeros(self.shape, numpy.float)
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result[inds[i]] += self.thisPtr.coeff(inds[i])   
            
        return result 

    def nonzero(self): 
        """
        Return a tuple of arrays corresponding to nonzero elements. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds = numpy.zeros(self.getnnz(), dtype=numpy.int64) 
        
        if self.getnnz() != 0:
            self.thisPtr.nonZeroInds(&inds[0])
        
        return (inds, )

    def __getitem__(self, ind):
        """
        Get a value or set of values from the array (denoted A). Currently 3 types of parameters 
        are supported. If i is an integer then the corresponding element of the array 
        is returned. If i is an arrays of ints then we return the corresponding 
        values of A[i[k]] (note: i,j must be sorted in ascending order). If i
         is a slice e.g. a[1:5] then we return the submatrix corresponding to 
        the slice. 
        """        
        if type(ind) == tuple: 
            if len(ind)==1: 
                ind = ind[0]
            else: 
                raise ValueError("Invalid input " + str(ind))
        
        
        if type(ind) == numpy.ndarray or type(ind) == slice:
            indList = []            
            if type(ind) == numpy.ndarray: 
                indList.append(ind) 
            elif type(ind) == slice: 
                if ind.start == None: 
                    start = 0
                else: 
                    start = ind.start
                if ind.stop == None: 
                    stop = self.shape[0]
                else:
                    stop = ind.stop  
                indArr = numpy.arange(start, stop)
                indList.append(indArr)
            
            return self.subArray(indList[0])
        else:
            #Deal with negative indices
            if ind < 0: 
                ind += self.thisPtr.rows()

            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid row index " + str(ind))       
            return self.thisPtr.coeff(ind)     

    def subArray(self, numpy.ndarray[numpy.int_t, ndim=1, mode="c"] inds): 
        """
        Explicitly perform an array slice to return a submatrix with the given
        indices. Only works with ascending ordered indices. This is similar 
        to using numpy.ix_. 
        """
        cdef numpy.ndarray[int, ndim=1, mode="c"] indsC 
        cdef csarray1d_long result = csarray1d_long(inds.shape[0])     
        
        indsC = numpy.ascontiguousarray(inds, dtype=numpy.int32) 
        
        if inds.shape[0] != 0: 
            self.thisPtr.slice(&indsC[0], indsC.shape[0], result.thisPtr) 
        return result 

    def sum(self): 
        """
        Sum all of the elements in this array. 
        """      
        return self.thisPtr.sumValues()

    def mean(self,): 
        """
        Find the mean value of this array. 
        """
        if self.thisPtr.size() != 0:
            return self.sum()/float(self.thisPtr.size())
        else: 
            return float("nan")

    def __abs__(self): 
        """
        Return a matrix whose elements are the absolute values of this array. 
        """
        cdef csarray1d_long result = csarray1d_long(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[long](self.thisPtr.abs())
        return result 

    def copy(self): 
        """
        Return a copied version of this array. 
        """
        cdef csarray1d_long result = csarray1d_long(self.shape)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[long](deref(self.thisPtr))
        return result 

         
    def __mul__(self, double x):
        """
        Return a new array multiplied by a scalar value x. 
        """
        cdef csarray1d_long result = self.copy() 
        result.thisPtr.scalarMultiply(x)
        return result 

    def min(self): 
        """
        Find the minimum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef long minVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            minVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) < minVal: 
                minVal = self.thisPtr.coeff(inds[i])
            
        return minVal 

    def max(self): 
        """
        Find the maximum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef long maxVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            maxVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) > maxVal: 
                maxVal = self.thisPtr.coeff(inds[i])
            
        return maxVal 

        
    def var(self): 
        """
        Return the variance of the elements of this array. 
        """
        cdef double mean = self.mean() 
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef double result = 0
        
        if self.size == 0: 
            return float("nan")
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result += (self.thisPtr.coeff(inds[i]) - mean)**2
        
        result += (self.size - self.getnnz())*mean**2
        result /= float(self.size)
        
        return result 
    
    def std(self): 
        """
        Return the standard deviation of the array elements. 
        """
        return numpy.sqrt(self.var())

    def __neg__(self): 
        """
        Return the negation of this array. 
        """
        cdef csarray1d_long result = csarray1d_long(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[long](self.thisPtr.negate())
        return result 

    def __add__(csarray1d_long self, csarray1d_long a): 
        """
        Add two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot add matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_long result = csarray1d_long(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[long](self.thisPtr.add(deref(a.thisPtr)))
        return result    
        
    def __sub__(csarray1d_long self, csarray1d_long a): 
        """
        Subtract two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot subtract matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_long result = csarray1d_long(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[long](self.thisPtr.subtract(deref(a.thisPtr)))
        return result 

    def ones(self): 
        """
        Fill the array with ones. 
        """
        self.thisPtr.fill(1)

    def hadamard(self, csarray1d_long a): 
        """
        Find the element-wise matrix (hadamard) product. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot elementwise multiply matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef csarray1d_long result = csarray1d_long(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[long](self.thisPtr.hadamard(deref(a.thisPtr)))
        return result 

    def dotCsarray1d(self, csarray1d_long a): 
        if self.shape != a.shape: 
            raise ValueError("Cannot compute dot product of matrices of shapes " + str(self.shape) + " and " + str(a.shape))
            
        return self.thisPtr.dot(deref(a.thisPtr))



    shape = property(__getShape)
    size = property(__getSize)
    ndim = property(__getNDim)        


cdef class csarray1d_float:
    def __cinit__(self, shape):
        """
        Create a new dynamic array.
        """
        cdef int shapeVal         
        
        if type(shape) == tuple and len(shape) == 1: 
            shapeVal = shape[0]
        elif type(shape) == int: 
            shapeVal = shape
        else: 
            raise ValueError("Invalid input: " + str(shape))
            
        self.thisPtr = new SparseVectorExt[float](shapeVal) 
            
    def __dealloc__(self): 
        """
        Deallocate the SparseVectorExt object.  
        """
        del self.thisPtr
        
    def __getNDim(self): 
        """
        Return the number of dimensions of this array. 
        """
        return 1 
        
    def __getShape(self):
        """
        Return the shape of this array
        """
        return (self.thisPtr.size(), )
        
    def __getSize(self): 
        """
        Return the size of this array, that is the number of elements. 
        """
        return self.thisPtr.size()   
        
    def getnnz(self): 
        """
        Return the number of non-zero elements in the array. 
        """
        return self.thisPtr.nonZeros()

    def __setitem__(self, ind, val):
        """
        Set elements of the array. If i is integers then the corresponding 
        value in the array is set. 
        """

        if type(ind) == numpy.ndarray : 
            self.put(val, ind)
        elif type(ind) == int:
            ind = int(ind) 
            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid index " + str(ind))   
            self.thisPtr.insertVal(ind, val) 
        else:
            raise ValueError("Invalid index " + str(ind))  
            
            

    def put(self, val, numpy.ndarray[numpy.int_t, ndim=1] inds not None): 
        """
        Insert a value or array of values into the array. 
        """
        cdef unsigned int ix 
        self.reserve(len(inds))
        
        if type(val) == numpy.ndarray: 
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val[ix])
        else:
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val)
                
    def reserve(self, int n): 
        """
        Reserve n nonzero entries  
        """
        self.thisPtr.reserve(n)
        
    def toarray(self): 
        """
        Convert this sparse matrix into a numpy array. 
        """
        cdef numpy.ndarray[double, ndim=1, mode="c"] result = numpy.zeros(self.shape, numpy.float)
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result[inds[i]] += self.thisPtr.coeff(inds[i])   
            
        return result 

    def nonzero(self): 
        """
        Return a tuple of arrays corresponding to nonzero elements. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds = numpy.zeros(self.getnnz(), dtype=numpy.int64) 
        
        if self.getnnz() != 0:
            self.thisPtr.nonZeroInds(&inds[0])
        
        return (inds, )

    def __getitem__(self, ind):
        """
        Get a value or set of values from the array (denoted A). Currently 3 types of parameters 
        are supported. If i is an integer then the corresponding element of the array 
        is returned. If i is an arrays of ints then we return the corresponding 
        values of A[i[k]] (note: i,j must be sorted in ascending order). If i
         is a slice e.g. a[1:5] then we return the submatrix corresponding to 
        the slice. 
        """        
        if type(ind) == tuple: 
            if len(ind)==1: 
                ind = ind[0]
            else: 
                raise ValueError("Invalid input " + str(ind))
        
        
        if type(ind) == numpy.ndarray or type(ind) == slice:
            indList = []            
            if type(ind) == numpy.ndarray: 
                indList.append(ind) 
            elif type(ind) == slice: 
                if ind.start == None: 
                    start = 0
                else: 
                    start = ind.start
                if ind.stop == None: 
                    stop = self.shape[0]
                else:
                    stop = ind.stop  
                indArr = numpy.arange(start, stop)
                indList.append(indArr)
            
            return self.subArray(indList[0])
        else:
            #Deal with negative indices
            if ind < 0: 
                ind += self.thisPtr.rows()

            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid row index " + str(ind))       
            return self.thisPtr.coeff(ind)     

    def subArray(self, numpy.ndarray[numpy.int_t, ndim=1, mode="c"] inds): 
        """
        Explicitly perform an array slice to return a submatrix with the given
        indices. Only works with ascending ordered indices. This is similar 
        to using numpy.ix_. 
        """
        cdef numpy.ndarray[int, ndim=1, mode="c"] indsC 
        cdef csarray1d_float result = csarray1d_float(inds.shape[0])     
        
        indsC = numpy.ascontiguousarray(inds, dtype=numpy.int32) 
        
        if inds.shape[0] != 0: 
            self.thisPtr.slice(&indsC[0], indsC.shape[0], result.thisPtr) 
        return result 

    def sum(self): 
        """
        Sum all of the elements in this array. 
        """      
        return self.thisPtr.sumValues()

    def mean(self,): 
        """
        Find the mean value of this array. 
        """
        if self.thisPtr.size() != 0:
            return self.sum()/float(self.thisPtr.size())
        else: 
            return float("nan")

    def __abs__(self): 
        """
        Return a matrix whose elements are the absolute values of this array. 
        """
        cdef csarray1d_float result = csarray1d_float(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[float](self.thisPtr.abs())
        return result 

    def copy(self): 
        """
        Return a copied version of this array. 
        """
        cdef csarray1d_float result = csarray1d_float(self.shape)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[float](deref(self.thisPtr))
        return result 

         
    def __mul__(self, double x):
        """
        Return a new array multiplied by a scalar value x. 
        """
        cdef csarray1d_float result = self.copy() 
        result.thisPtr.scalarMultiply(x)
        return result 

    def min(self): 
        """
        Find the minimum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef float minVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            minVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) < minVal: 
                minVal = self.thisPtr.coeff(inds[i])
            
        return minVal 

    def max(self): 
        """
        Find the maximum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef float maxVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            maxVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) > maxVal: 
                maxVal = self.thisPtr.coeff(inds[i])
            
        return maxVal 

        
    def var(self): 
        """
        Return the variance of the elements of this array. 
        """
        cdef double mean = self.mean() 
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef double result = 0
        
        if self.size == 0: 
            return float("nan")
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result += (self.thisPtr.coeff(inds[i]) - mean)**2
        
        result += (self.size - self.getnnz())*mean**2
        result /= float(self.size)
        
        return result 
    
    def std(self): 
        """
        Return the standard deviation of the array elements. 
        """
        return numpy.sqrt(self.var())

    def __neg__(self): 
        """
        Return the negation of this array. 
        """
        cdef csarray1d_float result = csarray1d_float(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[float](self.thisPtr.negate())
        return result 

    def __add__(csarray1d_float self, csarray1d_float a): 
        """
        Add two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot add matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_float result = csarray1d_float(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[float](self.thisPtr.add(deref(a.thisPtr)))
        return result    
        
    def __sub__(csarray1d_float self, csarray1d_float a): 
        """
        Subtract two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot subtract matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_float result = csarray1d_float(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[float](self.thisPtr.subtract(deref(a.thisPtr)))
        return result 

    def ones(self): 
        """
        Fill the array with ones. 
        """
        self.thisPtr.fill(1)

    def hadamard(self, csarray1d_float a): 
        """
        Find the element-wise matrix (hadamard) product. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot elementwise multiply matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef csarray1d_float result = csarray1d_float(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[float](self.thisPtr.hadamard(deref(a.thisPtr)))
        return result 

    def dotCsarray1d(self, csarray1d_float a): 
        if self.shape != a.shape: 
            raise ValueError("Cannot compute dot product of matrices of shapes " + str(self.shape) + " and " + str(a.shape))
            
        return self.thisPtr.dot(deref(a.thisPtr))



    shape = property(__getShape)
    size = property(__getSize)
    ndim = property(__getNDim)        


cdef class csarray1d_double:
    def __cinit__(self, shape):
        """
        Create a new dynamic array.
        """
        cdef int shapeVal         
        
        if type(shape) == tuple and len(shape) == 1: 
            shapeVal = shape[0]
        elif type(shape) == int: 
            shapeVal = shape
        else: 
            raise ValueError("Invalid input: " + str(shape))
            
        self.thisPtr = new SparseVectorExt[double](shapeVal) 
            
    def __dealloc__(self): 
        """
        Deallocate the SparseVectorExt object.  
        """
        del self.thisPtr
        
    def __getNDim(self): 
        """
        Return the number of dimensions of this array. 
        """
        return 1 
        
    def __getShape(self):
        """
        Return the shape of this array
        """
        return (self.thisPtr.size(), )
        
    def __getSize(self): 
        """
        Return the size of this array, that is the number of elements. 
        """
        return self.thisPtr.size()   
        
    def getnnz(self): 
        """
        Return the number of non-zero elements in the array. 
        """
        return self.thisPtr.nonZeros()

    def __setitem__(self, ind, val):
        """
        Set elements of the array. If i is integers then the corresponding 
        value in the array is set. 
        """

        if type(ind) == numpy.ndarray : 
            self.put(val, ind)
        elif type(ind) == int:
            ind = int(ind) 
            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid index " + str(ind))   
            self.thisPtr.insertVal(ind, val) 
        else:
            raise ValueError("Invalid index " + str(ind))  
            
            

    def put(self, val, numpy.ndarray[numpy.int_t, ndim=1] inds not None): 
        """
        Insert a value or array of values into the array. 
        """
        cdef unsigned int ix 
        self.reserve(len(inds))
        
        if type(val) == numpy.ndarray: 
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val[ix])
        else:
            for ix in range(len(inds)): 
                self.thisPtr.insertVal(inds[ix], val)
                
    def reserve(self, int n): 
        """
        Reserve n nonzero entries  
        """
        self.thisPtr.reserve(n)
        
    def toarray(self): 
        """
        Convert this sparse matrix into a numpy array. 
        """
        cdef numpy.ndarray[double, ndim=1, mode="c"] result = numpy.zeros(self.shape, numpy.float)
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result[inds[i]] += self.thisPtr.coeff(inds[i])   
            
        return result 

    def nonzero(self): 
        """
        Return a tuple of arrays corresponding to nonzero elements. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds = numpy.zeros(self.getnnz(), dtype=numpy.int64) 
        
        if self.getnnz() != 0:
            self.thisPtr.nonZeroInds(&inds[0])
        
        return (inds, )

    def __getitem__(self, ind):
        """
        Get a value or set of values from the array (denoted A). Currently 3 types of parameters 
        are supported. If i is an integer then the corresponding element of the array 
        is returned. If i is an arrays of ints then we return the corresponding 
        values of A[i[k]] (note: i,j must be sorted in ascending order). If i
         is a slice e.g. a[1:5] then we return the submatrix corresponding to 
        the slice. 
        """        
        if type(ind) == tuple: 
            if len(ind)==1: 
                ind = ind[0]
            else: 
                raise ValueError("Invalid input " + str(ind))
        
        
        if type(ind) == numpy.ndarray or type(ind) == slice:
            indList = []            
            if type(ind) == numpy.ndarray: 
                indList.append(ind) 
            elif type(ind) == slice: 
                if ind.start == None: 
                    start = 0
                else: 
                    start = ind.start
                if ind.stop == None: 
                    stop = self.shape[0]
                else:
                    stop = ind.stop  
                indArr = numpy.arange(start, stop)
                indList.append(indArr)
            
            return self.subArray(indList[0])
        else:
            #Deal with negative indices
            if ind < 0: 
                ind += self.thisPtr.rows()

            if ind < 0 or ind>=self.thisPtr.rows(): 
                raise ValueError("Invalid row index " + str(ind))       
            return self.thisPtr.coeff(ind)     

    def subArray(self, numpy.ndarray[numpy.int_t, ndim=1, mode="c"] inds): 
        """
        Explicitly perform an array slice to return a submatrix with the given
        indices. Only works with ascending ordered indices. This is similar 
        to using numpy.ix_. 
        """
        cdef numpy.ndarray[int, ndim=1, mode="c"] indsC 
        cdef csarray1d_double result = csarray1d_double(inds.shape[0])     
        
        indsC = numpy.ascontiguousarray(inds, dtype=numpy.int32) 
        
        if inds.shape[0] != 0: 
            self.thisPtr.slice(&indsC[0], indsC.shape[0], result.thisPtr) 
        return result 

    def sum(self): 
        """
        Sum all of the elements in this array. 
        """      
        return self.thisPtr.sumValues()

    def mean(self,): 
        """
        Find the mean value of this array. 
        """
        if self.thisPtr.size() != 0:
            return self.sum()/float(self.thisPtr.size())
        else: 
            return float("nan")

    def __abs__(self): 
        """
        Return a matrix whose elements are the absolute values of this array. 
        """
        cdef csarray1d_double result = csarray1d_double(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[double](self.thisPtr.abs())
        return result 

    def copy(self): 
        """
        Return a copied version of this array. 
        """
        cdef csarray1d_double result = csarray1d_double(self.shape)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[double](deref(self.thisPtr))
        return result 

         
    def __mul__(self, double x):
        """
        Return a new array multiplied by a scalar value x. 
        """
        cdef csarray1d_double result = self.copy() 
        result.thisPtr.scalarMultiply(x)
        return result 

    def min(self): 
        """
        Find the minimum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef double minVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            minVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) < minVal: 
                minVal = self.thisPtr.coeff(inds[i])
            
        return minVal 

    def max(self): 
        """
        Find the maximum element of this array. 
        """
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef double maxVal 
        
        if self.size == 0: 
            return float("nan")
        elif self.getnnz() != self.size: 
            maxVal = 0 
        
        (inds,) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            if self.thisPtr.coeff(inds[i]) > maxVal: 
                maxVal = self.thisPtr.coeff(inds[i])
            
        return maxVal 

        
    def var(self): 
        """
        Return the variance of the elements of this array. 
        """
        cdef double mean = self.mean() 
        cdef numpy.ndarray[long, ndim=1, mode="c"] inds
        cdef unsigned int i
        cdef double result = 0
        
        if self.size == 0: 
            return float("nan")
        
        (inds, ) = self.nonzero()
            
        for i in range(inds.shape[0]): 
            result += (self.thisPtr.coeff(inds[i]) - mean)**2
        
        result += (self.size - self.getnnz())*mean**2
        result /= float(self.size)
        
        return result 
    
    def std(self): 
        """
        Return the standard deviation of the array elements. 
        """
        return numpy.sqrt(self.var())

    def __neg__(self): 
        """
        Return the negation of this array. 
        """
        cdef csarray1d_double result = csarray1d_double(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[double](self.thisPtr.negate())
        return result 

    def __add__(csarray1d_double self, csarray1d_double a): 
        """
        Add two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot add matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_double result = csarray1d_double(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[double](self.thisPtr.add(deref(a.thisPtr)))
        return result    
        
    def __sub__(csarray1d_double self, csarray1d_double a): 
        """
        Subtract two matrices together. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot subtract matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef int shapeVal = self.shape[0]
        cdef csarray1d_double result = csarray1d_double(shapeVal)
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[double](self.thisPtr.subtract(deref(a.thisPtr)))
        return result 

    def ones(self): 
        """
        Fill the array with ones. 
        """
        self.thisPtr.fill(1)

    def hadamard(self, csarray1d_double a): 
        """
        Find the element-wise matrix (hadamard) product. 
        """
        if self.shape != a.shape: 
            raise ValueError("Cannot elementwise multiply matrices of shapes " + str(self.shape) + " and " + str(a.shape))
        
        cdef csarray1d_double result = csarray1d_double(self.shape[0])
        del result.thisPtr
        result.thisPtr = new SparseVectorExt[double](self.thisPtr.hadamard(deref(a.thisPtr)))
        return result 

    def dotCsarray1d(self, csarray1d_double a): 
        if self.shape != a.shape: 
            raise ValueError("Cannot compute dot product of matrices of shapes " + str(self.shape) + " and " + str(a.shape))
            
        return self.thisPtr.dot(deref(a.thisPtr))



    shape = property(__getShape)
    size = property(__getSize)
    ndim = property(__getNDim)        


