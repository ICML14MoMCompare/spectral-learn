

cdef extern from "include/SparseVectorExt.h":  
   cdef cppclass SparseVectorExt[T]:  
      SparseVectorExt() 
      SparseVectorExt(SparseVectorExt[T]) 
      SparseVectorExt(int)
#      double norm()
      int nonZeros()
      int rows()
      int size() 
      SparseVectorExt[T] abs()
      SparseVectorExt[T] add(SparseVectorExt[T]&)
      T dot(SparseVectorExt[T]&)
      SparseVectorExt[T] hadamard(SparseVectorExt[T]&)
      SparseVectorExt[T] negate()
      SparseVectorExt[T] subtract(SparseVectorExt[T]&)
      T coeff(int)
#      T sum()
      T sumValues()
      void insertVal(int, T) 
      void fill(T)
      void nonZeroInds(long*)
      void reserve(int)
      void scalarMultiply(double)
      void slice(int*, int, SparseVectorExt[T]*) 
      
cdef class csarray1d_signed_char:
    cdef SparseVectorExt[signed char] *thisPtr   
cdef class csarray1d_short:
    cdef SparseVectorExt[short] *thisPtr   
cdef class csarray1d_int:
    cdef SparseVectorExt[int] *thisPtr   
cdef class csarray1d_long:
    cdef SparseVectorExt[long] *thisPtr   
cdef class csarray1d_float:
    cdef SparseVectorExt[float] *thisPtr   
cdef class csarray1d_double:
    cdef SparseVectorExt[double] *thisPtr   
