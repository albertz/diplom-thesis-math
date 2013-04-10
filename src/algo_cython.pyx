
include "interrupt.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

from sage.matrix.matrix_integer_dense import \
	Matrix_integer_dense #cdef class


cdef extern from "algo_cpp.cpp":
	void test_algo()

def test():
	test_algo()



