
include "interrupt.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

from sage.rings.integer_ring import ZZ
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense


cdef extern from "algo_cpp.cpp":
	void test_algo()
	cdef cppclass M2T:
		int a,b,c
	cdef cppclass CurlS_Generator:
		M2T getNextS()
	cdef cppclass PrecisionF:
		int B
	cdef cppclass ReductionMatrices_Calc:
		ReductionMatrices_Calc()
		void init(int D, int HermWeight)
		PrecisionF curlF
		CurlS_Generator curlS
		size_t matrixRowCount, matrixColumnCount
		void calcMatrix()
		void getMatrix(mpz_t* out)

def test():
	test_algo()

cdef M2T_matrix(M2T m):
	M = MatrixSpace(ZZ, 2, 2)
	return M([m.a, m.b, m.b, m.c])

cdef class Calc:
	cdef ReductionMatrices_Calc calc
	cdef int D, HermWeight
	def init(self, int D, int HermWeight):
		self.D = D
		self.HermWeight = HermWeight
		self.calc.init(D, HermWeight)
		# start limit
		self.calc.curlF.B = 20
	def getNextS(self):
		cdef M2T m = self.calc.curlS.getNextS()
		return M2T_matrix(m)
	def calcMatrix(self):
		if self.D == 0: raise RuntimeError, "you have to call init first"
		self.calc.calcMatrix()
	def getMatrix(self):
		M = MatrixSpace(ZZ, self.calc.matrixRowCount, self.calc.matrixColumnCount)
		cdef Matrix_integer_dense m = M.zero_matrix().__copy__()
		self.calc.getMatrix(<mpz_t*>m._entries)
		return m

