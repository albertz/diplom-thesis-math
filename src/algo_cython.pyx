
include "interrupt.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

from sage.rings.integer_ring import ZZ
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense


cdef extern from "algo_cpp.cpp":
	void test_algo()
	cdef cppclass CurlS_Generator:
		void getNextS()
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

cdef class Calc:
	cdef ReductionMatrices_Calc calc
	def init(self, int D, int HermWeight):
		self.calc.init(D, HermWeight)
		# start limit
		self.calc.curlF.B = 20
	def getNextS(self):
		self.calc.curlS.getNextS()
	def calcMatrix(self):
		self.calc.calcMatrix()
	def getMatrix(self):
		M = MatrixSpace(ZZ, self.calc.matrixRowCount, self.calc.matrixColumnCount)
		cdef Matrix_integer_dense m = M.zero_matrix().__copy__()
		self.calc.getMatrix(<mpz_t*>m._entries)
		return m

def test_algo_cython():
	calc = Calc()
	calc.init(D = -4, HermWeight = 10)

	calc.getNextS()
	calc.getNextS()
	
	calc.calcMatrix()
	return calc.getMatrix()
