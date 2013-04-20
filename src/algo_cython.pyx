
include "interrupt.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

from sage.rings.integer_ring import ZZ, CC
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.matrix.matrix_integer_dense import Matrix_integer_dense


cdef extern from "algo_cpp.cpp":
	void test_algo()
	cdef cppclass M2T_O:
		int a,b1,b2,c
	cdef cppclass CurlS_Generator:
		M2T_O getNextS()
	cdef cppclass PrecisionF:
		int B
	cdef cppclass ReductionMatrices_Calc:
		ReductionMatrices_Calc()
		void init(int D, int HermWeight) except +
		PrecisionF curlF
		CurlS_Generator curlS
		size_t matrixRowCount, matrixColumnCount
		void calcMatrix()
		void getMatrix(mpz_t* out)

def test():
	test_algo()

cdef M2T_O_matrix(M2T_O m, int D):
	M = MatrixSpace(CC, 2, 2)
	b = m.b1 + m.b2 * (D + sqrt(D)) * 0.5
	return M([m.a, b, b.conjugate(), m.c])

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
		cdef M2T_O m = self.calc.curlS.getNextS()
		return M2T_O_matrix(m, self.D)
	def calcMatrix(self):
		if self.D == 0: raise RuntimeError, "you have to call init first"
		self.calc.calcMatrix()
	def getMatrix(self):
		"""
		:rtype : Matrix_integer_dense
		"""
		M = MatrixSpace(ZZ, self.calc.matrixRowCount, self.calc.matrixColumnCount)
		cdef Matrix_integer_dense m = M.zero_matrix().__copy__()
		self.calc.getMatrix(<mpz_t*>m._entries)
		return m

