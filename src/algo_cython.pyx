from sage.matrix.constructor import matrix

include "interrupt.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

from sage.functions.other import sqrt as ssqrt
from sage.rings.integer_ring import ZZ
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.matrix.matrix_integer_dense import Matrix_integer_dense


cdef extern from "algo_cpp.cpp":
	void test_algo()
	cdef cppclass ElemOfCurlO:
		int b1,b2
	cdef cppclass ElemOfCurlOdual:
		int b1,b2
	cdef cppclass M2_O:
		ElemOfCurlO a,b,c,d
	cdef cppclass M2_Odual:
		ElemOfCurlOdual a,b,c,d
	cdef cppclass M2T_O:
		int a,b1,b2,c
	cdef cppclass CurlS_Generator:
		M2T_O getNextS()
		void clearMatrices()
	cdef cppclass PrecisionF:
		int B
	cdef cppclass ReductionMatrices_Calc:
		ReductionMatrices_Calc()
		void init(int D, int HermWeight) except +
		PrecisionF curlF
		CurlS_Generator curlS
		void calcReducedCurlF() except +

		void calcMatrix() except +
		size_t matrixRowCount, matrixColumnCount
		void getMatrix(mpz_t* out)

		void calcMatrixTranslated(const M2_O& tS, const M2_O& tT, const int l) except +
		size_t matrixRowCountTrans, matrixColumnCountTrans, matrixCountTrans;
		size_t matrixRowDenomTrans;
		void getMatrixTrans(mpz_t* out, int matrixIndex)

def test():
	test_algo()

cdef M2T_O_fromC(M2T_O m, int D):
	"""
	:rtype : Matrix_symbolic_dense
	"""
	b = m.b1 + m.b2 * (D + ssqrt(D)) * 0.5
	return matrix(2, 2, [m.a, b, b.conjugate(), m.c])

def M2T_O_toC(m, int D):
	pass

cdef class Calc:
	# You need a recent Cython (e.g. >=0.19) for this.
	cdef ReductionMatrices_Calc calc
	cdef int D, HermWeight
	cdef public size_t matrixColumnCount
	def init(self, int D, int HermWeight, int B_cF=20):
		self.D = D
		self.HermWeight = HermWeight
		self.calc.init(D, HermWeight)
		# start limit
		# this is never changed at the moment
		self.calc.curlF.B = B_cF
	def getNextS(self):
		"""
		:rtype : Matrix_symbolic_dense
		"""
		cdef M2T_O m = self.calc.curlS.getNextS()
		return M2T_O_fromC(m, self.D)
	def curlS_clearMatrices(self):
		self.calc.curlS.clearMatrices()
	def calcReducedCurlF(self):
		if self.D == 0: raise RuntimeError, "you have to call init first"
		self.calc.calcReducedCurlF()
		self.matrixColumnCount = self.calc.matrixColumnCount
	def calcMatrix(self):
		"""
		:rtype : Matrix_integer_dense
		"""
		if self.D == 0: raise RuntimeError, "you have to call init first"
		self.calc.calcMatrix()
		M = MatrixSpace(ZZ, self.calc.matrixRowCount, self.calc.matrixColumnCount)
		cdef Matrix_integer_dense m = M.zero_matrix().__copy__()
		self.calc.getMatrix(m._entries)
		return m
	def calcMatrixTrans(self, tS, tT):
		cdef M2_O _tS
		cdef M2_O _tT
