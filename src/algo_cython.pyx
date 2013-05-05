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

		void calcMatrixTrans(const M2_O& tS, const M2_O& tT, int l) except +
		size_t matrixRowCountTrans, matrixColumnCountTrans, matrixCountTrans
		size_t matrixRowDenomTrans
		void getMatrixTrans(mpz_t* out, int matrixIndex) except +

def test():
	test_algo()

cdef M2T_O_fromC(M2T_O m, int D):
	"""
	:rtype : Matrix_symbolic_dense
	"""
	b = m.b1 + m.b2 * (D + ssqrt(D)) * 0.5
	return matrix(2, 2, [m.a, b, b.conjugate(), m.c])

cdef ElemOfCurlO O_toC(a, int D):
	cdef ElemOfCurlO b
	b.b2 = ZZ(a.imag() * 2 / ssqrt(-D))
	b.b1 = ZZ(a.real() - b.b2 * D / 2)
	return b

cdef M2_O M2_O_toC(m, int D):
	assert m.nrows() == 2 and m.ncols() == 2
	cdef M2_O _m
	_m.a = O_toC(m[0][0], D)
	_m.b = O_toC(m[0][1], D)
	_m.c = O_toC(m[1][0], D)
	_m.d = O_toC(m[1][1], D)
	return _m

cdef class Calc:
	# You need a recent Cython (e.g. >=0.19) for this.
	cdef ReductionMatrices_Calc calc
	cdef int D, HermWeight
	cdef public size_t matrixColumnCount
	cdef public size_t matrixRowCountTrans, matrixColumnCountTrans, matrixCountTrans, matrixRowDenomTrans

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

	def calcMatrixTrans(self, tS, tT, l):
		if self.D == 0: raise RuntimeError, "you have to call init first"
		cdef M2_O _tS
		cdef M2_O _tT
		_tS = M2_O_toC(tS, self.D)
		_tT = M2_O_toC(tT, self.D)
		self.calc.calcMatrixTrans(_tS, _tT, l)
		self.matrixRowCountTrans = self.calc.matrixRowCountTrans
		self.matrixColumnCountTrans = self.calc.matrixColumnCountTrans
		self.matrixCountTrans = self.calc.matrixCountTrans
		self.matrixRowDenomTrans = self.calc.matrixRowDenomTrans

		M = MatrixSpace(ZZ, self.calc.matrixRowCountTrans, self.calc.matrixColumnCountTrans)
		ms = [None] * self.matrixCountTrans
		for i in range(self.matrixCountTrans):
			m = M.zero_matrix().__copy__()
			self.calc.getMatrixTrans((<Matrix_integer_dense> m)._entries, i)
			ms[i] = m
		return ms
