# Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
# Copyright (c) 2013, Albert Zeyer, www.az2000.de
# This code is under the GPL v3 or later, see License.txt in the root directory of this project.

include "interrupt.pxi"
include "stdsage.pxi"
include "cdefs.pxi"

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.list cimport list

from sage.rings.integer import Integer
from sage.rings.number_field.number_field import QQ, ZZ, QuadraticField
from sage.matrix.constructor import matrix
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
	cdef cppclass M2T_Odual:
		int a,b1,b2,c
	cdef cppclass CurlS_Generator:
		list[M2T_O] matrices
		M2T_O getNextS()
		void clearMatrices()
		string getState() except +
		void setState(string state) except +
	cdef cppclass PrecisionF:
		int B
	cdef cppclass ReductionMatrices_Calc:
		ReductionMatrices_Calc()
		void init(int D, int HermWeight, string curlSiterType) except +
		PrecisionF curlF
		CurlS_Generator curlS
		vector[M2T_Odual] reducedCurlFList
		void calcReducedCurlF() except +

		void calcMatrix() except +
		size_t matrixRowCount, matrixColumnCount
		void getMatrix(mpz_t* out)
		void dumpMatrix() except +

		void calcMatrixTrans(const M2_O& tS, const M2_O& tT, int lS, int lT) except +
		size_t matrixRowCountTrans, matrixColumnCountTrans, matrixCountTrans
		size_t matrixRowDenomTrans
		void getMatrixTrans(mpz_t* out, int matrixIndex) except +
		void dumpMatrixTrans(int matrixIndex) except +

def test():
	test_algo()

cdef M2T_O_fromC(M2T_O m, int D):
	K = QuadraticField(D)
	Droot = K(D).sqrt()
	b = m.b1 + m.b2 * (D + Droot) * QQ(0.5)
	return matrix(K, 2, 2, [m.a, b, b.conjugate(), m.c])

cdef M2T_Odual_fromC(M2T_Odual m, int D):
	K = QuadraticField(D)
	Droot = K(D).sqrt()
	b = m.b1 / Droot + m.b2 * (1 + Droot) * QQ(0.5)
	return matrix(K, 2, 2, [m.a, b, b.conjugate(), m.c])

cdef ElemOfCurlO O_toC(a, int D) except *:
	cdef ElemOfCurlO b
	K = QuadraticField(D)
	Droot = K(D).sqrt()
	DrootHalf_coordinates_in_terms_of_powers = (Droot / 2).coordinates_in_terms_of_powers()
	try:
		real_part, b2 = DrootHalf_coordinates_in_terms_of_powers(a)
		b1 = real_part - b2 * D / 2
		b.b1 = ZZ(b1)
		b.b2 = ZZ(b2)
	except TypeError:
		print "cannot convert %r to CurlO(%i)" % (a, D)
		raise
	return b

cdef M2_O M2_O_toC(m, int D) except *:
	assert m.nrows() == 2 and m.ncols() == 2
	cdef M2_O _m
	_m.a = O_toC(m[0][0], D)
	_m.b = O_toC(m[0][1], D)
	_m.c = O_toC(m[1][0], D)
	_m.d = O_toC(m[1][1], D)
	return _m

cdef M2T_O M2T_O_toC(m, int D) except *:
	assert m.nrows() == 2 and m.ncols() == 2
	assert m[0,1] == m[1,0].conjugate()
	cdef M2T_O _m
	_m.a = m[0,0]
	_m.c = m[1,1]
	cdef ElemOfCurlO b = O_toC(m[0,1], D)
	_m.b1 = b.b1
	_m.b2 = b.b2
	return _m

cdef class Calc:
	"""
	This is mostly just a wrapper class to the C++ code in `algo_cpp.cpp`.
	All of the documentation can be found there. See also the text of the thesis.

	After you call `init()` on an instance of this class, you must also
	call `calcReducedCurlF()` to calculate the reduced index set of `\curlF`.

	Then, call `getNextS()` to get another `S`. You should call `curlS_clearMatrices()`
	every time before you do a further call to `getNextS()`.

	Then you can call `calcMatrix()` or `calcMatrixTrans()`.

	See the code in `algo.py` and `helpers.py` for a demonstration and usage.
	"""

	# You need a recent Cython (e.g. >=0.19) for this.
	cdef ReductionMatrices_Calc calc
	cdef public int D, HermWeight, B_cF
	cdef public object curlSiterType
	cdef public size_t matrixColumnCount
	cdef public size_t matrixRowCountTrans, matrixColumnCountTrans, matrixCountTrans, matrixRowDenomTrans
	cdef public object params
	cdef public object curlS

	def init(self, int D, int HermWeight, int B_cF=10, string curlSiterType="generic"):
		self.D = D
		self.HermWeight = HermWeight
		self.B_cF = B_cF
		self.curlSiterType = curlSiterType
		self.calc.init(D, HermWeight, curlSiterType)
		# start limit
		# this is never changed at the moment
		self.calc.curlF.B = B_cF
		self.params = (D,HermWeight,B_cF)
		self.curlS = ()

	def getNextS(self):
		"""
		:rtype : Matrix_symbolic_dense
		"""
		cdef M2T_O m = self.calc.curlS.getNextS()
		S = M2T_O_fromC(m, self.D)
		S.set_immutable()
		self.curlS += (S,)
		return S

	def curlS_clearMatrices(self):
		self.calc.curlS.clearMatrices()
		self.curlS = ()

	def calcReducedCurlF(self):
		if self.D == 0: raise RuntimeError, "you have to call init first"
		self.calc.calcReducedCurlF()
		self.matrixColumnCount = self.calc.matrixColumnCount

	def getReducedCurlF(self):
		size = self.calc.reducedCurlFList.size()
		l = [None] * size
		for i in range(size):
			l[i] = M2T_Odual_fromC(self.calc.reducedCurlFList[i], self.D)
			l[i].set_immutable()
		return tuple(l)

	def calcMatrix(self):
		"""
		:rtype : Matrix_integer_dense
		For documentation, see the C++ function `calcMatrix` in `algo_cpp.cpp`.
		"""
		if self.D == 0: raise RuntimeError, "you have to call init first"
		self.calc.calcMatrix()

		M = MatrixSpace(ZZ, self.calc.matrixRowCount, self.calc.matrixColumnCount)
		cdef Matrix_integer_dense m = M.zero_matrix().__copy__()
		self.calc.getMatrix(m._entries)
		return m

	def dumpMatrix(self):
		self.calc.dumpMatrix()

	def _getMatrixTrans(self, M, int i):
		cdef Matrix_integer_dense m = M.zero_matrix().__copy__()
		self.calc.getMatrixTrans(m._entries, i)
		return m

	def calcMatrixTrans(self, tS, tT, lS, lT):
		"""
		For documentation, see the C++ function `calcMatrixTrans` in `algo_cpp.cpp`.
		"""
		if self.D == 0: raise RuntimeError, "you have to call init first"
		cdef M2_O _tS
		cdef M2_O _tT
		_tS = M2_O_toC(tS, self.D)
		_tT = M2_O_toC(tT, self.D)
		self.calc.calcMatrixTrans(_tS, _tT, lS, lT)
		self.matrixRowCountTrans = self.calc.matrixRowCountTrans
		self.matrixColumnCountTrans = self.calc.matrixColumnCountTrans
		self.matrixCountTrans = self.calc.matrixCountTrans
		self.matrixRowDenomTrans = self.calc.matrixRowDenomTrans

		M = MatrixSpace(ZZ, self.calc.matrixRowCountTrans, self.calc.matrixColumnCountTrans)
		ms = [None] * self.matrixCountTrans
		for i in range(self.matrixCountTrans):
			ms[i] = self._getMatrixTrans(M, i)
		return ms

	def dumpMatrixTrans(self):
		for i in self.matrixCountTrans:
			self.calc.dumpMatrixTrans(i)

	def __getstate__(self):
		return ("prot2", self.params, self.curlSiterType, self.curlS, self.calc.curlS.getState())

	def __setstate__(self, state):
		protocol, params, curlSiterType, curlS, curlSiterState = state
		assert protocol == "prot2"
		D, HermWeight, B_cF = params
		self.init(D=D, HermWeight=HermWeight, B_cF=B_cF, curlSiterType=curlSiterType)
		self.calcReducedCurlF()
		self.curlS = curlS
		self.calc.curlS.matrices.clear()
		cdef M2T_O _S
		for S in curlS:
			_S = M2T_O_toC(S, D)
			self.calc.curlS.matrices.push_back(_S)
		self.calc.curlS.setState(curlSiterState)
