from sage.functions.other import floor
from sage.matrix.constructor import matrix, Matrix
from sage.modular.congroup import Gamma0
from sage.modular.modform.constructor import ModularForms
from sage.rings.arith import xgcd as orig_xgcd
from sage.rings.integer import Integer
from sage.rings.number_field.number_field import CyclotomicField, QQ, ZZ
from utils import *

# via Martin. while this is not in Sage:
import cusp_expansions


@persistent_cache(name="matrixTrans")
def _calcMatrixTrans(calc, tS, tT, lS, lT):
	try:
		ms = calc.calcMatrixTrans(tS, tT, lS, lT)
	except Exception:
		print (calc.params, calc.curlS, tS, tT, lS, lT)
		raise

	# Each matrix is for a zeta**i factor, where zeta is the n-th root of unity.
	# And n = calc.matrixCountTrans.
	assert len(ms) == calc.matrixCountTrans
	order = len(ms)

	K = CyclotomicField(order)
	zeta = K.gen()
	Kcoords = zeta.coordinates_in_terms_of_powers()

	assert len(K.power_basis()) == K.degree()
	new_ms = [matrix(QQ, ms[0].nrows(), ms[0].ncols()) for i in range(K.degree())]
	for l in range(order):
		coords = Kcoords(zeta**l)
		for i,m in enumerate(coords):
			new_ms[i] += ms[l] * m
	ms = new_ms

	denom = calc.matrixRowDenomTrans
	denom, ms = reduceNRow(denom=denom, mats=ms)

	return denom, order, ms


def calcMatrixTrans(calc, R):
	"""
	Returns a triple `(denom, order, ms)` where
	`ms` is a list of matrices. Each matrix is a factor to
	`zeta**i` where zeta is the order's root of unity,
	i.e. `zeta = CyclotomicField(order).gen()`.
	And we have `len(ms) == CyclotomicField(order).degree()`.

	The matrices represent the linear map of a Hermitian
	modular form `f` to the Elliptic modular form `(f|R)[S]`,
	both represented by Fourier expansions.
	(`R` is given as a parameter to this function. `S` is saved
	internally by the `calc` structure.)

	This uses the C++ function `Calc::calcMatrixTrans()` in  algo_cpp.cpp`.

	An Elliptic modular form given by vector `v` corresponds to the expansion

		\sum_{n \ge 0} v_n q^{n / denom}.

	INPUT:

	- `calc` -- The C++ Calc structure.

	- `R` -- The transformation matrix `R` in `(f|R)[S]`.

	OUTPUT:

	- The triple `(denom, order, ms)` as explained above.

	"""
	tS = R.submatrix(0,0,2,2)
	tT = R.submatrix(2,0,2,2)
	lS = _curlO_matrix_denom(tS, D=calc.D)
	lT = _curlO_matrix_denom(tT, D=calc.D)
	tS *= lS
	tT *= lT
	tS.set_immutable()
	tT.set_immutable()
	return _calcMatrixTrans(calc, tS, tT, lS, lT)



cuspExpansionsCache = PersistentCache(name="cuspExpansions")
def cuspExpansions(level, weight, prec):
	"""
	A cached version of `cusp_expansions.ModularFormsCuspExpansions._for_modular_forms()`.
	"""
	cacheIdx = (level, weight)
	if cacheIdx in cuspExpansionsCache:
		ce_prec,ce = cuspExpansionsCache[cacheIdx]
		if ce_prec >= prec: return ce
	verbose("calc ModularFormsCuspExpansions at level %i with weight %i and prec %i ..." % (level, weight, prec))
	ce = cusp_expansions.ModularFormsCuspExpansions._for_modular_forms(level, weight, prec)
	cuspExpansionsCache[cacheIdx] = prec, ce
	return ce


ellipBaseMatrixCache = PersistentCache(name="ellipBaseMatrix") # level,weight -> mat,prec
def getElliptModFormsBasisMatrix(level, weight, precision):
	"""
	Calculates the Echelon basis matrix of the Elliptic modular forms space.

	INPUT:

	- `level` -- The level for the modular group `\Gamma = \Gamma_0(level)`.

	- `weight` -- The weight of the Elliptic modular forms.

	- `precision` -- The precision of the Elliptic modular form Fourier expansions.

	OUTPUT:

	- A tuple `(dim,m)`, where `m` is a matrix which is the Echelon basis matrix of
	  the Elliptic modular forms over `\Gamma_0(level)` with weight `weight`
	  such that `m.ncols() == precision`.
	  `dim` is the dimension of the modular forms. You should check that `m.rank() == dim`,
	  otherwise it might not make sense to work with `m`.
	"""

	cacheIdx = (level, weight)
	if cacheIdx in ellipBaseMatrixCache and ellipBaseMatrixCache[cacheIdx][1] >= precision:
		fe_expansion_matrix_l = ellipBaseMatrixCache[cacheIdx][0]
		cut_matrix = fe_expansion_matrix_l[:,:precision]
		dim = fe_expansion_matrix_l.rank()
		return dim, cut_matrix
	#n = 2
	#while n < precision:
	#	n **= 2
	n = precision
	n = max(10, n)
	mf = ModularForms(Gamma0(level), weight)
	dim = mf.dimension()
	while True:
		fe_expansion_matrix_l = matrix(QQ, [b.qexp(n).padded_list(n) for b in mf.basis()])
		fe_expansion_matrix_l.echelonize()
		if fe_expansion_matrix_l.rank() == dim: break
		n += 10
	assert fe_expansion_matrix_l.rank() == dim, "%i != %i" % (fe_expansion_matrix_l.rank(), dim)
	ellipBaseMatrixCache[cacheIdx] = (fe_expansion_matrix_l, n)
	cut_matrix = fe_expansion_matrix_l[:,:precision]
	return dim, cut_matrix


@persistent_cache(name="restrMatrix")
def calcRestrictMatrix(calc):
	"""
	Calculates the matrix of the linear map `f \mapsto f[S]`
	where `f` is a Hermitian modular form.
	(`S` as well as all other parameters are saved internally in the `calc` structure.)

	This uses the C++ function `Calc::calcMatrix()` in `algo_cpp.cpp`.

	INPUT:

	- `calc` -- The C++ Calc structure.

	OUTPUT:

	- A matrix `m` which represents the linear map of the Fourier expansions.
	  The index set of the Fourier expansions of the Hermitian modular forms
	  can be returned by `herm_modform_indexset()`.

	"""
	return calc.calcMatrix()


def cusp_matrix(cusp):
	"""
	Returns a matrix M such that M * \infinity == cusp.
	"""
	if cusp == 0:
		return matrix(ZZ,2,2,[0,-1,1,0])
	else:
		a = cusp.numerator()
		c = cusp.denominator()
		_div, d, b = orig_xgcd(a, -c)
		assert _div == 1
		return matrix(ZZ,2,2,[a,b,c,d])


def calcPrecisionDimension(B_cF, S):
	"""
	Like the C++ implementation of `calcPrecisionDimension()`.
	"""
	assert S[0,1] == S[1,0].conjugate()
	s,t,u = S[0,0], S[0,1], S[1,1]
	s, u = QQ(s), QQ(u)
	precDim = B_cF * (s + u - 2 * abs(t))
	precDim = floor(precDim)
	return precDim


def _calcMatrix_py(D, HermWeight, S, B_cF):
	"""
	This is a Python implementation of the C++ `calcMatrix()` function.
	This is just for testing.
	"""

	assert S[0,1] == S[1,0].conjugate()
	s,t,u = S[0,0], S[0,1], S[1,1]
	from algo import herm_modform_indexset
	indexset = herm_modform_indexset(D=D, B_cF=B_cF)

	precDim = B_cF * (s + u - 2 * abs(t))
	precDim = floor(precDim)
	matrixRowCount = precDim
	matrixColumnCount = len(indexset)

	# ...
	raise NotImplementedError



def _simplify(a):
	try: return QQ(a)
	except: pass
	if hasattr(a, "simplify_full"):
		return a.simplify_full()
	return simplify(a)

class CurlO:
	# We set `b = b1 + b2 (D + \sqrt{D})/2`.
	def __init__(self, D):
		self.D = D
		assert (D*D - D) % 4 == 0
		self.field = QuadraticField(D)
		self.maxorder = self.field.maximal_order()
		self.Droot = self.field(D).sqrt()
		self.DrootHalf_coordinates_in_terms_of_powers = (self.Droot / 2).coordinates_in_terms_of_powers()
	def divmod(self, a, b):
		"""
		Returns q,r such that a = q*b + r.
		"""
		# Note that this implementation is quite naive!
		# Later, we can do better with QuadraticForm(...). (TODO)
		# Also, read here: http://www.fen.bilkent.edu.tr/~franz/publ/survey.pdf

		if b == 0: raise ZeroDivisionError
		a1,a2 = self.as_tuple_b(a)
		b1,b2 = self.as_tuple_b(b)

		#B = matrix([
		#	[b1, -b2 * (self.D**2 - self.D)/4],
		#	[b2, b1 + b2*self.D]
		#])
		Bdet = b1*b1 + b1*b2*self.D + b2*b2*(self.D**2 - self.D)/4
		Bdet = _simplify(Bdet)
		assert Bdet > 0
		qq1 = (a1*b1 + a1*b2*self.D + a2*b2*(self.D**2 - self.D)/4) / Bdet
		qq2 = (-a1*b2 + a2*b1) / Bdet
		assert _simplify(self.from_tuple_b(qq1,qq2) * b - a) == 0

		# not sure on this
		q1,q2 = int(round(qq1)), int(round(qq2))
		if (q1-qq1)*(q2-qq2) < 0:
			q1,q2 = int(floor(qq1)), int(floor(qq2))
		#print a1,a2,b1,b2,qq1,qq2,q1,q2
		q = self.from_tuple_b(q1,q2)
		# q * b + r == a
		r = _simplify(a - q * b)
		# Note that this works for -D < 5.27. See the text: 13.5.13,divmod
		euc_r, euc_b = map(self.euclidean_func, [r, b])
		assert euc_r < euc_b, "%r < %r; r=%r, b=%r, a=%r, b=%r" % (euc_r, euc_b, r, b, a, b)
		return q,r
	def euclidean_func(self, x):
		return self.field(x).abs()
	def divides(self, a, b):
		q,r = self.divmod(a, b)
		return r == 0
	def xgcd(self, a, b):
		if a == b: return a, 1, 0
		if a == -b: return a, 1, 0
		if a == 1: return 1, 1, 0
		if b == 1: return 1, 0, 1
		if a == 0: return b, 0, 1
		if b == 0: return a, 1, 0
		a1,a2 = self.as_tuple_b(a)
		b1,b2 = self.as_tuple_b(b)
		if a2 == b2 == 0:
			return orig_xgcd(a1, b1)
		if a1 == b1 == 0:
			d,s,t = orig_xgcd(a2, b2)
			B2 = (self.D + self.Droot) / 2
			return d * B2, s, t

		# a,b = map(self.field, [a,b])
		# gcd = 1
		# factors_a,factors_b = [dict(list(x.factor())) for x in [a,b]]
		# for q in factors_a.keys():
		# 	if q in factors_b:
		# 		e = min(factors_a[q], factors_b[q])
		# 		factors_a[q] -= e
		# 		factors_b[q] -= e
		# 		gcd *= q * e

		a,b = map(self.maxorder, [a,b])

		euc_a = self.euclidean_func(a)
		euc_b = self.euclidean_func(b)
		#assert abs_a != abs_b, "%r" % ((a,b,abs_a,abs_b),)
		if euc_a < euc_b:
			d,s,t = self.xgcd(b, a)
			return d,t,s
		# We have abs_b <= abs_a now.
		q,r = self.divmod(a, b)
		#assert q != 0
		#assert _simplify(abs(r)) < _simplify(abs(b))
		# q * b + r == a.
		assert q * b + r == a
		# => a - b*q == r
		d2,s2,t2 = self.xgcd(b, r)
		# d2 = b * s2 + r * t2
		assert d2 == b * s2 + r * t2
		# => d2 = b * s2 + (a - b*q) * t2
		# => d2 = a * t2 + b * (s2 - q * t2)
		d = d2
		s = _simplify(t2)
		t = _simplify(s2 - q * t2)
		assert d == a * s + b * t
		assert self.divides(a, d)
		assert self.divides(b, d)
		return d, s, t

	def gcd(self, a, b):
		d,_,_ = self.xgcd(a, b)
		return d
	def common_denom(self, *args):
		tupleargs = [None] * len(args) * 2
		for i in range(len(args)):
			tupleargs[2*i],tupleargs[2*i+1] = self.as_tuple_b(args[i])
		return matrix(QQ, 1,len(tupleargs), tupleargs).denominator()
	def as_tuple_b(self, a):
		real_part, b2 = self.DrootHalf_coordinates_in_terms_of_powers(a)
		b1 = real_part - b2 * self.D / 2
		return (b1, b2)
	def from_tuple_b(self, b1, b2):
		b1 = QQ(b1)
		b2 = QQ(b2)
		return self.field(b1 + b2 * (self.D + self.Droot) / 2)
	def __contains__(self, item):
		try:
			b1,b2 = self.as_tuple_b(item)
			b1,b2 = ZZ(b1), ZZ(b2)
		except TypeError: return False
		else: return True

def test_curlO():
	space = CurlO(-3)
	assert space.xgcd(1337,43) == orig_xgcd(1337,43)

class CurlOdual:
	def __init__(self, D):
		assert D < 0
		assert (D*D - D) % 4 == 0
		self.D = D
		self.field = QuadraticField(D)
		self.Droot = self.field(D).sqrt()
		self.DrootHalf_coordinates_in_terms_of_powers = (self.Droot / 2).coordinates_in_terms_of_powers()
	def from_tuple_b(self, b1, b2):
		b1, b2 = QQ(b1), QQ(b2)
		b = b1 / self.Droot + b2 * (1 + self.Droot) * QQ(0.5)
		return b

def M2T_Odual((a, b1, b2, c), D):
	space = CurlOdual(D)
	b = space.from_tuple_b(b1, b2)
	m = matrix(space.field, 2, 2, [a, b, b.conjugate(), c])
	m.set_immutable()
	return m


def solveR(M, S, space):
	"""
	Let M = [[a,b;c,d]] \in \SL_2(\ZZ).
	Let S \in \Her_2(\curlO) and S > 0.
	Define tM = [[a I, b S; c S^{-1}, d I]] \in \Sp_2(1/det(S) * \curlO).
	We find gamma \in \Sp_2(\curlO), R \in \Sp_2(\K)
	such that tM = gamma R.
	This algorithm is written in a way that it should work in other cases, too, though.
	"""
	assert isinstance(M, Matrix)
	assert isinstance(S, Matrix)
	assert S.nrows() == 2 and S.ncols() == 2
	assert M.nrows() == 2 and M.ncols() == 2
	assert _simplify(S[0][0]) > 0 and _simplify(S[1][1]) > 0 and _simplify(S.det()) > 0, "S is not positive definite"
	assert M.det() == 1, "M is not in \SL_2(\ZZ)"
	#Ring = S.base_ring()
	#print type(Ring), Ring
	#Ring = SymbolicRing()
	#Ring = I.base_ring() # Symbolic Ring
	Ring = space.field

	A1 = matrix(Ring, 2,2, M[0][0])
	B1 = matrix(Ring, 2,2, M[0][1] * S)
	C1 = matrix(Ring, 2,2, M[1][0] * S.inverse())
	D1 = matrix(Ring, 2,2, M[1][1])
	A1,B1,C1,D1 = [m.apply_map(_simplify) for m in [A1,B1,C1,D1]]
	def make4x4matrix(A1,B1,C1,D1):
		return matrix(Ring, 4,4,
			[A1[0][0],A1[0][1],B1[0][0],B1[0][1]] +
			[A1[1][0],A1[1][1],B1[1][0],B1[1][1]] +
			[C1[0][0],C1[0][1],D1[0][0],D1[0][1]] +
			[C1[1][0],C1[1][1],D1[1][0],D1[1][1]]
		)
	tM = tM1 = make4x4matrix(A1,B1,C1,D1)
	def make4x4matrix_embed(a1,a4,b1,b4,c1,c4,d1,d4):
		return matrix(Ring, 4,4,
			[a1,0,b1,0] +
			[0,a4,0,b4] +
			[c1,0,d1,0] +
			[0,c4,0,d4]
		)
	J = make4x4matrix_embed(0,0,-1,-1,1,1,0,0)
	I4 = make4x4matrix_embed(1,1,0,0,0,0,1,1)
	assert (tM1.conjugate_transpose() * J * tM1).apply_map(_simplify) == J
	l = space.common_denom(A1[0][0], C1[0][0])
	d,Ag11,Bg11 = space.xgcd(A1[0][0] * l, C1[0][0] * l)
	Cg11 = -C1[0][0] * l / d
	Dg11 = A1[0][0] * l / d
	assert Ag11 * Dg11 - Bg11 * Cg11 == 1, "{0}".format(tM1)
	assert all([x in space for x in [Ag11,Bg11,Cg11,Dg11]]), "%r" % ((A1[0][0],C1[0][0],l,d,(Ag11,Bg11,Cg11,Dg11),),)
	Dg14 = Ag14 = 0
	Bg14 = 1
	Cg14 = -1
	G1 = make4x4matrix_embed(Ag11,Ag14,Bg11,Bg14,Cg11,Cg14,Dg11,Dg14)
	tM2 = G1 * tM1
	assert tM2[2][0] == 0
	assert tM2[3][0] == 0
	c22,c24 = tM2[2][1],tM2[3][1]
	l = space.common_denom(c22,c24)
	d,Dg21,Dg22 = space.xgcd(c22 * l, c24 * l)
	if d == 0:
		G2 = I4
	else:
		Dg23 = -c24 * l / d
		Dg24 = c22 * l / d
		assert Dg21 * Dg24 - Dg22 * Dg23 == 1, "%r" % ([(c22,c24),(d,l),(Dg21, Dg22, Dg23, Dg24)],)
		assert all([x in space for x in [Dg21,Dg22,Dg23,Dg24]])
		Dg2 = matrix(Ring, 2,2, [Dg21,Dg22,Dg23,Dg24])
		assert Dg2.det() == 1
		Ag2 = Dg2.conjugate_transpose().inverse()
		G2 = make4x4matrix(Ag2,matrix(Ring,2,2,0),matrix(Ring,2,2,0),Dg2)
	tM3 = G2 * tM2
	assert tM3[2][0] == 0
	assert tM3[3][0] == 0
	assert tM3[3][1] == 0
	if tM3[2][1] == 0:
		G3 = I4
	else:
		assert tM3[0][0] == 0 # a_3,1 in our proof
		Cg34 = Bg34 = 0
		Ag34 = Dg34 = 1
		a32,c32 = tM3[0][1],tM3[2][1]
		l = space.common_denom(a32,c32)
		d,Ag31,Bg31 = space.xgcd(a32 * l, c32 * l)
		Cg31 = -c32 * l / d
		Dg31 = a32 * l / d
		assert Ag31 * Dg31 - Bg31 * Cg31 == 1, "%r" % ([(a32,c32),(d,l),(Ag31,Bg31,Cg31,Dg31)],)
		assert all([x in space for x in [Ag31,Bg31,Cg31,Dg31]])
		G3 = make4x4matrix_embed(Ag31,Ag34,Bg31,Bg34,Cg31,Cg34,Dg31,Dg34)
	tM4 = G3 * tM3
	assert tM4[2][0] == 0
	assert tM4[2][1] == 0
	assert tM4[3][0] == 0
	assert tM4[3][1] == 0

	R = tM4.apply_map(_simplify)
	gamma = (G1.inverse() * G2.inverse() * G3.inverse()).apply_map(_simplify) # G1,G2,G3 are in \Sp_2(\curlO).
	assert (tM - gamma * R).apply_map(_simplify) == 0, "\n%r ==\n%r *\n%r (=\n%r)" % (tM, gamma, R, (gamma * R).apply_map(_simplify))
	assert (gamma.conjugate_transpose() * J * gamma).apply_map(_simplify) == J
	assert (R.submatrix(0,0,2,2) * R.submatrix(2,2,2,2).conjugate_transpose()).apply_map(_simplify) == 1
	return gamma, R, tM

def test_solveR():
	space = CurlO(-3)
	a,b,c,d = 2,1,1,1
	s,t,u = 5,space.Droot,1
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 0,-1,1,0
	s,t,u = 1,0,2
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,2,1
	s,t,u = 1,0,4
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,3,1
	s,t,u = 1,2,16
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 0,-1,1,0
	s,t,u = 2, QQ(0.5) * space.Droot - QQ(0.5), 2
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,3,1
	s,t,u = 3, QQ(0.5) * space.Droot - QQ(1.5), 3
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,2,1
	s,t,u = 3, QQ(0.5) * space.Droot - QQ(1.5), 3
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 2,-1,3,-1
	s,t,u = 4, QQ(-0.5) * space.Droot + QQ(5/2.0), 4
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	return gamma,R,tM

def _curlO_matrix_denom(mat, D):
	space = CurlO(D)
	denom = space.common_denom(*mat.list())
	denom = int(ZZ(denom))
	for v in mat.list():
		assert v * denom in space, "%r (D=%r)" % (mat, D)
	return denom


def _toInt(a):
	a = _simplify(a)
	a = ZZ(a)
	a = int(a)
	return a


def takeEveryNRow(mat, n):
	"""
	INPUT:

	- `mat` -- a matrix with `mat.nrows()` being a multiple of `n`

	- `n` -- an integer

	OUTPUT:

	- A matrix which has only the rows [0,n,2*n,...] or the original
	  matrix `mat` in case every other rows are zero -
	  otherwise ``None``.

	"""

	assert mat.nrows() % n == 0, "%i, %i" % (mat.nrows(), n)
	newm = matrix(mat.base_ring(), mat.nrows() / n, mat.ncols())
	for i in range(mat.nrows()):
		if i % n == 0:
			newm[i / n] = mat[i]
		else:
			if mat[i] != 0:
				return None
	return newm


def toLowerCyclBase(ms, old_order, new_order):
	"""
	Let's

		K_old = CyclotomicField(old_order) ,
		K_new = CyclotomicField(new_order) .

	We transform the matrices in power base from `K_old` to `K_new`.

	The way this is implemented works only if `old_order`
	is a multiple of `new_order`.

	INPUT:

	- `ms` -- A list of matrices where every matrix is a factor to
	          `zeta_old**i` where `zeta_old = K_old.gen()`
	          and `len(ms) == K_old.degree()`.

	- `old_order` -- The order of `K_old`.

	- `new_order` -- The order of `K_new`.

	OUTPUT:

	- A list of matrices `new_ms` where every matrix is a factor to
	  `zeta_new**i` where `zeta_new = K_new.gen()` and
	  `len(new_ms) == K_new.degree()`.

	"""

	assert isinstance(ms, list) # list of matrices
	assert old_order % new_order == 0

	K_old = CyclotomicField(old_order)
	old_degree = int(ZZ(K_old.degree()))
	K_new = CyclotomicField(new_order)
	new_degree = int(ZZ(K_new.degree()))
	assert old_degree % new_degree == 0
	assert len(ms) == old_degree

	new_ms = [None] * new_degree
	for i in range(old_degree):
		i2,rem = divmod(i, old_degree / new_degree)
		if rem == 0:
			new_ms[i2] = ms[i]
		else:
			if ms[i] != 0:
				return None
	return new_ms


def toCyclPowerBase(M, order):
	"""
	Let's

		K = CyclotomicField(order).

	INPUT:

	- `M` -- A matrix over the cyclomotic field `K`.

	- `order` -- The order of `K`, the cyclomotic field.

	OUTPUT:

	- A list of matrices `ms` in power base where every matrix
	  is a factor to `zeta**i` where `zeta = K.gen()`
	  and `len(ms) == K.degree()`.

	"""

	K = CyclotomicField(order)
	zeta = K.gen()
	Kcoords = zeta.coordinates_in_terms_of_powers()

	assert len(K.power_basis()) == K.degree()
	ms = [matrix(QQ,M.nrows(),M.ncols()) for i in range(K.degree())]
	for y in range(M.nrows()):
		for x in range(M.ncols()):
			try:
				v_ = M[y,x]
				v = K(v_)
				coords = Kcoords(v)
			except TypeError:
				print "type of {1} ({2}) is not valid in Cyclomotic field of order {0}".format(order, M[y,x], type(M[y,x]))
				raise
			assert len(coords) == K.degree()
			for i in range(K.degree()):
				ms[i][y,x] = coords[i]
	return ms



def reduceNRow(denom, mats):
	"""
	INPUT:

	- `denom` -- an integer. it's a multiple of the nrows() of the matrices in `mats`.

	- `mats` -- a list of matrices of the same size.

	OUTPUT:

	- `(denom_new, mats_new)` where `denom_new` is a new integer which divides `denom`
	  and `mats_new` is a list of matrices of the same size, where

	      mats_new[0].nrows() == mats[0].nrows() / (denom / denom_new) .

	  It uses `takeEveryNRow()` for the reduction of rows.
	"""

	for (p,e) in Integer(denom).factor():
		assert e >= 0
		while e > 0:
			if mats[0].nrows() % p != 0: break
			mats_new = [takeEveryNRow(m, p) for m in mats]
			if all([m is not None for m in mats_new]):
				mats = mats_new
				denom = int(denom / p)
				e -= 1
			else:
				break
	return denom, mats


def addRows(mat, n):
	"""
	INPUT:

	- `mat` -- a matrix.
	- `n` -- an integer.

	OUTPUT:

	- A matrix with `mat.nrows() * n` rows. Every i*n-th row is the
	  i-th row of `mat`, every other row is zero.
	"""

	newm = matrix(mat.base_ring(), mat.nrows() * n, mat.ncols())
	for i in range(newm.nrows()):
		if i % n == 0:
			newm[i] = mat[i / n]
	return newm
