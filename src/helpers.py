from sage.functions.other import floor
from sage.matrix.constructor import matrix
from sage.modular.congroup import Gamma0
from sage.modular.modform.constructor import ModularForms
from sage.rings.arith import xgcd as orig_xgcd
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