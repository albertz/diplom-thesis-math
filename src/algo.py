from time import time
from sage.calculus.functional import simplify
from sage.functions.other import sqrt as ssqrt
from sage.functions.other import imag, real
from sage.matrix.constructor import matrix
from sage.matrix.matrix2 import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.modular.congroup import Gamma0
from sage.modular.modform.constructor import ModularForms
from sage.modules.free_module import FreeModule
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject
from sage.symbolic.all import I
from sage.rings.arith import xgcd as orig_xgcd, lcm
from sage.rings.number_field.number_field import QQ, ZZ, CyclotomicField
from sage.rings.power_series_ring import PowerSeriesRing
from sage.symbolic.ring import SymbolicRing
from sage.symbolic.expression import Expression
import algo_cython as C
from utils import *
from utils import _simplify, _curlO_matrix_denom # seems the above does not import "_"-prefixed symbols

# via Martin. while this is not in Sage:
import cusp_expansions

def reloadC():
	"""
	This is just for testing to reload the C (Cython) module
	after it was recompiled.
	Note that this code is very unsafe! It will likely crash
	when you have other references on the C module.
	This is only for debugging and development!
	"""
	global C
	import ctypes
	try:
		libdl = ctypes.CDLL("libdl.so")
	except Exception:
		# MacOSX:
		libdl = ctypes.CDLL("libdl.dylib")
	libdl.dlclose.argtypes = [ctypes.c_void_p]
	so = ctypes.PyDLL(C.__file__)
	assert(libdl.dlclose(so._handle) == 0)
	reload(C)

def test_algo_calcMatrix():
	calc = C.Calc()
	calc.init(D = -4, HermWeight = 10)

	calc.getNextS()
	calc.getNextS()

	calc.calcMatrix()
	return calc.getMatrix()


matrixTransCache = PersistentCache("matrixTrans.cache.sobj")
def calcMatrixTrans(calc, R):
	"""
	Returns a triple (denom, order, ms) where
	ms is a list of matrices. Each matrix is a factor to
	zeta**i where zeta is the order's root of unity,
	i.e. zeta = CyclotomicField(order).gen().
	And we have len(ms) == CyclotomicField(order).degree().

	The matrices represent the linear map of a Hermitian
	modular form f to the Elliptic modular form (f|R)[S].
	(R is given as a parameter to this function. S is saved
	internally by the calc structure.)

	This uses the C++ function `Calc::calcMatrixTrans()` in  algo_cpp.cpp`.
	"""
	tS = R.submatrix(0,0,2,2)
	tT = R.submatrix(2,0,2,2)
	lS = _curlO_matrix_denom(tS, D=calc.D)
	lT = _curlO_matrix_denom(tT, D=calc.D)
	tS *= lS
	tT *= lT
	tS.set_immutable()
	tT.set_immutable()
	cacheIdx = (calc.params, calc.curlS, tS, tT, lS, lT)
	if cacheIdx in matrixTransCache:
		return matrixTransCache[cacheIdx]

	t = time()
	try:
		ms = calc.calcMatrixTrans(tS, tT, lS, lT)
	except Exception:
		print (cacheIdx, R)
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

	if time() - t > 2.0:
		print "calculation of matrixTrans took %f secs" % (time() - t)
		matrixTransCache[cacheIdx] = calc.matrixRowDenomTrans, order, new_ms
	return calc.matrixRowDenomTrans, order, new_ms


cuspExpansionsCache = PersistentCache("cuspExpansions.cache.sobj")
def cuspExpansions(level, weight, prec):
	cacheIdx = (level, weight)
	if cacheIdx in cuspExpansionsCache:
		ce_prec,ce = cuspExpansionsCache[cacheIdx]
		if ce_prec >= prec: return ce
	verbose("calc ModularFormsCuspExpansions at level %i with weight %i and prec %i ..." % (level, weight, prec))
	ce = cusp_expansions.ModularFormsCuspExpansions._for_modular_forms(level, weight, prec)
	cuspExpansionsCache[cacheIdx] = prec, ce
	return ce

ellipBaseMatrixCache = PersistentCache("ellipBaseMatrix.cache.sobj") # level,weight -> mat,prec
def getElliptModFormsBasisMatrix(level, weight, precision):
	cacheIdx = (level, weight)
	if cacheIdx in ellipBaseMatrixCache and ellipBaseMatrixCache[cacheIdx][1] >= precision:
		return ellipBaseMatrixCache[cacheIdx][0][:,:precision]
	#n = 2
	#while n < precision:
	#	n **= 2
	n = precision
	mf = ModularForms(Gamma0(level), weight)
	fe_expansion_matrix_l = matrix(QQ, [b.qexp(n).padded_list(n) for b in mf.basis()])
	fe_expansion_matrix_l.echelonize()
	ellipBaseMatrixCache[cacheIdx] = (fe_expansion_matrix_l, n)
	return fe_expansion_matrix_l[:,:precision]


restrMatrixCache = PersistentCache("restrMatrix.cache.sobj") # by (calc.params,calc.curlS)
def calcRestrictMatrix(calc):
	cacheIdx = (calc.params, calc.curlS)
	if cacheIdx in restrMatrixCache:
		return restrMatrixCache[cacheIdx]
	t = time()
	mat = calc.calcMatrix()
	if time() - t > 2.0:
		print "calculation of restriction matrix took %f secs" % (time() - t)
		restrMatrixCache[cacheIdx] = mat
	return mat


def toLowerCyclBase(ms, old_order, new_order):
	# We expect to have ms in power_base.
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

def _takeEveryNRow(mat, n):
	assert mat.nrows() % n == 0
	newm = matrix(mat.base_ring(), mat.nrows() / n, mat.ncols())
	for i in range(mat.nrows()):
		if i % n == 0:
			newm[i / n] = mat[i]
		else:
			if mat[i] != 0:
				return None
	return newm

def _toInt(a):
	a = _simplify(a)
	a = ZZ(a)
	a = int(a)
	return a

def test_herm_modform_space(calc, herm_modform_space, used_curlS_denoms, testSCount = 3):
	"""
	It uses the C++ calc structure to search for additional S matrices
	which have other denominators than those in used_curlS_denoms.
	For testScount matrices with unique denominators, it calculates
	the restriction matrix via the C++ calc structure for f \mapsto f[S].
	When mapping herm_modform_space, we must only get Elliptic modular forms.
	We check whether they are in ModularForms(\Gamma_0(l), prec).
	"""

	HermWeight = calc.HermWeight
	curlS_denoms = set(used_curlS_denoms)

	while testSCount > 0:
		calc.curlS_clearMatrices()
		S = calc.getNextS()
		l = S.det()
		l = _toInt(l)
		if l in curlS_denoms: continue

		curlS_denoms.add(l)
		testSCount -= 1
		verbose("testing with S={0}, det={1}".format(S, l))

		verbose("calc restriction matrix...")
		M_S = calcRestrictMatrix(calc) # matrix over integer ring
		M_S = M_S.matrix_over_field() # matrix over rational field

		precLimit = M_S.nrows() # \cF(S)

		# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
		verbose("get elliptic modform space with precision %i ..." % precLimit)
		fe_expansion_matrix_l = getElliptModFormsBasisMatrix(l, 2*HermWeight, precLimit)
		ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()

		verbose("calc M_S * herm_modforms ...")
		m = M_S * herm_modform_space.basis_matrix().transpose()

		m_module = m.column_module()
		assert m_module.is_subspace(ell_modform_fe_expansions_l), \
			"%r not subspace of %r" % (m_module, ell_modform_fe_expansions_l)


def herm_modform_indexset(D, B_cF):
	"""
	This is the precision of the index set for the Fourier coefficients of the
	Hermitian modular forms (calculated by `herm_modform_space()`).

	The index set \Lambda is the set of positive definite 2*2 matrices over \curlO^#.
	\curlF \subset \Lambda is a precision such that for [a,b,c] \in \curlF,
	we have 0 \le a,c, \le B_cF.

	This function returns all the reduced matrices of \curlF. The reduction
	is given by `reduce_GL()` in `reduceGL.hpp`.

	The order of \curlF is predefined and fixed. It is such that
	[a,b,c] < [s,t,u]  <==>  max(a,c) < max(s,u).
	Thus, for B1 < B2, we have

		L1 = herm_modform_indexset(D, B1)
		L2 = herm_modform_indexset(D, B2)
		assert len(L1) < len(L2)
		assert L1 == L2[:len(L1)]
	"""

	HermWeight = 6 # Just some dummy value. The calc-init wants something valid.
	calc = C.Calc()
	calc.init(D=D, HermWeight=HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()
	return calc.getReducedCurlF()


def herm_modform_space_dim(D, HermWeight):
	"""
	Calculates and returns the dimension of the vector space
	of Hermitian modular forms of weight `HermWeight` over \Gamma,
	where \Gamma = \Sp_2(\curlO) and \curlO is the maximal order
	of \QQ(\sqrt{D}).
	"""

	if D == -3:
		# dern2003graded, Thm 7
		R = PowerSeriesRing(ZZ, name="t", default_prec = HermWeight + 1)
		t = R.gen()
		dims = (1 + t**45) / (1 - t**4 ) / ( 1 - t**6 ) / ( 1 - t**9 ) / ( 1 - t**10 ) / ( 1 - t**12 )
		return dims[HermWeight]
	#elif D == -4:
		# dern2003graded, Corollary 9 and Lemma 3
		# TODO...
		#R = PowerSeriesRing(ZZ, name="t", default_prec = HermWeight + 1)
		#t = R.an_element()
	else:
		raise NotImplementedError, "dimension calculation of Hermitian modular form with D = %i not implemented" % D


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


def _intersect_modform_cusp_info(calc, S, l, precLimit, herm_modform_fe_expannsion):
	assert l == S.det()
	assert list(calc.curlS) == [S]
	D = calc.D
	HermWeight = calc.HermWeight

	for cusp in Gamma0(l).cusps():
		if cusp == Infinity: continue
		M = cusp_matrix(cusp)

		try:
			gamma, R, tM = solveR(M, S, space=CurlO(D))
		except Exception:
			print (M, S)
			raise
		R.set_immutable() # for caching, we need it hashable

		herm_modforms = herm_modform_fe_expannsion.echelonized_basis_matrix().transpose()
		ell_R_denom, ell_R_order, M_R = calcMatrixTrans(calc, R)
		CycloDegree_R = CyclotomicField(ell_R_order).degree()
		print "M_R[0] nrows, ell_R_denom, ell_R_order, Cyclo degree:", \
			M_R[0].nrows(), ell_R_denom, ell_R_order, CycloDegree_R

		# The maximum precision we can use is M_R[0].nrows().
		# However, that can be quite huge (e.g. 600).
		ce_prec = min(precLimit, M_R[0].nrows())

		ce = cuspExpansions(level=l, weight=2*HermWeight, prec=ce_prec)
		ell_M_denom, ell_M = ce.expansion_at(SL2Z(M))
		ell_M_order = ell_R_order # not sure here. just try the one from R. toCyclPowerBase would fail if this doesn't work
		# CyclotomicField(l / prod(l.prime_divisors())) should also work.

		# Not sure if this is always the case but seems so.
		assert ell_R_denom >= ell_M_denom
		if ell_R_denom > ell_M_denom:
			assert ell_R_denom % ell_M_denom == 0
			M_R = [_takeEveryNRow(M_R_i, ell_R_denom / ell_M_denom) for M_R_i in M_R]
			assert all([M_R_i is not None for M_R_i in M_R])
			ell_R_denom = ell_M_denom
		assert ell_R_denom == ell_M_denom

		# ell_M rows are the elliptic FE. M_R[i] columns are the elliptic FE.
		# We expect that M_R gives a higher precision for the ell FE. I'm not sure
		# if this is always true but we expect it here (maybe not needed, though).
		print "precision of M_R[0], ell_M, wanted:", M_R[0].nrows(), ell_M.ncols(), ce_prec
		assert ell_M.ncols() >= ce_prec
		prec = min(M_R[0].nrows(), ell_M.ncols())
		# cut to have same precision
		M_R = [M_R_i[:prec,:] for M_R_i in M_R]
		ell_M = ell_M[:,:prec]
		assert ell_M.ncols() == M_R[0].nrows() == prec

		print "M_R[0] rank, herm rank, mult rank:", \
			M_R[0].rank(), herm_modforms.rank(), (M_R[0] * herm_modforms).rank()
		ell_R = [M_R_i * herm_modforms for M_R_i in M_R]

		# I'm not sure on this. Seems to be true and it simplifies things in the following.
		assert ell_M_order <= ell_R_order, "{0}".format((ell_M_order, ell_R_order))
		assert ell_R_order % ell_M_order == 0, "{0}".format((ell_M_order, ell_R_order))

		# Transform to same Cyclomotic Field in same power base.
		ell_M2 = toCyclPowerBase(ell_M, ell_M_order)
		ell_R2 = toLowerCyclBase(ell_R, ell_R_order, ell_M_order)
		# We must work with the matrix. maybe we should transform hf_M instead to a
		# higher order field instead, if this ever fails (I'm not sure).
		assert ell_R2 is not None
		assert len(ell_M2) == len(ell_R2) # They should have the same power base & same degree now.
		print "ell_M2[0], ell_R2[0] rank with order %i:" % ell_M_order, ell_M2[0].rank(), ell_R2[0].rank()

		assert len(M_R) == len(ell_M2)
		for i in range(len(ell_M2)):
			ell_M_space = ell_M2[i].row_space()
			ell_R_space = ell_R2[i].column_space()
			merged = ell_M_space.intersection(ell_R_space)

			herm_modform_fe_expannsion_Ci = M_R[i].solve_right( merged.basis_matrix().transpose() )
			herm_modform_fe_expannsion_Ci_module = herm_modform_fe_expannsion_Ci.column_module()
			herm_modform_fe_expannsion_Ci_module += M_R[i].right_kernel()

			herm_modform_fe_expannsion = herm_modform_fe_expannsion.intersection( herm_modform_fe_expannsion_Ci_module )
			print "power", i, merged.dimension(), herm_modform_fe_expannsion_Ci_module.dimension(), \
				current_dimension, herm_modform_fe_expannsion.dimension()
			current_dimension = herm_modform_fe_expannsion.dimension()


def herm_modform_space(D, HermWeight, B_cF=10):
	"""
	This calculates the vectorspace of Fourier expansions to
	Hermitian modular forms of weight `HermWeight` over \Gamma,
	where \Gamma = \Sp_2(\curlO) and \curlO is the maximal order
	of \QQ(\sqrt{D}).

	Each Fourier coefficient vector is indexed up to a precision
	\curlF which is given by `B_cF` such that for every
	[a,b,c] \in \curlF \subset \Lambda, we have 0 \le a,c \le B_cF.

	The function `herm_modform_indexset()` returns reduced matrices
	of that precision index set \curlF.
	"""

	if HermWeight % 3 != 0:
		raise TypeError, "the modulform is trivial/zero if HermWeight is not divisible by 3"

	calc = C.Calc()
	calc.init(D = D, HermWeight = HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()
	reducedCurlFSize = calc.matrixColumnCount

	# Calculate the dimension of Hermitian modular form space.
	dim = herm_modform_space_dim(D=D, HermWeight=HermWeight)

	herm_modform_fe_expannsion = FreeModule(QQ, reducedCurlFSize)
	curlS = [] # all matrices S we have visited so far
	curlS_denoms = set() # the denominators of the visited matrices S

	verbose("current dimension: %i, wanted: %i" % (herm_modform_fe_expannsion.dimension(), dim))
	if dim == 0:
		print "dim == 0 -> exit"
		return

	# Iterate S \in Mat_2^T(\curlO), S > 0.
	while True:
		# Get the next S.
		calc.curlS_clearMatrices() # In the C++ internal curlS, clear previous matrices.
		S = calc.getNextS()
		l = S.det()
		l = _toInt(l)
		curlS += [S]
		curlS_denoms.add(l)
		verbose("trying S={0}, det={1}".format(S, l))

		# Step 4. Calculate restriction matrix. Via calc.calcMatrix() (algo_cpp.cpp).
		# Note that calcMatrix() depends on the current internal calc.curlS set
		# and on the internal calc.curlF. curlF only depends on B_cF which is not changed here.
		verbose("calc restriction matrix...")
		M_S = calcRestrictMatrix(calc) # matrix over integer ring
		M_S = M_S.matrix_over_field() # matrix over rational field

		# The maximum precision of Elliptic modular forms is given in
		# the text by \cF(S). This is also the number of rows of M_S.
		precLimit = M_S.nrows()

		# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
		verbose("get elliptic modform space with precision %i ..." % precLimit)
		fe_expansion_matrix_l = getElliptModFormsBasisMatrix(l, 2*HermWeight, precLimit)
		ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()

		verbose("calc M_S_module...")
		M_S_module = M_S.column_module()
		verbose("dimension of M_S column module: %i" % M_S_module.dimension())
		restriction_fe_expansions = ell_modform_fe_expansions_l.intersection( M_S_module )
		verbose("dimension of restriction_fe_expansions: %i" % restriction_fe_expansions.dimension())
		herm_modform_fe_expannsion_S = M_S.solve_right( restriction_fe_expansions.basis_matrix().transpose() )
		herm_modform_fe_expannsion_S_module = herm_modform_fe_expannsion_S.column_module()
		verbose("dimension of herm column module: %i" % herm_modform_fe_expannsion_S_module.dimension())
		verbose("calc M_S_right_kernel...")
		M_S_right_kernel = M_S.right_kernel()
		verbose("dimension of M_S right kernel: %i" % M_S_right_kernel.dimension())
		herm_modform_fe_expannsion_S_module += M_S_right_kernel

		verbose("intersecting herm_modform_fe_expannsion...")
		herm_modform_fe_expannsion = herm_modform_fe_expannsion.intersection( herm_modform_fe_expannsion_S_module )
		current_dimension = herm_modform_fe_expannsion.dimension()
		verbose("current dimension: %i, wanted: %i" % (current_dimension, dim))
		assert current_dimension >= dim

		if dim == current_dimension:
			break

		_intersect_modform_cusp_info(
			calc=calc, S=S, l=l, precLimit=precLimit*2,
			herm_modform_fe_expannsion=herm_modform_fe_expannsion)

		if dim == current_dimension:
			break

	# Test for some other S with other not-yet-seen denominator.
	test_herm_modform_space(
		calc, herm_modform_space=herm_modform_fe_expannsion,
		used_curlS_denoms=curlS_denoms
		)

	return herm_modform_fe_expannsion
