from sage.matrix.constructor import matrix
from sage.matrix.matrix2 import Matrix
from sage.misc.misc import verbose
from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.modular.congroup import Gamma0
from sage.modular.modform.constructor import ModularForms
from sage.modules.free_module import FreeModule
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.symbolic.all import I
from sage.rings.arith import xgcd as orig_xgcd
from sage.rings.number_field.number_field import QQ, ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.symbolic.ring import SymbolicRing
from sage.symbolic.expression import Expression
import algo_cython as C

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

def xgcd(a,b):
	if a.imag() == 0 and b.imag() == 0:
		return orig_xgcd(a,b)
	if a.imag() != 0:
		if (I*a).imag() != 0: raise NotImplementedError
		d,p,q = xgcd(I*a,b)
		return d,I*p,q
	if b.imag() != 0:
		if (I*b).imag() != 0: raise NotImplementedError
		d,p,q = xgcd(a,I*b)
		return d,p,I*q
	assert False

def gcd(a,b):
	d,_,_ = xgcd(a, b)
	return d

def solveR(M, S):
	"""
	Let M = [[a,b;c,d]] \in \SL_2(\ZZ).
	Let S \in \Her_2(\curlO) and S > 0.
	Define tM = [[a I, b S; c S^{-1}, d I]] \in \Sp_2(\K).
	We find gamma \in \Sp_2(\curlO), R \in \Sp_2(\K)
	such that tM = gamma R.
	"""
	assert isinstance(M, Matrix)
	assert isinstance(S, Matrix)
	assert S.nrows() == 2 and S.ncols() == 2
	assert M.nrows() == 2 and M.ncols() == 2
	assert S[0][0] > 0 and S[1][1] > 0 and S.det() > 0, "S is not positive definite"
	assert M.det() == 1, "M is not in \SL_2(\ZZ)"
	#Ring = S.base_ring()
	#print type(Ring), Ring
	#Ring = SymbolicRing()
	Ring = I.base_ring() # Symbolic Ring

	A1 = matrix(Ring, 2,2, M[0][0])
	B1 = M[0][1] * S
	C1 = M[1][0] * S.inverse()
	D1 = matrix(Ring, 2,2, M[1][1])
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
	assert tM1.conjugate_transpose() * J * tM1 == J
	l = C1.denominator()
	d = gcd(A1[0][0] * l, C1[0][0] * l)
	if not d: d = 1
	Cg11 = C1[0][0] * l / d
	Dg11 = -A1[0][0] * l / d
	d,Ag11,Bg11 = xgcd(Dg11, -Cg11)
	assert d == 1, "{0}".format(tM1)
	Dg14 = Ag14 = 0
	Bg14 = 1
	Cg14 = -1
	G1 = make4x4matrix_embed(Ag11,Ag14,Bg11,Bg14,Cg11,Cg14,Dg11,Dg14)
	tM2 = G1 * tM1
	assert tM2[2][0] == 0
	assert tM2[3][0] == 0
	c22,c24 = tM2[2][1],tM2[3][1]
	d = gcd(c22,c24)
	if not d: d = 1
	Dg23 = c24 / d
	Dg24 = -c22 / d
	d,Dg21,Dg22 = xgcd(Dg24, -Dg23)
	if d == 0:
		G2 = I4
	else:
		assert d == 1, "{0}".format(tM2)
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
		d = gcd(a32, c32)
		if not d: d = 1
		Cg31 = c32 / d
		Dg31 = -a32 / d
		d,Ag31,Bg31 = xgcd(Dg31, -Cg31)
		assert d == 1
		G3 = make4x4matrix_embed(Ag31,Ag34,Bg31,Bg34,Cg31,Cg34,Dg31,Dg34)
	tM4 = G3 * tM3
	assert tM4[2][0] == 0
	assert tM4[2][1] == 0
	assert tM4[3][0] == 0
	assert tM4[3][1] == 0

	R = tM4
	gamma = G1.inverse() * G2.inverse() * G3.inverse()
	assert tM == gamma * R
	assert gamma.conjugate_transpose() * J * gamma == J
	assert R.submatrix(0,0,2,2) * R.submatrix(2,2,2,2).conjugate_transpose() == 1
	return gamma, R, tM

def test_solveR():
	a=2
	b=c=d=1
	s=5
	t=I
	u=1
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S)

	a=0
	b=-1
	c=1
	d=0
	s=1
	t=0
	u=2
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S)

	a,b,c,d=1,0,2,1
	s,t,u = 1,0,4
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S)

	return gamma,R,tM

# our own verbose function because I just want our msgs, not other stuff
def verbose(msg):
	print msg

def modform(D, HermWeight, B_cF=10):
	"Main algo"

	if HermWeight % 3 != 0:
		raise TypeError, "the modulform is trivial/zero if HermWeight is not divisible by 3"

	calc = C.Calc()
	calc.init(D = D, HermWeight = HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()
	reducedCurlFSize = calc.matrixColumnCount

	# Step 1. Iterate through square-free numbers l, starting at 1.
	# Init curlS = {}.
	curlS = []

	herm_modform_fe_expannsion = FreeModule(QQ, reducedCurlFSize)

	# Calculate the dimension of Hermitian modular form space.
	if D == -3:
		# dern2003graded, Thm 7
		R = PowerSeriesRing(ZZ, name="t", default_prec = HermWeight + 1)
		t = R.gen()
		dims = (1 + t**45) / (1 - t**4 ) / ( 1 - t**6 ) / ( 1 - t**9 ) / ( 1 - t**10 ) / ( 1 - t**12 )
		dim = dims[HermWeight]
	#elif D == -4:
		# dern2003graded, Corollary 9 and Lemma 3
		# TODO...
		#R = PowerSeriesRing(ZZ, name="t", default_prec = HermWeight + 1)
		#t = R.an_element()
	else:
		raise NotImplementedError, "dimension calculation of Hermitian modular form with D = %i not implemented" % D

	verbose("current dimension: %i, wanted: %i" % (herm_modform_fe_expannsion.dimension(), dim))

	while True:
		# Step 3. Iterate S \in Mat_2^T(Z). Add to curlS. iterate by denominator.
		# S_11 and S_22 (diagonal entries) are positive.
		# S positive definite.
		# S can be changed arbitrarily by GL(2, \ZZ).
		calc.curlS_clearMatrices()
		S = calc.getNextS()
		curlS += [S]
		verbose("trying S={0}, det={1}".format(S, S.det()))

		# Step 4. Calculate restriction matrix. Via calc.calcMatrix() (algo_cpp.cpp).
		# Note that calcMatrix() depends on the current internal calc.curlS set
		# and on the internal calc.curlF. curlF only depends on B_cF which is not changed here.
		verbose("calc restriction matrix...")
		M_S = calc.calcMatrix() # matrix over integer ring
		M_S = M_S.matrix_over_field() # matrix over rational field

		precLimit = M_S.nrows() # \cF(S)

		# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
		l = S.det()
		l = ZZ(l)
		mf = ModularForms(Gamma0(l), 2 * HermWeight)
		fe_expansion_matrix_l = matrix(QQ, [b.qexp(precLimit).padded_list(precLimit) for b in mf.basis()])
		fe_expansion_matrix_l.echelonize()
		#assert fe_expansion_matrix_l.rank() == precLimit, "{0} != {1}".format(fe_expansion_matrix_l.rank(), precLimit)

		# or:  fe_expansion_matrix[:n2,:].row_module()
		ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()

		verbose("calc M_S_module...")
		M_S_module = M_S.column_module()
		restriction_fe_expansions = ell_modform_fe_expansions_l.intersection( M_S_module )
		herm_modform_fe_expannsion_S = M_S.solve_right( restriction_fe_expansions.basis_matrix().transpose() )
		herm_modform_fe_expannsion_S_module = herm_modform_fe_expannsion_S.column_module()
		verbose("calc M_S_right_kernel...")
		M_S_right_kernel = M_S.right_kernel()
		herm_modform_fe_expannsion_S_module += M_S_right_kernel

		verbose("intersecting herm_modform_fe_expannsion...")
		herm_modform_fe_expannsion = herm_modform_fe_expannsion.intersection( herm_modform_fe_expannsion_S_module )
		current_dimension = herm_modform_fe_expannsion.dimension()
		verbose("current dimension: %i, wanted: %i" % (current_dimension, dim))
		assert current_dimension >= dim

		# Step 5. dimension check
		if dim == current_dimension:
			break

		# cusp info:
		ce = cusp_expansions.ModularFormsCuspExpansions._for_modular_forms(l, HermWeight*2)
		for cusp in Gamma0(l).cusps():
			if cusp == Infinity: continue
			if cusp == 0:
				M = matrix(ZZ,2,2,[0,-1,1,0])
			else:
				a = cusp.numerator()
				c = cusp.denominator()
				_div, d, b = xgcd(a, -c)
				assert _div == 1
				M = matrix(ZZ,2,2,[a,b,c,d])
				del a,b,c,d

			try:
				gamma, R, tM = solveR(M, S)
			except:
				print repr((M, S))
				raise
			tS = tM.submatrix(0,0,2,2)
			tT = tM.submatrix(2,0,2,2)
			ms = calc.calcMatrixTrans(tS * l, tT * l, l)
			print ms

			M = SL2Z(M)
			for f in herm_modform_fe_expannsion_S_module.gens():
				g = M_S * f
				if g == 0: continue
				g_inbase = fe_expansion_matrix_l.solve_left(g)
				# g in ModularForms(l, 2 k), M = [[a,b;c,d]] the cusp representation with c = M \infty.
				# this calculates f|M
				f_M = ce.expansion_at(M, g_inbase)

				print f_M



		# Step 5. dimension check
		if dim == current_dimension:
			break


	# TODO:
	# Otherwise, reconstruct fourier expansion.
	return herm_modform_fe_expannsion
