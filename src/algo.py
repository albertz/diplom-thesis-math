import sage
from sage.calculus.functional import simplify
from sage.functions.other import sqrt as ssqrt
from sage.functions.other import imag, real
from sage.matrix.constructor import matrix
from sage.matrix.matrix2 import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.misc import verbose
from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.modular.congroup import Gamma0
from sage.modular.modform.constructor import ModularForms
from sage.modules.free_module import FreeModule
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject
from sage.structure.sequence import Sequence_generic, Sequence
from sage.symbolic.all import I
from sage.rings.arith import xgcd as orig_xgcd, lcm
from sage.rings.number_field.number_field import QQ, ZZ, CyclotomicField
from sage.rings.power_series_ring import PowerSeriesRing
from sage.symbolic.ring import SymbolicRing
from sage.symbolic.expression import Expression
import algo_cython as C

# via Martin. while this is not in Sage:
import cusp_expansions

# for debugging
import better_exchook


# It seems that Sage load/save uses the standard pickle module.
# The standard pickle cannot save the Sage Expression objects for some reason.
# We extend the standard pickler.
import pickle, types, marshal, sys
CellType = type((lambda x: lambda: x)(0).func_closure[0])
def makeCell(value): return (lambda: value).func_closure[0]
def getModuleDict(modname): return __import__(modname).__dict__
class Pickler(pickle.Pickler):
	def __init__(self, *args, **kwargs):
		if not "protocol" in kwargs:
			kwargs["protocol"] = pickle.HIGHEST_PROTOCOL
		pickle.Pickler.__init__(self, *args, **kwargs)
	dispatch = pickle.Pickler.dispatch.copy()

	def save_func(self, obj):
		try:
			self.save_global(obj)
			return
		except pickle.PicklingError:
			pass
		assert type(obj) is types.FunctionType
		self.save(types.FunctionType)
		self.save((
			obj.func_code,
			obj.func_globals,
			obj.func_name,
			obj.func_defaults,
			obj.func_closure,
			))
		self.write(pickle.REDUCE)
		self.memoize(obj)
	dispatch[types.FunctionType] = save_func

	def save_code(self, obj):
		assert type(obj) is types.CodeType
		self.save(marshal.loads)
		self.save((marshal.dumps(obj),))
		self.write(pickle.REDUCE)
		self.memoize(obj)
	dispatch[types.CodeType] = save_code

	def save_cell(self, obj):
		assert type(obj) is CellType
		self.save(makeCell)
		self.save((obj.cell_contents,))
		self.write(pickle.REDUCE)
		self.memoize(obj)
	dispatch[CellType] = save_cell

	# We also search for module dicts and reference them.
	def intellisave_dict(self, obj):
		if len(obj) <= 5: # fastpath
			self.save_dict(obj)
			return
		for modname, mod in sys.modules.iteritems():
			if not mod: continue
			moddict = mod.__dict__
			if obj is moddict:
				self.save(getModuleDict)
				self.save((modname,))
				self.write(pickle.REDUCE)
				self.memoize(obj)
				return
		self.save_dict(obj)
	dispatch[types.DictionaryType] = intellisave_dict

	# Some types in the types modules are not correctly referenced,
	# such as types.FunctionType. This is fixed here.
	def fixedsave_type(self, obj):
		try:
			self.save_global(obj)
			return
		except pickle.PicklingError:
			pass
		for modname in ["types"]:
			moddict = sys.modules[modname].__dict__
			for modobjname,modobj in moddict.iteritems():
				if modobj is obj:
					self.write(pickle.GLOBAL + modname + '\n' + modobjname + '\n')
					self.memoize(obj)
					return
		self.save_global(obj)
	dispatch[types.TypeType] = fixedsave_type

	# Wrap _batch_setitems (e.g. for dicts) so that our representations stays fixed
	# (the order of dict.keys() can be different at each run).
	orig_batch_setitems = pickle.Pickler._batch_setitems
	def _batch_setitems(self, items):
		items = sorted(items)
		self.orig_batch_setitems(iter(items))

	# Wrap save_reduce so that we can catch a few cases (e.g. set)
	# to fix up the representation so that it stays fixed (as for dicts).
	orig_save_reduce = pickle.Pickler.save_reduce
	def save_reduce(self, func, args, state=None, listitems=None, dictitems=None, obj=None):
		if func is set:
			assert len(args) == 1
			args = (sorted(args[0]),)
		self.orig_save_reduce(func=func, args=args, state=state, listitems=listitems, dictitems=dictitems, obj=obj)

	# avoid pickling instances of ourself. this mostly doesn't make sense and leads to trouble.
	# however, also doesn't break. it mostly makes sense to just ignore.
	def __getstate__(self): return None
	def __setstate__(self, state): pass

# The extended Pickler above writes pickle-data such that the standard unpickler
# should be able to unpickle it.
# However, Sage has problems to generate some objects, so we catch those and automatically
# create equivalent structures.
class Unpickler(pickle.Unpickler):
	dispatch = dict(pickle.Unpickler.dispatch)

	def create(self, cls, args):
		if cls is Sequence_generic:
			return Sequence(list(*args))
		return cls.__new__(cls, *args)

	def wrapped_load_newobj(self):
		args = self.stack.pop()
		cls = self.stack[-1]
		obj = self.create(cls, args)
		self.stack[-1] = obj
	dispatch[pickle.NEWOBJ] = wrapped_load_newobj



# We dont use these:
#from sage.misc.db import save, load
# But these:
def save(obj, filename):
	f = open(filename, "w")
	pickler = Pickler(f)
	pickler.dump(obj)
def load(filename):
	f = open(filename)
	unpickler = Unpickler(f)
	obj = unpickler.load()
	return obj

class PersistentCache:
	def __init__(self, name):
		self.name = name
		self.dict = {}
		self.load()
	def __getitem__(self, item):
		return self.dict[item]
	def __setitem__(self, key, value):
		self.dict[key] = value
		self.save()
	def __contains__(self, item):
		return item in self.dict
	def load(self):
		try:
			self.dict = load(self.name)
		except IOError:
			pass
		except Exception:
			better_exchook.better_exchook(*sys.exc_info())
			raise
		else:
			assert isinstance(self.dict, dict)
	def save(self):
		try:
			save(self.dict, self.name)
		except Exception:
			print self.name, self.dict
			raise

# our own verbose function because I just want our msgs, not other stuff
def verbose(msg):
	print msg

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
	del l, d
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
	l = matrix(1,2,(c22,c24)).denominator()
	d = gcd(c22 * l, c24 * l)
	if not d: d = 1
	Dg23 = c24 * l / d
	Dg24 = -c22 * l / d
	del l, d
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
		l = matrix(1,2,(a32,c32)).denominator()
		d = gcd(a32 * l, c32 * l)
		if not d: d = 1
		Cg31 = c32 * l / d
		Dg31 = -a32 * l / d
		del l, d
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
	a,b,c,d = 2,1,1,1
	s,t,u = 5,I,1
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S)

	a,b,c,d = 0,-1,1,0
	s,t,u = 1,0,2
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S)

	a,b,c,d = 1,0,2,1
	s,t,u = 1,0,4
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S)

	a,b,c,d = 1,0,3,1
	s,t,u = 1,2,16
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S)

	return gamma,R,tM

def _curlO_matrix_denom(mat, D):
	assert D < 0
	mat_real = mat.apply_map(real)
	mat_imag = mat.apply_map(lambda x: simplify(imag(x) * 2 / ssqrt(-D)))
	denom_mat_real = ZZ(mat_real.denominator())
	denom_mat_imag = ZZ(mat_imag.denominator())
	assert denom_mat_real * mat_real in MatrixSpace(ZZ,mat.nrows(),mat.ncols())
	assert denom_mat_imag * mat_imag in MatrixSpace(ZZ,mat.nrows(),mat.ncols())
	denom = lcm(denom_mat_real, denom_mat_imag)
	return int(ZZ(denom))

matrixTransCache = PersistentCache("matrixTrans.cache.sobj")
def calcMatrixTrans(calc, R):
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

	ms = calc.calcMatrixTrans(tS, tT, lS, lT)

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

	matrixTransCache[cacheIdx] = calc.matrixRowDenomTrans, order, new_ms
	return calc.matrixRowDenomTrans, order, new_ms

def calcElliptViaReduct(calc, f, R, l):
	tS = R.submatrix(0,0,2,2)
	tT = R.submatrix(2,0,2,2)
	try:
		ms = calc.calcMatrixTrans(tS * l, tT * l, l)
	except Exception:
		print (calc.params, calc.curlS, R * l, l)
		raise
	matrixCountTrans = calc.matrixCountTrans
	denom = calc.matrixRowDenomTrans

	g = 0
	K = CyclotomicField(matrixCountTrans)
	zeta, = K.gens()
	for i in range(matrixCountTrans):
		g += ms[i] * f * (zeta ** i)
	return denom, g


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
def getElliptModule(level, weight, precision):
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
	mat = calc.calcMatrix()
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
	if dim == 0:
		print "dim == 0 -> exit"
		return

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
		M_S = calcRestrictMatrix(calc) # matrix over integer ring
		M_S = M_S.matrix_over_field() # matrix over rational field
		#print M_S

		precLimit = M_S.nrows() # \cF(S)

		l = S.det()
		l = ZZ(l)
		# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
		verbose("get elliptic modform space with precision %i ..." % precLimit)
		fe_expansion_matrix_l = getElliptModule(l, 2*HermWeight, precLimit)
		ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()
		#print ell_modform_fe_expansions_l

		verbose("calc M_S_module...")
		M_S_module = M_S.column_module()
		verbose("dimension of M_S column module: %i" % M_S_module.dimension())
		#print M_S_module
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

		# Step 5. dimension check
		if dim == current_dimension:
			break

		# cusp info:
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
			except Exception:
				print (M, S)
				raise
			R.set_immutable() # for caching, we need it hashable

			herm_modforms = herm_modform_fe_expannsion.echelonized_basis_matrix().transpose()
			ell_R_denom, ell_R_order, M_R = calcMatrixTrans(calc, R)
			CycloDegree_R = CyclotomicField(ell_R_order).degree()
			print "M_R[0] nrows, ell_R_denom, ell_R_order, Cyclo degree:", \
				M_R[0].nrows(), ell_R_denom, ell_R_order, CycloDegree_R

			ce = cuspExpansions(level=l, weight=2*HermWeight, prec=M_R[0].nrows())
			ell_M_denom, ell_M = ce.expansion_at(SL2Z(M))
			ell_M_order = ell_M_denom # we expect that a CyclomoticField of the order of the denom can represent all entries

			print "1) M_R[0] rank, herm rank, mult rank:", \
				M_R[0].rank(), herm_modforms.rank(), (M_R[0] * herm_modforms).rank()

			# Not sure if this is always the case but seems so.
			assert ell_R_denom >= ell_M_denom
			if ell_R_denom > ell_M_denom:
				assert ell_R_denom % ell_M_denom == 0
				M_R = [_takeEveryNRow(M_R_i, ell_R_denom / ell_M_denom) for M_R_i in M_R]
				assert all([M_R_i is not None for M_R_i in M_R])
				ell_R_denom = ell_M_denom
			assert ell_R_denom == ell_M_denom

			print "2) M_R[0] rank, herm rank, mult rank:", \
				M_R[0].rank(), herm_modforms.rank(), (M_R[0] * herm_modforms).rank()

			# ell_M rows are the elliptic FE. M_R[i] columns are the elliptic FE.
			# We expect that M_R gives a higher precision for the ell FE. I'm not sure
			# if this is always true but we expect it here (maybe not needed, though).
			print "precision of M_R[0], ell_M:", M_R[0].nrows(), ell_M.ncols()
			assert ell_M.ncols() >= M_R[0].nrows()
			prec = min(M_R[0].nrows(), ell_M.ncols())
			# cut to have same precision
			M_R = [M_R_i[:prec,:] for M_R_i in M_R]
			ell_M = ell_M[:,:prec]
			assert ell_M.ncols() == M_R[0].nrows() == prec

			print "3) M_R[0] rank, herm rank, mult rank:", \
				M_R[0].rank(), herm_modforms.rank(), (M_R[0] * herm_modforms).rank()
			ell_R = [M_R_i * herm_modforms for M_R_i in M_R]

			# I'm not sure on this. Seems to be true and it simplifies things in the following.
			assert ell_R_order % ell_M_order == 0, "{0}".format((ell_M_order, ell_R_order))

			# Transform to same Cyclomotic Field in same power base.
			ell_M2 = toCyclPowerBase(ell_M, ell_M_order)
			ell_R2 = toLowerCyclBase(ell_R, ell_R_order, ell_M_order)
			# We must work with the matrix. maybe we should transform hf_M instead to a
			# higher order field instead, if this ever fails (I'm not sure).
			assert ell_R2 is not None
			assert len(ell_M2) == len(ell_R2) # They should have the same power base & same degree now.
			print "ell_M2[0], ell_R2[0] rank with order %i:" % ell_M_order, ell_M2[0].rank(), ell_R2[0].rank()

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

		if dim == current_dimension:
			break

	# TODO:
	# Otherwise, reconstruct fourier expansion.
	return herm_modform_fe_expannsion
