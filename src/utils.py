from sage.calculus.functional import simplify
from sage.functions.other import sqrt as ssqrt
from sage.functions.other import imag, real
from sage.matrix.constructor import matrix
from sage.matrix.matrix2 import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.structure.sequence import Sequence_generic, Sequence
from sage.symbolic.all import I
from sage.rings.arith import xgcd as orig_xgcd, lcm
from sage.rings.arith import gcd as orig_gcd
from sage.rings.number_field.number_field import QQ, ZZ
from sage.symbolic.ring import SymbolicRing
from sage.symbolic.expression import Expression


# for debugging
import better_exchook


# our own verbose function because I just want our msgs, not other stuff
def verbose(msg):
	print msg


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


def _simplify(a):
	if hasattr(a, "simplify_full"):
		return a.simplify_full()
	return simplify(a)

class CurlO:
	# We set `b = b1 + b2 (D + \sqrt{D})/2`.
	def __init__(self, D):
		self.D = D
		assert (D*D - D) % 4 == 0
	def divmod(self, a, b):
		"Returns q,r such that a = q*b + r."
		a1,a2 = self.as_tuple_b(a)
		b1,b2 = self.as_tuple_b(b)
		B = matrix([
			[b1, -b2 * (self.D**2 - self.D)/4],
			[b2, b1 + b2*self.D]
		])
		qq1,qq2 = B.solve_right(vector((a1,a2)))
		# We want r2 = 0, i.e. imag(r) = 0.
		# r1 = a1 - (q*b)_1 = a1 - q1*b1 + q2*b2 * (D*D-D)/4.
		# r2 = a2 - (q*b)_2 = a2 - q1*b2 - q2*b1 - D*q2*b2 = 0.
		# => q1 * b2 + q2 * (b1 + D*b2) = a2.
		q1,q2 = int(qq1), int(qq2) # not sure on this
		print qq1,qq2,q1,q2
		q = self.from_tuple_b(q1,q2)
		# q * b + r == a
		r = _simplify(a - q * b)
		return q,r
	def xgcd(self, a, b):
		a1,a2 = self.as_tuple_b(a)
		b1,b2 = self.as_tuple_b(b)
		if a2 == b2 == 0:
			return orig_xgcd(a1, a2)
		B2 = self.D + ssqrt(self.D)
		if a1 == b1 == 0:
			d,p,q = orig_xgcd(a2, b2)
			return d * B2, p, q
		if a1 == b2 == 0:
			d,p,q = orig_xgcd(a2, b1)
			return d * B2, p, q * B2
		if a2 == b1 == 0:
			d,p,q = orig_xgcd(a1, b2)
			return d * B2, p * B2, q

		# imag(gcd) = 0 =>
		# a2*p1 + a1*p2 + a2*D*p2 + b2*q1 + b1*q2 + b2*D*q2 == 0
		p1 = a1 + a2 * self.D
		p2 = -a2
		try:
			p1,p2 = int(ZZ(p1)), int(ZZ(p2))
		except Exception:
			print p1,p2
			raise
		pgcd = orig_gcd(p1,p2)
		p1,p2 = p1/pgcd, p2/pgcd
		p = self.from_tuple_b(p1,p2)

		q1 = b1 + b2 * self.D
		q2 = -b2
		q1,q2 = int(ZZ(q1)), int(ZZ(q2))
		qgcd = orig_gcd(q1,q2)
		q1,q2 = q1/qgcd, q2/qgcd
		q = self.from_tuple_b(q1,q2)

		d = a * p + b * q
		d = _simplify(d)
		return d, p, q

	def gcd(self, a, b):
		d,_,_ = self.xgcd(a, b)
		return d
	def common_denom(self, a, b):
		a1,a2 = self.as_tuple_b(a)
		b1,b2 = self.as_tuple_b(b)
		return matrix(1,4,[a1,a2,b1,b2]).denominator()
	def as_tuple_b(self, a):
		b2 = imag(a) * 2 / ssqrt(-self.D)
		b1 = real(a) - b2 * self.D / 2
		return (b1, b2)
	def from_tuple_b(self, b1, b2):
		return b1 + b2 * (self.D + ssqrt(self.D)) / 2
	def __contains__(self, item):
		try: self.as_tuple_b(item)
		except TypeError: return False
		else: return True

def test_curlO():
	space = CurlO(-3)
	assert space.xgcd(1337,43) == orig_xgcd(1337,43)

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
	Ring = I.base_ring() # Symbolic Ring

	A1 = matrix(Ring, 2,2, M[0][0])
	B1 = M[0][1] * S
	C1 = M[1][0] * S.inverse()
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
	d = space.gcd(A1[0][0] * l, C1[0][0] * l)
	if not d: d = 1
	Cg11 = C1[0][0] * l / d
	Dg11 = -A1[0][0] * l / d
	del l, d
	d,Ag11,Bg11 = space.xgcd(Dg11, -Cg11)
	assert d == 1, "{0}".format(tM1)
	Dg14 = Ag14 = 0
	Bg14 = 1
	Cg14 = -1
	G1 = make4x4matrix_embed(Ag11,Ag14,Bg11,Bg14,Cg11,Cg14,Dg11,Dg14)
	tM2 = G1 * tM1
	assert tM2[2][0] == 0
	assert tM2[3][0] == 0
	c22,c24 = tM2[2][1],tM2[3][1]
	l = space.common_denom(c22,c24)
	d = space.gcd(c22 * l, c24 * l)
	if not d: d = 1
	Dg23 = c24 * l / d
	Dg24 = -c22 * l / d
	del l, d
	d,Dg21,Dg22 = space.xgcd(Dg24, -Dg23)
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
		l = space.common_denom(a32,c32)
		d = space.gcd(a32 * l, c32 * l)
		if not d: d = 1
		Cg31 = c32 * l / d
		Dg31 = -a32 * l / d
		del l, d
		d,Ag31,Bg31 = space.xgcd(Dg31, -Cg31)
		assert d == 1
		G3 = make4x4matrix_embed(Ag31,Ag34,Bg31,Bg34,Cg31,Cg34,Dg31,Dg34)
	tM4 = G3 * tM3
	assert tM4[2][0] == 0
	assert tM4[2][1] == 0
	assert tM4[3][0] == 0
	assert tM4[3][1] == 0

	R = tM4.apply_map(_simplify)
	gamma = (G1.inverse() * G2.inverse() * G3.inverse()).apply_map(_simplify) # G1,G2,G3 are in \Sp_2(\curlO).
	assert tM == gamma * R
	assert gamma.conjugate_transpose() * J * gamma == J
	assert (R.submatrix(0,0,2,2) * R.submatrix(2,2,2,2).conjugate_transpose()).apply_map(_simplify) == 1
	return gamma, R, tM

def test_solveR():
	space = CurlO(-3)
	a,b,c,d = 2,1,1,1
	s,t,u = 5,ssqrt(-3),1
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
	s,t,u = 2, QQ(0.5) * ssqrt(-3) - QQ(0.5), 2
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,3,1
	s,t,u = 3, QQ(0.5) * ssqrt(-3) - QQ(1.5), 3
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	return gamma,R,tM

def _curlO_matrix_denom(mat, D):
	assert D < 0
	mat_real = mat.apply_map(real * 2)
	mat_imag = mat.apply_map(lambda x: simplify(imag(x) * 2 / ssqrt(-D)))
	denom_mat_real = ZZ(mat_real.denominator())
	denom_mat_imag = ZZ(mat_imag.denominator())
	assert denom_mat_real * mat_real in MatrixSpace(ZZ,mat.nrows(),mat.ncols())
	assert denom_mat_imag * mat_imag in MatrixSpace(ZZ,mat.nrows(),mat.ncols())
	denom = lcm(denom_mat_real, denom_mat_imag)
	return int(ZZ(denom))