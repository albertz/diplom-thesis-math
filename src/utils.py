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
from sage.rings.number_field.number_field import QQ, ZZ, QuadraticField
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
		self.field = QuadraticField(D)
		self.Droot = self.field(D).sqrt()
	def divmod(self, a, b):
		"""
		Returns q,r such that a = q*b + r.
		"""
		# Note that this implementation is quite naive!
		# Later, we can do better with QuadraticForm(...). (TODO)

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

		q1,q2 = int(round(qq1)), int(round(qq2)) # not sure on this
		#print a1,a2,b1,b2,qq1,qq2,q1,q2
		q = self.from_tuple_b(q1,q2)
		# q * b + r == a
		r = _simplify(a - q * b)
		# Note that this works for -D < 5.27. See the text.
		assert _simplify(abs(r)) < _simplify(abs(b))
		return q,r
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
			B2 = (self.D + ssqrt(self.D)) / 2
			return d * B2, s, t

		abs_a = _simplify(abs(a))
		abs_b = _simplify(abs(b))
		#assert abs_a != abs_b, "%r" % ((a,b,abs_a,abs_b),)
		if abs_a < abs_b:
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
		return matrix(1,len(tupleargs),tupleargs).denominator()
	def as_tuple_b(self, a):
		b2 = _simplify(imag(a) * 2 / ssqrt(-self.D))
		b1 = _simplify(real(a) - b2 * self.D / 2)
		return (b1, b2)
	def from_tuple_b(self, b1, b2):
		return b1 + b2 * (self.D + self.Droot) / 2
	def __contains__(self, item):
		try:
			b1,b2 = self.as_tuple_b(item)
			b1,b2 = ZZ(b1), ZZ(b2)
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
	#Ring = I.base_ring() # Symbolic Ring
	Ring = space.field

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

	return gamma,R,tM

def _curlO_matrix_denom(mat, D):
	space = CurlO(D)
	denom = space.common_denom(*mat.list())
	denom = int(ZZ(denom))
	for v in mat.list():
		assert v * denom in space, "%r (D=%r)" % (mat, D)
	return denom


# Hack for reload handling
def reimportMeIntoAlgoModule():
	import sys
	if "algo" in sys.modules:
		mod = sys.modules["algo"]
		for attr in globals().keys():
			if attr.startswith("__"): continue
			if hasattr(mod, attr):
				setattr(mod, attr, globals()[attr])
reimportMeIntoAlgoModule()
