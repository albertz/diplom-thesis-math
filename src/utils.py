from threading import currentThread
import time, os
from sage.matrix.matrix2 import Matrix
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.structure.sequence import Sequence_generic, Sequence
from sage.symbolic.ring import SymbolicRing
from sage.symbolic.expression import Expression


# for debugging
import better_exchook


# our own verbose function because I just want our msgs, not other stuff
def verbose(msg):
	print msg


MyDir = os.path.dirname(__file__) or os.getcwd()
MyDir = os.path.abspath(MyDir)


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
			if modname == "__main__": continue
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
	def _key_repr(self, key):
		from StringIO import StringIO
		key_sstr = StringIO()
		Pickler(key_sstr).dump(key)
		key_str = key_sstr.getvalue()
		import hashlib
		m = hashlib.md5()
		m.update(key_str)
		return m.hexdigest()
	def _filename_for_key(self, key):
		return MyDir + "/cache/" + self.name + "_" + self._key_repr(key) + ".sobj"
	def _filename_pattern(self):
		return MyDir + "/cache/" + self.name + "_*.sobj"
	def __getitem__(self, key):
		try:
			cache_key, value = load(self._filename_for_key(key))
		except IOError:
			raise KeyError, "cache file not found"
		except Exception as exc:
			better_exchook.better_exchook(*sys.exc_info())
			if isinstance(exc, KeyError): raise Exception # no KeyError exception fall-through
			raise
		else:
			if cache_key == key: return value
			raise KeyError, "key repr collidation"
	def __setitem__(self, key, value):
		try:
			save((key, value), self._filename_for_key(key))
		except Exception:
			print self.name, key
			raise
	def __contains__(self, key):
		try: self.__getitem__(key)
		except KeyError: return False
		else: return True
	def items(self):
		from glob import glob
		for fn in glob(self._filename_pattern()):
			try:
				key, value = load(fn)
			except Exception:
				print "Exception on loading cache file %s" % fn
				better_exchook.better_exchook(*sys.exc_info())
				raise
			yield key, value
	def keys(self):
		for key,value in self.items():
			yield key
	def __iter__(self):
		return self.keys()


def persistent_cache(name, index=None, timeLimit=2, ignoreNone=True):
	def decorator(func):
		from functools import wraps
		import algo_cython as C
		import inspect
		funcargs = tuple(inspect.getargspec(func).args)
		funcargdefaults = tuple(inspect.getargspec(func).defaults or [])
		assert len(funcargdefaults) <= len(funcargs)
		NotSpecifiedFlag = object()
		cache = PersistentCache(name=name)
		@wraps(func)
		def cached_function(*args, **kwargs):
			kwargs = kwargs.copy()
			if len(args) > len(funcargs): raise TypeError, "too many args"
			for i,arg in enumerate(args):
				key = funcargs[i]
				if key in kwargs: raise TypeError, "%r specified in both args and kwargs" % arg
				kwargs[key] = arg
			args = list(funcargdefaults) + [NotSpecifiedFlag] * (len(funcargs) - len(funcargdefaults))
			for key,value in kwargs.items():
				if key not in funcargs:
					raise TypeError, "kwarg %r is unknown" % key
				i = funcargs.index(key)
				args[i] = value
			for key,value in zip(funcargs,args):
				if value is NotSpecifiedFlag:
					raise TypeError, "arg %r is not specified" % key
			cacheidx = ()
			if index is not None:
				cacheidx = index(*args)
			else:
				for arg in args:
					if isinstance(arg, C.Calc):
						cacheidx += (arg.params, arg.curlS,)
					else:
						cacheidx += (arg,)
			if cacheidx in cache:
				return cache[cacheidx]
			t = time.time()
			res = func(*args)
			if res is None and ignoreNone:
				return None
			if timeLimit is None or time.time() - t > timeLimit:
				print "calculation of %r took %f secs" % (func, time.time() - t)
				cache[cacheidx] = res
			return res
		return cached_function
	return decorator


def convert_old_cache(name):
	old_cache_fn = MyDir + "/" + name + ".cache.sobj"
	assert os.path.exists(old_cache_fn)
	d = load(old_cache_fn)
	assert isinstance(d, dict)
	new_cache = PersistentCache(name=name)
	for key,value in d.items():
		new_cache[key] = value



# The following code is partly from [MusicPlayer](https://github.com/albertz/music-player/)
# but all written by me, thus I transfer this part to GPL here.

def attrChain(base, *attribs, **kwargs):
	default = kwargs.get("default", None)
	obj = base
	for attr in attribs:
		if obj is None: return default
		obj = getattr(obj, attr, None)
	if obj is None: return default
	return obj

# This is needed in some cases to avoid pickling problems with bounded funcs.
def funcCall(attrChainArgs, args=()):
	f = attrChain(*attrChainArgs)
	return f(*args)

class ExecingProcess:
	def __init__(self, target, args, name):
		self.target = target
		self.args = args
		self.name = name
		self.daemon = True
		self.pid = None
	def start(self):
		assert self.pid is None
		def pipeOpen():
			readend,writeend = os.pipe()
			readend = os.fdopen(readend, "r")
			writeend = os.fdopen(writeend, "w")
			return readend,writeend
		self.pipe_c2p = pipeOpen()
		self.pipe_p2c = pipeOpen()
		pid = os.fork()
		if pid == 0: # child
			self.pipe_c2p[0].close()
			self.pipe_p2c[1].close()
			if "/local/bin/" in sys.argv[0]: #'/Applications/sage-5.9/local/bin/sage-ipython'
				SageBin = os.path.normpath(os.path.dirname(sys.argv[0]) + "/../../sage")
				assert os.path.exists(SageBin), "%r" % SageBin
			else:
				assert False, "add code for: %r" % sys.argv
			args = [
				SageBin,
				"run.py",
				"--forkExecProc",
				str(self.pipe_c2p[1].fileno()),
				str(self.pipe_p2c[0].fileno())]
			print "parent pid: %r, pid: %r, args: %r" % (os.getppid(), os.getpid(), args)
			os.execv(args[0], args)
		else: # parent
			self.pipe_c2p[1].close()
			self.pipe_p2c[0].close()
			self.pid = pid
			self.writeend = self.pipe_p2c[1]
			self.readend = self.pipe_c2p[0]
			self.pickler = Pickler(self.writeend)
			self.pickler.dump(self.name)
			self.pickler.dump(self.target)
			self.pickler.dump(self.args)
			self.writeend.flush()
			self.unpickler = Unpickler(self.readend)

	@staticmethod
	def checkExec():
		if "--forkExecProc" in sys.argv:
			argidx = sys.argv.index("--forkExecProc")
			writeFileNo = int(sys.argv[argidx + 1])
			readFileNo = int(sys.argv[argidx + 2])
			readend = os.fdopen(readFileNo, "r")
			writeend = os.fdopen(writeFileNo, "w")
			unpickler = Unpickler(readend)
			name = unpickler.load()
			print "ExecingProcess child %s (pid %i)" % (name, os.getpid())
			try:
				target = unpickler.load()
				args = unpickler.load()
			except EOFError:
				print "Error: unpickle incomplete"
				raise SystemExit
			pickler = Pickler(writeend)
			target(*(args + (readend, unpickler, writeend, pickler)))
			print "ExecingProcess child %s (pid %i) finished" % (name, os.getpid())
			raise SystemExit

class ExecingProcess_ConnectionWrapper(object):
	def __init__(self, readend, unpickler, writeend, pickler):
		self.readend = readend
		self.unpickler = unpickler
		self.writeend = writeend
		self.pickler = pickler
	def send(self, value):
		self.pickler.dump(value)
		self.writeend.flush()
	def recv(self):
		return self.unpickler.load()
	def poll(self, timeout=None):
		import select
		rlist, wlist, xlist = select.select([self.readend.fileno()], [], [], timeout)
		return bool(rlist)
	def close(self):
		self.readend.close()
		self.writeend.close()

class AsyncTask:
	def __init__(self, func, name=None):
		self.name = name or repr(func)
		self.func = func
		self.parent_pid = os.getpid()
		self.proc = ExecingProcess(
			target = self._asyncCall,
			args = (self,),
			name = self.name + " worker process")
		self.proc.daemon = True
		self.proc.start()
		self.child_pid = self.proc.pid
		assert self.child_pid
		self.conn = ExecingProcess_ConnectionWrapper(
			readend=self.proc.readend,
			unpickler=self.proc.unpickler,
			writeend=self.proc.writeend,
			pickler=self.proc.pickler)

	@staticmethod
	def _asyncCall(self, readend, unpickler, writeend, pickler):
		assert self.isChild
		self.conn = ExecingProcess_ConnectionWrapper(
			readend=readend,
			unpickler=unpickler,
			writeend=writeend,
			pickler=pickler)
		try:
			self.func(self)
		except KeyboardInterrupt:
			print "Exception in AsyncTask", self.name, ": KeyboardInterrupt"
		except:
			print "Exception in AsyncTask", self.name
			sys.excepthook(*sys.exc_info())
		finally:
			self.conn.close()

	def put(self, value):
		self.conn.send(value)

	def get(self):
		thread = currentThread()
		try:
			thread.waitQueue = self
			res = self.conn.recv()
		except EOFError: # this happens when the child died
			raise ForwardedKeyboardInterrupt()
		except Exception:
			raise
		finally:
			thread.waitQueue = None
		return res

	def poll(self, **kwargs):
		return self.conn.poll(**kwargs)

	@property
	def isParent(self):
		return self.parent_pid == os.getpid()

	@property
	def isChild(self):
		if self.isParent: return False
		# No check. The Sage wrapper binary itself forks again, so this is wrong.
		#assert self.parent_pid == os.getppid(), "%i, %i, %i" % (self.parent_pid, os.getppid(), os.getpid())
		return True

	# This might be called from the module code.
	# See OnRequestQueue which implements the same interface.
	def setCancel(self):
		self.conn.close()
		if self.isParent and self.child_pid:
			import signal
			os.kill(self.child_pid, signal.SIGINT)
			self.child_pid = None

class ForwardedKeyboardInterrupt(Exception): pass

def asyncCall(func, name=None):
	def doCall(queue):
		res = None
		try:
			res = func()
			queue.put((None,res))
		except KeyboardInterrupt as exc:
			print "Exception in asyncCall", name, ": KeyboardInterrupt"
			queue.put((ForwardedKeyboardInterrupt(exc),None))
		except BaseException as exc:
			print "Exception in asyncCall", name
			sys.excepthook(*sys.exc_info())
			queue.put((exc,None))
	task = AsyncTask(func=doCall, name=name)
	# If there is an unhandled exception in doCall or the process got killed/segfaulted or so,
	# this will raise an EOFError here.
	# However, normally, we should catch all exceptions and just reraise them here.
	exc,res = task.get()
	if exc is not None:
		raise exc
	return res

# END of the part of MusicPlayer code.


class Parallelization_Worker:
	def __init__(self, _id):
		self.id = _id
		self.jobidcounter = 0
		self.task = AsyncTask(func=self._work, name="Parallel worker")
	def _handle_job(self, queue, jobid, func, name):
		try:
			res = func()
			queue.put((self.id, jobid, func, None, res))
		except KeyboardInterrupt as exc:
			print "Exception in asyncCall", name, ": KeyboardInterrupt"
			queue.put((self.id, jobid, func, ForwardedKeyboardInterrupt(exc), None))
		except BaseException as exc:
			print "Exception in asyncCall", name
			sys.excepthook(*sys.exc_info())
			queue.put((self.id, jobid, func, exc, None))
	def _work(self, queue):
		while True:
			jobid, func, name = queue.get()
			self._handle_job(queue=queue, jobid=jobid, func=func, name=name)
	def put_job(self, func, name):
		self.jobidcounter += 1
		jobid = self.jobidcounter
		self.task.put((jobid, func, name))
	def is_ready(self):
		return self.jobidcounter == 0
	def get_result(self, block=False, timeout=None):
		if not block or timeout is not None:
			poll_kwargs = {}
			if timeout is not None: poll_kwargs["timeout"] = timeout
			if not self.task.poll(**poll_kwargs):
				from Queue import Empty
				raise Empty
		selfid, jobid, func, exc, res = self.task.get()
		assert selfid == self.id
		if jobid == self.jobidcounter:
			self.jobidcounter = 0
		return func, exc, res


class Parallelization:
	def __init__(self, task_limit=4):
		import multiprocessing
		self.task_count = 0
		self.task_limit = task_limit
		self.task_iter = None
		self.task_queue = multiprocessing.Queue()
		self.workers = [Parallelization_Worker(i) for i in range(self.task_limit)]

	def get_next_result(self):
		from Queue import Empty
		while True:
			self.maybe_queue_tasks()
			for w in self.workers:
				if w.is_ready(): continue
				try: restuple = w.get_result(timeout=0.1)
				except Empty: pass
				else:
					self.task_count -= 1
					self.maybe_queue_tasks()
					return restuple

	def get_all_ready_results(self):
		results = []
		from Queue import Empty
		for w in self.workers:
			if w.is_ready(): continue
			try: restuple = w.get_result(timeout=0)
			except Empty: pass
			results += [restuple]
		return results

	def maybe_queue_tasks(self):
		while self.task_count < self.task_limit:
			next_task = next(self.task_iter)
			self.exec_task(func=next_task)

	def exec_task(self, func, name=None):
		jobidcounter, w = min([(w.jobidcounter, w) for w in self.workers])
		if name is None: name=repr(func)
		self.task_count += 1
		w.put_job(func=func, name=name)

def reloadC():
	"""
	This is just for testing to reload the C (Cython) module
	after it was recompiled.
	Note that this code is very unsafe! It will likely crash
	when you have other references on the C module.
	This is only for debugging and development!
	"""
	import algo
	C = algo.C
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


