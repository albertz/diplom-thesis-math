# -*- coding: utf-8 -*-
# Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
# Copyright (c) 2013, Albert Zeyer, www.az2000.de
# This code is under the GPL v3 or later, see License.txt in the root directory of this project.

from threading import currentThread
import time, os
from sage.matrix.matrix2 import Matrix
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.structure.sequence import Sequence_generic, Sequence
from sage.symbolic.ring import SymbolicRing
from sage.symbolic.expression import Expression


# For debugging. It's a drop-in replacement for sys.excepthook.
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
	from StringIO import StringIO
	s = StringIO()
	pickler = Pickler(s)
	pickler.dump(obj)
	try:
		f = open(filename, "w")
		f.write(s.getvalue())
	except BaseException:
		# try again
		f = open(filename, "w")
		f.write(s.getvalue())
		# but reraise
		raise

def load(filename):
	f = open(filename)
	unpickler = Unpickler(f)
	obj = unpickler.load()
	return obj


CacheEnabled = True
if not int(os.environ.get("CACHE", 1)):
	print "Cache disabled"
	CacheEnabled = False


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
		if not CacheEnabled: raise KeyError
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
		if not CacheEnabled: return
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
		if not CacheEnabled: return
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
		if not CacheEnabled: return
		for key,value in self.items():
			yield key
	def __iter__(self):
		return self.keys()

def pickle_dumps(obj):
	from StringIO import StringIO
	s = StringIO()
	pickler = Pickler(s)
	pickler.dump(obj)
	return s.getvalue()

def pickle_loads(s):
	from StringIO import StringIO
	ss = StringIO(s)
	unpickler = Unpickler(ss)
	return unpickler.load()


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


# This was needed to convert some cache files which were created with an earlier version
# of this code. This probably can be deleted at some later point. (2013-06-12)
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
		self.isChild = False

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
			self.isChild = True
			self.pipe_c2p[0].close()
			self.pipe_p2c[1].close()
			if "/local/bin/" in sys.executable: #'/Applications/sage-5.9/local/bin/python'
				SageBin = os.path.normpath(os.path.dirname(sys.executable) + "/../../sage")
				assert os.path.exists(SageBin), "%r" % ([SageBin, sys.executable, sys.argv, sys.exec_prefix])
			else:
				assert False, "%r" % ([sys.executable, sys.argv, sys.exec_prefix])
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

	def isFinished(self):
		assert not self.isChild
		pid, status = os.waitpid(self.pid, os.WNOHANG)
		if pid != 0: assert pid == self.pid
		return pid != 0

	def __del__(self):
		if not self.isChild and self.pid:
			import signal
			os.kill(self.pid, signal.SIGINT)
			os.waitpid(self.pid, 0)

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
		# We exec() instead of just fork() because of problems with non-fork-safe libs in Sage.
		# See: http://ask.sagemath.org/question/2627/sigill-in-forked-process
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
		except Exception:
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
	"""
	See documentation of the `Parallelization` class.
	"""
	def __init__(self, _id):
		self.id = _id
		self.joblist = []
		self._init_task()
	def _init_task(self):
		self.task = AsyncTask(func=self._work, name="Parallel worker")
	def _handle_job(self, queue, func, name):
		try:
			try:
				res = func()
			except KeyboardInterrupt as exc:
				print "Exception in asyncCall", name, ": KeyboardInterrupt"
				queue.put((self.id, func, ForwardedKeyboardInterrupt(exc), None))
			except BaseException as exc:
				print "Exception in asyncCall", name, "<pid %i>" % os.getpid(), ": %r" % exc
				sys.excepthook(*sys.exc_info())
				queue.put((self.id, func, exc, None))
			else:
				queue.put((self.id, func, None, res))
		except IOError as exc:
			# This is probably a broken pipe or so.
			print "parallel worker <pid %i> IOError: %r" % (os.getpid(), exc)
			raise SystemExit
	def _work(self, queue):
		while True:
			func, name = queue.get()
			self._handle_job(queue=queue, func=func, name=name)
	def put_job(self, func, name):
		job = (func, name)
		self.joblist += [job]
		self.task.put(job)
	def is_ready(self):
		return len(self.joblist) == 0
	def get_result(self, block=False, timeout=None):
		if not block or timeout is not None:
			poll_kwargs = {}
			if timeout is not None: poll_kwargs["timeout"] = timeout
			if not self.task.poll(**poll_kwargs):
				from Queue import Empty
				raise Empty
		selfid, func, exc, res = self.task.get()
		assert selfid == self.id
		assert len(self.joblist) > 0
		self.joblist.pop(0)
		return func, exc, res
	def fixup_broken_proc(self):
		# We expect that the proc has terminated. If it hasn't,
		# this function shouldn't be called and this is a bug.
		assert self.task.proc.isFinished()
		self._init_task() # reinit new proc
		# push all incompleted jobs
		for job in self.joblist:
			self.task.put(job)


class Parallelization:
	"""
	This class is the base class to parallize some calculation. It creates multiple
	worker processes via the `Parallelization_Worker` class. These are own independent
	processes and we do the communication via serialization (via pickling) on pipes.

	You queue any number of tasks which will be executed on any worker. You can get
	any available results via `get_next_result()` which is blocking or `get_all_ready_results()`.

	Queueing the tasks works by setting the `task_iter` attribute to any iterator.
	That iterator must yield tasks, which are callable objects.

	The managing/parent process can call `maybe_queue_tasks()` to start the calculation
	of new tasks, if we are idling. This is also called automatically in `get_next_result()`.

	You can also call `exec_task()` manually to queue a task.
	"""

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
				try:
					restuple = self._stable_proc_communicate(
						w, lambda: w.get_result(timeout=0.1))
				except Empty: pass
				else:
					self.task_count -= 1
					self.maybe_queue_tasks()
					return restuple

	def _stable_proc_communicate(self, worker, action):
		while True:
			try:
				return action()
			except IOError as e: # broken pipe or so
				# Expect that the proc has been terminated.
				# (If this is wrong, this need some more checking code what happened
				#  so that we can recover accordingly.)
				print "%r proc exception: %r" % (action, e)
				print "Reinit proc ..."
				worker.fixup_broken_proc()
				print "And retry."

	def get_all_ready_results(self):
		results = []
		from Queue import Empty
		for w in self.workers:
			while True:
				# Try to get all results and break if there aren't anymore.
				if w.is_ready(): break
				try:
					restuple = self._stable_proc_communicate(
						w, lambda: w.get_result(timeout=0))
				except Empty:
					break
				else:
					self.task_count -= 1
					results += [restuple]

		return results

	def maybe_queue_tasks(self):
		count = 0
		while self.task_count < self.task_limit:
			next_task = next(self.task_iter)
			self.exec_task(func=next_task)
			count += 1
		return count

	def exec_task(self, func, name=None):
		_, w = min([(len(w.joblist), w) for w in self.workers])
		if name is None: name=repr(func)
		self.task_count += 1
		self._stable_proc_communicate(
			w, lambda: w.put_job(func=func, name=name))

	def get_pending_tasks(self):
		joblist = []
		for w in self.workers:
			joblist += w.joblist # list of (func,name)
		return joblist


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


