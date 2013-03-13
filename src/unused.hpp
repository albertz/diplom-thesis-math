template<typename T>
struct BaseComplex {
	T real;
	T imag;
	BaseComplex(T _real = 0, T _imag = 0) : real(_real), imag(_imag) {}
	BaseComplex conjugate() const { return BaseComplex(real, -imag); }
};

typedef BaseComplex<Int> Complex;


template<typename T>
struct LinearGenerator {
	T elem;
	LinearGenerator(T _start = 0) : elem(_start) {}
	void next() { ++elem; }
	void prev() { --elem; }
	T get() { return elem; }
};

template<typename T, typename TGen>
struct PlusMinusGenerator {
	TGen gen;
	bool plus;
	PlusMinusGenerator(TGen _gen = 0) : gen(_gen) {}
	void next() {
		if(plus && gen.get() != -gen.get()) {
			plus = false;
		} else {
			plus = true;
			gen.next();
		}
	}
	void prev() {
		if(!plus) {
			plus = true;
		} else {
			plus = false;
			gen.prev();
			if(gen.get() == -gen.get()) plus = true;
		}
	}
	T get() { if(plus) return gen.get(); else -gen.get(); }
};

template<typename T1, typename T2>
struct Pair {
	T1 a; T2 b;
	Pair(T1 _a = 0, T2 _b = 0) : a(_a), b(_b) {}
};

template<typename T1, typename T1Gen, typename T2, typename T2Gen>
struct CombinedGenerator {
	T1Gen gen1, gen1_base, gen1_outer;
	T2Gen gen2, gen2_base, gen2_outer;
	
	CombinedGenerator(T1Gen _gen1 = 0, T2Gen _gen2 = 0)
	: gen1_base(_gen1), gen2_base(_gen2) {
		gen1 = gen1_outer = gen1_base;
		gen2 = gen2_outer = gen2_base;
	}
	
	void next() {
		if(gen1.get() == gen1_base.get()) {
			gen1_outer.next();
			gen2_outer.next();
			gen1 = gen1_outer;
			gen2 = gen2_base;
		}
		else {
			gen1.prev();
			gen2.next();
		}
	}
	
	void prev() {
		if(gen1.get() == gen1_outer.get()) {
			gen1_outer.prev();
			gen2_outer.prev();
			gen1 = gen1_base;
			gen2 = gen2_outer;
		}
		else {
			gen1.next();
			gen2.prev();
		}
	}
	
	typedef typename Pair<T1,T2> Type;
	Type get() {
		return Pair<T1,T2>(gen1.get(), gen2.get());
	}
};

template<typename T>
struct Matrix2Generator {
	typedef typename LinearGenerator<T> LinGen;
	typedef typename PlusMinusGenerator<T,LinGen> BaseGen;
	typedef typename CombinedGenerator<T,BaseGen,T,BaseGen> Vector2Gen;
	typedef typename Vector2Gen::Type Vector2Type;
	typedef typename CombinedGenerator<Vector2Type,Vector2Gen,Vector2Type,Vector2Gen> Vector4Gen;
	Vector4Gen vec4gen;
	
	void next() { vec4gen.next(); }
	void prev() { vec4gen.prev(); }
	
	typedef typename Matrix2<T> Matrix2Type;
	typedef Matrix2Type Type;
	Type get() {
		Matrix2Type m(
					  vec4gen.get().a.get().a,
					  vec4gen.get().a.get().b,
					  vec4gen.get().b.get().a,
					  vec4gen.get().b.get().b
					  );
		return m;
	}
};

template<typename T>
struct GL2Iterator {
	Matrix2Generator<T> mgen;
	
	GL2Iterator() { while(!isValid()) mgen.next(); }
	bool isValid() { return get().det() > 0; }
	void next() { do { mgen.next(); } while(!isValid()); }
	void prev() { do { mgen.prev(); } while(!isValid()); }
	
	typedef typename Matrix2Generator<T>::Type Type;
	Type get() { return mgen.get(); }
};
