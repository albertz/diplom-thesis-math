
#include "reduceGL.hpp"
#include "structs.hpp"


// In many cases, we wont use this variable and hardcode
// it for degree 2. However, we define it here,
// so we can point it out in some cases.
static const int HermDegree = 2;


template<typename T>
struct Matrix2 {
	T a,b,c,d; // [[a,b],[c,d]]
	Matrix2(T _a = 0, T _b = 0, T _c = 0, T _d = 0)
	: a(_a), b(_b), c(_c), d(_d) {}
	T det() { return a*c - b*d; }
};

struct CurlS_Generator {
	
	M2T getNextS();
};

struct ElemOfS {
	int calcPrecisionDimension() {}
	M2T matrix;
};

struct Odual {
	
};

struct PrecisionF {
	
};

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

struct ReductionMatrices_Calc {
	int HermWeight; // k in the paper. usually <20
	//Odual oDual; // not sure if i need to specify it explicitely..., probably not
	int D; // discriminant. usually D in {-2,-3,-4}
	
	//CurlO curlO;
	//Gamma gamma;
	//Character nu;
	
	CurlS_Generator curlS;
	PrecisionF curlF;
	
	int wantedDimension; // dim M_k^H
	int calcDimension() {
		
	}
	
	int calcOutDim(ElemOfS elemOfS) {
		// TODO...
		// this is \calcF(S) in the text
		return 10;
	}
	
	bool getA(ElemOfF aChar, Matrix2<Complex> m, Complex& res) {
		// aChar is characterizing our function a:\F->\Q like:
		//   a(aChar) = 1
		//   a(aChar[U]) = det(U)^-k for all U \in GL_2(\curlO)
		//   a(x) = 0 elsewhere.
		// Thus, we try to find some U with
		//   m = aChar[U] .
		
		// we are trying to solve this.
		// WolframAlpha input:
		// solve {{a1,a2},{a3,a4}} == conjugatetranspose({{u1,u2},{u3,u4}}) * {{b1,b2},{b3,b4}} * {{u1,u2},{u3,u4}}
		// gives us {{a1 -> b1 u1 Conjugate[u1] + b2 u3 Conjugate[u1], a2 -> b1 u2 Conjugate[u3] + b2 u4 Conjugate[u3], a3 -> b3 u1 Conjugate[u2] + b4 u3 Conjugate[u2], a4 -> b3 u2 Conjugate[u4] + b4 u4 Conjugate[u4]}}
		
		res = 0;
		
		auto a1 = m.a;
		auto a2 = m.b;
		auto a3 = m.c;
		auto a4 = m.d;
		auto b1 = aChar.a;
		auto b2 = aChar.b;
		auto b3 = aChar.c;
		auto b4 = aChar.d;
		
		GL2Iterator<CurlO> uGen;
		while(true) {
			auto u = uGen.get();
			uGen.next();
			auto u1 = u.a;
			auto u2 = u.b;
			auto u3 = u.c;
			auto u4 = u.d;
			
			if(a1 != b1 * u1 * u1.conjugate() + b2 * u3 * u1.conjugate())
				continue;
			if(a2 != b1 * u2 * u3.conjugate() + b2 * u4 * u3.conjugate())
				continue;
			if(a3 != b3 * u1 * u2.conjugate() + b4 * u3 * u2.conjugate())
				continue;
			if(a4 != b3 * u2 * u4.conjugate() + b4 * u4 * u4.conjugate())
				continue;
			
			// found it!
			res = pow(u.det(), -HermWeight);
			return true;
		}
		
		return false;
	}
		
	void calcReducedMatrix(const ElemOfF& T, ElemOfF& reducedMatrix, Matrix2<ElemOfO>& U) {
		// ...
		
	}
	
	std::vector<ElemOfF> curlFList;
	std::map<ElemOfF,size_t> reducedCurlFMap; // reducedMatrix(\cF) -> index in list
	std::vector<ElemOfF> reducedCurlFList; // reducedMatrix(\cF)
	void calcReducedCurlF() {
		curlFList.clear();
		reducedCurlFList.clear();
		for(ElemOfF T : curlF) {
			curlFList.push_back(T);
			ElemOfF F; // reduced matrix
			Matrix2<ElemOfO> U; // F[U] = T
			calcReducedMatrix(T, F, U);
			if(reducedCurlFMap.find(U) != reducedCurlFMap.end()) {
				reducedCurlFMap[U] = reducedCurlFList.size();
				reducedCurlFList.push_back(F);
			}
		}
	}
	
	// This defines a_F(T) for all F \in reducedMatrix(\cF), T \in \cF.
	// F is indexed by reducedCurlFList and T is indexed by curlFList.
	typedef ElemOfO ValueOfA;
	std::vector<ValueOfA> A;
	void calcA() {
		for(ElemOfF T : curlFList) {
			ElemOfF F; // reduced matrix
			Matrix2<ElemOfO> U; // F[U] = T
			calcReducedMatrix(T, F, U);
			// ...
		}
	}
	
	// a_F(T)
	ValueOfA evalA(ElemOfF aRepr, ElemOfF T) {
		struct hermitian_form_with_character_evaluation reduced;
		reduceGL(T, D, reduced);
		const int sign = 0; // 0 or 1
		const int nu_exp = 0; // 0 or 1
		if(aRepr == reduced.matrix)
			return
				Pow(reduced.character.transposition, sign) *
				Pow(reduced.character.nu, nu_exp) *
				Pow(reduced.character.determinant, -HermWeight);
		return 0;
	}
	
	// a_F[S](n)
	ValueOfA evalA_S_n(ElemOfF aRepr, ElemOfS S, int n) {
		// = \sum_{T \in \cF, tr(ST) = n} a_F(T)
		ValueOfA sum = 0;
		for(ElemOfF T : curlFList) {
			if(trace(S,T) == n) {
				sum += evalA(aRepr, T);
			}
		}
		return sum;
	}
	
	void calcOneColumn(ElemOfF elemOfF, vector<ValueOfA>& outVec) {
		int outVecIndex = 0;
		for(ElemOfS S : curlS) {
			for(int i = 0; i < S.calcPrecisionDimension(); ++i, ++outVecIndex) {
				auto& out = outVec[outVecIndex];
				out = evalA_S_n(elemOfF, S, i);
			}
		}
	}

	std::vector<ValueOfA> matrix; // flat. column \times row
	void calcMainMatrix() {
		calcA();
				
		size_t rowCount = 0;
		for(ElemOfS S : curlS) {
			rowCount += S.calcPrecisionDimension();
		}
		
		size_t column = 0;
		for(ElemOfF F : reducedCurlFList) {
			calcOneColumn( F, &matrix[rowCount * column], &matrix[rowCount * (column + 1)] );
			++column;
		}
	}
	
	void loop() {
		while(calcDimension() < wantedDimension) {
			M2T S = curlS.getNextS();
			
		}
	}
};



