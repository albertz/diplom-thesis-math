
#include "reduceGL.hpp"
#include "structs.hpp"
#include <map>
#include <vector>
#include <list>


// In many cases, we wont use this variable and hardcode
// it for degree 2. However, we define it here,
// so we can point it out in some cases.
static const int HermDegree = 2;


typedef M2T ElemOfS;

struct CurlS_Generator {
	std::list<ElemOfS> matrices;
	std::list<ElemOfS>::iterator begin() { return matrices.begin(); }
	std::list<ElemOfS>::iterator end() { return matrices.end(); }
	void getNextS() {
		// TODO... (or in Python?)
		matrices.push_back(M2T(2,1,1));
	}
};

// calculates trace(S * T)
// is always an integer
Int trace(M2T S, M2T_O T) {
	// = S.a * T.a + 2 * Re(S.b * \overline{T.b}) + S.c * T.c
	// = S.a * T.a + S.c * T.c + 2 * S.b * Re(T.b)
	// = S.a * T.a + S.c * T.c + S.b * T.b2
	return S.a * T.a + S.c * T.c + S.b * T.b2;
}

struct PrecisionF {
	int D; // fundamental discriminant
	Int B; // limit
	PrecisionF(int _D = 0, Int _B = 0) : D(_D), B(_B) {}
	
	struct Iter {
		const PrecisionF& F;
		M2T_O cur;
		bool hitEnd;
		Iter(const PrecisionF& _F, bool _end = false) : F(_F), hitEnd(_end) {}
		
		bool isValid() {
			if(cur.det4D(F.D) < 0) return false;
			if(cur.a < 0 || cur.a >= F.B) return false;
			if(cur.c < 0 || cur.c >= F.B) return false;
			return true;
		}
		bool _hardLimitCheck() {
			auto &a = cur.a, &b1 = cur.b1, &b2 = cur.b2, &c = cur.c;
			// det4D >= 0 <=> 4(-D)ac >= 4b1^2 - 4(-D)b1b2 + (-D)(1-D)b2^2
			// Thus, when we have 4(-D)ac >= 4b1^2 - 4(-D)|b1b2| + (-D)(1-D)b2^2,
			// we are always safe that we don't miss any values. Of course,
			// we must still check for validity because we will get invalid values.
			// 4b1^2 - 4(-D)|b1b2| = 4|b1| (|b1| - (-D)|b2|).
			// This is always strongly monotonly increasing and sign-independent -> thus we can iterate.
			return 4*(-F.D)*a*c >= 4*b1*b1 - 4*(-F.D)*abs(b1*b2) + (-F.D)*(1-F.D)*b2*b2;
		}
		void next() {
			auto &a = cur.a, &b1 = cur.b1, &b2 = cur.b2, &c = cur.c;
			if(b2 > 0) { b2 *= -1; return; }
			b2 = -b2 + 1;
			if(!_hardLimitCheck()) {
				b2 = 0;
				if(b1 > 0) { b1 *= -1; return; }
				b1 = -b1 + 1;
			}
			if(!_hardLimitCheck()) {
				b1 = b2 = 0;
				c ++;
			}
			if(c >= F.B) {
				c = b1 = b2 = 0;
				a ++;
			}
			if(a >= F.B)
				hitEnd = true;
		}
		Iter& operator++() {
			do {
				next();
			} while(!isValid() && !hitEnd);
			return *this;
		}
		M2T_O operator*() const { return cur; }
		bool operator==(const Iter& other) const {
			if(hitEnd && other.hitEnd) return true;
			if(!hitEnd && other.hitEnd) return false;
			if(hitEnd && !other.hitEnd) return false;
			return cur == other.cur;
		}
		bool operator!=(const Iter& other) const { return !(*this == other); }
	};
	Iter begin() { return Iter(*this); }
	Iter end() { return Iter(*this, true); }
};


int calcPrecisionDimension(const PrecisionF& F, ElemOfS S) {
	// tr([s,t,u] * [a,b,c]) >= max(a,c) * (s + u - 2|t|)
	// T=[a,b,c] \in \cF - \Lambda => max(a,c) >= B
	assert(S.a + S.c - 2 * abs(S.b) > 0); // this is always the case if S is positive definite
	assert(F.B > 0);
	return F.B * (S.a + S.c - 2 * abs(S.b));
}

typedef M2T_O ElemOfF;
typedef Int ValueOfA;

struct ReductionMatrices_Calc {
	int HermWeight; // k in the paper. usually <20
	int D; // discriminant. usually D in {-2,-3,-4}
	
	ReductionMatrices_Calc() {
		HermWeight = 0;
		D = 0;
	}
	
	//CurlO curlO;
	//Gamma gamma;
	//Character nu;
	
	CurlS_Generator curlS;
	PrecisionF curlF;
	
	std::map<ElemOfF,size_t> reducedCurlFMap; // reducedMatrix(\cF) -> index in list
	std::vector<ElemOfF> reducedCurlFList; // reducedMatrix(\cF)
	void calcReducedCurlF() {
		reducedCurlFMap.clear();
		reducedCurlFList.clear();
		for(ElemOfF T : curlF) {
			struct hermitian_form_with_character_evaluation reduced;
			reduce_GL(T, D, reduced);
			if(reducedCurlFMap.find(reduced.matrix) == reducedCurlFMap.end()) {
				reducedCurlFMap[reduced.matrix] = reducedCurlFList.size();
				reducedCurlFList.push_back(reduced.matrix);
			}
		}
	}
		
	// a_F(T)
	ValueOfA evalA(ElemOfF aRepr, ElemOfF T) {
		struct hermitian_form_with_character_evaluation reduced;
		reduce_GL(T, D, reduced);
		const int sign = 0; // 0 or 1
		const int nu_exp = 0; // 0 or 1
		if(aRepr == reduced.matrix) {			
			ValueOfA result = Pow(reduced.character.detValue(D), -HermWeight);
			if(sign) result *= reduced.character.transposition;
			if(nu_exp) result *= reduced.character.nu;
			return result;
		}
		return 0;
	}
	
	// a_F[S](n)
	ValueOfA evalA_S_n(ElemOfF aRepr, ElemOfS S, int n) {
		// = \sum_{T \in \cF, tr(ST) = n} a_F(T)
		ValueOfA sum = 0;
		for(ElemOfF T : curlF) {
			if(trace(S,T) == n) {
				sum += evalA(aRepr, T);
			}
		}
		return sum;
	}
	
	template<typename TIter>
	void calcOneColumn(ElemOfF elemOfF, TIter outVec, TIter outVecEnd) {
		for(ElemOfS S : curlS) {
			for(int i = 0; i < calcPrecisionDimension(curlF, S); ++i, ++outVec) {
				auto& out = *outVec;
				out = evalA_S_n(elemOfF, S, i);
			}
		}
		assert(outVec == outVecEnd);
	}

	std::vector<ValueOfA> matrix; // flat. column \times row
	void calcMainMatrix() {				
		size_t rowCount = 0;
		for(ElemOfS S : curlS) {
			rowCount += calcPrecisionDimension(curlF, S);
		}
		
		calcReducedCurlF();
		matrix.resize(rowCount * reducedCurlFList.size());
		
		// TODO: reorder loops for performance.
		// reduce_GL is expensive, thus we should iterate
		// through curlF only once.
		
		size_t column = 0;
		for(ElemOfF F : reducedCurlFList) {
			calcOneColumn( F, &matrix[rowCount * column], &matrix[rowCount * (column + 1)] );
			++column;
		}
	}
	
};



void test_algo_PrecisionF() {
	using namespace std;
	PrecisionF curlF;
	curlF.D = -4;
	curlF.B = 20;
	size_t c = 0;
	for(ElemOfF T : curlF) {
		cout << T << endl;
		++c;
	}
	cout << "count: " << c << endl;
}

void test_algo_calcReducedCurlF() {
	using namespace std;
	ReductionMatrices_Calc calc;
	calc.HermWeight = 10;
	calc.D = calc.curlF.D = -4;
	calc.curlF.B = 10;
	calc.calcReducedCurlF();
	cout << "size of reducedMatrix(curlF): " << calc.reducedCurlFList.size() << endl;	
}

void test_algo() {
	using namespace std;
	ReductionMatrices_Calc calc;
	calc.HermWeight = 10;
	calc.D = calc.curlF.D = -4;
	calc.curlF.B = 10;
	calc.curlS.getNextS();
	calc.calcMainMatrix();
	cout << "size of reducedMatrix(curlF): " << calc.reducedCurlFList.size() << endl;
	cout << "size of matrix: " << calc.matrix.size() << endl;
}

