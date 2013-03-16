
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
		matrices.push_back(M2T(1,1,0,1));
	}
};

int calcPrecisionDimension(ElemOfS S) {
	// TODO...
	return 10;
}

Int trace(M2T m1, M2T m2) {
	// calculates trace(m1 * m2)
	return m1.a * m2.a + 2 * m1.b1 * m2.b1 + 2 * m1.b2 * m2.b2 + m1.c * m2.c;
}

struct PrecisionF {
	Int B; // limit
	PrecisionF() : B(0) {}
	
	struct Iter {
		const PrecisionF& F;
		M2T cur;
		bool hitEnd;
		Iter(const PrecisionF& _F, bool _end = false) : F(_F), hitEnd(_end) {}
		
		bool isValid() {
			if(cur.det() < 0) return false;
			if(cur.a < 0 || cur.a >= F.B) return false;
			if(cur.c < 0 || cur.c >= F.B) return false;
			return true;
		}
		void next() {
			if(cur.b2 > 0) { cur.b2 *= -1; return; }
			cur.b2 = -cur.b2 + 1;
			if(cur.det() < 0) {
				cur.b2 = 0;
				if(cur.b1 > 0) { cur.b1 *= -1; return; }
				cur.b1 = -cur.b1 + 1;
			}
			if(cur.det() < 0) {
				cur.b1 = cur.b2 = 0;
				cur.c ++;
			}
			if(cur.c >= F.B) {
				cur.c = cur.b1 = cur.b2 = 0;
				cur.a ++;
			}
			if(cur.a >= F.B)
				hitEnd = true;
		}
		Iter& operator++() {
			do {
				next();
			} while(!isValid() && !hitEnd);
			return *this;
		}
		M2T operator*() const { return cur; }
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


typedef M2T ElemOfF;
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
			if(reducedCurlFMap.find(reduced.matrix) != reducedCurlFMap.end()) {
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
			for(int i = 0; i < calcPrecisionDimension(S); ++i, ++outVec) {
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
			rowCount += calcPrecisionDimension(S);
		}
		
		calcReducedCurlF();
		matrix.resize(rowCount * reducedCurlFList.size());
		
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
	curlF.B = 20;
	size_t c = 0;
	for(ElemOfF T : curlF) {
		cout << T << endl;
		++c;
	}
	cout << "count: " << c << endl;
}

void test_algo() {
	ReductionMatrices_Calc calc;
	calc.HermWeight = 10;
	calc.D = -2;
	calc.curlF.B = 10;
	calc.curlS.getNextS();
	calc.calcMainMatrix();
}

