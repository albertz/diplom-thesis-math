
#include "reduceGL.hpp"
#include "structs.hpp"
#include <map>
#include <vector>
#include <list>
#include <gmp.h>


// In many cases, we wont use this variable and hardcode
// it for degree 2. However, we define it here,
// so we can point it out in some cases.
static const int HermDegree = 2;


typedef M2T_O ElemOfS;

struct CurlS_Generator {
	std::list<ElemOfS> matrices;
	std::list<ElemOfS>::iterator begin() { return matrices.begin(); }
	std::list<ElemOfS>::iterator end() { return matrices.end(); }

	/*
	Just now, we iterare reduced matrices in Her_2(\Z).
	This might not be correct, we might need all reduced
	matrices in Her_2(\Z) here.
	It is not possible to iterate them by det(S)!
	This could be a problem. Maybe it is enough to limit Im(S).
	*/

	ElemOfS cur; // matrix S
	Int curDenom; // denom of S^-1 = det S = ac - b^2

	CurlS_Generator() : curDenom(1) {}

	bool isValid() {
		if(cur.a > cur.c) return false;
		if(cur.det() != curDenom) return false;
		if(gcd(cur.a, cur.b, cur.c) > 1) return false;
		return true;
	}
	void next() {
		auto &a = cur.a, &b = cur.b, &c = cur.c;

		c ++;
		if(c <= curDenom + b*b) {
			a = (curDenom + b*b) / c;
			return;
		}

		c = 0;
		if(b > 0) { b *= -1; return; }
		b = -b + 1;
		
		// It must hold: det >= 3b^2.
		// Thus, when b becomes too huge, we have all matrices with that det.
		if(3*b*b > curDenom) {
			a = b = c = 0;
			curDenom ++;
		}
	}
	CurlS_Generator& operator++() {
		do {
			next();
		} while(!isValid());
		return *this;
	}

	M2T getNextS() {
		Int oldDenom = curDenom;
		++(*this); // the very first ([0,0,0]) is not valid
		if(curDenom != oldDenom)
			matrices.clear();
		matrices.push_back(cur);
		return cur;
	}
};

// calculates trace(S * T)
// is always an integer
Int trace(M2T_O S, M2T_Odual T) {
	// = (S.a * T.a + S.b * \overline{T.b}) + (\overline{S.b} * T.b + S.b * T.b)
	// = S.a * T.a + 2 * Re(S.b * \overline{T.b}) + S.c * T.c
	// = S.a * T.a + S.c * T.c + (S.b1 * T.b2 - S.b2 * T.b1)
	return S.a * T.a + S.c * T.c + S.b1 * T.b2 - S.b2 * T.b1;
}

struct PrecisionF {
	int D; // fundamental discriminant
	Int B; // limit
	PrecisionF(int _D = 0, Int _B = 0) : D(_D), B(_B) {}
	
	struct Iter {
		const PrecisionF& F;
		M2T_Odual cur;
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
		M2T_Odual operator*() const { return cur; }
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
	DOMAIN_CHECK(S.a + S.c - 2 * abs(S.b) > 0); // this is always the case if S is positive definite
	DOMAIN_CHECK(F.B > 0);
	return F.B * (S.a + S.c - 2 * abs(S.b));
}

typedef M2T_Odual ElemOfF;
typedef Int ValueOfA;

struct ReductionMatrices_Calc {
	int HermWeight; // k in the paper. usually <20
	int D; // discriminant. usually D in {-2,-3,-4}
	
	ReductionMatrices_Calc() {
		HermWeight = 0;
		D = 0;
		matrixRowCount = matrixColumnCount = 0;
	}
		
	void init(int _D, int _HermWeight) {
		D = curlF.D = _D;
		
		DOMAIN_CHECK(D < 0);
		// Fundamental discriminant properties on D:
		DOMAIN_CHECK(
			   (Mod(D, 4) == 1) /* TODO: and check square-free */
			   || (Mod(D, 4) == 0 && (Mod(Div(D, 4), 4) == 2 || Mod(Div(D, 4), 4) == 3)));
		DOMAIN_CHECK(Mod(D*D - D, 4) == 0); // is implied by the above, but check anyway...

		HermWeight = _HermWeight;
		
		DOMAIN_CHECK(HermWeight > 0);
		// By Dern, Satz 1.10, to imply that we don't have the trivial cases.
		if(D == -4)
			DOMAIN_CHECK(Mod(HermWeight, 2) == 0);
		if(D == -3) {
			// nu = det^j
			// then: Mod(HermWeight, 3) == j
			// TODO...
		}
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
		if(aRepr == reduced.matrix)
			return reduced.character.value(D, -HermWeight);
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
		LOGIC_CHECK(outVec == outVecEnd);
	}

	std::vector<ValueOfA> matrix; // flat. format: [[0]*ColumnCount]*RowCount
	size_t matrixRowCount, matrixColumnCount;
	void calcMatrix() {
		matrixRowCount = 0;
		for(ElemOfS S : curlS) {
			matrixRowCount += calcPrecisionDimension(curlF, S);
		}
		
		calcReducedCurlF();
		matrixColumnCount = reducedCurlFList.size();
		
		matrix.clear();
		matrix.resize(matrixRowCount * matrixColumnCount);

		// reduce_GL is expensive, thus we iterate through curlF only once.
		for(ElemOfF T : curlF) {
			struct hermitian_form_with_character_evaluation reduced;
			reduce_GL(T, D, reduced);
			size_t column = reducedCurlFMap[reduced.matrix];
			for(ElemOfS S : curlS) {
				int traceNum = trace(S,T);
				if(traceNum >= calcPrecisionDimension(curlF, S)) continue;
				size_t row = traceNum;
				size_t matrixIndex = row * matrixColumnCount + column;
				matrix[matrixIndex] += reduced.character.value(D, -HermWeight);
			}
		}
	}
	
	void getMatrix(mpz_t* out) {
		for(size_t i = 0; i < matrixColumnCount * matrixRowCount; ++i) {
			mpz_set_si((mpz_ptr)out, matrix[i]);
			++out;
		}
	}
	
	void dumpMatrix() {
		using namespace std;
		LOGIC_CHECK(matrix.size() == matrixRowCount * matrixColumnCount);
		cout << "matrix " << matrixRowCount << "*" << matrixColumnCount << endl;
		cout << "[" << endl;
		for(size_t row = 0; row < matrixRowCount; ++row) {
			if(row > 0) cout << "," << endl;
			cout << " [";
			for(size_t column = 0; column < matrixColumnCount; ++column) {
				if(column > 0) cout << ",";
				cout << matrix[row * matrixColumnCount + column];
			}
			cout << "]";
		}
		cout << "]" << endl;
	}
};



void test_algo_CurlSGen() {
	using namespace std;
	CurlS_Generator curlS;
	size_t c = 0;
	size_t denomLimit = 10;
	while(true) {
		++curlS;
		if(curlS.curDenom > denomLimit) break;
		cout << curlS.curDenom << ", " << curlS.cur << endl;
		++c;
	}
	cout << "count: " << c << endl;
}

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
	calc.init(-4, 10);
	calc.curlF.B = 10;
	calc.calcReducedCurlF();
	cout << "size of reducedMatrix(curlF): " << calc.reducedCurlFList.size() << endl;	
}

void test_algo() {
	using namespace std;
	ReductionMatrices_Calc calc;
	calc.init(-4, 10);
	calc.curlF.B = 20;
	calc.curlS.getNextS();
	calc.curlS.getNextS();
	
	{
		Timer timer("calcMatrix");
		calc.calcMatrix();
	}
	{
		size_t c = 0;
		for(ElemOfF T : calc.curlF) { ++c; (void)T; }
		cout << "size of curlF: " << c << endl;
	}
	cout << "size of reducedMatrix(curlF): " << calc.reducedCurlFList.size() << endl;
	cout << "size of matrix: " << calc.matrix.size() << endl;
	//calc.dumpMatrix();
}

