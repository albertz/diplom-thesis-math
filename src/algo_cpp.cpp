
#include "reduceGL.hpp"
#include "structs.hpp"
#include <map>
#include <vector>
#include <list>
#include <gmp.h>
#include <iostream>

// In many cases, we wont use this variable and hardcode
// it for degree 2. However, we define it here,
// so we can point it out in some cases.
static const int HermDegree = 2;


typedef M2T_O ElemOfS;

template<typename T>
struct InfiniteIterIntf {
	typedef InfiniteIterIntf Self;
	virtual ~InfiniteIterIntf() {}
	virtual bool isValid() const = 0;
	virtual void next() = 0;
	virtual T get() const = 0;
	Self& operator++() {
		do {
			next();
		} while(!isValid());
		return *this;
	}
	virtual T operator*() const { return get(); }
	virtual bool operator==(const Self& other) const {
		return (**this) == (*other);
	}
	bool operator!=(const Self& other) const { return !(*this == other); }
};

struct _InfIterM2T_O : InfiniteIterIntf<M2T_O> {
	int D;
	_InfIterM2T_O(int _D) : D(_D) {}
};

struct M2T_O_PosDefSortedGeneric_Iterator : _InfIterM2T_O {
	M2T_O cur;
	M2T_O_PosDefSortedGeneric_Iterator(int _D) : _InfIterM2T_O(_D) {}
	M2T_O get() const { return cur; }
	bool isValid() const {
		if(cur.det(D) < 0) return false;
		if(cur.a < 0) return false;
		if(cur.c < 0) return false;
		return true;
	}
	bool _hardLimitCheck() {
		auto &a = cur.a, &b1 = cur.b1, &b2 = cur.b2, &c = cur.c;
		// det >= 0 <=> a*c - b1*b1 - D*b2 - Div(D*D-D, 4) - b2*b2 >= 0.
		// Thus, when we have ac >= b1*b1 - (-D)*|b2| + Div(D*D-D, 4) + b2*b2,
		// we are always safe that we don't miss any values. Of course,
		// we must still check for validity because we will get invalid values.
		return a*c >= b1*b1 - (-D) * abs(b2) + Div(D*D-D, 4) + b2*b2;
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
		if(c >= a) {
			c = b1 = b2 = 0;
			a ++;
		}
	}
};

struct M2T_O_PosDefSortedZZ_Iterator : _InfIterM2T_O {
	ElemOfS cur; // matrix S
	Int curDenom; // denom of S^-1 = det S = ac - b^2
	
	M2T_O_PosDefSortedZZ_Iterator(int _D) : _InfIterM2T_O(_D), curDenom(1) {}

	/*
	 Just now, we iterare reduced matrices in Her_2(\Z).
	 This might not be correct, we might need all reduced
	 matrices in Her_2(\Z) here.
	 It is not possible to iterate them by det(S)!
	 This could be a problem. Maybe it is enough to limit Im(S).
	 */
	
	M2T_O get() const { return cur; }
	bool isValid() const {
		if(cur.a > cur.c) return false;
		if(cur.det(D) != curDenom) return false;
		if(cur.gcd() > 1) return false;
		return true;
	}
	void next() {
		auto &a = cur.a, &b = cur.b1, &c = cur.c;
		
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
};

struct CurlS_Generator {
	int D;
	std::list<ElemOfS> matrices;
	std::list<ElemOfS>::iterator begin() { return matrices.begin(); }
	std::list<ElemOfS>::iterator end() { return matrices.end(); }
	size_t size() { return matrices.size(); }
	
	std::string iterType;
	std::auto_ptr<_InfIterM2T_O> iter;
	
	CurlS_Generator() : D(0) {}
	void init(int _D, const std::string& _iterType) {
		D = _D;
		_InfIterM2T_O* _iter = NULL;
		iterType = _iterType;
		if(iterType == "ZZ")
			_iter = new M2T_O_PosDefSortedZZ_Iterator(D);
		else if(iterType == "generic")
			_iter = new M2T_O_PosDefSortedGeneric_Iterator(D);
		else {
			std::cerr << "CurlS::init: iterType " << iterType << " is invalid" << std::endl;
			LOGIC_CHECK(false);
		}
		iter = std::auto_ptr<_InfIterM2T_O>(_iter);
	}
	
	CurlS_Generator& operator++() {
		do {
			iter->next();
		} while(!iter->isValid());
		return *this;
	}

	M2T_O getNextS() {
		Int oldDenom = ElemOfS(**iter).det(D);
		++(*this); // the very first ([0,0,0]) is not valid
		ElemOfS cur = ElemOfS(**iter);
		Int curDenom = cur.det(D);
		if(curDenom != oldDenom)
			matrices.clear();
		matrices.push_back(cur);
		return cur;
	}
	
	void clearMatrices() {
		matrices.clear();
	}
};


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
	DOMAIN_CHECK(S.det(F.D) > 0);
	DOMAIN_CHECK(F.B > 0);
	// Note: B(s + u) is \in \Z. B(s+u) >= 2 B |t|.
	// 2 B |t| might be in \R-\Z, but we still have
	// Gauss_upper(2 B |t|) <= B(s+u).
	// Those can be equal though, thus this can be zero.
	// With increasing B, this is always increasing, though.
	return (F.B * S.a + F.B * S.c - S.absBMupper(F.D, 2 * F.B));
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
		
	void init(int _D, int _HermWeight, const std::string& curlSiterType) {
		D = curlF.D = _D;
		curlS.init(D, curlSiterType);
		
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

		DOMAIN_CHECK(Mod(HermWeight, 3) == 0); // the modulform is trivial/zero if HermWeight is not divisible by 3
	}
	
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
		matrixColumnCount = reducedCurlFList.size();
	}
		
	// this calcs the matrix for the map f \mapsto f[S]
	std::vector<ValueOfA> matrix; // flat. format: [[0]*ColumnCount]*RowCount
	size_t matrixRowCount, matrixColumnCount;
	void calcMatrix() {
		matrixRowCount = 0;
		LOGIC_CHECK(curlS.size() > 0);
		for(ElemOfS S : curlS) {
			matrixRowCount += calcPrecisionDimension(curlF, S);
		}
		
		LOGIC_CHECK(reducedCurlFList.size() > 0);
		LOGIC_CHECK(reducedCurlFList.size() == matrixColumnCount);
		
		matrix.clear();
		matrix.resize(matrixRowCount * matrixColumnCount);

		// reduce_GL is expensive, thus we iterate through curlF only once.
		for(ElemOfF T : curlF) {
			struct hermitian_form_with_character_evaluation reduced;
			reduce_GL(T, D, reduced);
			size_t column = reducedCurlFMap[reduced.matrix];
			size_t rowStart = 0;
			for(ElemOfS S : curlS) {
				int traceNum = trace(S,T);
				size_t row = rowStart + traceNum;
				rowStart += calcPrecisionDimension(curlF, S);
				if(row >= rowStart) continue;
				size_t matrixIndex = row * matrixColumnCount + column;
				auto a_T = reduced.character.value(D, -HermWeight);
				matrix[matrixIndex] += a_T;
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
	
	// This calcs the matrix for the map f \mapsto (f|R)[S] up to a factor,
	// where R = [[tS,tT;0,tS^{-1}^{*}]] and tS,tT \in Mat_2(\K).
	// We have
	//   (a|R_i)[S](n) =
	//     det((tS)^{-1}^{*})^{-k}
	//     * \sum_{T, tr(T * tS * S * tS^{*}) = n} a(T) * e^{2 \pi i tr(T * tT * tS^{*})}
	// The function just calculates the function without the det(..) factor.
	// The tS,tT parameters given to this function might be multiplied with some factors lS,lT
	// so that the matrices are in \curlO.
	// Note that the right part are elements of the CyclomoticField(lS*lT).
	// The FE represent quotients of the form q^{n/l} where l = matrixRowDenomTrans = lS * lS.
	std::vector<ValueOfA> matrixTrans; // flat. format: [[[0]*ColumnCount]*RowCount]*ZetaOrder
	size_t matrixRowCountTrans, matrixColumnCountTrans, matrixCountTrans;
	size_t matrixRowDenomTrans;
	void calcMatrixTrans(const M2_O& tS, const M2_O& tT, const Int lS, const Int lT) {
		using namespace std;
		cout << "calcMatrixTrans: tS = " << tS << ", tT = " << tT << ", tS^* = " << tS.conjugate_transpose(D) << ", lS = " << lS << ", lT = " << lT << endl;
		
		matrixColumnCountTrans = matrixColumnCount;
		LOGIC_CHECK(matrixColumnCountTrans > 0);
		
		matrixRowCountTrans = 0;
		LOGIC_CHECK(curlS.size() > 0);
		for(ElemOfS S : curlS) {
			// TODO: is calcPrecisionDimension correct here?
			matrixRowCountTrans += calcPrecisionDimension(curlF, S);
		}
	
		// TODO: WARNING: All the stuff about matrixRowDenomTrans and trace (below) has to be studied and worked out.
		matrixRowDenomTrans = lS * lS;
		matrixRowCountTrans *= matrixRowDenomTrans;
		
		matrixTrans.clear();
		matrixCountTrans = lS * lT;
		matrixTrans.resize(matrixRowCountTrans * matrixColumnCountTrans * matrixCountTrans);
		
		// reduce_GL is expensive, thus we iterate through curlF only once.
		for(ElemOfF T : curlF) {
			M2_Odual T_M2 = M2_Odual_from_M2T_Odual(T, D);
			struct hermitian_form_with_character_evaluation reduced;
			reduce_GL(T, D, reduced);
			size_t column = reducedCurlFMap[reduced.matrix];
			size_t rowStart = 0;
			for(ElemOfS S : curlS) {
				M2_O S_M2 = M2_O_from_M2T_O(S, D);
				// n = tr(T tS S tS^*)
				/*
				auto t1 = T_M2.mulMat(tS, D);
				auto t2 = t1.mulMat(S_M2, D);
				auto _t3 = tS.conjugate_transpose(D);
				auto t3 = t2.mulMat(_t3, D);
				auto t4 = t3.trace();
				*/
				Int traceNum =
					T_M2
					.mulMat(tS, D)
					.mulMat(S_M2, D)
					.mulMat(tS.conjugate_transpose(D), D)
					.trace()
					.asInt(D);
				//cout << "traceNum=" << traceNum << ", T_M2=" << T_M2 << ", S_M2 = " << S_M2 << endl;
				size_t row = rowStart + traceNum;
				rowStart += calcPrecisionDimension(curlF, S) * matrixRowDenomTrans;
				//cout << "traceNum: " << traceNum << ", next rowStart: " << rowStart << endl;
				if(row >= rowStart) continue;
				size_t matrixIndex = row * matrixColumnCountTrans + column;
				auto a_T = reduced.character.value(D, -HermWeight);
				// factor = tr(T tT tS^*)
				auto factor_exp =
					T_M2
					.mulMat(tT, D)
					.mulMat(tS.conjugate_transpose(D), D)
					.trace()
					.asInt(D);
				factor_exp = Mod(factor_exp, lS*lT);
				matrixIndex += factor_exp * matrixRowCountTrans * matrixColumnCountTrans;
				matrixTrans[matrixIndex] += a_T;
			}
		}
	}
	
	void getMatrixTrans(mpz_t* out, int matrixIndex) {
		LOGIC_CHECK(matrixIndex >= 0 && matrixIndex < matrixCountTrans);
		for(size_t i = 0; i < matrixColumnCountTrans * matrixRowCountTrans; ++i) {
			mpz_set_si((mpz_ptr)out, matrixTrans[matrixIndex * matrixRowCountTrans * matrixColumnCountTrans + i]);
			++out;
		}
	}
	
	void dumpMatrixTrans(int matrixIndex) {
		using namespace std;
		LOGIC_CHECK(matrixIndex >= 0 && matrixIndex < matrixCountTrans);
		LOGIC_CHECK(matrixTrans.size() == matrixColumnCountTrans * matrixRowCountTrans * matrixCountTrans);
		cout << "matrix " << matrixRowCountTrans << "*" << matrixColumnCountTrans << endl;
		cout << "[" << endl;
		for(size_t row = 0; row < matrixRowCountTrans; ++row) {
			if(row > 0) cout << "," << endl;
			cout << " [";
			for(size_t column = 0; column < matrixColumnCountTrans; ++column) {
				if(column > 0) cout << ",";
				cout << matrixTrans[matrixIndex * matrixRowCountTrans * matrixColumnCountTrans + row * matrixColumnCountTrans + column];
			}
			cout << "]";
		}
		cout << "]" << endl;
	}

};



void test_algo_CurlSGen() {
	using namespace std;
	CurlS_Generator curlS;
	curlS.D = -4;
	size_t c = 0;
	size_t denomLimit = 10;
	while(true) {
		++curlS;
		Int curDenom = (**curlS.iter).det(curlS.D);
		if(curDenom > denomLimit) break;
		cout << curDenom << ", " << (**curlS.iter) << endl;
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
	calc.init(-4, 10, "ZZ");
	calc.curlF.B = 10;
	calc.calcReducedCurlF();
	cout << "size of reducedMatrix(curlF): " << calc.reducedCurlFList.size() << endl;	
}

void test_algo() {
	using namespace std;
	
	{
		ReductionMatrices_Calc calc;
		calc.init(-4, 18, "ZZ");
		calc.curlF.B = 20;
		calc.calcReducedCurlF();
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

	{

		ReductionMatrices_Calc calc;
		const int D = -3;
		calc.init(D, 6, "ZZ");
		calc.curlF.B = 10;
		calc.calcReducedCurlF();
		calc.curlS.matrices.push_back(M2T_O(1,0,0,4));

		{
			Timer timer("calcMatrixTrans");
			calc.calcMatrixTrans(M2_O_from_M2T_O(M2T_O(4,0,0,4), D), M2_O_from_M2T_O(M2T_O(8,0,0,2), D), 4, 4);
		}

		{
			Timer timer("calcMatrixTrans");
			calc.calcMatrixTrans(M2_O_from_M2T_O(M2T_O(-2,0,0,1), D), M2_O_from_M2T_O(M2T_O(0,0,0,0), D), 2, 2);
		}
	}
}

