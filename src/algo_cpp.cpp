// Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
// Copyright (c) 2013, Albert Zeyer, www.az2000.de
// This code is under the GPL v3 or later, see License.txt in the root directory of this project.

#include "reduceGL.hpp"
#include "structs.hpp"
#include <map>
#include <vector>
#include <list>
#include <gmp.h>
#include <iostream>
#include <memory>
#include <string>

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
	virtual std::string getState() const = 0;
	virtual void setState(const std::string& state) = 0;
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
	std::string getState() const {
		return cur.getState();
	}
	void setState(const std::string& state) {
		DOMAIN_CHECK(state.size() == 4 * sizeof(Int));
		cur.setState(state);
	}
	bool isValid() const {
		if(cur.det(D) <= 0) return false;
		if(cur.a <= 0) return false;
		if(cur.c <= 0) return false;
		if(cur.gcd() > 1) return false;
		return true;
	}
	bool _hardLimitCheckB1() {
		auto &a = cur.a, &b1 = cur.b1, &c = cur.c;
		return a*c*(1-D) >= b1*b1;
	}
	bool _hardLimitCheckB2() {
		auto &a = cur.a, &b2 = cur.b2, &c = cur.c;
		return a*c*4 >= b2*b2*(-D);
	}
	void next() {
		auto &a = cur.a, &b1 = cur.b1, &b2 = cur.b2, &c = cur.c;
		if(b2 > 0) { b2 *= -1; return; }
		b2 = -b2 + 1;
		if(!_hardLimitCheckB2()) {
			b2 = 0;
			if(b1 > 0) { b1 *= -1; return; }
			b1 = -b1 + 1;
		}
		if(!_hardLimitCheckB1()) {
			b1 = b2 = 0;
			c ++;
		}
		if(c > a) {
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
	
	std::string getState() const {
		return cur.getState() + int_to_bin(curDenom);
	}
	void setState(const std::string& state) {
		DOMAIN_CHECK(state.size() == 5 * sizeof(Int));
		cur.setState(state.substr(0, 4 * sizeof(Int)));
		curDenom = bin_to_int<Int>(&state[sizeof(Int)*4]);
	}
	
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
		DOMAIN_CHECK(D < 0);
		_InfIterM2T_O* _iter = NULL;
		iterType = _iterType;
		if(iterType == "ZZ")
			_iter = new M2T_O_PosDefSortedZZ_Iterator(D);
		else if(iterType == "generic")
			_iter = new M2T_O_PosDefSortedGeneric_Iterator(D);
		else {
			std::cerr << "CurlS::init: iterType " << iterType << " is invalid" << std::endl;
			DOMAIN_CHECK(false);
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
	
	// state without list of current matrices (this->matrices).
	std::string getState() const {
		LOGIC_CHECK(D < 0);
		LOGIC_CHECK(iterType != "");
		LOGIC_CHECK(iter.get() != NULL);
		std::string res;
		res += int_to_bin<Int>(D);
		res += string_to_bin(iterType);
		res += string_to_bin(iter->getState());
		return res;
	}
	
	void setState(const std::string& state) {
		DOMAIN_CHECK(state.size() >= 4);
		const char* s = &state[0];
		const char* end = s + state.size();
		uint32_t slen = 0;
		int _D = bin_to_int<Int>(s); s += sizeof(Int);
		std::string _iterType = bin_to_string(s, end, slen); s += slen; slen = 0;
		std::string _iterState = bin_to_string(s, end, slen); s += slen; slen = 0;
		DOMAIN_CHECK(s == end);
		init(_D, _iterType);
		LOGIC_CHECK(iter.get() != NULL);
		iter->setState(_iterState);
	}
};


struct PrecisionF {
	int D; // fundamental discriminant
	Int B; // limit
	PrecisionF(int _D = 0, Int _B = 0) : D(_D), B(_B) {}
	
	struct Iter {
		const PrecisionF& F;
		M2T_Odual cur;
		Int curB;
		bool hitEnd;
		Iter(const PrecisionF& _F, bool _end = false) : F(_F), curB(0), hitEnd(_end) {}
		
		bool isValid() {
			if(cur.detD(F.D) < 0) return false;
			if(cur.a < 0 || cur.a >= F.B) return false;
			if(cur.c < 0 || cur.c >= F.B) return false;
			return true;
		}
		bool _hardLimitCheckB1() {
			auto &a = cur.a, &b1 = cur.b1, &c = cur.c;
			return a*c*(F.D*F.D-F.D) >= b1*b1;
		}
		bool _hardLimitCheckB2() {
			auto &a = cur.a, &b2 = cur.b2, &c = cur.c;
			return a*c*4 >= b2*b2;
		}
		void next() {
			auto &a = cur.a, &b1 = cur.b1, &b2 = cur.b2, &c = cur.c;
			if(b2 > 0) { b2 *= -1; return; }
			b2 = -b2 + 1;
			if(!_hardLimitCheckB2()) {
				b2 = 0;
				if(b1 > 0) { b1 *= -1; return; }
				b1 = -b1 + 1;
			}
			if(!_hardLimitCheckB1()) {
				b1 = b2 = 0;
				c ++;
			}
			if(c > curB) {
				b1 = b2 = 0;
				if(a == curB) {
					curB ++;
					a = 0;
					c = curB;
				}
				else {
					a ++;
					if(a < curB)
						c = curB;
					else
						c = 0;
				}
			}
			if(curB >= F.B)
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
	// See text about \curlF(S).
	// This is the precision limit of the restriction Elliptic modular forms.
	DOMAIN_CHECK(S.det(F.D) > 0);
	DOMAIN_CHECK(F.B > 0);
	return std::max(0, F.B * std::min(S.a, S.c) - S.absBMupper(F.D, 2 * F.B));
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
		for(auto _T = curlF.begin(); _T != curlF.end(); ++_T) {
			ElemOfF T = *_T;
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
		for(auto S = curlS.begin(); S != curlS.end(); ++S) {
			matrixRowCount += calcPrecisionDimension(curlF, (ElemOfS)*S);
		}
		
		LOGIC_CHECK(reducedCurlFList.size() > 0);
		LOGIC_CHECK(reducedCurlFList.size() == matrixColumnCount);
		
		matrix.clear();
		matrix.resize(matrixRowCount * matrixColumnCount);

		// reduce_GL is expensive, thus we iterate through curlF only once.
		for(auto _T = curlF.begin(); _T != curlF.end(); ++_T) {
			ElemOfF T = *_T;
			struct hermitian_form_with_character_evaluation reduced;
			reduce_GL(T, D, reduced);
			size_t column = reducedCurlFMap[reduced.matrix];
			size_t rowStart = 0;
			for(auto _S = curlS.begin(); _S != curlS.end(); ++_S) {
				ElemOfS S = *_S;
				int traceNum = trace(S,T, D);
				LOGIC_CHECK(traceNum >= 0);
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
	// The function just calculates the function without the det(\overline{tS})^k factor.
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
		for(auto S = curlS.begin(); S != curlS.end(); ++S) {
			// TODO: is calcPrecisionDimension correct here?
			matrixRowCountTrans += calcPrecisionDimension(curlF, (ElemOfS)*S);
		}
	
		matrixRowDenomTrans = lS * lS;
		
		matrixTrans.clear();
		matrixCountTrans = lS * lT;
		matrixTrans.resize(matrixRowCountTrans * matrixColumnCountTrans * matrixCountTrans);
		
		// reduce_GL is expensive, thus we iterate through curlF only once.
		for(auto _T = curlF.begin(); _T != curlF.end(); ++_T) {
			ElemOfF T = *_T;
			M2_Odual T_M2 = M2_Odual_from_M2T_Odual(T, D);
			struct hermitian_form_with_character_evaluation reduced;
			reduce_GL(T, D, reduced);
			size_t column = reducedCurlFMap[reduced.matrix];
			size_t rowStart = 0;
			for(auto _S = curlS.begin(); _S != curlS.end(); ++_S) {
				ElemOfS S = *_S;
				M2_O S_M2 = M2_O_from_M2T_O(S, D);
				// n = tr(T tS S tS^*)
				Int traceNum =
					T_M2
					.mulMat(tS, D)
					.mulMat(S_M2, D)
					.mulMat(tS.conjugate_transpose(D), D)
					.trace()
					.asInt(D);
				//cout << "traceNum=" << traceNum << ", T_M2=" << T_M2 << ", S_M2 = " << S_M2 << endl;
				size_t row = rowStart + traceNum;
				rowStart += calcPrecisionDimension(curlF, S);
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
				// Without the det(\overline{tS})^k factor, see the header comment.
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



void test_algo_CurlSGen(int D, int B_cF, const std::string& iterType, ssize_t denomLimit, ssize_t countLimit) {
	using namespace std;
	PrecisionF curlF;
	curlF.D = D;
	curlF.B = B_cF;
	CurlS_Generator curlS;
	curlS.init(D, iterType);
	size_t c = 0;
	while(true) {
		++curlS;
		auto S = **curlS.iter;
		Int curDenom = S.det(curlS.D);
		if(denomLimit >= 0 && curDenom > denomLimit) break;
		cout << "i=" << c << ", det=" << curDenom << ", S=" << S << ", prec=" << calcPrecisionDimension(curlF, S) << endl;
		++c;
		if(countLimit >= 0 && c > countLimit) break;
	}
	cout << "count: " << c << endl;
}

void test_algo_CurlSGen_ZZ() {
	test_algo_CurlSGen(-4, 7, "ZZ", 10, -1);
}

void test_algo_CurlSGen_generic() {
	test_algo_CurlSGen(-3, 7, "generic", -1, 1000);
}

void test_algo_PrecisionF() {
#ifndef OLDGCC
	using namespace std;
	PrecisionF curlF;
	curlF.D = -3;
	curlF.B = 6;
	size_t c = 0;
	for(ElemOfF T : curlF) {
		cout << T << endl;
		++c;
	}
	cout << "count: " << c << endl;
#endif
}

void test_algo_calcReducedCurlF() {
	using namespace std;
	ReductionMatrices_Calc calc;
	calc.init(-3, 6, "generic");
	calc.curlF.B = 10;
#ifndef OLDGCC
	size_t c = 0;
	for(ElemOfF T : calc.curlF) (void)T, ++c;
	cout << "size of curlF: " << c << endl;
#endif
	calc.calcReducedCurlF();
	cout << "size of reducedMatrix(curlF): " << calc.reducedCurlFList.size() << endl;	

	ReductionMatrices_Calc calc2;
	calc2.init(-3, 6, "generic");
	calc2.curlF.B = 20;
#ifndef OLDGCC
	for(ElemOfF T : calc2.curlF) (void)T, ++c;
	cout << "size of curlF2: " << c << endl;
#endif
	calc2.calcReducedCurlF();
	cout << "size of reducedMatrix(curlF2): " << calc2.reducedCurlFList.size() << endl;
	
	LOGIC_CHECK(calc.reducedCurlFList.size() <= calc2.reducedCurlFList.size());
	for(size_t i = 0; i < calc.reducedCurlFList.size(); ++i) {
		LOGIC_CHECK(calc.reducedCurlFList[i] == calc2.reducedCurlFList[i]);
	}

#ifndef OLDGCC
	// another test whether all reductions are reducing to themselves.
	for(auto Tr : calc.reducedCurlFMap) {
		struct hermitian_form_with_character_evaluation reduced;
		reduce_GL(Tr.first, calc.D, reduced);
		LOGIC_CHECK(reduced.matrix == Tr.first);
	}
	cout << "all reduced matrices reduce to themselves" << endl;

	// are all reduced matrices also in \curlF?
	for(auto Tr : calc.reducedCurlFMap) {
		bool found = false;
		for(auto T : calc.curlF) {
			if(T == Tr.first) {
				found = true;
				break;
			}
		}
		LOGIC_CHECK(found);
	}
	cout << "all reduced matrices are in curlF" << endl;
#endif
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
#ifndef OLDGCC
		{
			size_t c = 0;
			for(ElemOfF T : calc.curlF) { ++c; (void)T; }
			cout << "size of curlF: " << c << endl;
		}
#endif
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

