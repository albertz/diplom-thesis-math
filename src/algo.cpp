
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
	void getNextS();
};

int calcPrecisionDimension(ElemOfS S) {
	// TODO...
	return 10;
}

Int trace(M2T m1, M2T m2) {
	// calculates trace(m1 * m2)
	// TODO..
	return 0;
}

struct Odual {
	
};

struct PrecisionF {
	struct Iter {
		
		M2T operator*() const { return M2T(); }
		Iter& operator++() { return *this; }
		bool operator!=(const Iter&) const { return true; }
	};
	Iter begin() { return Iter(); }
	Iter end() { return Iter(); }
};


typedef M2T ElemOfF;
typedef Int ValueOfA;

struct ReductionMatrices_Calc {
	int HermWeight; // k in the paper. usually <20
	//Odual oDual; // not sure if i need to specify it explicitely..., probably not
	int D; // discriminant. usually D in {-2,-3,-4}
	
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
		
		size_t column = 0;
		for(ElemOfF F : reducedCurlFList) {
			calcOneColumn( F, &matrix[rowCount * column], &matrix[rowCount * (column + 1)] );
			++column;
		}
	}
	
};



void test_algo() {
	
}

