
#include "reduceGL.hpp"
#include "structs.hpp"


// In many cases, we wont use this variable and hardcode
// it for degree 2. However, we define it here,
// so we can point it out in some cases.
static const int HermDegree = 2;


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
	
	std::map<ElemOfF,size_t> reducedCurlFMap; // reducedMatrix(\cF) -> index in list
	std::vector<ElemOfF> reducedCurlFList; // reducedMatrix(\cF)
	void calcReducedCurlF() {
		reducedCurlFMap.clear();
		reducedCurlFList.clear();
		for(ElemOfF T : curlF) {
			curlFList.push_back(T);
			struct hermitian_form_with_character_evaluation reduced;
			reduceGL(T, D, reduced);
			if(reducedCurlFMap.find(reduced.matrix) != reducedCurlFMap.end()) {
				reducedCurlFMap[reduced.matrix] = reducedCurlFList.size();
				reducedCurlFList.push_back(reduced.matrix);
			}
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
		for(ElemOfF T : curlF) {
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



