
typedef int Int;
//typedef ... M2; // Matrix 2x2

// In many cases, we wont use this variable and hardcode
// it for degree 2. However, we define it here,
// so we can point it out in some cases.
static const int HermDegree = 2;

struct M2T {
	// [[a,b],[\bar b, c]]
	Int a, b, c;
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


struct ReductionMatrices_Calc {
	int HermWeight; // k in the paper. usually <20
	//Odual oDual; // not sure if i need to specify it explicitely..., probably not
	
	//CurlO curlO;
	//Gamma gamma;
	//Character nu;
	
	CurlS_Generator curlS;
	PrecisionF curlF;
	
	int wantedDimension; // dim M_k^H
	int calcDimension() {
		
	}
	
	void calcOneColumn(ElemOfF elemOfF, vector& outVec) {
		int outVecIndex = 0;
		for(ElemOfS S : curlS) {
			for(int i = 0; i < S.calcPrecisionDimension(); ++i, ++outVecIndex) {
				auto& out = outVec[outVecIndex];
				
				// we have "a" given by the base which is identified by elemOfF.
				
				// calculate a[S]
				
				
				// {{a1,a2},{a3,a4}} = conjugatetranspose({{u1,u2},{u3,u4}}) * {{b1,b2},{b3,b4}} * {{u1,u2},{u3,u4}}
				// -> {{a1, a2}, {a3, a4}} == {{(b1 u1 + b2 u3) Conjugate[u1], (b1 u2 + b2 u4) Conjugate[u3]}, {(b3 u1 + b4 u3) Conjugate[u2], (b3 u2 + b4 u4) Conjugate[u4]}}
			}
		}
	}
	
	
	void loop() {
		while(calcDimension() < wantedDimension) {
			M2T S = curlS.getNextS();
			
		}
	}
};



