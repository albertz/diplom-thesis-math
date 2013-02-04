
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
	
	void calcOneColumn(ElemOfF elemOfF, vector& out) {
		
	}
	
	
	void loop() {
		while(calcDimension() < wantedDimension) {
			M2T S = curlS.getNextS();
			
		}
	}
};



