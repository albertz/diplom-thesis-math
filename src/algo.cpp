
typedef int Int;
//typedef ... M2; // Matrix 2x2

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
	CurlS_Generator
};

