
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
	
	
	CurlS_Generator curlS;
	
	
	
};

