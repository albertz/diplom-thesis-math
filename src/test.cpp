// Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
// Copyright (c) 2013, Albert Zeyer, www.az2000.de
// This code is under the GPL v3 or later, see License.txt in the root directory of this project.

// Compile:
// clang++ -std=c++11 test.cpp reduceGL.cpp algo.cpp

void test_reduceGL();
void test_algo_CurlSGen_ZZ();
void test_algo_CurlSGen_generic();
void test_algo_PrecisionF();
void test_algo_calcReducedCurlF();
void test_algo();

int main() {
//	test_reduceGL();
//	test_algo_CurlSGen_ZZ();
	test_algo_CurlSGen_generic();
//	test_algo_PrecisionF();
//	test_algo_calcReducedCurlF();
//	test_algo();
}
