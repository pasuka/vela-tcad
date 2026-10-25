// Minimal error-injection runtime for the adapter contract, not a numerical solver.
#include <cstdlib>
#ifdef _WIN32
#define EXPORT extern "C" __declspec(dllexport)
#else
#define EXPORT extern "C"
#endif
EXPORT int MKL_Set_Interface_Layer(int) { return 0; }
EXPORT int MKL_Set_Num_Threads_Local(int value) {
    static thread_local int threads = 1;
    const int previous = threads; threads = value; return previous;
}
EXPORT void pardisoinit(void**, const int*, int* iparm) {
    for (int i = 0; i < 64; ++i) iparm[i] = 0;
}
EXPORT void pardiso(void**, const int*, const int*, const int*, const int* phase,
    const int* n, const void* a, const int*, const int*, int*, const int*, int*,
    const int*, void* b, void* x, int* error) {
    const char* fail = std::getenv("VELA_TEST_MKL_FAIL_PHASE");
    *error = fail && std::atoi(fail) == *phase ? -4 : 0;
    if (*error || *phase != 33) return;
    // Tests use a diagonal system only for successful calls to this fixture.
    for (int i = 0; i < *n; ++i)
        static_cast<double*>(x)[i] = static_cast<double*>(b)[i] / static_cast<const double*>(a)[i];
}
