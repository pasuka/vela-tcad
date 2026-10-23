// Loadable negative fixture: deliberately exports no MKL symbols.
#ifdef _WIN32
extern "C" __declspec(dllexport) int vela_mkl_negative_fixture() { return 0; }
#else
extern "C" int vela_mkl_negative_fixture() { return 0; }
#endif
