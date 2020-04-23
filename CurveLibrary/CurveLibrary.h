// CurveLibrary.h : Declares all DLL functions
#define EXPORT comment(linker, "/EXPORT:" __FUNCTION__ "=" __FUNCDNAME__)

LONG __stdcall squareXL(int, int&);

BOOL __stdcall ir_measurement_multiXL(BOOL const);
BOOL __stdcall run_all_multiXL(BOOL const);

