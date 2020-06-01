// CurveLibrary.h : Declares all DLL functions
#pragma once
#define EXPORT comment(linker, "/EXPORT:" __FUNCTION__ "=" __FUNCDNAME__)

LONG __stdcall squareXL(int x, int& y);

BOOL __stdcall ir_measurement_multiXL(BOOL const save_curves);
BOOL __stdcall run_all_multiXL(BOOL const compute_curves, double& norm_test);
