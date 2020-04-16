// CurveLibrary.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "CurveLibrary.h"



int WINAPI testFN()
{
#pragma EXPORT
	return 10;
}

// ?
// https://stackoverflow.com/questions/538134/exporting-functions-from-a-dll-with-dllexport?rq=1