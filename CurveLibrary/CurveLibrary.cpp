// CurveLibrary.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "CurveLibrary.h"


// This is an example of an exported variable
CURVELIBRARY_API int nCurveLibrary=0;

// This is an example of an exported function.
CURVELIBRARY_API int fnCurveLibrary(void)
{
    return 0;
}

// This is the constructor of a class that has been exported.
CCurveLibrary::CCurveLibrary()
{
    return;
}
