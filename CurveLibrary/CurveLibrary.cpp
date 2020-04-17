// CurveLibrary.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "CurveLibrary.h"


int __stdcall squareXL(int x, int &y) {
	#pragma EXPORT
	y += 100;
	return x * x;
}
