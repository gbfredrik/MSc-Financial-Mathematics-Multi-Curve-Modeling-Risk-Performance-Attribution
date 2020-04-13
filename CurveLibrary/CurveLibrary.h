// The following ifdef block is the standard way of creating macros which make exporting
// from a DLL simpler. All files within this DLL are compiled with the CURVELIBRARY_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see
// CURVELIBRARY_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef CURVELIBRARY_EXPORTS
#define CURVELIBRARY_API __declspec(dllexport)
#else
#define CURVELIBRARY_API __declspec(dllimport)
#endif

// This class is exported from the dll
class CURVELIBRARY_API CCurveLibrary {
public:
	CCurveLibrary(void);
	// TODO: add your methods here.
};

extern CURVELIBRARY_API int nCurveLibrary;

CURVELIBRARY_API int fnCurveLibrary(void);
