# dotSimplex
A Network Simplex implementation for Discrete Optimal Transport.

These files requires a recent Microsoft Visual Studio compiler (MS VS 2019), and must be compiled in release mode (turn on all the optimization flags), and to enable the "OpenMP" support.


The files are:

- TestSolver.cpp : basic example to show how to solve instances on R^2 and on R^300 (such as word mover distances
- DOT_Vars.h : basic support for handling variables and measuring memory consumption
- DOT_Simplex :  our incremental version of the Network Simplex algorithm that support constrant addition and removal
- DOT_R2_Solver.h : the column generation algorithm for problems defined on R^2
- DOT_Rk_Solver.h : the column generation algorithm for problems defined on R^300 (but it is valid for any k)

Data files:
- The random instances on R^2 are created with given seed at runtime
- The Classical Images of the DOTmark benchmark are available online, as noted in the paper.
- For generating the random instances on R^300 we used the Glove word embedding.


All the MS compiler options are:

/GS /Qpar /GL /W3 /Gy /Zc:wchar_t /Zi /Gm- /O2 /sdl /Fd"x64\Release\vc142.pdb" /Zc:inline /fp:fast /D "_MBCS" /errorReport:prompt /WX- /Zc:forScope 
/arch:AVX2 /Gd /Oi /MD /openmp /std:c++latest /FC /Fa"x64\Release\" /EHsc /nologo /Fo"x64\Release\" /Ot /Fp"x64\Release\dotSimplex.pch" /diagnostics:classic


All the linker options are:

/OUT:"dotSimplex\msvc\dotSimplex\x64\Release\dotSimplex.exe" 
/MANIFEST /LTCG:incremental /NXCOMPAT /PDB:"dotSimplex\msvc\dotSimplex\x64\Release\dotSimplex.pdb" 
/DYNAMICBASE "kernel32.lib" "user32.lib" "gdi32.lib" "winspool.lib" "comdlg32.lib" "advapi32.lib" "shell32.lib" "ole32.lib" "oleaut32.lib" "uuid.lib" "odbc32.lib" "odbccp32.lib" 
/DEBUG:FULL /MACHINE:X64 /OPT:REF /PGD:"dotSimplex\msvc\dotSimplex\x64\Release\dotSimplex.pgd" 
/SUBSYSTEM:CONSOLE /MANIFESTUAC:"level='asInvoker' uiAccess='false'" /ManifestFile:"x64\Release\dotSimplex.exe.intermediate.manifest" /OPT:ICF /ERRORREPORT:PROMPT /NOLOGO /TLBID:1 

