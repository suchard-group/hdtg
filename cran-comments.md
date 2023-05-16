## Initial submission of patch
* fixed undefined behavior (uninitialized value) in 'see2neon.h' that
  caused a WARNING on MacOS

## R CMD check results
* There were no ERRORs and WARNINGs.
* There are (occasionally) 2 NOTEs: 
  - "  Compilation used the following non-portable flag(s):
    '-mavx' '-mavx2' '-mfma' '-msse' '-msse3' '-msse4.1' '-msse4.2'
    '-mssse3'"

  But this is after configure file checks for system availability.
  
  - checking installed package size ... NOTE
    installed size is 10.3Mb
    sub-directories of 1Mb or more:
      libs 10.2Mb

  This occurs on systems (like 'r-devel-linux-x86_64-fedora-clang') that include debug
  symbols in their compilation; 'hdtg' performance is heavily dependent on many template
  instantiations that generate a large library when including debug symbols.  Future
  availability of C++17 'if (constexpr ...)' should decrease library size substantially.  
  
## Test environments
* local OS X install, R 4.1.2
* CentOS Linux 7, R 4.1.0, gcc 10.2.0
* win-builder (devel, release, oldrel)
* mac-builder (R 4.2.0, macOS 11.5.2, Apple M1)  
  
## Downstream dependencies
There are currently no downstream dependencies.
