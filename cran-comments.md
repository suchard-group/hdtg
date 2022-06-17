## Test environments
* local OS X install, R 4.1.2
* CentOS Linux 7, R 4.1.0, gcc 10.2.0
* win-builder (devel, release, oldrel)
* mac-builder (R 4.2.0, macOS 11.5.2, Apple M1)

## R CMD check results
* There were no ERRORs and WARNINGs.
* There are 2 NOTEs:
  - New submission
  - "  Compilation used the following non-portable flag(s):
    '-mavx' '-mavx2' '-mfma' '-msse' '-msse3' '-msse4.1' '-msse4.2'
    '-mssse3'"

    But this is after configure file checks for system availability.

## Downstream dependencies
There are currently no downstream dependencies.