## Resubmission of archived package

This is a resubmission of the `hdtg` package, which was archived on CRAN on 2023-05-01 due to unresolved issues. 
All issues have now been addressed.

## Initial submission of patch
* fixed undefined behavior (uninitialized value) in 'see2neon.h' that
  caused a WARNING on MacOS

## R CMD check results
* No ERRORs or WARNINGs.
* One NOTE remains, expected for resubmission of an archived package:
  - "New submission: Package was archived on CRAN"
  
## Test environments
* local OS X install, R 4.1.2
* CentOS Linux 7, R 4.1.0, gcc 10.2.0
* win-builder (devel, release, oldrel)
* mac-builder (R 4.2.0, macOS 11.5.2, Apple M1)  
  
## Downstream dependencies
There are currently no downstream dependencies.
