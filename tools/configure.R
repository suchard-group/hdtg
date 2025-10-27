# Gets run when compiling

makevars_in <- file.path("src", "Makevars.in")
#makevars_win_in <- file.path("src", "Makevars.win.in")

makevars_out <- file.path("src", "Makevars")
#makevars_win_out <- file.path("src", "Makevars.win")

txt <- readLines(makevars_in)
#txt_win <- readLines(makevars_win_in)

if (getRversion() < "4.2") {
    if (!any(grepl("^CXX_STD", txt))) {
        txt <- c("CXX_STD = CXX14", txt)
    }

#    if (!any(grepl("^CXX_STD", txt_win))) {
#        txt_win <- c("CXX_STD = CXX14", txt_win)
#    }
}

platform <- R.version$platform
on_cran_debian <- grepl("linux-gnu", platform) && grepl("x86_64", platform)

flags <- "PKG_CXXFLAGS = -I. -I./xsimd"
# if (!on_cran_debian) {
#   message("SIMD optimizations (AVX/SSE) ENABLED during compilation.")
#   if (RcppXsimd::supportsSSE()) {
#   	flags <- paste(flags, "-DUSE_SSE", RcppXsimd::getSSEFlags())
#   	message("SSE")
#   }
#   if (RcppXsimd::supportsAVX()) {
#   	flags <- paste(flags, "-DUSE_AVX -mfma", RcppXsimd::getAVXFlags())
#   	message("AVX")
#   }
# } else {
#   message("SIMD optimizations DISABLED for Debian/CRAN.")
# }

# Always write Makevars
txt <- c(flags, txt)
cat(txt, file = makevars_out, sep = "\n")

# 
# if (.Platform$OS.type == "unix") {
# 	cat(txt, file = makevars_out, sep = "\n")
# } else {
# #	cat(txt_win, file = makevars_win_out, sep = "\n")
# }

