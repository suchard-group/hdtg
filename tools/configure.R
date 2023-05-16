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

flags <- "PKG_CXXFLAGS = -I."
if (RcppXsimd::supportsSSE()) {
	flags <- paste(flags, "-DUSE_SSE", RcppXsimd::getSSEFlags())
}

if (RcppXsimd::supportsAVX()) {
	flags <- paste(flags, "-DUSE_AVX -mfma", RcppXsimd::getAVXFlags())
}

flags_win <- "PKG_CXXFLAGS = -I. -I../inst/include -DUSE_SIMD -DUSE_SSE"

txt <- c(flags, txt)
#txt_win <- c(flags_win, txt_win)


if (.Platform$OS.type == "unix") {
	cat(txt, file = makevars_out, sep = "\n")
} else {
#	cat(txt_win, file = makevars_win_out, sep = "\n")
}

