Thank you for your ongoing help!  Please see responses to previous comments as bullets below:

Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_phylosem.html>.

Please correct before 2023-09-30 to safely retain your package on CRAN.

Packages in Suggests should be used conditionally: see 'Writing R Extensions'.
This needs to be corrected even if the missing package(s) become available.
It can be tested by checking with _R_CHECK_DEPENDS_ONLY_=true.

* We have removed `phylosignal` (which was removed from CRAN) from SUGGESTS.

* We have also changed the examples section to use `if(require(semPlot)){...}`, i.e., using SUGGESTS packages conditionally.  The other packages in SUGGESTS packages are used when building vignettes.  We have added `VignetteDepends` statements to declare those conditional dependencies on SUGGESTS packages.