## Resubmission

This is a resubmission to fix '1 warning' on CRAN checks as requested.

### Comments from Victoria Wimmer (thanks for these!) with responses as bullets

If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

* We have added references (including an in-revisions description of the package itself)

In addition , we see:
Non-FOSS package license (file LICENSE)
--> Please use a FOSS license

* We now list in the DESCRIPTION: `License: GPL-3`

Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags in up to 15 .Rd files, e.g.:
      AIC.phylosem.Rd: \value
      as_fitted_DAG.Rd: \value
      as_phylo4d.Rd: \value
      as_sem.Rd: \value
      average.compare_phylosem.Rd: \value
      average.Rd: \value
      ...

* We have added a \value output for all package functions

\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or create
additionally small toy examples to allow automatic testing.
(You could also replace \dontrun{} with \donttest, if it takes longer
than 5 sec to be executed, but it would be preferable to have automatic
checks for functions. Otherwise, you can also write some tests.)

* We have a `testthtat` check that runs on GitHub, as automatic check on `phylosem`.
* We have removed the `\dontrun{}` from the `phylosem` example as requested.

Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. e.g.:
...
oldwd <- getwd()           # code line i
on.exit(setwd(oldwd))        # code line i+1
...
setwd(...)            # somewhere after
...
e.g.: R/phylosem.R
If you're not familiar with the function, please check ?on.exit. This
function makes it possible to restore options before exiting a function
even if the function breaks. Therefore it needs to be called immediately
after the option change within a function.

* This use of `setwd()` has now been removed
