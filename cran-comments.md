## Overview

Plotting parameters have been restored with on.exit() in the plotting() function (see R/plotting.R lines 31-32).
A generic print function has been added (print.missmi) to address the problem with the information that was printed to the console (see R/missmi.R lines 60-88).

── R CMD check results ──── GPAbin 1.0.6 ────
Duration: 2m 9.4s

❯ checking for future file timestamps ... NOTE
  unable to verify current time

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## Test environment

-   local R installation: R version 4.5.0 (2025-04-11 ucrt) -- "How About a Twenty-Six"
-   other platforms checked using `check_rhub()`: `linux`, `macos`, `windows`
