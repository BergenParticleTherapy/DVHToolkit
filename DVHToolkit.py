# Changelog version 1.1 (Released 2019-02-15)
# - Added LLH + parameter bound check for ParametricBootstrap
# - Added change in basin hopping parameter names for D% logit model
# - Fixed patient.calculateNTCP (wrong formula!!), used for ParametricBooststrap
# - For plotting of individual BSd models in the sigmoid, now using CI limits on parameters
#        for model selection (draw/nodraw) instead of LLH limit -- the latter cannot be
#        trusted since its limit depend on the cohort in question (and the cohort is BSd)
#   -> More trust in CI method now. Best to use for tox asymmetry is profile likelihood.
# - Fixed ECLIPSE loading
# - Added "plan" options for choosing Eclipse patients
# - Switched to "python" engine for reading CSV files, a bit more stable...

# Changelog version 1.11 (Released 2019-02-15)
# - Added support for multiple structures in DVH set for the same patient
# - Fixed extrapolated DVH file I/O

# Changelog version 1.12 (Released 2019-02-18)
# - Autoselect ECLIPSE input with support for two columns (Dose, Volume)
# - Changed DVH calculation to simple interpolation for increased stability
# - Added Mean Dose readout from ECLIPSE files

# Changelog version 1.2
# - Finalized and verified confidence interval calculations, using bootstrap subtractions
#   - Removed "pivotal CI" option, substituted with "bootstrap subtraction correction"
#     which acts on the optimal identified parameter set as well

# Changelog version 1.21
# - Added DVH plotting options

# Changelog version 1.3
# - Full support for python 3
# - Added 2-sample t-test CI tooltip
# - The CI calculation now acts directly on TD50 which is the variable of interest
#       (in logit this is -a/b)
# - More DVH plotting options (compare median values of tox / non tox vs cohorts)

# Changelog version 1.31
# - Fixed ECLIPSE input in python 3
# - More DVH plotting options

# Changelog version 1.32
# - Cleaned the output a bit
# - Still needs a bit of work. Shows conflicting TD50 values in Logit!!

# Changelog version 1.33
# - Added DVH plotting options: Mean value of all cohort, compare between struct + plan combos

# Changelog version 1.34
# - Added specific NTCP output

# Changelog version 1.35
# - Added more flexibility in the DVH calculations (both ways + comma separated list)
# - Big update in the flexibility of the DVH plotting (grouping, multiple windows, storing plots)

# Changelog version 1.4 (Released 2020-01-02)
# - Removed manual entry for "Number of Plan: " lines in start of ECLIPSE files, added automatic detection
# - Full restructuring of program code

# Changelog version 1.41 (Released 2020-01-06)
# - More flexibility when plotting aggregated DVHs

# Changelog version 1.42 (Released 2020-01-09)
# - Added dialog for renaming/combining plans and structure

# Changelog version 1.43 (Released 2020-01-10)
# - Allow for empty structures in ECLIPSE file (now: skip,

# Changelog version 1.44 (Released 2020-02-03)
# - Now using pandas for CSV / XLSX conversion from DVH (fr
# - More flexibility in DVH plotting
# - RayStation DVH Input

# Changelog version 1.45 (Released 2020-02-26)
# - Added customimization of plot window for aggrevated DVH

# Changelog version 1.5 (Released 2020-10-05)
# - Changed GUI layout for smaller footprint
#   - Moved NTCP and Bootstrapping options to their own dialog windows
# - Lots of small bugfixes
#   - zeropadding when calculating median / mean of several distinct cohort
#   - Perform aggregate calculations on all cohorts for tox / notox

# Changelog version 1.6 (Released 2020-10-15)
# - Rewrite of the bootstrap code
# - New object: ParameterSpace, to handle the model paramaters from bootstrapping: Plot, analyse etc.
# - With distinct cohorts (patient groups in the left panel), NTCP calculations and bootstrapping will be performed seperately on them.

# Changelog version 1.61 (Released 2020-10-16)
# - Bugfix in calculation of TD50 during bootstrapping

import sys
from tkinter import *
from modules import MainMenu
print(sys.executable)
print(sys.version)

root = Tk()
mainmenu = MainMenu(root)
root.mainloop()
