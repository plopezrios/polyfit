POLYFIT
=======
Toolbox for performing and manipulating polynomial fits to data with
statistical uncertainties.

POLYFIT uses Monte Carlo resampling to obtain uncertainties on
quantities of interest, such as values and derivatives of the fit
function, intersections between fits to different datasets, etc.
Early versions of POLYFIT have been used in e.g.,
https://doi.org/10.1103/PhysRevB.98.161105 and
https://doi.org/10.1103/PhysRevLett.120.177701 .

POLYFIT is a stand-alone Fortran utility which uses a LAPACK library.
The code can be compiled with:

```
gfortran -o polyfit polyfit.f90 -llapack
```

The compile.sh script is provided for reference only.

POLYFIT can be run with, e.g.,

```
./polyfit
./polyfit < input > output
```

POLYFIT uses a command-line interface which can be scripted easily,
and the main output entities are easily grep-able (notice the lines
starting FIT, EVAL, and INTR in the example below).  Type 'help' at
the prompt for instructions.

Usage example
=============
The following, which uses data reported in Phys. Rev. B 98, 161105
(2018), summarizes the main capabilities of the code:


```
$ cat ta_dmc_ecorr_fs.dat
#rs    N  ---------- E +/- dE ----------
0.5    7  0.156770148194  0.000016980203
0.5    9  0.105752542851  0.000005619661
0.5   11  0.074467669241  0.000007081585
0.5   15  0.039049986116  0.000006263365
0.5   19  0.019155648089  0.000009561691
0.5   27 -0.001492801089  0.000009103566
0.5   33 -0.010065166634  0.000008560291
0.5   40 -0.016201333524  0.000094154615
0.5   57 -0.024517042686  0.000064299643
0.5   81 -0.029691063830  0.000035556828
0.5   93 -0.031258004106  0.000009228552
0.5  123 -0.033580685974  0.000033209748
0.5  147 -0.034639462939  0.000020005378
0.5  171 -0.035424228589  0.000011498968
0.5  179 -0.035622941555  0.000011718580
0.5  203 -0.036037163776  0.000015034583
0.5  251 -0.036750741281  0.000014417195
0.5  305 -0.037212342213  0.000013281108
0.5  515 -0.037989783855  0.000005866217
1.0    7  0.023447832983  0.000009345429
1.0    9  0.008535880330  0.000005445480
1.0   11 -0.000433666914  0.000005857720
1.0   15 -0.010092235908  0.000002589109
1.0   19 -0.015649492865  0.000007801581
1.0   27 -0.021209709060  0.000008768095
1.0   33 -0.023684653303  0.000033929297
1.0   40 -0.024988618589  0.000030250542
1.0   57 -0.027222014513  0.000047314379
1.0   81 -0.028463844564  0.000022308815
1.0   93 -0.028895942329  0.000013667873
1.0  123 -0.029421765996  0.000017225948
1.0  147 -0.029641279502  0.000033247655
1.0  171 -0.029884594743  0.000006840416
1.0  179 -0.029960018358  0.000014933905
1.0  203 -0.030002124681  0.000010741907
1.0  251 -0.030159975612  0.000018393002
1.0  305 -0.030294851628  0.000006599289
1.0  515 -0.030481575436  0.000002406513
1.0 1021 -0.030577406470  0.000002554204
5.0    7 -0.013453644024  0.000016020829
5.0   19 -0.014861414289  0.000005634764
5.0   33 -0.015179755760  0.000006038090
5.0   57 -0.015231683300  0.000006467927
5.0   93 -0.015245851200  0.000004085820
5.0  147 -0.015248988200  0.000004125453
5.0  203 -0.015256079278  0.000003169073
5.0  251 -0.015265327018  0.000006072534
5.0  305 -0.015261594895  0.000002068126
5.0  515 -0.015265406209  0.000000490339
$ cat input
# Enable input echo.
set echo
# Load datafile using x/y/dy from columns 2/3/4, weights from column 2
# (see wexp later), splitting into different datasets for different values
# of column 1.
load ta_dmc_ecorr_fs.dat type xydyw using 2 3 4 2 by 1
# Define X (independent variable in fit) as X=1/x, where x is the "x" value
# in the dataset.
set xscale reciprocal
# Set the exponent of the weights to 2 -- i.e., the fit weights are column
# 2 (as per the above) squared.
set wexp 2
# Set the fit exponents.
set fit 0 4/3 5/3 6/3
# Restrict the data range for the fit.
set range x>=15
# Perform the fit and report parameters.
fit
# Evaluate the value of the fit at X=0.
evaluate f at X=0
# Just to double-check, evaluate the first derivative of the fit at X=0,
# which must be zero.
evaluate f' at X=0
# Plot the fit (including uncertainties) to a file.
plot f at X=0:0.25:1000 to plot_ta_dmc_ecorr_fs.dat
# Find the intersection between each pair of datasets between X=0.01 and 0.1.
intersect between 0.01 0.1
$ ./polyfit < input > out
$ cat out
====================================
POLYFIT - polynomial fitting toolbox
====================================

Type "help" for a list of commands.

POLYFIT> 
Enabled input echo.

POLYFIT> # Load datafile using x/y/dy from columns 2/3/4, weights from column 2
POLYFIT> # (see wexp later), splitting into different datasets for different values
POLYFIT> # of column 1.
POLYFIT> load ta_dmc_ecorr_fs.dat type xydyw using 2 3 4 2 by 1

Loaded data from "ta_dmc_ecorr_fs.dat" as dataset #1, type xydyw, 19 data.
Loaded data from "ta_dmc_ecorr_fs.dat" as dataset #2, type xydyw, 20 data.
Loaded data from "ta_dmc_ecorr_fs.dat" as dataset #3, type xydyw, 10 data.

POLYFIT> # Define X (independent variable in fit) as X=1/x, where x is the "x" value
POLYFIT> # in the dataset.
POLYFIT> set xscale reciprocal

Set xscale to reciprocal for all sets.

POLYFIT> # Set the exponent of the weights to 2 -- i.e., the fit weights are column
POLYFIT> # 2 (as per the above) squared.
POLYFIT> set wexp 2

Set wexp for all sets.

POLYFIT> # Set the fit exponents.
POLYFIT> set fit 0 4/3 5/3 6/3

Fit form set to:
  Y = k1 + k2*X^1.3333333333333333 + k3*X^1.6666666666666667 + k4*X^2
  Shared coefficients reset to: none

POLYFIT> # Restrict the data range for the fit.
POLYFIT> set range x>=15

Range set for all sets.

POLYFIT> # Perform the fit and report parameters.
POLYFIT> fit

Fit parameters:

       Set      i             ki                  dki       
     -------------------------------------------------------
FIT      1      1   -3.877698393408E-02   9.979537652320E-06
FIT      1      2    3.240753252141E+00   6.340463599246E-02
FIT      1      3    3.820102694279E-01   3.949454437262E-01
FIT      1      4   -3.158386058865E+00   6.042788345717E-01
FIT      2      1   -3.065007079828E-02   3.595943822482E-06
FIT      2      2    7.009000856591E-01   4.067981234242E-02
FIT      2      3    3.308420296878E-01   2.562282710891E-01
FIT      2      4   -4.594152624760E-01   3.944106277702E-01
FIT      3      1   -1.527297002434E-02   2.108904165945E-06
FIT      3      2    8.035107968844E-02   2.449426893071E-02
FIT      3      3   -5.363036341001E-01   1.583437516559E-01
FIT      3      4    1.005441487685E+00   2.510078282147E-01
     -------------------------------------------------------

Fit assessment:

        Measure            Value                Stderr      
     -------------------------------------------------------
FIT     chi^2        6.849303445188E-04   1.432051256552E-04
FIT     chi^2/Ndf    2.283101148396E-05   4.773504188506E-06
FIT     RMS(Y-f)     1.502843092500E-02   1.567362253790E-03
     -------------------------------------------------------

Fit in XMGRACE format:
  Set #1: y=-3.8776983934078940E-002+3.2407532521406792*x^1.3333333333333333+0.38201026942791011*x^1.6666666666666667-3.1583860588649100*x^2
  Set #2: y=-3.0650070798276103E-002+0.70090008565909034*x^1.3333333333333333+0.33084202968781667*x^1.6666666666666667-0.45941526247603687*x^2
  Set #3: y=-1.5272970024344511E-002+8.0351079688436772E-002*x^1.3333333333333333-0.53630363410012438*x^1.6666666666666667+1.0054414876847220*x^2

POLYFIT> # Evaluate the value of the fit at X=0.
POLYFIT> evaluate f at X=0

       set              X                    f                   df       
      --------------------------------------------------------------------
EVAL     1    0.000000000000E+00  -3.877697920277E-02   1.012313420160E-05
EVAL     2    0.000000000000E+00  -3.065011372705E-02   3.616442398085E-06
EVAL     3    0.000000000000E+00  -1.527299882541E-02   2.119293647237E-06
      --------------------------------------------------------------------

POLYFIT> # Just to double-check, evaluate the first derivative of the fit at X=0,
POLYFIT> # which must be zero.
POLYFIT> evaluate f' at X=0

       set              X                    f                   df       
      --------------------------------------------------------------------
EVAL     1    0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
EVAL     2    0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
EVAL     3    0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
      --------------------------------------------------------------------

POLYFIT> # Plot the fit (including uncertainties) to a file.
POLYFIT> plot f at X=0:0.25:1000 to plot_ta_dmc_ecorr_fs.dat

Plot saved to "plot_ta_dmc_ecorr_fs.dat".

POLYFIT> # Find the intersection between each pair of datasets between X=0.01 and 0.1.
POLYFIT> intersect between 0.01 0.1

Intersections:
       Sets          X0                   DX0                  Y0                   DY0         
INTR   1   2   1.405982396176E-02   2.305068495826E-05  -2.809133828798E-02   1.380193045185E-05
INTR   1   3   2.589765766944E-02   2.631303203173E-05  -1.519874461894E-02   6.952810779026E-06
INTR   2   3   5.481057313347E-02   4.141980752275E-05  -1.482100287717E-02   1.064268386707E-05
POLYFIT> 
Quitting.
```
