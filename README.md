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
The following, which uses data reported in
https://doi.org/10.1103/PhysRevB.98.161105, summarizes the main
capabilities of the code.  This example can be found under
`examples/heg_demo`.


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
load ta_dmc_ecorr_fs.dat type xydyw using 2,3,4,2 by $1
# Define X (independent variable in fit) as X=1/x, where x is the "x" value
# in the dataset.
set X 1/x
# Set the exponent of the weights to 2 -- i.e., the fit weights are column
# 2 (as per the above) squared.
set wexp 2
# Set the fit exponents.
set fit 0,4/3,5/3,6/3
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
plot f at X=0:0.25:1001 to plot_ta_dmc_ecorr_fs.dat
# Find the intersection between each pair of datasets between X=0.01 and 0.1.
intersect between 0.01 0.1
$ ../../polyfit < input > out
$ cat out
====================================
POLYFIT - polynomial fitting toolbox
====================================

Type "help" for a list of commands.

POLYFIT> POLYFIT> 
Enabled input echo.

POLYFIT> # Load datafile using x/y/dy from columns 2/3/4, weights from column 2
POLYFIT> # (see wexp later), splitting into different datasets for different values
POLYFIT> # of column 1.
POLYFIT> load ta_dmc_ecorr_fs.dat type xydyw using 2,3,4,2 by $1

Loaded data from "ta_dmc_ecorr_fs.dat" as dataset #1, type xydyw, 19 data.
Loaded data from "ta_dmc_ecorr_fs.dat" as dataset #2, type xydyw, 20 data.
Loaded data from "ta_dmc_ecorr_fs.dat" as dataset #3, type xydyw, 10 data.

POLYFIT> # Define X (independent variable in fit) as X=1/x, where x is the "x" value
POLYFIT> # in the dataset.
POLYFIT> set X 1/x

Set X=1/x for set #1.
Set X=1/x for set #2.
Set X=1/x for set #3.

POLYFIT> # Set the exponent of the weights to 2 -- i.e., the fit weights are column
POLYFIT> # 2 (as per the above) squared.
POLYFIT> set wexp 2

Set wexp for set #1.
Set wexp for set #2.
Set wexp for set #3.

POLYFIT> # Set the fit exponents.
POLYFIT> set fit 0,4/3,5/3,6/3

Fit #1 set to:
  Y = k1 + k2*X^1.3333333333333333 + k3*X^1.6666666666666667 + k4*X^2
  Shared coefficients reset to: none

POLYFIT> # Restrict the data range for the fit.
POLYFIT> set range x>=15

Range set.

POLYFIT> # Perform the fit and report parameters.
POLYFIT> fit

Fit parameters:

       Set      i             ki                  dki       
     -------------------------------------------------------
FIT      1      1   -3.877665448852E-02   9.997444091167E-06
FIT      1      2    3.238534447804E+00   6.249673728661E-02
FIT      1      3    3.956825777627E-01   3.891929731896E-01
FIT      1      4   -3.178939483208E+00   5.957812738965E-01
FIT      2      1   -3.065004306850E-02   3.522825931985E-06
FIT      2      2    7.009220545903E-01   3.972702931493E-02
FIT      2      3    3.314873972041E-01   2.503286825419E-01
FIT      2      4   -4.612722504276E-01   3.855576990937E-01
FIT      3      1   -1.527302011212E-02   2.137060593291E-06
FIT      3      2    8.082765122309E-02   2.490596186053E-02
FIT      3      3   -5.391506362269E-01   1.611772026377E-01
FIT      3      4    1.009653032164E+00   2.556818746999E-01
     -------------------------------------------------------

Fit assessment:

        Measure            Value                Stderr      
     -------------------------------------------------------
FIT     chi^2        6.871658180043E-04   1.443024871481E-04
FIT     chi^2/Ndf    2.290552726681E-05   4.810082904936E-06
FIT     RMS(Y-f)     1.505220249673E-02   1.577013018238E-03
     -------------------------------------------------------

Fit in XMGRACE format:
  Set #1: y=-3.8776654488518227E-002+3.2385344478042342*x^1.3333333333333333+0.39568257776266547*x^1.6666666666666667-3.1789394832081102*x^2
  Set #2: y=-3.0650043068504659E-002+0.70092205459025703*x^1.3333333333333333+0.33148739720408921*x^1.6666666666666667-0.46127225042757053*x^2
  Set #3: y=-1.5273020112118715E-002+8.0827651223087382E-002*x^1.3333333333333333-0.53915063622694726*x^1.6666666666666667+1.0096530321643367*x^2

POLYFIT> # Evaluate the value of the fit at X=0.
POLYFIT> evaluate f at X=0

       set              X                    f                   df       
      --------------------------------------------------------------------
EVAL     1    0.000000000000E+00  -3.877700446503E-02   1.005932917176E-05
EVAL     2    0.000000000000E+00  -3.065006796243E-02   3.604887247724E-06
EVAL     3    0.000000000000E+00  -1.527302407867E-02   2.125590973695E-06
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
POLYFIT> plot f at X=0:0.25:1001 to plot_ta_dmc_ecorr_fs.dat

Plot saved to "plot_ta_dmc_ecorr_fs.dat".

POLYFIT> # Find the intersection between each pair of datasets between X=0.01 and 0.1.
POLYFIT> intersect between 0.01 0.1

Intersections:
       Sets          X0                   DX0                  Y0                   DY0                missfrac      
INTR   1   2   1.405926632414E-02   2.259229508922E-05  -2.809157279243E-02   1.354188631307E-05   0.000000000000E+00
INTR   1   3   2.589723259051E-02   2.641208564420E-05  -1.519881315775E-02   7.017734351594E-06   0.000000000000E+00
INTR   2   3   5.481088932129E-02   4.155159859833E-05  -1.482108079002E-02   1.081200999301E-05   0.000000000000E+00

POLYFIT> 
Quitting.
```
