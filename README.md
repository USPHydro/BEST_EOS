'''

Copyright (c) 2021.

This file was created by: 

Paolo Parotto, Department of Physics, University of Wuppertal

Debora Mroczek, Department  of  Physics,  UIUC

Jamie Stafford, Department of Physics, University of Houston

'''

READ EoS.readme

Local lattice-only options added to this version:

./EoS -L  <param_file>  generates the lattice-only Taylor-expanded EoS
./EoS -LF <param_file>  applies Gaussian smoothing to the lattice-only pressure
./EoS -LH <param_file>  smoothly matches the lattice-only EoS to HRG
./EoS -LFH <param_file> combines smoothing and HRG matching

The -S option can be combined with the commands above to use the
strangeness-neutral input and append _SN to the output directory.

For lattice-only output without HRG matching, the generated directories are
LATonly, LATonly_smooth, LATonly_SN, or LATonly_smooth_SN.

For lattice-only output with HRG matching, the generated directories are
CROSSOVER_T<T0>, CROSSOVER_smooth_T<T0>, CROSSOVER_T<T0>_SN, or
CROSSOVER_smooth_T<T0>_SN. Here <T0> is the transition temperature at
muB = 0, rounded to the nearest integer; for T0 = 155 MeV this gives T155.
The HRG matching uses the same shifted pseudocritical curve as the full EoS,
with T_match(muB) = T0 + kappa/T0 * muB^2 - 23 MeV and DeltaT = 17 MeV.


To invert the EOS tables from TmuB to enB you should type 

python ./utilities/mapEOS_TmuB_to_enB.py [Folder with EOS Tables WITHOUT / in the end]
