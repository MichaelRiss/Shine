Shine can simulate NMR (Nuclear Magnetic Resonance) NOESY experiments.
As a specialty it simulates 3D experiments with two INEPT steps for 
experiments such as CCH, CNH, NNH. 

For more information about these experiments see:

Tammo Diercks, Murray Coles & Horst Kessler "An efficient strategy for assignment of cross-peaks in 3D heteronuclear NOESY experiments", Journal of Biomolecular NMR, 15, 177-180, 1999

The simulation code is based on the fortran program "spirit". See:

Leiming Zhu, H. Jane Dyson and Peter E. Wright "A NOESY-HSQC simulation program, SPIRIT", Journal of Biomolecular NMR, 11(1), 17-29, 1998

The original "spirit" is able to handle NOESY experiments with one 
INEPT step like HCH and HNH.
It has been ported to C/C++ and extended to incorporate double
INEPT experiments. To avoid a name clash with "spirit" the 
C++ LL parser framework of the Boost project the program has a new 
name: "Shine".

Currently Shine is in the state of a working prototype.
The simulation works with common parameter settings and the resulting 
spectra look sound. Some more exotic parameter settings have been 
implemented but are not tested yet.

### Usage - HOWTO:

To conduct a simulation you need:
- a file with simulation parameters
- a file with the sequence of the protein
- a file with the chemical shifts of the protein
- a file with the atom coordinates of the protein (PDB-File)


#### The file with simulation parameters:
The parameters are key-value pairs with the following format:
`#<key>: <value>`

```
#freq:	<freq1> <freq2> <freq3>
```
 The frequencys of the spectrometer for the different dimensions. 
 <freq1> is the aquisition frequency and therefore always a
 proton dimension. <freq3> is always a dimension with a hetero atom.
 <freq2> is either a proton dimension or a dimension with a 
 hetero atom depending on the type of experiment. In experiments
 with one INEPT step it's a proton dimension, in experiments with 
 two INEPT steps it's the second hetero dimension.

```
#spwd:  <spwd1> <spwd2> <spwd3>
```
 The sweepwidth of the spectrometer for the different dimensions.
 For a explanation of the dimensions see #freq.

```
#size:  <size1> <size2> <size3>
```
 The size of the spectra for the different dimensions.
 For a explanation of the dimensions see #freq.

```
#rfpt:  <rfpt1> <rfpt2> <rfpt3>
```
 The reference point expressed in spectra coordinates.
 The spectra is simply a 3D array of numbers. The array coordinates
 start at 0,0,0 and end at size1-1, size2-1, size3-1 (ANSI C array notation).
 But you usually don't want to see the array coordinates in your spectra 
 viewer, you want to see the corresponding chemical shifts. The scaling between
 shifts and array can be calculated from the spectrometer frequency/sweepwidth
 and the array size. What's left is the "pinning" of a spectra coordinate to 
 a chemical shift or vice versa. This is done by this parameter and the 
 following parameter.
 NOTE: Although the implementation uses a C/C++ data structures and therefore
       the arrays start at 0.0, you have to provide fortran style coordinates
       here and let the arrays start at 1.0. This is a compatibility decision 
       in order to stay in sync with the nmrview parameterfile. 
       This might change in the future.
       Confused? 
       Just put a 1.0 in here if the reference point is at the start of 
       the array. :)

```
#rfpm:  <rfpm1> <rfpm2> <rfpm3>
```
 The reference point expressed in chemical shifts.
 See also #rfpt.

```
#lwfc: <lwfc1> <lwfc2> <lwfc3>
```
 The line width factor correction for the different dimensions.
 With this parameter it's possible to control the general line width
 for one dimension of the spectra. The default value is 1.0, greater
 values result in broader line shapes.

```
#seqfile: <filename>
```
 The name of the file with the sequence.
 Format see below.

```
#shiftfile: <filename>
```
 The name of the file with the chemical shifts.
 Format see below.

```
#pdbfile: <filename>
```
 The name of the PDB-file.

```
#label: <value>
```
 Parameter controlling which atoms are labeled. Possible <value>s are:
 "C13", "N15", "double".
 Until now only the "double" option has been tested.

```
#tumbling: <value>
```
 Controls if the protein is tumbling isotropic (around all axis) or
 if it has a prefered axis of rotation. Currently this parameter
 has no effect.

```
#symatom1_resID: <value>
#symatom1_Name: <value>
#symatom2_resID: <value>
#symatom2_Name: <value>
```
 When set these parameters define the prefered axis of tumbling.
 These parameters have not been tested yet. Only isotropic tumbling
 has been used.

```
#isocorrel_ns: <value>
```
 The isocorrelation time in ns, e.g. 10.0

```
#relax_delay_sec: <value>
```
 The relaxation delay in seconds, e.g. 1.0

```
#Mixing_times_sec: <nrmix> <mix1> [<mix2> ...]
```
 It's possible to define a whole set of mixing times. The experiments
 for all these mixing times should then be simulated in one go.
 But currently the simulation supports only one mixing time.
 <nrmix> specifies the number of mixing times (set it to 1).
 The following values specify each mixing time (use only one).

```
#noemin: <value>
```
 The minimum intensity for peaks. Peaks with lower intensities are ignored.
 You will have to play around with this setting. Experiments with C13 INEPT
 steps tend to produce weak peak intensities and you will set this parameter 
 lower (e.g. 0.000002), experiments with N15 produce stronger peaks so you
 should set this parameter higher (e.g. 0.002).

```
#seq_start: <value>
```
 The number of the first residue in the chain. Sometimes the chains you want
 to analyze are fractions of larger proteins and therefore the residueID 
 numbers don't start with 1. In such cases you have to specify the start 
 residueIDs here. This is necessary to synchronize sequence, shifts and PDB 
 information.

The parameters so far are more or less common to all experiments. The 
following parameters describe special properties of experiments or experiment
groups, so they only make sense in a certain context.

The most important parameter is:
```
#heteromode: single|double
```
 This parameter controls whether a experiment has one INEPT step ("single") 
 or two INEPT steps ("double").
 For "single" mode these parameters are used:
 ```
 #heterobound:
 #hetbound_sensitivity_enhanced:
 #delay_proton:
 #delay_hetbound:
 #inept_hetbound:
 #inept_hetbound_rev:
 #inept_hetbound_MQ:
 ```

 For "double" mode these parameters are used:
 ```
 #heterobound:
 #hetbound_sensitivity_enhanced:
 #heterofree:
 #hetfree_sensitivity_enhanced:
 #delay_hetbound:
 #delay_hetfree:
 #inept_hetbound:
 #inept_hetbound_rev:
 #inept_hetbound_MQ:
 #inept_hetfree:
 #inept_hetfree_rev:
 #inept_hetfree_MQ:
 ```

Some parameters appear in both modes and they usually have a similar meaning.

As explained at #freq: the first dimension corresponds to
the aquisition dimension and is therefore always a proton dimension.
The third dimension is always a hetero dimension (C13 or N15).
First and third dimension are "bound", that means the proton of the
first dimension is bound to the hetero atom of the third dimension.
Therefore the third dimension is called heterobound or short "hetbound".
The second dimension is either a proton or a hetero dimension, depending
on the experiment type. The atoms of the second dimension are "free", 
that means they are not bound to atoms of the first or third dimension.
If used as hetero dimension it's called "hetfree".

The meaning of parameters in "single" mode:

```
#heterobound: C13|N15
```
 Specifies the type of the hetero atom of the third dimension. This basically
 chooses between HCH and HNH experiments.

```
#hetbound_sensitivity_enhanced: true|false
```
 Whether sensitivity enhancement is used. While implemented this code path
 is not yet tested, so set "false" here.
 If you want to use sensitivity enhancement you also have to set 
 #inept_hetbound_MQ:.

```
#delay_proton: true|false
```
 Whether a initial delay of 1/(2*sw) should be used in the second dimension
 (free proton). This has an effect on folded spectra, folded peaks will have
 a intensity with negative sign.
 Set to "false", "true" hasn't been tested.

```
#delay_hetbound: true|false
```
 Whether a initial delay of 1/(2*sw) should be used in the third dimension
 (heterobound dimension). This has an effect on folded spectra, folded peaks 
 will have a intensity with negative sign.
 Set to "false", "true" hasn't been tested.

```
#inept_hetbound: <value>
```
 First INEPT delay in ns. Typical values: 1.7 for C13 dimensions and 2.5 for
 N15 dimensions.

```
#inept_hetbound_rev: <value>
```
 Second INEPT delay (reverse INEPT) in ns. 
 Typical values again: 1.7 for C13 dimensions and 2.5 for N15 dimensions.

```
#inept_hetbound_MQ: <value>
```
 Third INEPT delay in ns; used for sensitivity enhanced experiments. 
 Only necessary when #hetbound_sensitivity_enhanced: is set to true.
 This option has not been tested.

The meaning of parameters in "double" mode:

```
#heterobound: C13|N15
```
 Specifies the type of hetero atom of the third dimension (hetero atom bound
 to the proton of the first dimension).

```
#heterofree: C13|N15
```
 Specifies the type of hetero atom in the second dimension (free hetero atom).

These two parameters specify the type of the "double" INEPT experiment.
```
 #heterobound: C13
 #heterofree: C13   => CCH
```

```
 #heterobound: N15
 #heterofree: N15   => NNH
```

```
 #heterobound: N15
 #heterofree: C13   => CNH (H bound to N15)
```

```
#hetbound_sensitivity_enhanced: true|false
```
 Whether sensitivity enhancement is used. While implemented this code path
 is not yet tested, so set "false" here.
 If you want to use sensitivity enhancement you also have to set 
 #inept_hetbound_MQ:.

```
#hetfree_sensitivity_enhanced: true|false
```
 Whether sensitivity enhancement is used. While implemented this code path
 is not yet tested, so set "false" here.
 If you want to use sensitivity enhancement you also have to set 
 #inept_hetfree_MQ:.

```
#delay_hetbound: true|false
```
 Whether a initial delay of 1/(2*sw) should be used in the third dimension
 (heterobound dimension). This has an effect on folded spectra, folded peaks 
 will have a intensity with negative sign.
 Set to "false", "true" hasn't been tested.

```
#delay_hetfree: true|false
```
 Whether a initial delay of 1/(2*sw) should be used in the second dimension
 (free hetero dimension). This has an effect on folded spectra, folded peaks 
 will have a intensity with negative sign.
 Set to "false", "true" hasn't been tested.

```
#inept_hetbound: <value>
```
 First INEPT delay for the bound hetero atom (third dimension) in ns. 
 Typical values: 1.7 for C13 and 2.5 for N15.

```
#inept_hetbound_rev: <value>
```
 Second INEPT delay (reverse INEPT) for the bound hetero atom (third dimension)
 in ns. 
 Typical values: 1.7 for C13 and 2.5 for N15.

```
#inept_hetbound_MQ: <value>
```
 Third INEPT delay for the bound hetero atom (third dimension) in ns; used 
 for sensitivity enhanced experiments. 
 Only necessary when #hetbound_sensitivity_enhanced: is set to true.
 This option has not been tested.

```
#inept_hetfree: <value>
```
 First INEPT delay for the free hetero atom (second dimension) in ns. 
 Typical values: 1.7 for C13 and 2.5 for N15.

```
#inept_hetfree_rev: <value>
```
 Second INEPT delay (reverse INEPT) for the free hetero atom (second dimension)
 in ns. 
 Typical values: 1.7 for C13 and 2.5 for N15.

```
#inept_hetfree_MQ: <value>
```
 Third INEPT delay for the free hetero atom (second dimension) in ns; used 
 for sensitivity enhanced experiments. 
 Only necessary when #hetfree_sensitivity_enhanced: is set to true.
 This option has not been tested.


#### The file with the sequence of the protein:
This file contains a sequence of letters. Each letter describes a
amino acid according to the usual encoding. 
A == Alanine
R == Arginine
...
There must not be gaps between the letters.
The residueID of the first amino acid has to be set as #seq_start:
in the simulation parameter file.


#### The file with the chemical shifts:
This file contains four columns. In the first column there is the 
residueID (the number of the amino acid in the chain).
In the second column there is the three-letter code of the amino acid.
ALA == Alanine
ARG == Arginine
...
In the third column there is the name of the atom in PDB encoding.
The fourth column contains the chemical shift of the atom.
Missing entries are automatically set to NaN which are ignored during the
course of simulation.


#### The file with the atom coordinates is a standard PDB-File.


To run a simulation copy all four files in a directory and type:
`Shine <parameterfile>`

<parameterfile> has to point to the file with the simulation parameters.
A file called "output.txt" will be generated. 
It basically contains a peaklist with peak coordinates relative to the
spectra coordinates (and not expressed in chemical shifts).

To generate a spectra:
`make3D_C++ <peaklistfile> <spectrafile>`
This will generate a plain spectra file, i.e. without submatrix structure.

As most NMR viewers need spectra files with submatrix structure there is
one conversion tool included to convert the plain spectra file to 
nmrview format. It's called "plain2nv" and needs the plain spectra and
the ".par" file with the spectra description in nmrview format.

### Examples:
In the examples folder you will find parameter files to simulate
the CCH-, CNH-, HCH-, HNH- and NNH-experiment on the kdp protein.
Type:
```
cd examples/kdp/CCH
../../../Shine simulationparameters
../../../make3D_C++ output.txt cchnoesy.plain
../../../plain2nv/plain2nv cchnoesy cchnoesy.mat
cd ../CNH
../../../Shine simulationparameters
../../../make3D_C++ output.txt cnhnoesy.plain
../../../plain2nv/plain2nv cnhnoesy cnhnoesy.mat
cd ../HCH
../../../Shine simulationparameters
../../../make3D_C++ output.txt hchnoesy.plain
../../../plain2nv/plain2nv hchnoesy hchnoesy.mat
cd ../HNH
../../../Shine simulationparameters
../../../make3D_C++ output.txt hnhnoesy.plain
../../../plain2nv/plain2nv hnhnoesy hnhnoesy.mat
cd ../NNH
../../../Shine simulationparameters
../../../make3D_C++ output.txt nnhnoesy.plain
../../../plain2nv/plain2nv nnhnoesy nnhnoesy.mat
```

This should generate nmrview spectra in the corresponding folders.
Try to start with these examples and go on from there in small steps.


### Hacking the Code

The code is unfortunately in a messy state, you just have to look at the
numerous untested options above. :-\
It's a fresh prototype which just succeeded simulating the first 
experiments with sound results. So you should always be cautious
with the results.

The program was coded in three epochs. The first epoch was the port of 
the fortran code to C++, the second epoch introduced different data 
structures and the third epoch added double INEPT capability.
Maybe this helps to understand why e.g. there are C++ data structures with
pointers to model a protein and on the other side there are fortran style 
matrices for the simulation. 

Some words to the data structures:
During development I struggled with incomplete data sets, shifts were 
missing, atoms were named in different nomenclatures etc.
To solve this problem I separated the construction in two steps:
The first step reads in the sequence (which has to be complete!)
and creates an "empty" protein data structure. "Empty" means that all
fields for chemical shift or geometric atom coordinates are initialized
with "NaN" (Not a Number) to indicate missing data.
The second step reads the file with the chemical shifts and the PDB file 
and as far as data is available it is filled into the data structure.


Things to do (if [I|you] have some spare time):
- Implement a plain2sparky converter.
  The ucsf-format used by sparky looks like a normal submatrix format,
  so it shouldn't be a big problem.

- Optimization. The combination of a complex C++ pointer structure and
  the use of matrices is far from efficient. Either we switch back to a
  pure matrix approach or we have to come up with a clever approach to
  combine these two paradigms.

- Those who have experience with sensitivity enhancement or even
  the opportunity to conduct real experiments might have a look at
  the sensitivity enhancement options 
  ```
  #hetbound_sensitivity_enhanced:
  #hetfree_sensitivity_enhanced:
  #inept_hetbound_MQ:
  #inept_hetfree_MQ:
  ```
  I ported them from the original spirit, but I don't have the 
  means/example spectra to test them.
  Chances are high that it works out of the box or with little bug
  squashing.

- Other features of that kind are the delay - parameters
  ```
  #delay_proton:
  #delay_hetbound:
  #delay_hetfree:
  ```
  But as they are used in areas of code which have seen substantial
  change during my modifications it's less likely they will work correctly
  out of the box. Sorry.

- Leak parameters: spirit uses for the approximation of the INEPT step
  "leak parameters". There are some issues with that:
  - One leak parameter was not declared correctly (leak2 in makepp.f) 
    so that the initialization did not result in 0.75 but in 0.0.
    Right now Shine behaves the same and is bug(?) compatible.
    This might change in the next version. See the LEAK2 definition in
    NMRdata.cc if want to change it right now.
  - For some atoms the leak parameters are missing. These atoms wont
    produce peaks in the simulation. If possible these leak parameters
    should be found and inserted. Have a look at the commented regions 
    in BuildSequence.cc. 
  - The CA atom of Proline gets the leak parameter of the N atom.
    This might result in the correct value for that atom, but might also
    be a typo. Again, the correct leak parameter should be found and
    this value should be tested for correctness.

- At some point during the implementation of the double INEPT step I 
  dropped support for running a simulation for several mixing times at
  once. This performance enhancing feature complicated the access to the
  data structures and I didn't have a pressing application for that.
  But if someone really needs that feature he might add it.

- The same goes for the feature to approximate internal protein movement 
  by supplying several PDB files and defining a population. While this
  is definitely a useful feature I dropped it because I had not the 
  opportunity to test it. This feature should be reintroduced.


### Outlook (things that would be nice to have):

- A standard open source library for NMR experiment simulation with an
  active comunity contributing new experiments.
  Additionally: standard file formats, standard nomenclature, ... ;)

- Operator formalism:
  I heard some talks about operator formalism and it's application
  for the description of NMR experiments. Unfortunately I don't have
  a strong background in nuclear spin physics, but this toolkit looks
  like it could be the base for a generic NMR experiment simulator. 
  Just my 5 c. ;)


Michael Riss
