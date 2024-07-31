# RelEPro
RelEPro (renamed from RelED) is a utility program to calculate one-electron properties based on
the eXact 2-Component (X2C) all-electron quasi-relativistic Hamiltonian and its local approximations.

## Recent Changes
Version 1.2.1 (June 04, 2024).

1. Several incompatibilities with gfortran 10+ have been fixed.

Version 1.2.0 (May 14, 2024).

1. The first released version.

## Features

Hamiltonians:
* non-relativistic
* spin-free quasi-relativistic X2C and its local approximations X2C-DLXR, X2C-AXR, X2C-DLU, and X2C-AU
* two-component quasi-relativistic X2C and its local approximations (FCHK file only)

Properties:
* electric dipole moment (DIP)
* contact density (CD) <sup>(1,2)</sup> with the picture-change-error correction, being used to study the Mössbauer isomer shift (IS)
* effective contact density (ED) <sup>(1,2)</sup>, being used to study the Mössbauer isomer shift (IS)
* electric field gradient (EFG) <sup>(2)</sup> as well as the related nuclear quadrupole coupling constant (NQCC) and Mössbauer nuclear quadrupole splitting (NQS) for <sup>57</sup>Fe

<font size=2>

1. The finite nuclear charge distribution must be used in the quantum chemistry calculations. See the following Table for the supported quantum chemistry programs.

2. It is necessary to modify the basis functions of the interested heavy atom to obtain reliable results. For CD and ED, the s-functions should be decontracted and supplemented with 2 to 4 additional very steep primitive s-functions. Optionally, do the same things for p-functions. For EFG, the p-functions should be decontracted and supplemented with 2 to 4 very steep primitive p-functions. Do the same things for d- or f-functions if there are also *d* or *f* valence orbitals.

</font>

## Methodology
According to the linear response approach, the general formula of a first-order one-electron property can be represented as

∂E/∂μ = tr(**P**<sup>QR</sup> ∂**H**<sup>X2C</sup>/∂μ)

where **P**<sup>QR</sup> is the quasi-relativistic (QR) density matrix constructed using occupied molecular orbitals (MOs) and ∂**H**<sup>X2C</sup>/∂μ is the first-order derivative of the X2C Hamiltonian with respect to an external perturbation μ depending on the property to be calculated.

The quasi-relativistic canonical, pseudo-canonical, or natural MOs may be calculated at the X2C (= Dirac-exact NESC), BSS (= IODKH = IOTC), or DKH<sub>n</sub> (scalar-)relativistic level of theory (e.g. X2C/DFT, DKH2/CCSD, and so on) using a third-party quantum chemistry program (see the following table),
and saved in a [MOLDEN](https://www.theochem.ru.nl/molden/) or [FCHK](https://gaussian.com) data file.
However, the MOLDEN file is not strictly defined, and therefore it has to
be converted to an "official" form with the help of
[Molden2AIM](https://github.com/zorkzou/Molden2AIM) or [Multiwfn](http://sobereva.com/multiwfn/).

**TABLE 1. Quasi-relativistic Hamiltonians with the finite nuclear model implemented in different quantum chemistry programs**

| Program      | Version | QR in **P**<sup>QR</sup> | Data format                  |
|--------------|---------|--------------------------|------------------------------|
| BAGEL        |         | sf-DKH2                  | MOLDEN                       |
| BDF          | ≥ 2018  | sf-X2C                   | MOLDEN                       |
| CFour        | ≥ 2.1   | sf-X2C                   | MOLDEN                       |
| Cologne      | ≥ 2010  | sf-X2C                   | FCHK                         |
|              |         | X2C                      | FCHK                         |
| Columbus     |         | sf-X2C, sf-BSS, sf-DKHn  | MOLDEN; through (Open)Molcas |
| Dalton       | ≥ 2015  | sf-DKH2                  | MOLDEN                       |
| Gaussian     | ≥ 09    | sf-DKH2, sf-DKH4         | FCHK                         |
|              |         | DKH4                     | FCHK                         |
| (Open)Molcas |         | sf-X2C, sf-BSS, sf-DKHn  | MOLDEN                       |
| Molpro       | ≥ 2019  | sf-X2C, sf-DKHn          | MOLDEN                       |
| MRCC         |         | sf-X2C                   | MOLDEN; through CFour        |
| NWChem       | ≥ 7.0   | sf-DKH2, sf-DKH3, sf-X2C | MOLDEN                       |
| ORCA         | ≥ 6.0   | sf-DKH2, sf-X2C          | MOLDEN                       |
| PySCF        |         | sf-X2C                   | MOLDEN                       |
| Turbomole    | ≥ 7.3   | sf-X2C, sf-BSS, sf-DKHn  | MOLDEN                       |

<font size=2>

  *Some other quantum chemistry programs with quasi-relativistic Hamiltonians do not support the finite nuclear model or cannot save MOLDEN/FCHK files, and therefore are not listed here.*

</font>

RelEPro can calculate the one-electron integrals in ∂**H**<sup>X2C</sup>/∂μ and combine them together with **P**<sup>QR</sup>
to get the one-electron property ∂E/∂μ. For theoretical details, please refer to our papers.
* ED by X2C (called Dirac-exact NESC at that time): [Ref.1](#references)
* ED by X2C and its local approximations: [Ref.2](#references)
* EFG by X2C (called Dirac-exact NESC at that time): [Ref.3](#references)
* EFG by X2C and its local approximations: [Ref.4](#references)

The local X2C Hamiltonians have been defiend in the following papers.
If there is no HA-HA bond in a molecule (HA: 5*d* or heavier atoms),
the most efficient approximations DLU and AU are usually accurate enough.
* AXR (atomic **X** with full **R**), originally called FATM (from atoms to molecules): [Ref.5](#references)
* DLXR (diagonal local **X** with full **R**) and DLU (diagonal
local unitary transformation): [Ref.6](#references)
* AU (atomic unitary transformation): [Ref.7](#references)

## Compilation

Run `make` in the RelEPro/source directory, where `gfortran` is used by default. In Makefile you may also specify other Fortran 90 compilers like `nvf90` (`pgf90`). If the compilation is successful you will see the binary program in RelEPro/bin .

For `ifort` + `mkl`, run `make -f Makefile-mkl` instead.

## Running the program

1. (Ignore this step for FCHK files.) Save a new MOLDEN file using [Molden2AIM](https://github.com/zorkzou/Molden2AIM) or [Multiwfn](http://sobereva.com/multiwfn/).
2. Edit the input file `job.inp` (see below for the descriptions) and the batch script `run1.bat` (Windows) or `run1.sh` (Linux/MacOS) in RelEPro/work as needed.
3. (a) Windows. Double-click `run1.bat` and check the results in `job.out`.

   (b) Linux/MacOS. Run the following commands in the terminal and
   check the results in `job.out`.
```bash
chmod +x run1.sh
./run1.sh
```

## Description of input file

The available keywords and options are grouped by namelists.

### `$Contrl` group

The keywords and their default options are
```plaintext
$contrl
  prop=' '  iham=0  minz=0  popu=0
$end
```

The options of `prop` are `+dip`, `+cd`, `+ed`, and `+efg`,
which perform DIP, CD, ED, and EFG calculations, respectively.
These options may also be combined together. For example, `prop='+efg+ed'` will calculate both EFG and ED.

`iham` specifies the Hamiltonian. Its options are 0 (=4), 1 (non-relativistic), 2 (X2C), 3 (X2C-DLXR), 4 (X2C-AXR), 5 (X2C-DLU), 6 (X2C-AU).

`minz` specifies the atoms to calculate CD, ED, or EFG. Ir can be 0 (all the atoms), < 0 (heavy atoms with Z = |minz|), or > 0 (heavy atoms with Z ≥ minz).

`popu = 1` does population analysis for ED and EFG. For EFG, however, population analysis of spinor orbitals has not been implemented yet.

### `$QCData` group

There is only one keyword `fnam` to specify the name of the MOLDEN/FCHK data file. For example,
```plaintext
$qcdata
  fnam="fef6.molden"
$end
```

### `$PChar` group (optional for DIP and EFG)

External point charges needed in the DIP and EFG calculations are not included in the MOLDEN/FCHK data file.
They have to be provided after the `$PChar` group explicitly.
The unit of Cartesian coordinates can be Bohr (`unit=0`; default) or Angstrom (`unit=1`).
```plaintext
$pchar unit=0 $end
Q1  x1  y1  z1
Q2  x2  y2  z2
...
```

### `$NQMDat` group (optional for EFG)

In the NQCC and Mössbauer NQS calculations, the nuclear quadrupole moment (NQM) parameters of common nuclides are taken from [Ref.8](#references).
The NQM parameters of the other nuclides may be found in [Refs.9-11](#references) and provided after the `$NQMDat` group (in millibarn).
```plaintext
$nqmdat $end
i_Atom  NQM_value_i
j_Atom  NQM_value_j
...
```

## Limitations and known issues

* The *L*- or *sp*-shells in the Pople basis sets are not supported.
They have to be split into *s*- and *p*-shells and reordered.

<a id="references"></a>
## References
1. M. Filatov, W. Zou, and D. Cremer. *Analytic calculation of contact densities and Mössbauer isomer shifts using the normalized elimination of the small-component formalism*.
J. Chem. Theory Comput., 8:875–882, 2012. [DOI: 10.1021/ct2008632](https://doi.org/10.1021/ct2008632)
2. H. Zhu, C. Gao, M. Filatov, and W. Zou. *Mössbauer isomer shifts and effective contact densities obtained by the exact two-component (X2C) relativistic method and its local variants*.
Phys. Chem. Chem. Phys., 22:26776–26786, 2020. [DOI: 10.1039/d0cp04549g](https://doi.org/10.1039/D0CP04549G)
3. M. Filatov, W. Zou, and D. Cremer. *Relativistically Corrected Electric Field Gradients Calculated with the Normalized Elimination of the Small Component Formalism*.
J. Chem. Phys., 137:054113, 2012. [DOI: 10.1063/1.4742175](https://dx.doi.org/10.1063/1.4742175)
4. W. Li, M. Filatov, and W. Zou. *Calculation of electric field gradients with the exact two-component (X2C) quasi-relativistic method and its local approximations*.
Phys. Chem. Chem. Phys., 26:18333–18342, 2024. [DOI: 10.1039/d4cp01567c](https://doi.org/10.1039/D4CP01567C)
5. D. Peng, W. Liu, Y. Xiao, and L. Cheng. *Making four- and two-component relativistic density functional methods fully equivalent based on the idea of "from atoms to molecule"*.
J. Chem. Phys., 127:104106, 2007. [DOI: 10.1063/1.2772856](https://doi.org/10.1063/1.2772856)
6. D. Peng and M. Reiher. *Local relativistic exact decoupling*. J. Chem. Phys., 136:244108, 2012. [DOI: 10.1063/1.4729788](https://doi.org/10.1063/1.4729788)
7. W. Zou, G. Guo, B. Suo, and W. Liu. *Analytic Energy Gradients and Hessians of Exact Two-Component Relativistic Methods: Efficient Implementation and Extensive Applications*. J. Chem. Theory Comput., 16:1541-1554, 2020. [DOI: 10.1021/acs.jctc.9b01120](https://doi.org/10.1021/acs.jctc.9b01120)
8. P. Pyykkö. *Year-2017 nuclear quadrupole moments*. Mol. Phys., 116:1328–1338, 2018. [DOI: 10.1080/00268976.2018.1426131](https://doi.org/10.1080/00268976.2018.1426131)
9. W. Hüttner, editor. *Dipole Moments, Quadrupole Coupling Constants, Hindered Rotation and Magnetic Interaction Constants of Diamagnetic Molecules*. Springer, Berlin, 2002.
10. N. J. Stone. Table of nuclear electric quadrupole moments. Atomic Data and Nuclear Data Tables,
111-112:1–28, 2016. [DOI: 10.1016/j.adt.2015.12.002](http://dx.doi.org/10.1016/j.adt.2015.12.002)
11. N. J. Stone. *Table of nuclear electric quadrupole moments*. [Technical Report No. INDC(NDS)-0833](https://www-nds.iaea.org/publications/indc/indc-nds-0833/), International Atomic
Energy Agency, INDC International Nuclear Data Committee, Austria, October 2021.
