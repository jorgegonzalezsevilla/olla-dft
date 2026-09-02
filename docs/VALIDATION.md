# Validation

Every number in this document comes from a real Quantum ESPRESSO 6.6 run,
compiled from source, that was prepared, executed and post-processed with
Olla-DFT itself: generate the inputs, run QE, read the results. Nothing here
is a mock-up. The raw data live in the repository: the `examples/` folders
carry the inputs, the `.dat` files and the figures of each case together with
the exact `olla-dft` commands that produced them, and `tests/datos/` holds the
real QE outputs (bands, DOS, phonons, `epsilon.x`, XANES, TDDFPT, `pwcond.x`,
electron-phonon, MD) that the test suite reads. Reference values that were
validated once against experiment or literature are frozen in
`tests/referencias.py`, where they act as regression detectors. The inputs
were also checked with the Quantum ESPRESSO parser bundled in ASE.

Three systems were chosen to cover the difficult cases: a semiconductor
(Si), a metal (Al) and a spin-polarised metal (bcc Fe); the remaining
sections extend the check to spectra, phonons, surfaces, Wannier functions,
polarisation, thermal transport and figures.

## Silicon (semiconductor): bands and DOS

`scf → nscf → dos.x → projwfc.x → bands → bands.x` ran without editing a single
generated file (`examples/demo_Si/`). Olla-DFT detects an indirect gap with
the VBM at Γ and the CBM along Γ→X (the characteristic Δ minimum of Si) and a
direct gap at Γ; both agree with the known LDA values. Integrating the DOS up
to E_F recovers the eight valence electrons.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Indirect gap (VBM Γ, CBM on Γ→X) | 0.524 eV | LDA value | known LDA result for Si |
| Direct gap at Γ | 2.56 eV | LDA value | known LDA result for Si |
| Valence electrons from ∫DOS up to E_F | 7.98 of 8 | 8 | electron count (0.2 % mesh error) |

## Aluminium (metal)

The analysis classifies Al as metallic from bands crossing E_F, with the √E
shape of a nearly free electron gas. The energy origin shifts to the Fermi
energy only, as it should for a metal.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Classification | metal (bands cross E_F) | metal | — |
| DOS(E_F) | 0.44 states/eV | √E free-electron shape | nearly-free-electron model |

## bcc iron (spin-polarised metal)

With `--mag Fe=0.7` the calculation converges to 2.28 μB per cell
(`examples/demo_Fe/`). Integrating the spin-resolved DOS, Olla-DFT recovers
the moment independently, which confirms that both spin channels are read and
separated correctly.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Magnetic moment from pw.x | 2.28 μB/cell | 2.22 μB | experiment |
| Moment from integrating the spin-resolved DOS | 2.27 μB | 2.28 μB | the pw.x value above |

## Optics (Si)

The `epsilon.x` spectrum (`examples/demo_propiedades/`) passes two independent
tests: the f-sum rule, ∫E·ε₂(E)dE = (π/2)ħω_p², holds with factor 1.000
against the plasmon frequency the code itself reports, and the Kramers-Kronig
transform of ε₂ reproduces the ε₁ written by `epsilon.x`. With a scissor of
0.65 eV (experimental gap minus LDA gap) the ε₂ peak falls exactly on the E₂
critical point of silicon.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| f-sum rule factor | 1.000 | 1 | exact sum rule |
| Plasmon frequency | 16.95 eV | 16.7 eV | experiment |
| Kramers-Kronig ε₁ vs. `epsilon.x` ε₁ | 0.3 % difference | 0 | analytic KK relation |
| ε₂ peak with scissor 0.65 eV | 4.30 eV | 4.30 eV | E₂ critical point of Si |

## Phonons (Si)

DFPT on a 2×2×2 q mesh (`examples/demo_propiedades/`). Frequencies at Γ, X,
L and W lie within 1–6 % of the inelastic neutron data, with the correct
degeneracies and no imaginary frequencies; C_v at 300 K agrees with
experiment to 1 %.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Γ TO/LO | 508.9 cm⁻¹ | 517 cm⁻¹ | inelastic neutron scattering |
| X TA | 140.8 cm⁻¹ | 150 cm⁻¹ | inelastic neutron scattering |
| X LA/LO | 406.5 cm⁻¹ | 410 cm⁻¹ | inelastic neutron scattering |
| X TO | 455.9 cm⁻¹ | 463 cm⁻¹ | inelastic neutron scattering |
| L TA | 107.7 cm⁻¹ | 114 cm⁻¹ | inelastic neutron scattering |
| L LA | 372.3 cm⁻¹ | 378 cm⁻¹ | inelastic neutron scattering |
| L LO | 408.5 cm⁻¹ | 417 cm⁻¹ | inelastic neutron scattering |
| L TO | 484.8 cm⁻¹ | 490 cm⁻¹ | inelastic neutron scattering |
| C_v at 300 K (2-atom cell) | 0.411 meV/K | 0.415 meV/K (20 J/(mol·K)) | experiment |

## Work function (graphene)

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Φ (graphene) | 4.54 eV | 4.6 eV | experiment |

## X-ray diffraction

Checked against the PDF cards: Si (111)/(220)/(311)/(400) and NaCl
(200)/(220)/(222)/(400)/(111)/(311) with Δ2θ < 0.05°. Indices are given in the
conventional cell even when the input is the primitive one, which is what
makes the hkl comparable with the literature (`examples/demo_laminar/`).

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Si (111), (220), (311), (400) positions | Δ2θ < 0.05° | PDF 27-1402 | ICDD PDF card |
| NaCl (200), (220), (222), (400), (111), (311) positions | Δ2θ < 0.05° | PDF 05-0628 | ICDD PDF card |

## Amorphous SiO₂

A melt-and-quench of 24 atoms at 2.2 g/cm³ with MACE (1.8 ps) gives the
continuous random network of corner-sharing SiO₄ tetrahedra that is the
structure of silica glass, without a single O–O bond. Coordination is measured
with a per-pair cutoff taken from the covalent radii, not with a single global
radius: with a global 3 Å cutoff the oxygens would appear "bonded" to each
other, which is the classic error of this analysis.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Si–O coordination | 4.00 | 4 | SiO₄ tetrahedra |
| O–Si coordination | 2.00 | 2 | corner-sharing network |
| O–O bonds | 0 | 0 | silica glass structure |
| Si–O distance | 1.690 Å | 1.61 Å | experiment |

## Wannier functions (Si)

The four valence bands of silicon are wannierised onto s orbitals centred on
the bonds. The centres come out 0.6788 Å from each atom along the three
directions, i.e. at 1.1756 Å = √3·a/8: the bond midpoint, to four decimals and
without anything imposing it — the Berry phase puts it there. The spread
converges with the mesh and the 6³ value matches Marzari and Vanderbilt. The
interpolation reproduces the DFT bands at points that were NOT in the mesh;
transforming the eigenvalues directly, without gauge, is 3.5 to 5.9 times
worse. Ω splits into Ω_I + Ω_D + Ω_OD exactly, Ω_I moves by no more than
10⁻¹² Å² during minimisation (it is gauge invariant, and moving would signal a
bug), and the gradient of the spread functional was checked against its own
numerical derivative: that is where a missing 1/N_k factor showed up, which
with a 4×4×4 mesh makes the step 64 times too long.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Wannier centre distance from the atom | 1.1756 Å | 1.1756 Å (√3·a/8) | bond midpoint, geometry |
| Spread per function, 4³ / 6³ / 8³ k mesh | 1.605 / 1.901 / 2.047 Å² | 1.93 Å² | Marzari and Vanderbilt |
| Max. interpolation error with gauge, 4³ / 6³ / 8³ | 275 / 108 / 48 meV | DFT bands off-mesh | direct pw.x eigenvalues |
| Max. error without gauge, 4³ / 6³ / 8³ | 962 / 473 / 281 meV | — | same |
| Change of Ω_I during minimisation | ≤ 10⁻¹² Å² | 0 | gauge invariance |

## Berry-phase polarisation (Si and cubic BN)

The ionic part is checked against its exact formula, Σ Z_a·f_a, and agrees
with what pw.x writes to the last printed decimal (5·10⁻⁶) in both systems,
including the ion-by-ion folding modulo 1 that Quantum ESPRESSO applies when
some valence charge is odd, which is also what halves the polarisation
quantum. Displacing a silicon atom by 0.16 Å moves the ionic phase by 0.204 and
the electronic phase by −0.204: the Born effective charge comes out as a
cancellation, not from nothing happening. The electronic phase of the
distorted silicon computed with `lberry` matches the one given by the Wannier
centres of the same system: two routes that share no line of code and reach
the same number, because they are the same Berry phase.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Ionic phase vs. Σ Z_a·f_a | agrees to 5·10⁻⁶ | pw.x output | last decimal printed by pw.x |
| Z* of Si (atom displaced 0.16 Å) | 0 to 10⁻¹⁴ e | 0 | exact in a homopolar crystal |
| Z* of cubic BN, 6×6 mesh, 11 points per string, 60 Ry | 1.94 e (2.01 coarse mesh) | 1.92 e | literature |
| Electronic phase, `lberry` vs. Wannier centres | 0.0884 vs. 0.0892 | same phase | two independent routes |

## Lattice thermal conductivity (Si)

The 57 displaced configurations of a 2×2×2 supercell were computed with pw.x
and the phonon Boltzmann equation in the RTA gives κ at 300 K below the
measured value by the amount expected from the RTA (10–15 % under the exact
solution) plus the fc3 cutoff. The temperature dependence is the T⁻¹ of
Umklapp processes, and half of κ is carried by phonons with mean free path
above 1.0 µm, which is exactly what mean-free-path spectroscopy measures in
silicon and the number that explains why nanostructuring silicon works so
well for thermoelectrics. The same calculation with MACE forces instead of
DFT takes 8 seconds instead of 40 minutes, reproduces the exponent and misses
the absolute value by a factor of 2: that is why the report says so every
time the forces do not come from DFT. Supercell convergence was checked
through that cheap route (2×2×2 → 3×3×3 moves κ from 50.1 to 50.8 W/m·K, i.e.
nothing), which is exactly what it is for.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| κ at 300 K, DFT forces, RTA | 101 W/m·K (96 with natural isotopes) | ~140 W/m·K | experiment |
| Temperature exponent, DFT forces | κ ∝ T⁻¹·¹⁶ | T⁻¹ | Umklapp scattering |
| Mean free path carrying half of κ | 1.0 µm | ~1 µm | mean-free-path spectroscopy |
| κ at 300 K, MACE forces | 51 W/m·K | 101 W/m·K (DFT) | this work |
| Temperature exponent, MACE forces | κ ∝ T⁻¹·⁰⁶ | T⁻¹ | Umklapp scattering |
| Supercell convergence (MACE), 2×2×2 → 3×3×3 | 50.1 → 50.8 W/m·K | converged | this work |

## Surfaces with ESM (Al(111))

With `esm_bc='bc1'` the vacuum level is zero by construction, so the work
function is directly −E_F without fitting any plateau, and it stops depending
on the vacuum — with 8, 12 and 16 Å the energy changes by 6·10⁻⁶ Ry and E_F by
0.4 meV — so one can work with half a cell. The module refuses to run `bc1`
with a net charge, which is an ill-posed problem that pw.x computes anyway (it
gave −379 and −677 Ry for the same slab with two different vacua), and it
centres the slab at z = 0 before writing anything, because ESM measures z from
the centre of the cell and a slab left where ASE puts it lands on the cell
boundary. With `bc3` and charge it checks that Φ(q) is a straight line before
reporting a capacitance: in Al(111) with ±0.04 e it is not yet (16 %
deviation), and in that case it says so instead of printing the number.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Φ of Al(111), `bc1` | 4.24 eV | 4.24–4.26 eV | experiment |
| Energy change, vacuum 8/12/16 Å | 6·10⁻⁶ Ry | 0 | vacuum independence of ESM |
| E_F change, vacuum 8/12/16 Å | 0.4 meV | 0 | vacuum independence of ESM |
| Linearity of Φ(q) at ±0.04 e | 16 % deviation (reported, no C given) | straight line | plane-capacitor electrostatics |

## Wannier disentanglement (Si)

With eight sp³ functions extracted from twelve DFT bands, projection alone
leaves the valence bands far from DFT; choosing the subspace by the
Souza-Marzari-Vanderbilt method lowers the error by a factor of 4.3. Ω_I drops
and then does NOT move during the spread minimisation, which is what has to
happen: it is gauge invariant. Bands inside the frozen window are reproduced
exactly, as the method promises.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Valence-band error, projection only | 899 meV | — | DFT bands |
| Valence-band error, after disentanglement | 208 meV | factor 4.3 better | DFT bands |
| Ω_I before → after disentanglement | 12.37 → 10.37 Å² | decreases | Souza-Marzari-Vanderbilt |
| Change of Ω_I during minimisation | 5·10⁻¹⁴ Å² | 0 | gauge invariance |
| Error inside the frozen window | 3·10⁻¹³ eV | 0 | exact by construction |

## Charged ESM capacitance (Al(111))

The capacitance was validated without any literature number, with
electrostatics: for a plane capacitor 1/C = d/ε₀, and that slope depends
neither on the material nor on the functional. Measuring C at four different
separations between 4 and 11 Å gives a straight line with slope 1/ε₀. This
validates the formula, the area and the unit conversion at once, and with them
the charged branch of the module.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Slope of 1/C vs. d (4 separations, 4–11 Å) | 1/ε₀ within 0.4 % | 1/ε₀ | plane-capacitor electrostatics |
| R² of the linear fit | 0.99998 | 1 | — |

## Figures

The figures were checked by measuring the resulting PDF (178.0 mm exactly for
double column, 86.0 mm for single column), by inspecting that the only
embedded fonts are those of the chosen family and that they are embedded as
TrueType, and by reviewing legibility at real print size rather than enlarged.
The palette was validated with a colour-separation checker that simulates
protanopia and deuteranopia.

| Quantity | Olla-DFT | Reference | Source of reference |
|---|---|---|---|
| Double-column figure width | 178.0 mm | 178 mm | journal specification |
| Single-column figure width | 86.0 mm | 86 mm | journal specification |
| Embedded fonts | chosen family only, TrueType | — | PDF inspection |
| Colour separation of the first four series (OKLab ΔE) | ≥ 11 | ≥ 8 (safe threshold) | protanopia/deuteranopia simulation |

## Built-in checks: `olla-dft selftest`

pytest checks that the code does what the code says; `selftest` compares
Olla-DFT with the world. Each check computes a quantity that somebody has
measured or derived, and contrasts it with that value and its source. `olla-dft
selftest` runs the quick checks (no Quantum ESPRESSO, seconds); `--full` adds
the ones that run pw.x on small systems (about ten minutes, needs a working
`pw.x` and `--pseudo-dir`); `--mlip` adds, separately, the check that needs
MACE; `--list` prints the table below without running anything. Tolerances
are relative except where the reference is zero, where they are absolute.

| Key | Check | Reference | Tolerance | Source | Needs |
|---|---|---|---|---|---|
| `madelung` | Madelung constant, simple cubic lattice, α_M | 2.8372974 | 1·10⁻⁵ | classical Ewald-sum value for a point charge in a neutralising background | — |
| `lorenz` | Lorenz number of a free-electron gas, L/L₀ | 1.0 | 12 % | Sommerfeld limit, L₀ = (π²/3)(k_B/e)² = 2.44·10⁻⁸ W·Ω/K² | — |
| `npw` | Plane waves of Si at 30 Ry, N_PW | 725 | 6 % | what pw.x reports for the Si primitive cell (V = 39.5 Å³) at 30 Ry | — |
| `sackur` | Translational entropy of N₂ at 298 K | 150.4 J/(mol·K) | 1 % | Sackur-Tetrode at 1 bar; NIST-JANAF tables | — |
| `allen_dynes` | Allen-Dynes T_c of aluminium | 1.18 K | 12 % | experimental T_c of Al with λ = 0.44, ω_log = 270 K (Allen-Dynes 1975) and µ* = 0.12; µ* is a fitted parameter, with 0.10 the same formula gives 1.9 K | — |
| `allen_dynes_mu` | Sensitivity of Allen-Dynes to µ*, T_c(0.10)/T_c(0.12) | 1.56 | 5 % | the formula is exponential in µ*: raising it from 0.10 to 0.12 lowers T_c to two thirds | — |
| `born2d` | Sheet moduli of an isotropic sheet, Y_2D | 341.8 N/m | 1 % | C11 = 352, C12 = 60 N/m (graphene, DFT), Y = C11 − C12²/C11 | — |
| `gap_invariante` | Band alignment removes the arbitrary zero, ΔE_v of a material with itself | 0 eV | 1·10⁻⁹ | exact identity | — |
| `ewald_escala` | Madelung constant is scale independent, \|α(L=3) − α(L=30)\| | 0 | 1·10⁻⁶ | exact invariance of the Ewald sum under a change of units | — |
| `chern_qwz` | Chern number of the Qi-Wu-Zhang insulator (lower band, m = −1) | −1 | 1·10⁻¹⁰ | Qi, Wu and Zhang, Phys. Rev. B 74, 085308 (2006) | — |
| `her_pt` | HER: platinum sits at the top of the volcano, ΔG_H* | −0.09 eV | 5 % | Nørskov et al. 2005, Pt(111) | — |
| `oer_ruo2` | OER overpotential of RuO₂(110), η | 0.48 V | 10 % | Man et al. 2011 (ChemCatChem), ΔG(OH) = 0.77, ΔG(O) = 2.16, ΔG(OOH) = 3.87 eV | — |
| `escala_oer` | OER scaling relation, ΔG(OOH) − ΔG(OH) | 3.2 eV | 10 % | universal scaling relation, 3.2 ± 0.2 eV on nearly every oxide surface | — |
| `escala_eta_min` | Scaling limit of the OER overpotential, η_min | 0.37 V | 2 % | Man et al. 2011 | — |
| `umklapp` | κ_L of silicon decays as 1/T, exponent n | 1.0 | 25 % | above the Debye temperature Umklapp scattering gives κ ∝ 1/T; the exponent, not κ, is checked | MACE (`--mlip`, ~25 s) |
| `fonon_si` | Optical mode of Si at Γ, ω(Γ) | 520 cm⁻¹ | 10 % | experimental Raman of silicon, 520.7 cm⁻¹ at 300 K | QE (~20 s) |
| `wannier_si` | Wannier centre of the Si–Si bond, \|r̄\| | 1.17563 Å | 2 % | bond centre of the diamond structure at √3·a/8 with a = 5.43 Å | QE (~30 s) |
| `condensador` | Charged ESM: slope of 1/C vs. distance over 1/ε₀ | 1.0 | 6 % | plane-capacitor electrostatics, independent of material, pseudopotential and functional | QE (~90 s) |
| `born_si` | Born effective charge of silicon, Z* | 0 e | 0.05 | exactly zero in a homopolar crystal by the acoustic sum rule | QE (~60 s) |
| `gamma_al` | Surface energy of Al(111), γ | 1.10 J/m² | 25 % | full-potential LDA (Vitos et al. 1998) gives 1.20 J/m²; polycrystalline experiment, 1.14 | QE (~60 s) |
| `bulk_si` | Bulk modulus of Si by strain, B | 95 GPa | 15 % | LDA gives 93–97 GPa (Nielsen and Martin 1985); experiment, 98 | QE (~50 s) |
| `sitio_h_al` | H on Al(111): the hollow site beats the top site, E_ads(top) − E_ads(hollow) | 5.6 eV | 60 % | hydrogen chemisorbs in the fcc(111) hollow; the order hollow < bridge < top is textbook | QE (~60 s) |

A check that comes out wrong is not always a code failure: a tight
tolerance, a different pseudopotential or a low cutoff also move it. What it
does mean is that the number changed and the reason has to be found.

## Test suite

`tests/` contains 977 pytest tests that run without Quantum ESPRESSO
(`python -m pytest -q`, under a minute). They read the real QE outputs in
`tests/datos/`, compare against the frozen references in
`tests/referencias.py`, validate every `olla-dft` command quoted in the
`examples/` READMEs and in the recipes against the actual argparse tree, and
run the whole program with the output forced to cp1252 to make sure no report
dies on a legacy Windows console. `tests/barrido_cli.sh` is the complementary
command-level regression sweep over already-computed QE outputs, with the
expected exit code declared on every line.
