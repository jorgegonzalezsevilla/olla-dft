# The physics behind Olla-DFT

## Electronic structure

This chapter documents the physics that Olla-DFT actually implements — not what the Quantum ESPRESSO manual promises — in the commands that prepare an electronic-structure calculation (`gen`, `kpath`, `info`, `prim`, `conv`, `supercell`) and in those that read its results (`bands`, `gap`, `dos`, `plot`, `effmass`, `fermi`, `unfold`, `wannier`, `topology`, `berry`, `hubbard`, `align`). Each section states which question the command answers, which formulas the code uses (citing the Python function that contains them), which QE file or physical constant every number comes from, and in which cases the result is not valid. Constants and defaults were read from the source code of version 0.35.0 (`qekit/config.py`, `qekit/core/qeout.py`, and each module).

---

### `olla-dft gen` — generate the pw.x and post-processing inputs

**What it answers.** It translates a structure (CIF, POSCAR, XYZ with cell, pw.x input or output) into a consistent set of input files for pw.x, dos.x, projwfc.x and bands.x, choosing for you the cutoffs, the k-point mesh, the high-symmetry path, the number of bands and the treatment of occupations.

**Background for non-experts.** A plane-wave DFT calculation needs four numerical decisions before it starts: (1) how many plane waves to use to describe the orbitals (the *cutoff* `ecutwfc`, a maximum kinetic energy in Rydberg: the higher, the finer the description and the more expensive the run), (2) how many for the density (`ecutrho`, which for ultrasoft or PAW pseudopotentials must be considerably larger), (3) how many k-points to sample in the Brillouin zone (the "resolution" with which one integrates over the infinite crystal), and (4) how to distribute the electrons among bands near the Fermi level (in an insulator the occupation is fixed; in a metal it is "smeared" so that the sum over k converges). Olla-DFT takes (1) and (2) from the header of the pseudopotential itself, (3) from a spacing in reciprocal space (the same idea as VASP's `KSPACING`) and (4) from whether the user declares the system to be an insulator.

A pseudopotential is a "smoothed version" of an atom: it replaces the nucleus and the core electrons by an effective potential, and only the valence electrons are solved explicitly. Every UPF file carries recommended cutoffs and the number of valence electrons, and Olla-DFT reads them. For the band structure one also needs a *high-symmetry path*: a route of straight segments between special points of the Brillouin zone (Γ, X, L…). That path is decided by the seekpath library with the convention of Hinuma et al., and it refers to a standardized primitive cell, which is why `gen` uses that cell whenever the preset includes bands.

**Formulas.**

k-point mesh from a spacing (`qekit/core/kpoints.py: kgrid_from_spacing`):

$$
n_i = \max\!\left(1,\ \left\lceil \frac{|\mathbf{b}_i|}{\Delta k} \right\rceil\right), \qquad \mathbf{b}_i = 2\pi\,(\mathbf{A}^{-1})^{\mathsf T}_{i}
$$

- $n_i$: number of mesh points along reciprocal vector $i$ (dimensionless).
- $\mathbf{b}_i$: reciprocal-lattice vector **including the factor $2\pi$**, in Å⁻¹; $\mathbf{A}$ is the cell matrix with the vectors as rows (Å).
- $\Delta k$: requested spacing in Å⁻¹. Defaults: `kspacing = 0.20` (scf) and `kspacing_nscf = 0.12` (nscf/DOS), read from `qekit/config.py: DEFAULTS`. The `--klevel` levels are `coarse 0.30`, `medium 0.20`, `fine 0.15`, `very-fine 0.10` and `gamma` (Γ only).

If along some axis the widest gap between atoms exceeds `VACIO_MINIMO = 8.0` Å (`kpoints.direcciones_con_vacio`), that $n_i$ is forced to 1. The mesh is written without shift (`0 0 0`, Γ-centred); if it is $1\times1\times1$ the card is `K_POINTS gamma`.

Cell thickness and vacuum gap (`qekit/modules/inputgen.py: espesor_celda`, `hueco_vacio`):

$$
d_i = \frac{V}{|\mathbf{a}_j \times \mathbf{a}_k|}, \qquad h_{\text{Å}} = d_i \cdot \max_m \left(f^{(i)}_{m+1} - f^{(i)}_m\right)
$$

- $d_i$: height of the cell along the normal to the plane of the other two vectors (Å); $V$ is the volume (Å³).
- $h_{\text{Å}}$: the largest gap between sorted fractional coordinates $f^{(i)}$ along axis $i$ (including the gap that wraps around the periodic boundary), converted to Å.

Recommended cutoffs (`qekit/core/pseudo.py: recommend_cutoffs`):

$$
E_{\text{wfc}} = \max_s E^{\text{UPF}}_{\text{wfc},s}, \qquad
E_{\rho} = \max\!\left(\max_s E^{\text{UPF}}_{\rho,s},\ 4\,E_{\text{wfc}}\right)
$$

- $E^{\text{UPF}}_{\text{wfc},s}$, $E^{\text{UPF}}_{\rho,s}$: suggested cutoffs in the UPF header of species $s$ (Ry), read by `pseudo.suggested_cutoffs` (attributes `wfc_cutoff`/`rho_cutoff` in UPF v2, or the text "Suggested minimum cutoff for wavefunctions/charge density" in UPF v1). Values $\le 1$ are ignored.
- If no UPF declares a cutoff, `ecutwfc = 60.0` Ry and `ecutrho = dual × ecutwfc` with `dual = 8` (config) are used. The floor $4E_{\text{wfc}}$ is the physical minimum for plane waves.

Estimated number of bands for nscf/bands (`inputgen._estimate_nbnd`):

$$
n_{\text{bnd}} = \left\lceil 1.25\cdot\frac{N_{\text{el}}}{2} + 4 \right\rceil, \qquad N_{\text{el}} = \sum_{\text{atoms}} Z^{\text{UPF}}_{\text{val}}
$$

With `--nspin 2` it is enlarged to $\lfloor 1.2\,n_{\text{bnd}}\rfloor + 2$. If any UPF does not declare `z_valence`, `nbnd` is not written and pw.x uses its own default.

MD time step (`inputgen.build_pw_input`): $\mathrm{dt}_{\text{Ry}} = \mathrm{dt}_{\text{fs}} / 0.048378$, because pw.x asks for `dt` in Rydberg atomic units (`_FS_POR_UA = 4.8378e-2` fs).

Dipole correction (`inputgen._region_vacio`): the maximum of the sawtooth is placed at the centre of the vacuum gap, `emaxpos = centre`, and its decrease occupies `eopreg = clip(gap/3, 0.02, 0.20)` (fractions of the axis). $h_{\text{Å}} \ge 5$ Å is required.

Estimated cost of a hybrid (`inputgen.generate`), measured on 2-atom silicon: $\text{factor} \approx 3 + 2.6\,n_q$, with $n_q = n_{q1}n_{q2}n_{q3}$ the exact-exchange mesh. This is an empirical rule, not a formula.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_gen` reads the structure with `qekit/core/structure.py: load` (ASE; for `POSCAR/CONTCAR` it forces `format="vasp"`, and if the file carries several images it keeps the last one). It combines `--klevel`/`--kspacing`, `--mag` (which switches on `nspin=2`), `--hubbard`, `--soc`, `--functional`, `--exx-grid` and the MD options into a `GenOptions`.
2. `qekit/modules/inputgen.py: generate` decides the working cell: if the preset is `bands` or `all`, it calls `kpoints.get_kpath` (seekpath) and uses the **standardized primitive cell** it returns; with `--primitive` it uses `structure.primitive` (spglib); otherwise the cell as given.
3. `qekit/core/pseudo.py: resolve` looks for one UPF per element in `pseudo_dir` (a file whose name starts with the symbol followed by a non-alphabetic character, extension `.upf` case-insensitive). It honours `--pseudo El=file`, and `_coherencia_de_funcional` re-selects whatever is needed so that all pseudopotentials share the same functional (preference `PBE > PBESOL > REVPBE > PZ > BLYP`).
4. `pseudo.recommend_cutoffs` sets `ecutwfc`/`ecutrho`; the user can override them with `--ecutwfc`/`--ecutrho`.
5. `kpoints.kgrid_from_spacing` produces the scf and nscf meshes; `_estimate_nbnd` the number of bands.
6. `inputgen.build_pw_input` writes `&CONTROL` (with `tprnfor`, `tstress`, `outdir='./out'`), `&SYSTEM` (`ibrav=0`, cutoffs, occupations, spin, SOC, Hubbard, `tot_charge`, dipole, hybrid, `nosym`/`noinv`), `&ELECTRONS` (`conv_thr=1e-8`, `mixing_beta=0.7` with fixed occupations and `0.4` with smearing, `electron_maxstep=200`), `&IONS`/`&CELL` according to the preset (BFGS for relax, Verlet for MD with the requested thermostat, `press_conv_thr=0.05`), and the cards `ATOMIC_SPECIES` (masses from `ase.data.atomic_masses`), `ATOMIC_POSITIONS crystal`, `CELL_PARAMETERS angstrom`, the `HUBBARD` card if `--hubbard-style card`, and `K_POINTS`.
7. For bands, `kpoints.kpath_card` writes `K_POINTS crystal_b` with `band_points` points per segment (20 by default) and `KPATH.txt` with the labels; `build_bandsx_input` writes `bands_pp.in` (`lsym=.true.`). For DOS, `build_dos_input` and `build_projwfc_input` write `dos.in` and `projwfc.in` with `DeltaE = 0.02` eV.
8. `build_run_script` and `build_run_python_script` write `run.sh` (with `set -e -o pipefail` and `mpirun -np $NP`) and `run.py` in the order `pw.x → dos.x/projwfc.x → pw.x (bands) → bands.x`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Cell and positions | user file (CIF/POSCAR/…) | `structure.load` via `ase.io.read` |
| Primitive cell and k-path | seekpath library | `kpoints.get_kpath` with `symprec = 1e-4` Å (`structure.SYMPREC`) |
| Suggested cutoffs | UPF header (`wfc_cutoff`, `rho_cutoff`) | `pseudo.suggested_cutoffs`; reads only the first 20 000 characters |
| Valence electrons | `z_valence` in the UPF | `pseudo.z_valence`, used in `_estimate_nbnd` |
| Pseudopotential type and relativity | `pseudo_type`, `relativistic` in the UPF | `pseudo.pseudo_type`, `pseudo.relativistic` |
| `ecutwfc`, `dual`, `kspacing`, `degauss`, `smearing`, `nproc` | `~/.config/qekit/config.ini` or `config.DEFAULTS` | 60 Ry, 8, 0.20 Å⁻¹, 0.01 Ry, `cold`, 4 |
| Atomic masses | `ase.data.atomic_masses` | `ATOMIC_SPECIES` card |
| Hybrid parameters | table `inputgen.HIBRIDOS` | HSE: `exx_fraction 0.25`, `screening_parameter 0.106` bohr⁻¹; PBE0 0.25; B3LYP 0.20; Gau-PBE 0.24 |
| Hubbard orbital (card) | atomic number (`inputgen._orbital_hubbard`) | 3d (Z 21–30), 4d (39–48), 5d (72–80), 4f (57–71), 5f (89–103), `2p` otherwise |
| fs → Rydberg a.u. conversion | constant `inputgen._FS_POR_UA` | 4.8378e-2 fs |

**Limits and pitfalls.**

- The "automatic" cutoffs are those **suggested by the UPF**, not a convergence: the report says "(automático)". If the UPF does not declare them, the code falls back to 60 Ry / 480 Ry with no warning beyond the report.
- `--soc` writes `noncolin=.true.` and `lspinorb=.true.` only if every pseudopotential declares `relativistic="full"`: `inputgen.generate` calls `sweep.check_soc_pseudos` before writing anything and, if any of them is scalar-relativistic or does not declare it, it stops with "el acoplamiento espín-órbita necesita pseudopotenciales TOTALMENTE RELATIVISTAS (relativistic='full'), y estos no lo son". The reason, quoted in the error itself: with scalar pseudopotentials lspinorb "devuelve un desdoblamiento espín-órbita de cero que parece un resultado válido". `--soc` and `--nspin 2` are rejected together.
- `--hubbard-style legacy` (the default) writes `lda_plus_u` and `Hubbard_U(i)`, a syntax removed in QE ≥ 7.1; for those versions `--hubbard-style card` is needed, which writes the `HUBBARD (ortho-atomic)` card. The orbital in the card is deduced only from the atomic number.
- With hybrids, `nqx` **must divide** the k mesh; otherwise the command stops with "la malla de intercambio exacto tiene que DIVIDIR la de k". The report also warns that with `1x1x1` "el resultado va a salir claramente sobrestimado" and that pw.x cannot do a `calculation='bands'` with EXX.
- Without `--mag`, `--nspin 2` starts with zero magnetization and the report warns: "sin magnetización inicial el cálculo suele converger a la solución no magnética".
- The dipole correction requires a vacuum gap ≥ 5 Å; otherwise the error reads "la corrección dipolar necesita vacío en la dirección …".
- The scf mesh of a `bands` preset is computed on the seekpath primitive cell, which may not be the one the user supplied; the report says so ("AVISO: se usó la celda primitiva estandarizada").
- For MD, `nosym` is forced, and a warning is issued for fewer than 20 atoms or `dt > 2` fs.
- `tot_charge` is compensated by a uniform background; the report reminds that the energy of a charged cell is not comparable to the neutral one.
- The mesh is always uniform and Γ-centred (no shift), as the `kpoints.py` docstring states: it is not a shifted Monkhorst-Pack mesh, and with even $n$ it contains Γ where an MP mesh would not.

**References.**

- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* **128**, 140 (2017) — seekpath k-path convention. DOI 10.1016/j.commatsci.2016.10.015.
- A. Togo, I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv:1808.01590 (2018).
- N. Marzari, D. Vanderbilt, A. De Vita, M. C. Payne, *Phys. Rev. Lett.* **82**, 3296 (1999) — "cold" smearing. DOI 10.1103/PhysRevLett.82.3296.
- J. Heyd, G. E. Scuseria, M. Ernzerhof, *J. Chem. Phys.* **118**, 8207 (2003) — HSE.
- L. Bengtsson, *Phys. Rev. B* **59**, 12301 (1999) — dipole correction for slabs.
- G. Prandini et al., *npj Comput. Mater.* **4**, 72 (2018) — the SSSP pseudopotential library and cutoffs.

---

### `olla-dft kpath` — high-symmetry path

**What it answers.** Which standard path through the Brillouin zone (and in which cell it is expressed) must be used to draw the band structure of this structure.

**Background for non-experts.** The Brillouin zone is the "unit cell" of wave-vector space; bands $E(\mathbf k)$ are drawn along a route through its highest-symmetry points, because that is where bands touch, cross or have their extrema. Which points and in which order depends on the space group and the cell shape, and there are several incompatible conventions in the literature. Olla-DFT delegates the choice to seekpath (the convention of Hinuma et al., the same used by Materials Cloud), which also **standardizes the cell**: the coordinates of the special points are only valid in that primitive cell, not necessarily in the one in the user's CIF.

**Formulas.** No formulas of its own: the command calls `seekpath.get_path` and transcribes its result. The only arithmetic is the "same cell" criterion (`qekit/core/kpoints.py: get_kpath`):

$$
\text{cell\_changed} = \neg\left[N_{\text{prim}} = N_{\text{in}} \ \wedge\ \max_{ij} |A^{\text{prim}}_{ij} - A^{\text{in}}_{ij}| \le 10^{-5}\ \text{Å}\right]
$$

- $A^{\text{prim}}$, $A^{\text{in}}$: cell matrices (Å) of the seekpath primitive and of the input; $N$: number of atoms.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_kpath` loads the structure with `structure.load`.
2. `kpoints.get_kpath` converts to the spglib tuple (`structure.to_spglib_cell`), calls `seekpath.get_path(..., symprec=1e-4)` and rebuilds the primitive cell (`structure.from_spglib_cell`) from `primitive_lattice`, `primitive_positions`, `primitive_types`.
3. `kpoints.kpath_text` prints the space group (`spacegroup_international`, `spacegroup_number`), the compacted path (`Γ — X — U | K — Γ …`), the fractional coordinates of every point (`point_coords`) with the "pretty" labels of `pretty_label` (GAMMA→Γ, DELTA_0→Δ0), and a warning if the cell changed.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Path and coordinates | seekpath library (`get_path`) | keys `path`, `point_coords` |
| Space group | seekpath (spglib underneath) | `spacegroup_international`, `spacegroup_number` |
| Symmetry tolerance | constant `structure.SYMPREC` | 1e-4 Å |
| Primitive cell | seekpath | `primitive_lattice/positions/types` |

**Limits and pitfalls.**

- The coordinates are **in the standardized primitive cell**. If `cell_changed` is true, the text warns: "el k-path está referido a la celda primitiva estandarizada, que difiere de la celda de entrada. Usa esa celda primitiva en el cálculo de bandas". Using those coordinates in the original cell gives a wrong path with no error whatsoever.
- With `symprec = 1e-4` Å, a relaxed structure with numerical noise can lose symmetry and receive the path of a lower space group; there is no command-line option to change the tolerance.
- seekpath ignores magnetism and spin-orbit coupling when choosing the space group.

**References.**

- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Comput. Mater. Sci.* **128**, 140 (2017). DOI 10.1016/j.commatsci.2016.10.015.
- W. Setyawan, S. Curtarolo, *Comput. Mater. Sci.* **49**, 299 (2010) — the alternative convention that seekpath does **not** use.

---

### `olla-dft info` — structure and symmetry

**What it answers.** What the structure file contains: formula, lattice parameters, volume, space group, point group, Wyckoff positions and how many atoms the primitive cell would have.

**Background for non-experts.** Before computing anything it is worth knowing whether the cell is the smallest possible one (the primitive) or a larger one (conventional or supercell), because the cost of the calculation grows with the number of atoms, and whether the structure has the symmetry one believes. spglib compares every atom with the candidate symmetry operations within a tolerance and returns the space group in international notation (for example `Fd-3m`, No. 227 for silicon), the Hall symbol, the point group and the Wyckoff letter of each site (a label for the kind of symmetry position each atom occupies).

**Formulas.** No formulas beyond the cell geometry that ASE computes (`atoms.cell.cellpar()` returns $a, b, c, \alpha, \beta, \gamma$ and `atoms.cell.volume` the volume $V = |\det \mathbf A|$ in Å³).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_info` → `structure.load`.
2. `qekit/core/structure.py: info_text` calls `symmetry_dataset` (`spglib.get_symmetry_dataset` with `symprec = 1e-4`) and `primitive` (`spglib.standardize_cell(to_primitive=True)`) to count the atoms of the primitive cell.
3. It prints formula, composition, number of atoms, volume, lattice parameters, space group (`international`, `number`), Hall symbol, point group, atoms in the primitive cell, Wyckoff positions (sorted set of `ds.wyckoffs`) and the cell vectors.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Space group, Hall, point group, Wyckoff | spglib library | `structure.symmetry_dataset` |
| Lattice parameters and volume | ASE (`Cell.cellpar`, `Cell.volume`) | a, b, c in Å; angles in degrees |
| Atoms in the primitive cell | spglib `standardize_cell` | `structure.primitive` |
| Tolerance | `structure.SYMPREC` | 1e-4 Å |

**Limits and pitfalls.**

- If spglib cannot determine the symmetry, the command fails with `RuntimeError("spglib no pudo determinar la simetría de la estructura")`.
- The tolerance is fixed (1e-4 Å); there is no `--symprec`.
- Only the distinct Wyckoff letters are listed, not how many atoms sit in each.

**References.**

- A. Togo, I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv:1808.01590 (2018).
- International Tables for Crystallography, Vol. A (IUCr) — notation of space groups and Wyckoff positions.

---

### `olla-dft prim` — standardized primitive cell

**What it answers.** Which is the smallest cell that reproduces the crystal by translation, written in spglib's standard orientation.

**Background for non-experts.** Many CIFs come in the conventional cell (the one that makes the symmetry visible, e.g. the 8-atom cube of silicon), but the calculation only needs the primitive cell (2 atoms in silicon). Reducing it saves a factor equal to the ratio of atom counts in cost, without changing the physics. "Standardized" means that spglib reorients it and expresses it with the choice of vectors fixed by the International Tables convention, so that two equivalent inputs give the same output.

**Formulas.** No arithmetic of its own: it is `spglib.standardize_cell(cell, to_primitive=True, symprec=1e-4)`.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_prim` → `structure.load`.
2. `structure.primitive` → spglib → `structure.from_spglib_cell` (an `Atoms` with `pbc=True`).
3. `structure.convert` writes the result according to the extension of `-o` (default `primitive.cif`): CIF, POSCAR/`.vasp` (with `direct=True, sort=True`) or any format ASE can infer.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Primitive cell | spglib `standardize_cell(to_primitive=True)` | `structure.primitive` |
| Tolerance | `structure.SYMPREC` | 1e-4 Å |
| Output format | file extension | `structure.convert` |

**Limits and pitfalls.**

- spglib's primitive cell is **not** necessarily the same as the seekpath one used by `gen -p bands` (seekpath applies its own additional standardization); for bands, let `gen` choose.
- When writing POSCAR the atoms are **reordered** by species (`sort=True`); if the user relied on a specific order (for example for a `--displace`), it is lost.
- If spglib fails: `RuntimeError("spglib no pudo estandarizar la celda")`.

**References.**

- A. Togo, I. Tanaka, arXiv:1808.01590 (2018).

---

### `olla-dft conv` — standardized conventional cell

**What it answers.** Which is the conventional cell (the one showing the full symmetry of the crystal system) of the structure.

**Background for non-experts.** It is the inverse operation to `prim`: start from any cell and obtain the "textbook" cell (cubic for silicon, hexagonal for graphite), useful for building surfaces, supercells or for comparison with diffraction data, even though it has more atoms than strictly needed for the calculation.

**Formulas.** No arithmetic of its own: `spglib.standardize_cell(cell, to_primitive=False, symprec=1e-4)`.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_conv` → `structure.load`.
2. `structure.conventional` → spglib → `from_spglib_cell`.
3. `structure.convert` writes `-o` (default `conventional.cif`).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Conventional cell | spglib `standardize_cell(to_primitive=False)` | `structure.conventional` |
| Tolerance | `structure.SYMPREC` | 1e-4 Å |

**Limits and pitfalls.**

- Same as `prim`: fixed tolerance, atom reordering in POSCAR, error if spglib cannot standardize.
- For a low-symmetry (triclinic) cell the "conventional" cell coincides with the primitive one and the command changes nothing.

**References.**

- A. Togo, I. Tanaka, arXiv:1808.01590 (2018).

---

### `olla-dft supercell` — build a supercell

**What it answers.** The structure repeated $n_x \times n_y \times n_z$ times along its three cell vectors.

**Background for non-experts.** A supercell is several cells glued together and treated as a single periodic unit. It is needed to place a defect, a dopant or an adsorbed molecule at low concentration, for molecular dynamics, or to compute phonons by finite displacements. The price is that the Brillouin zone shrinks by the same factor and the bands "fold" (see `unfold`).

**Formulas.** ASE's `atoms.repeat((nx, ny, nz))`: the new cell is $\mathbf a'_i = n_i \mathbf a_i$ and every atom is copied into the $n_x n_y n_z$ translations $\sum_i m_i \mathbf a_i$ with $0 \le m_i < n_i$.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_supercell` → `structure.load`.
2. `structure.supercell` checks that the three factors are ≥ 1 (otherwise `ErrorDeUso("los factores de la supercelda deben ser >= 1")`) and calls `Atoms.repeat`.
3. `structure.convert` writes `-o` (default `supercell.cif`).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Factors $n_x, n_y, n_z$ | user parameters (positional) | integers ≥ 1 |
| Repetition | ASE `Atoms.repeat` | diagonal multiples only |

**Limits and pitfalls.**

- Only **diagonal** supercells (multiples of each vector); general matrices $\mathbf A' = \mathbf M \mathbf a$ such as those `unfold` can recognize cannot be built.
- No symmetry reduction and no check that the supercell is "reasonable" (e.g. cubic).

**References.**

- A. H. Larsen et al., "The atomic simulation environment—a Python library for working with atoms", *J. Phys.: Condens. Matter* **29**, 273002 (2017). DOI 10.1088/1361-648X/aa680e.

---
### `olla-dft bands` — band structure and band gap

**What it answers.** How the energy of every electronic state varies along the high-symmetry path, whether the system has a gap or is metallic, where the valence-band maximum (VBM) and the conduction-band minimum (CBM) are, whether the gap is direct or indirect, and — with `--fat` — which atomic orbital each band "is made of".

**Background for non-experts.** In a crystal electrons do not have discrete energies but *bands*: for every wave vector $\mathbf k$ (a "direction and wavelength" of the electron) there is a list of allowed energies $\varepsilon_n(\mathbf k)$. Drawing them along a path in the Brillouin zone gives the classic "spaghetti" plot. The *band gap* is the distance between the highest occupied band and the lowest empty one. If the maximum of one and the minimum of the other are at the same $\mathbf k$ the gap is *direct* (the material absorbs and emits light efficiently); otherwise it is *indirect*. If some band crosses the Fermi level (the energy up to which states are filled), it is a metal and there is no gap.

*Fat bands* answer a different question: what weight each atomic orbital (nickel $d$, oxygen $p$) has in each state. projwfc.x projects every wavefunction onto atomic orbitals and writes the weights; Olla-DFT draws them as dots whose size is proportional to the weight on top of the bands.

**Formulas.**

Unit conversion of the pw.x XML (`qekit/core/qeout.py: read_xml`):

$$
E_{\text{eV}} = E_{\text{Ha}} \cdot 27.211386245988, \qquad
\mathbf{k}_{\text{Å}^{-1}} = \mathbf{k}_{2\pi/a} \cdot \frac{2\pi}{a_{\text{bohr}}\cdot 0.529177210903}, \qquad
\mathbf{k}_{\text{frac}} = \mathbf{k}_{\text{Å}^{-1}}\, \mathbf{B}^{-1}
$$

- $a_{\text{bohr}}$: `alat` from the XML (bohr). $\mathbf B$: matrix with the reciprocal vectors $\mathbf b_i = 2\pi(\mathbf A^{-1})^{\mathsf T}_i$ as rows (Å⁻¹).

Cumulative distance on the x axis (`qekit/modules/bands.py: _build_kdist`):

$$
x_0 = 0, \qquad x_i = x_{i-1} + \begin{cases} 0 & i \in \text{breaks} \\ |\mathbf{k}_i - \mathbf{k}_{i-1}| & \text{otherwise} \end{cases}
$$

- `breaks`: indices where two special points appear consecutively (a `U|K` discontinuity of the path), detected by `_detect_breaks`.

Metal / insulator classification (`bands.analyze_gap`), with `CROSS_TOL = 1e-6` eV and reference $E_{\text{ref}}$:

$$
\text{crosses}_n = \left[\min_{\mathbf k}\varepsilon_n < E_{\text{ref}} - \delta\right] \wedge \left[\max_{\mathbf k}\varepsilon_n > E_{\text{ref}} + \delta\right]
$$

If any band crosses, it is a metal. Otherwise $n_v = \max\{n : \max_{\mathbf k}\varepsilon_n \le E_{\text{ref}}+\delta\}$ and $n_c = \min\{n : \min_{\mathbf k}\varepsilon_n > E_{\text{ref}}-\delta\}$ (forcing $n_c \ge n_v+1$), and

$$
E_{\text{VBM}} = \max_{\mathbf k}\varepsilon_{n_v}(\mathbf k), \quad
E_{\text{CBM}} = \min_{\mathbf k}\varepsilon_{n_c}(\mathbf k), \quad
E_g = E_{\text{CBM}} - E_{\text{VBM}}, \quad
E_g^{\text{dir}} = \min_{\mathbf k}\left[\varepsilon_{n_c}(\mathbf k) - \varepsilon_{n_v}(\mathbf k)\right]
$$

The gap is direct if $\arg\max\varepsilon_{n_v} = \arg\min\varepsilon_{n_c}$ (same k-point index, not same coordinate).

Reference $E_{\text{ref}}$ (in this order): `<fermi_energy>` from the XML; if absent, `<highestOccupiedLevel>`; if also absent and `nspin = 1`, the midpoint $\tfrac12[\max_{\mathbf k}\varepsilon_{n_{occ}-1} + \min_{\mathbf k}\varepsilon_{n_{occ}}]$ with $n_{occ} = \mathrm{round}(N_{\text{el}}/2)$; as a last resort, the median of all energies.

Energy zero in the plot and in the exported data (`bands.reference_energy`): `--ref auto` uses the VBM when there is a gap and $E_F$ for a metal; `fermi`, `vbm`, `none` as their names say.

Weight of a selector in the fat bands (`bands.peso_de`):

$$
w_{n}(\mathbf k) = \sum_{i \in \text{selector}} |\langle \phi_i | \psi_{n\mathbf k}\rangle|^2
$$

- $|\langle \phi_i | \psi_{n\mathbf k}\rangle|^2$: the coefficients `psi = 0.498*[# 1] + …` from the text output of projwfc.x for atomic state $i$, read by `bands.leer_proyecciones`. **Not normalized**: the part missing up to 1 is the fraction of the wavefunction that falls in no atomic sphere, which `report_fat` quantifies as $1 - \langle\sum_i w_i\rangle$.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_bands` → `bands.load(path, prefix)`.
2. `qeout.find_xml` locates the XML (`./out/*.xml`, `./*.xml` or `*.save/data-file-schema.xml`, checking that it contains "espresso"); `qeout.read_xml` reads `<atomic_structure>` (cell in bohr), `<band_structure>` (`nbnd`, `nelec`, `lsda`, `noncolin`, `fermi_energy`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`), and every `<ks_energies>` (`k_point` with `weight`, `eigenvalues`, `occupations`). With `lsda`, the list of eigenvalues per k is split into two halves (up/down).
3. Labels come from `KPATH.txt` (`qeout.read_kpath_labels`) or, if it does not exist, from the `K_POINTS crystal_b` card of `bands.in` with `! G` comments (`qeout.read_crystal_b_card`). `qeout.match_labels_to_kpoints` assigns them to indices with tolerance `1e-3` in fractional coordinates, always moving forward and tolerating reciprocal-lattice translations.
4. `bands.analyze_gap` per spin channel; `bands.gap_report` prints the summary and the reminder that "los funcionales GGA/LDA subestiman el gap sistemáticamente (típicamente 30–50 %)".
5. With `--fat`: `bands.leer_proyecciones` reads `projwfc.out` (or `proj.out`, `projwfc_bands.out`) from the same bands calculation; `comprobar_compatibilidad` requires the same number of k-points; `peso_de` sums the states that match the selector (`Ni`, `Ni-d`, `d`, `atomo:3`).
6. `bands.export` writes `BAND.dat` (or `BAND_up.dat`/`BAND_dw.dat`), `KLABELS.dat` and `BAND_GAP.txt`; `bands.plot` draws with matplotlib (bands in ink, spin ↓ dashed, VBM as a circle and CBM as a square, fat bands as `scatter` with `s = w · fat_scale`, `fat_scale = 55`).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues, k-points, weights | `prefix.xml` from pw.x (`<ks_energies>`) | `qeout.read_xml`; Ha → eV |
| Fermi energy | `<fermi_energy>` in `prefix.xml` | only if the scf used smearing |
| HOMO / LUMO | `<highestOccupiedLevel>` / `<lowestUnoccupiedLevel>` | fixed occupations |
| Number of electrons, bands, spin | `<nelec>`, `<nbnd>`, `<lsda>`, `<noncolin>` | `nbnd` is recomputed from the length of the eigenvalue list |
| Cell and `alat` | `<atomic_structure alat=…>` and `<cell>` | bohr → Å with 0.529177210903 |
| High-symmetry labels | `KPATH.txt` from `olla-dft gen` or `bands.in` | matching tolerance 1e-3 |
| Orbital weights (fat bands) | `projwfc.out` (`psi = …` blocks) | `bands.leer_proyecciones`; `state #` gives atom, element and $l$ |
| Hartree in eV, bohr in Å | constants `qeout.HARTREE_EV`, `qeout.BOHR_ANG` | 27.211386245988 eV; 0.529177210903 Å (CODATA 2018) |

**Limits and pitfalls.**

- In a `bands` calculation with fixed occupations the XML may not carry `<fermi_energy>`; the reference is then `<highestOccupiedLevel>`, which QE inherits from the scf. If the scf used smearing, $E_F$ may fall in the middle of the gap or inside a flat band; a band touching $E_F \pm 10^{-6}$ eV at a single point is classified as a metal.
- If `nbnd` only covers the occupied bands, the report says "No hay bandas de conducción en el cálculo (aumenta nbnd para obtener el gap)".
- With `nspin = 2` each channel is analysed separately with the **same** reference; the report does not compute the global gap between different channels (e.g. VBM up and CBM down).
- The plot and `--ref auto` always use the analysis of channel 0 (spin up) to decide the zero and mark the extrema.
- "Direct" is decided by comparing k-point **indices**; two symmetry-equivalent k-points at different indices count as indirect.
- Fat bands: if projwfc.x was run on the DOS nscf rather than on the bands, `comprobar_compatibilidad` stops: "las bandas tienen N puntos k y las proyecciones M. No son del mismo cálculo". If more than 10 % of the mean weight falls outside atomic spheres, `report_fat` warns: "De media, un X % de cada función de onda NO cae dentro de ninguna esfera atómica".
- Weights of states with $l>3$ are labelled `l4`, `l5`… and cannot be selected by letter.
- SOC (`noncolin`) calculations are read as a single channel; projwfc weights with $j$ (`p_j1.5`) are grouped only by the orbital letter.

**References.**

- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017) — Quantum ESPRESSO (XML format, projwfc.x). DOI 10.1088/1361-648X/aa8f79.
- J. P. Perdew, M. Levy, *Phys. Rev. Lett.* **51**, 1884 (1983) and L. J. Sham, M. Schlüter, *Phys. Rev. Lett.* **51**, 1888 (1983) — why DFT underestimates the gap.
- CODATA 2018, E. Tiesinga et al., *Rev. Mod. Phys.* **93**, 025010 (2021) — constants.

---

### `olla-dft gap` — band-gap report only

**What it answers.** The same as the analysis part of `bands` — metal or not, VBM, CBM, fundamental gap and minimum direct gap per spin channel — without exporting data or plotting.

**Background for non-experts.** It is the most frequent question asked of a band calculation ("how big is the gap?") decoupled from the figure. It works equally on an `scf`, an `nscf` mesh or a `bands` path: it reads any pw.x XML with eigenvalues. On a mesh, the gap obtained is that of the sampled points, which can be larger than the true one if the extremum falls between points.

**Formulas.** Exactly those of `bands.analyze_gap` described under `olla-dft bands` (classification with `CROSS_TOL = 1e-6` eV, $E_g = E_{\text{CBM}} - E_{\text{VBM}}$, $E_g^{\text{dir}} = \min_{\mathbf k}[\varepsilon_{n_c} - \varepsilon_{n_v}]$).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_gap` → `bands.load(path, prefix)` (reads the XML and, if present, `KPATH.txt` or `bands.in` to label the points).
2. `bands.gap_report` loops over `range(nspin)` calling `analyze_gap`, prints the XML path, calculation type, `nbnd`, `nk`, `nelec`, $E_F$ if present, and the result per channel.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues and k | `prefix.xml` from pw.x | `qeout.read_xml` |
| Reference | `<fermi_energy>` → `<highestOccupiedLevel>` → electron count → median | `bands.analyze_gap` |
| Label of the k-point of the extremum | `KPATH.txt` / `bands.in`, or the fractional coordinates | `bands._label_for` |

**Limits and pitfalls.**

- On a symmetry-reduced scf/nscf mesh, "direct" can only be detected if VBM and CBM fall at the **same index** of the irreducible-point list.
- An XML without `<output>` (unfinished calculation) yields `FaltanDatos("… no contiene una sección <output>")`.
- `gap_report` calls `analyze_gap` twice per channel (once for the report and once to decide whether to print the GGA reminder); this is only cost, it does not change the result.

**References.**

- The same as for `olla-dft bands`.

---

### `olla-dft dos` — total and projected density of states

**What it answers.** How many electronic states there are per unit energy (DOS), how they are distributed among elements and orbitals (PDOS), how large the DOS is at the Fermi level, and — with `--dband` — the centre, width and filling of a projected band (the "d-band centre" of catalysis).

**Background for non-experts.** The DOS is the histogram of the energies of all states in the Brillouin zone: where it is high there are many states, where it is zero there is a gap. dos.x computes it from the eigenvalues of the nscf (with the tetrahedron method that `gen -p dos` requests); projwfc.x decomposes every state into atomic orbitals and gives a PDOS per atom and orbital. Summing the files by element and by orbital letter ($s, p, d, f$) gives the chemical decomposition that gets published.

The *d-band centre* is the mean energy of the $d$ PDOS of a transition metal relative to the Fermi level. It is an empirical descriptor: the closer to the Fermi level (less negative), the more strongly the surface adsorbs (Hammer–Nørskov model).

**Formulas.**

File columns (`qekit/modules/dos.py: read_dos_file`, `read_pdos_file`): from `<prefix>.dos` one takes $E$, DOS (one or two columns depending on spin) and the integrated DOS; from `pdos_atm#N(El)_wfc#M(l)` one takes the `ldos` column (already summed over $m$), or `ldosup`/`ldosdw` with spin, according to the number of columns $1 + n_s(1 + (2l+1))$.

Aggregated PDOS (`dos.load`): $\rho_{\text{El},l}(E) = \sum_{\text{atoms } a \in \text{El}} \sum_{\text{wfc with } l} \text{ldos}_{a,l}(E)$; if dos.x and projwfc.x use different energy meshes, the projections are interpolated linearly (`np.interp`, zero outside the range) onto the total-DOS mesh. Without `<prefix>.dos`, the total DOS is defined as $\sum_{\text{El},l}\rho_{\text{El},l}$.

DOS at the Fermi level (`dos.report`): $\rho(E_F) = \sum_s \rho_s(E_{i^*})$ with $i^* = \arg\min_i |E_i - E_F|$; it is called "compatible with a gap" when $\rho(E_F) < 10^{-3}$ states/eV.

Moments of a projected band (`dos.momentos`), with $e = E - E_F$ and $\rho$ the PDOS of the selector (per spin channel, integrated with the trapezoidal rule `np.trapezoid`):

$$
N = \int \rho(e)\,de, \qquad
\varepsilon_c = \frac{1}{N}\int e\,\rho(e)\,de, \qquad
W = \sqrt{\frac{1}{N}\int (e-\varepsilon_c)^2\rho(e)\,de}, \qquad
f = \frac{1}{N}\int_{e \le 0}\rho(e)\,de
$$

- $N$: integrated states (dimensionless); $\varepsilon_c$: centre (eV relative to $E_F$); $W$: rms width (eV); $f$: filling (fraction).
- With two channels, the reported value is the average of each quantity weighted by $N_s$; the "exchange splitting" is $\varepsilon_c^{\uparrow} - \varepsilon_c^{\downarrow}$.
- Relative tail: $\max(\rho)$ over the last $\max(3, n/50)$ points divided by the global $\max(\rho)$; if it exceeds 0.05 the band is cut off at the top.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_dos` → `dos.load(path, prefix)`.
2. `dos.load` tries to read the XML with `qeout.read_xml` to get $E_F$ and the prefix; it looks for `<prefix>.dos` (or `*.dos`) and all `*pdos_atm#*`; it parses the name with the expression `pdos_atm#(\d+)\(([A-Za-z]+)\)_wfc#(\d+)\(([A-Za-z])…\)` to obtain element and orbital letter (also with SOC's `p_j1.5`).
3. If the XML gives no $E_F$, it takes it from the `EFermi = …` comment in the header of the `.dos`.
4. It orders the projections by element (order of appearance) and orbital $s,p,d,f$; `by_element` sums orbitals.
5. `dos.report` prints the energy range, $E_F$, the origin of the zero (`reference_energy`: Fermi unless `--ref none`), channels, projections and $\rho(E_F)$.
6. `dos.export` writes `DOS.dat` (E, DOS[_up/_dw], DOS_integrada) and `PDOS.dat` (per element_orbital and totals per element); `dos.plot` / `dos.draw` draw the total with a fill, the PDOS by `--mode orbital|element|total`, spin ↓ mirrored downwards.
7. With `--dband El[-orb]`: `dos.momentos` requires $E_F$ (otherwise `ErrorDeUso("no se encontró la energía de Fermi…")`) and the key (El, orb); `report_momentos` prints centre, width, filling and warnings.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Total, integrated DOS | `<prefix>.dos` from dos.x | columns E, dos[, dosup, dosdw], int |
| PDOS per atom/orbital | `<prefix>.pdos.pdos_atm#N(El)_wfc#M(l)` from projwfc.x | column `ldos` (or `ldosup`, `ldosdw`) |
| Fermi energy | `<fermi_energy>` in `prefix.xml`; else `EFermi` in the `.dos` header | `dos.load` |
| Energy mesh | that of dos.x (`DeltaE = 0.02` eV in `gen`) | PDOS are interpolated onto it |
| "Gap" threshold | constant in `dos.report` | 1e-3 states/eV |
| Tail threshold | constant in `dos.momentos`/`report_momentos` | 5 % of the peak |

**Limits and pitfalls.**

- A `pdos_atm#` file with a different number of points from the first one is skipped, and `dos.load` records it in `DOSData.avisos`; the report prints it: "se han SALTADO N archivo(s) de projwfc.x cuya malla de energía no coincide con la del primero (… puntos), así que la PDOS está incompleta… Casi siempre es que hay archivos de dos corridas de projwfc.x mezclados en la misma carpeta". The exported data still lack that orbital: move the old files away and reload.
- The total DOS defined as a sum of PDOS (when `.dos` is missing) omits the part of the wavefunctions that falls outside atomic spheres: it will be smaller than the real one.
- `momentos` integrates the **whole** available range except for `--dband-emax`; if the PDOS has not decayed, it warns: "al final del rango todavía queda un X % del pico de PDOS. La banda está CORTADA por arriba, así que el centro sale más bajo de lo que debería". The text recommends "Vuelve a correr projwfc.x con un Emax mayor".
- The d-band centre "es una correlación empírica dentro de una misma familia de metales, no una ley" (report text).
- With SOC, the orbital letter is taken from `p_j1.5` → `p`; the $j$ components are summed.
- `--ref vbm` does not exist for the DOS: `reference_energy` only distinguishes `none` from the rest (always Fermi).

**References.**

- P. E. Blöchl, O. Jepsen, O. K. Andersen, *Phys. Rev. B* **49**, 16223 (1994) — tetrahedron method (dos.x, `tetrahedra_opt`).
- B. Hammer, J. K. Nørskov, *Surf. Sci.* **343**, 211 (1995); *Adv. Catal.* **45**, 71 (2000) — d-band centre model.
- P. Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017) — projwfc.x.

---

### `olla-dft plot` — combined bands + DOS figure

**What it answers.** It produces the standard figure of an electronic-structure paper: bands on the left and the DOS rotated on the right, sharing the energy axis and the same zero.

**Background for non-experts.** Bands say *where* in k-space the states are; the DOS says *how many* there are at each energy. Placing them side by side with the same zero (the VBM if there is a gap, the Fermi level if metallic) lets one read at a glance which orbitals form each band.

**Formulas.** None of its own: `qekit/modules/combined.py: plot` only draws; the zero is taken from the band analysis (`bands.reference_energy`) and applied to both panels.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_plot` loads `bands.load` and `dos.load` on the same folder.
2. It prints `bands.gap_report`.
3. `combined.plot` creates two axes with ratio `ratio = 2.6`, calls the band-drawing logic and `dos.draw(vertical=True)` with the same energy shift.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Bands and gap | `prefix.xml` (see `bands`) | `bands.load` |
| DOS/PDOS | `.dos` and `pdos_atm#` (see `dos`) | `dos.load` |
| Energy zero | `bands.reference_energy(bs, ref)` | the same for both panels |

**Limits and pitfalls.**

- Bands and DOS usually come from different calculations (path vs. mesh) with **the same scf**; if they come from different scfs, their $E_F$ do not coincide and the right panel is shifted without any warning.
- The zero of the DOS is that of the bands (VBM in `auto` with a gap), whereas `olla-dft dos` on its own would use the Fermi level.

**References.**

- Those of `bands` and `dos`.

---
### `olla-dft effmass` — effective mass by parabolic fit

**What it answers.** How much an electron at the bottom of the conduction band and a hole at the top of the valence band "weigh": the effective mass $m^*/m_e$ in each direction, which governs mobilities, effective densities of states and exciton levels.

**Background for non-experts.** Near an extremum a band looks like a parabola, just like the kinetic energy of a free particle $E = \hbar^2 k^2/2m$. The curvature of that parabola defines a mass: a strongly curved ("open") band corresponds to a light, fast carrier; a flat band to a heavy one. The mass may depend on direction (in silicon the electron has a longitudinal mass of ~0.92 and a transverse one of ~0.19), so parabolas must be fitted along specific straight lines in k-space.

Olla-DFT works in two stages. First it fits on the bands you already have (fast, but with few points and only along the path directions). Then it writes a dedicated `bands` calculation with very fine lines crossing the VBM and the CBM in three directions (for a valley away from Γ: the radial Γ→k₀ direction, "longitudinal", and two perpendicular ones, "transverse"; at Γ: [100], [110], [111]) and, when it finishes, fits one parabola per line.

**Formulas.**

Quadratic fit and mass (`qekit/modules/effmass.py: from_bands`, `collect_fine`, `_mass_from_quadratic`):

$$
E(k) \approx a\,k^2 + b\,k + c, \qquad \frac{m^*}{m_e} = \frac{\hbar^2/m_e}{2a}, \qquad \frac{\hbar^2}{m_e} = 7.6199682\ \text{eV·Å}^2
$$

- $k$: signed distance to the extremum along the line (Å⁻¹); $a$ in eV·Å²; the fit is `np.polyfit(x, y, 2)`.
- The sign is kept: $a<0$ (downward curvature) gives $m^*<0$, which is the report's convention for a hole.
- The linear term $b$ is fitted but **not** used: the extremum is assumed at $k=0$.

Fit quality (`effmass._r2`): $R^2 = 1 - \sum(y - \hat y)^2 / \sum (y - \bar y)^2$.

Valley directions (`effmass.valley_directions`): if $|\mathbf k_0| < 10^{-6}$ Å⁻¹ (extremum at Γ), $\{[100], [110]/\sqrt2, [111]/\sqrt3\}$; otherwise $\hat e_1 = \mathbf k_0/|\mathbf k_0|$ and two perpendiculars built with cross products.

Points of the fine line (`effmass.prepare`): $\mathbf k_j = \mathbf k_0 + t_j\,\hat e$, $t_j \in [-h, h]$ with `half_width = 0.06` Å⁻¹ and `npts = 21` (forced odd), converted to fractional coordinates with $\mathbf k_{\text{frac}} = \mathbf k\,\mathbf B^{-1}$.

Identification of the valence band in the fine calculation (`collect_fine`): $n_v = \mathrm{round}(N_{\text{el}}/2) - 1$ (0-based index) if `nspin = 1`; otherwise the last band whose maximum lies below $E_F$ (or the HOMO, or the median).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_effmass` requires `--bands-dir` (a folder with a finished bands calculation) unless `--collect`. It loads the structure and `bands.load(bands_dir)`.
2. `effmass.from_bands`: `bands.analyze_gap` gives VBM/CBM; for holes it takes the VBM band and those lying within `DEGEN_TOL = 0.05` eV below it at that k; for electrons the CBM band and those degenerate above it. For each band it relocates the extremum (`argmax`/`argmin`), delimits the segment without crossing special points or discontinuities (`_segment_bounds`), collects the points with $|k - k_0| \le$ `--window` (a half-width; by default `WINDOW_DEFAULT = PARABOLIC_MAX/2 = 0.06` Å⁻¹, widening up to `--min-points` = 7) on both sides if the extremum is interior or on one side if it sits on a special point (`_collect_window`), and fits.
3. It warns if there are fewer than 5 points ("solo N puntos: el ajuste no es confiable; haz el cálculo dedicado (effmass sin --collect y luego --collect)") or if the **total fitted span** ($k_{\max} - k_{\min}$, `MassFit.window`) exceeds `PARABOLIC_MAX = 0.12` Å⁻¹ ("tramo ajustado de X Å⁻¹ (límite parabólico 0.12): el camino no tiene puntos más finos; haz el cálculo dedicado"). With the default half-width a fit centred on the extremum spans exactly 0.12, so the warning fires only when the window had to be widened for lack of points.
4. `effmass.prepare` reduces to the primitive cell (`structure.primitive`), resolves pseudopotentials and cutoffs with `sweep.prepare_common`, writes `masa.in` (`calculation='bands'`, `K_POINTS crystal` with the 6 lines, `nbnd` equal to that of the bands calculation, `occupations='fixed'` because `insulator=True`) and `scf.in` (mesh `sweep.default_grid`), and saves `masa_meta.json` with the description of every line.
5. With `--run`, `runner.run_all` launches pw.x on `scf.in` and `masa.in`; with `--collect`, `effmass.collect_fine` reads `out/*.xml`, slices the k list into chunks of `npts`, computes $t_j = \pm|\mathbf k_j - \mathbf k_c|$ and fits all bands degenerate with the extremum.
6. `effmass.report` prints the table (carrier, band, $m^*/m_e$, $R^2$, points, Δk, direction) and `export` writes `MASA_EFECTIVA.dat`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues and Cartesian k | `prefix.xml` from pw.x | `qeout.read_xml` (previous bands and fine calculation) |
| VBM, CBM and their k | `bands.analyze_gap` | see `olla-dft bands` |
| $\hbar^2/m_e$ | constant `effmass.HBAR2_OVER_ME` | 7.6199682 eV·Å² |
| Number of electrons | `<nelec>` from the XML | to identify the valence in `collect_fine` |
| Window, minimum points, half-width, points per line | user parameters | `--window` (default `effmass.WINDOW_DEFAULT` = 0.06 Å⁻¹), `--min-points 7`, `--half-width 0.06`, `--points 21` |
| Parabolic limit | `effmass.PARABOLIC_MAX` | 0.12 Å⁻¹ of total span (slack `_TOL_VENTANA = 1e-6`) |
| Degeneracy tolerance | `effmass.DEGEN_TOL` | 0.05 eV |
| Description of the lines | `masa_meta.json` written by `prepare` | `effmass.load_meta` |

**Limits and pitfalls.**

- There is no "quick fit only" mode: with `--bands-dir` the command always fits on the path and then prepares the fine calculation in `--outdir` (or runs it with `--run`); with `--collect` it reads the fine one. This is what the module docstring describes.
- For a metal the command stops: "El sistema es metálico: no hay un extremo de banda aislado que ajustar".
- The report warns that the calculation **does not include spin-orbit coupling**: "cerca de Γ hay un triplete degenerado, no el par hueco pesado / hueco ligero del modelo de Luttinger".
- An $R^2$ of 1.0000 with 3 or 4 points "no dice nada — una parábola pasa exacta por tres puntos cualesquiera" (report text).
- `--window` is a **half-width**: the fitted span is up to twice it, and it is that span that is compared with `PARABOLIC_MAX`. Raising `--window` above 0.06 Å⁻¹ triggers the non-parabolic warning even if the fit looks "good" in $R^2$.
- `--collect` uses the **first** XML in `out/`; with several prefixes it may read the wrong one.
- The fine calculation is always written with `occupations='fixed'` and without spin; for magnetic systems `masa.in` must be edited by hand.
- The transverse lines for a valley away from Γ are chosen with an arbitrary cross product: they are not necessarily crystallographic axes, and in a non-ellipsoidal valley the two transverse masses may differ.

**References.**

- N. W. Ashcroft, N. D. Mermin, *Solid State Physics* (1976), ch. 12 — definition of effective mass.
- J. M. Luttinger, W. Kohn, *Phys. Rev.* **97**, 869 (1955) — k·p model of degenerate valence bands.
- Silicon reference values: M. Cardona, F. H. Pollak, *Phys. Rev.* **142**, 530 (1966).

---

### `olla-dft fermi` — Fermi surface in BXSF format

**What it answers.** Which bands cross the Fermi level and what the surface $\varepsilon_n(\mathbf k) = E_F$ of each one looks like, written to a BXSF file that XCrySDen or FermiSurfer render in 3D.

**Background for non-experts.** In a metal the states fill up to an energy $E_F$; the set of k-points whose energy is exactly $E_F$ forms a surface in the Brillouin zone, the Fermi surface. Its shape determines conductivity, quantum oscillations and many instabilities (charge-density waves, superconductivity). To draw it one needs $\varepsilon_n(\mathbf k)$ on a **complete and uniform** mesh of the Brillouin zone, which is what `olla-dft transport` produces (an nscf with `nosym`, `noinv`).

**Formulas.**

Bands crossing $E_F$ (`qekit/modules/transport.py: crossing_bands`), with `tol = 1e-6` eV:

$$
\min_{\mathbf k}\varepsilon_n(\mathbf k) < E_F - \delta \quad\wedge\quad \max_{\mathbf k}\varepsilon_n(\mathbf k) > E_F + \delta
$$

Mesh reconstruction (`transport.load`): fractional coordinates are brought to $[0,1)$ with $f \leftarrow f - \lfloor f + 10^{-6}\rfloor$, rounded to 6 decimals, and $n_i$ is the number of distinct values along each axis; $n_1 n_2 n_3 = N_k$ is required.

BXSF grid (`transport.export_bxsf`): $(n_i + 1)$ points per axis are written, repeating the first plane at the end (`np.pad(..., mode="wrap")`), in C order (last index fastest), with the reciprocal vectors $\mathbf b_i = 2\pi(\mathbf A^{-1})^{\mathsf T}_i$ in Å⁻¹ and the energies in eV.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_fermi` looks for `out/*.xml` inside `--outdir` (default `transporte`).
2. `transport.load` reads the XML with `qeout.read_xml`; it rejects an `scf`-type XML ("es de un cálculo SCF, no del nscf de malla densa"); it reconstructs the mesh and reorders the energies with `np.lexsort`. It also computes band velocities by finite differences (unused here) and warns if the mesh is smaller than 24×24×24 or has fewer than 12 000 points.
3. `transport.crossing_bands` lists the metallic bands; if none, it prints "Ninguna banda cruza E_F: el sistema no es metálico y no tiene superficie de Fermi".
4. `transport.export_bxsf` writes `superficie_fermi.bxsf` with `Fermi Energy`, the grid and one `BAND:` block per band.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues on the mesh | `prefix.xml` of the `olla-dft transport` nscf | `qeout.read_xml`; spin channel 0 only |
| Fermi energy | `<fermi_energy>` from the XML | `run.fermi`; without it, `ErrorDeUso("no hay nivel de Fermi…")` |
| Cell | `<atomic_structure>` from the XML | reciprocal vectors with $2\pi$ |
| Crossing tolerance | argument `tol` of `crossing_bands` | 1e-6 eV |

**Limits and pitfalls.**

- It only works on the `olla-dft transport` folder (same full-mesh nscf); a band path or a symmetry-reduced mesh fails with "los N puntos k no forman una malla uniforme".
- Only spin channel 0 is exported (`transport.load(spin=0)`); a ferromagnetic metal would need two files and the command does not produce them.
- $E_F$ is taken as-is from the nscf XML; with `occupations='fixed'` it does not exist and the command fails.
- The Fermi level of a dense nscf is not recomputed: it is the one inherited from the scf (coarser mesh).
- The reciprocal vectors are written in Å⁻¹ with the factor $2\pi$; the viewer must interpret them in those units.

**References.**

- A. Kokalj, *Comput. Mater. Sci.* **28**, 155 (2003) — XCrySDen and the BXSF format. DOI 10.1016/S0927-0256(03)00104-6.
- M. Kawamura, *Comput. Phys. Commun.* **239**, 197 (2019) — FermiSurfer. DOI 10.1016/j.cpc.2019.01.017.

---

### `olla-dft unfold` — band unfolding of a supercell

**What it answers.** What fraction of every supercell state "belongs" to each k-point of the primitive cell: the spectral weight that lets one see the band of the original material (and how much a defect, a dopant or disorder blurs it) from a supercell calculation.

**Background for non-experts.** A supercell of $N$ primitive cells has a Brillouin zone $N$ times smaller, so its bands come out *folded*: where the primitive cell had one band, there are $N$ branches piled up. Every supercell state is a sum of plane waves $e^{i(\mathbf K + \mathbf G)\cdot\mathbf r}$, and every plane wave has a well-defined wave vector. Asking "how much of this state lives at point $\mathbf k$ of the primitive cell?" has an exact answer: the sum of $|C(\mathbf G)|^2$ over the plane waves whose $\mathbf K + \mathbf G$ coincides with $\mathbf k$ modulo the primitive reciprocal lattice. If the supercell is perfect, every state has weight 1 at a single $\mathbf k$ and the primitive band is recovered; if there is a defect, the weight spreads and the band looks blurred. That blurring is the physical result.

**Formulas.**

Supercell matrix (`qekit/modules/unfold.py: matriz_supercelda`): $\mathbf M = \mathbf A_{\text{sc}}\,\mathbf a_{\text{prim}}^{-1}$, rounded to integers; accepted if $\max|\mathbf M - \mathrm{round}(\mathbf M)| \le 10^{-3}$. If it fails because of orientation, `_m_por_metricas` searches for an integer $\mathbf M$ such that $\mathbf G_{\text{sc}} = \mathbf M\,\mathbf G_p\,\mathbf M^{\mathsf T}$ with $\mathbf G = \mathbf X\mathbf X^{\mathsf T}$ the metric tensor (rotation-invariant), row by row among integer vectors of the right length. The primitive cell is then **re-derived** as $\mathbf a = \mathbf M^{-1}\mathbf A_{\text{sc}}$ so that both share axes. $N = |\det\mathbf M|$.

Coordinates: since $\mathbf B_{\text{sc}} = \mathbf M^{-\mathsf T}\mathbf b_{\text{prim}}$, a vector with coordinates $\mathbf c_{\text{sc}}$ in the supercell reciprocal basis has coordinates $\mathbf c_p = \mathbf c_{\text{sc}}\mathbf M^{-\mathsf T}$ in the primitive one, and $\mathbf k_{\text{prim}} = \mathbf k_{\text{sc}}\mathbf M^{-\mathsf T}$ (`desdoblar`).

Spectral weight (`unfold.pesos_de_k`), with $\mathbf m_0 = \mathbf k_{\text{prim}}\mathbf M^{\mathsf T} - \mathbf k_{\text{sc}}$ (must be integer to `TOL_ENTERO = 1e-4`; otherwise the weight is 0 because that $\mathbf k$ does not fold onto this $\mathbf K$):

$$
P_{n}(\mathbf k) = \frac{\sum_{\mathbf G \in S}\ \sum_{\sigma}|C_{n\sigma}(\mathbf G)|^2}{\sum_{\mathbf G}\ \sum_{\sigma}|C_{n\sigma}(\mathbf G)|^2}, \qquad
S = \left\{\mathbf G : (\mathbf G - \mathbf m_0)\,\mathbf M^{-\mathsf T} \in \mathbb{Z}^3\right\}
$$

- $C_{n\sigma}(\mathbf G)$: plane-wave coefficients of band $n$ (spinor component $\sigma$ if `npol = 2`) read from `wfc<N>.dat`; $\mathbf G$ given by its Miller indices in the supercell reciprocal basis.
- The denominator normalizes in case the coefficients are not normalized; $P_n \in [0,1]$.

Distance on the x axis (`unfold._distancias`): sum of $|\Delta\mathbf k|$ with $\mathbf k = \mathbf k_{\text{frac}}\,\mathbf b_{\text{prim}}$; a jump larger than 5 times the median of the non-zero steps counts as zero (branch change).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_unfold` loads the primitive structure and calls `unfold.desdoblar(path, primitive_cell, bandas=range(--bands), spin=--spin)`; `--spin` is `up` (default) or `dw`.
2. `desdoblar` reads the XML (`qeout.read_xml`), locates the `.save` folder (`_carpeta_save`) and the wavefunction files of the requested channel (`qekit/core/wfc.py: buscar_wfc(save, spin)`, sorted by k number): if `wfc.es_lsda` detects `wfcup*`/`wfcdw*`, it returns only the `wfc{up|dw}<N>.dat` of that channel; otherwise the `wfc<N>.dat` of a spin-unpolarized run. Without them: "El cálculo no guardó las funciones de onda: eso pasa con disk_io='nowf' o 'low'" (or, with lsda, "falta el canal '…'").
3. `matriz_supercelda` obtains $\mathbf M$; the k-points of the calculation are converted to primitive coordinates.
4. For every k-point, `wfc.leer_wfc` reads the unformatted Fortran file: record 1 (`ik`, `xk`, `ispin`, `gamma_only`, `scalef`), record 2 (`ngw`, `igwx`, `npol`, `nbnd`), record 3 (`b1,b2,b3`), record 4 (Miller indices) and one record per band with `npol·igwx` complex numbers (only the requested bands are materialized).
5. `pesos_de_k` computes $P_n(\mathbf k)$; the energies come from the XML, from the same spin channel as the wavefunctions (`res.eigenvalues[0]` for `up`, `[1]` for `dw`).
6. `unfold.report` prints $N$, $\mathbf M$, the weight distribution (mean, fraction > 0.9, fraction < 0.1) and warnings; `export` writes `UNFOLD.dat` (distance, $E - E_F$, weight) and `UNFOLD.txt`; `plot` draws a `scatter` with size $= 60\,P$ for weights > 0.005.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Coefficients $C(\mathbf G)$, Miller indices, `npol` | `out/<prefix>.save/wfc<N>.dat` from pw.x | `wfc.leer_wfc` (sequential Fortran format, little-endian) |
| Eigenvalues, fractional k, $E_F$, supercell cell | `prefix.xml` | `qeout.read_xml` |
| Spin channel | `--spin up|dw` (user) | `wfc.buscar_wfc`, `wfc.es_lsda`; in a spin-unpolarized run it changes nothing |
| Primitive cell | user file | re-derived as $\mathbf M^{-1}\mathbf A_{\text{sc}}$ |
| Integer tolerance | `unfold.TOL_ENTERO` | 1e-4 (and 1e-3 to accept $\mathbf M$) |

**Limits and pitfalls.**

- `disk_io='medium'` or `'high'` is needed in the supercell bands calculation; `olla-dft gen` does not set it by default.
- **One spin channel per run** is unfolded. For an `lsda` calculation the report warns: "el cálculo es de espín polarizado (lsda) y aquí solo se ha desdoblado el canal 'up' (wfcup<N>.dat y sus energías). El otro canal no se mezcla ni se suma: para verlo repite el desdoblamiento con --spin dw". The two channels are never combined in a single figure.
- If the supercell is relaxed and the primitive is not (or vice versa), $\mathbf M$ is not integer: "la celda de la supercelda no es un múltiplo entero de la primitiva (error …)".
- If almost all weights are 1, the report warns: "la supercelda parece PERFECTA (sin defecto ni desorden). En ese caso el desdoblamiento reproduce exactamente las bandas primitivas — que es la comprobación de que funciona, pero no un resultado nuevo".
- The k-points of the calculation are interpreted as supercell k-points and converted to the primitive; no primitive path is generated and no check is made that the supercell k-points are the correct folds of the desired path.
- Ultrasoft/PAW wavefunctions are handled only through their plane-wave part (the augmentation term $S$ is not included); the weight is that of the smooth part.

**References.**

- V. Popescu, A. Zunger, *Phys. Rev. B* **85**, 085201 (2012) — unfolding spectral weight. DOI 10.1103/PhysRevB.85.085201.
- P. V. C. Medeiros, S. Stafström, J. Björk, *Phys. Rev. B* **89**, 041407(R) (2014) — plane-wave unfolding (BandUP). DOI 10.1103/PhysRevB.89.041407.
- W. Ku, T. Berlijn, C.-C. Lee, *Phys. Rev. Lett.* **104**, 216401 (2010).

---
### `olla-dft wannier` — Wannier functions and band interpolation

**What it answers.** From a DFT calculation on a coarse k-point mesh it builds a small model $H_{mn}(\mathbf R)$ in a basis of localized functions (Wannier functions) with which the bands can be evaluated at **any** k-point without re-running pw.x; it also gives where each function is centred, how far it extends (its spread $\Omega$), and how closely the interpolated band matches the DFT one.

**Background for non-experts.** Bloch states $\psi_{n\mathbf k}$ are delocalized over the whole crystal. Their Fourier transform in k gives functions $|\mathbf R n\rangle$ localized around a cell $\mathbf R$: the Wannier functions. In that basis the Hamiltonian is a small matrix $H(\mathbf R)$ that decays with $|\mathbf R|$, and transforming back to k gives the band at any point (a "Fourier interpolation" that is exact at the starting points). The difficulty is that every $\psi_{n\mathbf k}$ is defined only up to a phase (and, with degenerate bands, up to a unitary rotation among them): that freedom is called the *gauge*. With an arbitrary gauge the Wannier functions are not localized and the interpolation is garbage. Marzari and Vanderbilt proposed choosing the gauge that minimizes the total spread $\Omega$ (the sum of the squared "widths"); a good starting point is to project onto trial atomic orbitals and orthonormalize.

When the bands of interest cross others (metals, conduction bands) there is no isolated group to transform: one must choose at every k a subspace of $J$ states that "connects smoothly" with that of its neighbours. That is the *disentanglement* of Souza, Marzari and Vanderbilt, with an *outer* window (where one may choose from) and optionally a *frozen* one (states that must be kept exactly). Olla-DFT implements both in Python, using only the overlaps and projections computed by `pw2wannier90.x` (shipped with QE), without needing wannier90; if the user has wannier90, it also reads its `seedname_hr.dat`.

**Formulas.**

Full mesh and $\mathbf b$ vectors (`qekit/modules/wannier.py: malla_completa`, `capas_b`, `residuo_completitud`): $\mathbf k_{ijk} = (i/n_1, j/n_2, k/n_3)$ in QE's order (last index fastest). Neighbour shells $\mathbf b = (h_1/n_1, h_2/n_2, h_3/n_3)\,\mathbf B$ are added by distance until, by least squares over the 6 independent components,

$$
\sum_{\mathbf b} w_{\mathbf b}\, b_\alpha b_\beta = \delta_{\alpha\beta}, \qquad \text{residual } = \left\|\textstyle\sum_{\mathbf b} w_{\mathbf b}\,\mathbf b\otimes\mathbf b - \mathbf 1\right\|_\infty < 10^{-5}
$$

- $w_{\mathbf b}$: weight of each shell (Å²); shells that add no rank (SVD) or with $|w| < 10^{-8}$ are discarded.

Projection gauge (`wannier.gauge_proyeccion`), with $A_{mn}(\mathbf k) = \langle\psi_{m\mathbf k}|g_n\rangle$ from the `.amn` and the SVD $A = u\,s\,v^\dagger$:

$$
U(\mathbf k) = A\,(A^\dagger A)^{-1/2} = u\,v^\dagger
$$

- $U$: $N_b\times J$ matrix with orthonormal columns (Löwdin); the smallest singular value $s_{\min}$ is reported (warning if $< 0.2$).

Gauge-invariant spread and disentanglement (`wannier.omega_I`, `gauge_desenredo`), with $M^{\mathbf k,\mathbf b}_{mn} = \langle u_{m\mathbf k}|u_{n,\mathbf k+\mathbf b}\rangle$ from the `.mmn`:

$$
\Omega_I = \frac{1}{N_k}\sum_{\mathbf k}\sum_{\mathbf b} w_{\mathbf b}\left[J - \left\|U^\dagger(\mathbf k)\,M^{\mathbf k,\mathbf b}\,U(\mathbf k+\mathbf b)\right\|_F^2\right]
$$

$$
Z(\mathbf k) = \sum_{\mathbf b} w_{\mathbf b}\, M^{\mathbf k,\mathbf b}\,U(\mathbf k+\mathbf b)\,U^\dagger(\mathbf k+\mathbf b)\,M^{\mathbf k,\mathbf b\,\dagger}
$$

- At every iteration $Z$ is mixed with the previous one ($Z \leftarrow \mu Z_{\text{new}} + (1-\mu)Z_{\text{old}}$, $\mu = 0.5$ initially, halved if $\Omega_I$ goes up), restricted to the bands in the outer window, the frozen ones are projected out ($Q Z Q$ with $Q = 1 - P_{\text{frozen}}$) and the $J - N_{\text{frozen}}$ eigenvectors of largest eigenvalue are taken. At most 200 steps, tolerance $10^{-10}$ Å². At the end it re-projects onto the trial orbitals inside the subspace ($U \leftarrow U\,\mathrm{polar}(U^\dagger A)$) to obtain a smooth starting gauge.

Real-space Hamiltonian and interpolation (`wannier.hamiltoniano_k`, `a_reales`, `interpolar`, `celda_wigner_seitz`):

$$
H(\mathbf k) = U^\dagger(\mathbf k)\,\mathrm{diag}\big(\varepsilon_n(\mathbf k)\big)\,U(\mathbf k), \qquad
H(\mathbf R) = \frac{1}{N_k}\sum_{\mathbf k} e^{-2\pi i\,\mathbf k\cdot\mathbf R}\,H(\mathbf k), \qquad
H^{\text{int}}(\mathbf k) = \sum_{\mathbf R}\frac{e^{2\pi i\,\mathbf k\cdot\mathbf R}}{\deg(\mathbf R)}\,H(\mathbf R)
$$

- $\mathbf k$ and $\mathbf R$ in fractional coordinates (hence the explicit $2\pi$). $\mathbf R$ runs over the vectors of the Wigner-Seitz cell of the $n_1\times n_2\times n_3$ superlattice; $\deg(\mathbf R)$ is the number of equidistant images (tolerance $10^{-5}$ Å²). The bands are the eigenvalues of $\tfrac12(H^{\text{int}} + H^{\text{int}\dagger})$.

Centres and spread (`wannier.dispersion`), equations 31 and 34–36 of Marzari-Vanderbilt, with $M^W = U^\dagger(\mathbf k) M^{\mathbf k,\mathbf b} U(\mathbf k+\mathbf b)$ and $\phi_n = \operatorname{Im}\ln M^W_{nn}$:

$$
\bar{\mathbf r}_n = -\frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\,\mathbf b\,\phi_n, \qquad
\Omega_n = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\left[\left(1-|M^W_{nn}|^2\right) + \phi_n^2\right] - |\bar{\mathbf r}_n|^2
$$

$$
\Omega_I = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\Big[J - \sum_{mn}|M^W_{mn}|^2\Big], \quad
\Omega_{OD} = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\sum_{m\ne n}|M^W_{mn}|^2, \quad
\Omega_D = \frac{1}{N_k}\sum_{\mathbf k,\mathbf b} w_{\mathbf b}\sum_n\left(\phi_n + \mathbf b\cdot\bar{\mathbf r}_n\right)^2
$$

- $\bar{\mathbf r}_n$ in Å, $\Omega$ in Å²; $\Omega = \sum_n\Omega_n = \Omega_I + \Omega_D + \Omega_{OD}$ (the report prints the sum as a check).

Minimization (`wannier._gradiente`, `_rotar`, `minimizar`), eqs. 52–57 of Marzari-Vanderbilt, with $R_{mn} = M_{mn}M^*_{nn}$, $T_{mn} = (M_{mn}/M_{nn})\,q_n$, $q_n = \phi_n + \mathbf b\cdot\bar{\mathbf r}_n$, $\mathcal A(B) = (B - B^\dagger)/2$, $\mathcal S(B) = (B + B^\dagger)/2i$:

$$
G(\mathbf k) = -\frac{4}{N_k}\sum_{\mathbf b} w_{\mathbf b}\left[\mathcal A(R^{\mathbf k,\mathbf b}) - \mathcal S(T^{\mathbf k,\mathbf b})\right], \qquad
U(\mathbf k) \leftarrow U(\mathbf k)\,\exp\!\left(-\Delta t\,G(\mathbf k)\right), \qquad
\Delta t_0 = \frac{\alpha}{4\sum_{\mathbf b} w_{\mathbf b}},\ \alpha = 2
$$

- If the step raises $\Omega$ it is halved up to 12 times; at most 500 steps (`--iterations`); stop when $|\Delta\Omega| < 10^{-10}$. It checks that $\Omega_I$ does not change (`deriva_I`).

Interpolated DOS (`wannier.dos_interpolada`): $\rho(E) = \frac{1}{N_k\,\sigma\sqrt{2\pi}}\sum_{\mathbf k,n} e^{-(E-\varepsilon_n(\mathbf k))^2/2\sigma^2}$ on an $N^3$ mesh (`--dos N`), $\sigma$ = `--sigma` 0.05 eV; it integrates to $J$ states per cell, **without** the spin factor 2. The header of `WANNIER_dos.dat` declares it via `wannier.DOS_UNIDADES`: "estados/eV/celda, sin factor de espín: integra a num_wann (x2 para comparar con dos.x sin espín)".

**How Olla-DFT computes it.**

1. *Prepare* (`qekit/cli.py: _cmd_wannier` → `wannier.prepare`): it translates `--projections` (`Si:sp3`, `O:p;Ti:d`, `f=0.25,0.25,0.25:s`, or `auto` = $s$ and $p$ on every atom) into $(l, m_r)$ orbitals from the `ORBITALES` table (wannier90 convention); it writes `1_scf.in`, `2_nscf.in` (full `--grid` mesh, 4×4×4 by default, `K_POINTS crystal`, `nosym`, `noinv`, `conv_thr 1e-10`, `nbnd = --bands` or $J$ + excluded), the `<prefix>.nnkp` (`escribir_nnkp`: real and reciprocal lattice, k-points, projections, `nnkpts` neighbours with their $\mathbf G$, `exclude_bands`), `3_pw2wan.in` (`write_amn`, `write_mmn`, `write_unk=.false.`), `<prefix>.win` (in case one prefers wannier90) and `4_bands.in` (DFT bands along the seekpath path, 30 points per segment, `outdir='./out_bandas'`).
2. *Run* (`--run` → `wannier.correr`): `pw.x` (scf), `pw.x` (nscf), `pw2wannier90.x`, `pw.x` (bands), in that order, stopping at the first failure.
3. *Collect* (`--collect` → `wannier.collect`): reads `.eig` (`leer_eig`), `.amn` (`leer_amn`), `.mmn` (`leer_mmn`, with $m$ running fastest → `reshape(order="F")`); recomputes shells and neighbours from the `.nnkp` (`_leer_nnkp`); if there are more bands than functions or `--window`/`--frozen` were given, `gauge_desenredo`; otherwise `gauge_proyeccion`. `dispersion` before and after `minimizar` (unless `--no-minimize`). `celda_wigner_seitz`, `hamiltoniano_k`, `a_reales`; it checks that `interpolar` reproduces the mesh (`error_malla` < `TOL_EXACTA = 1e-6` eV) and, if DFT bands exist (`out_bandas`, `--dft-bands`), compares at points that were not in the mesh. As a negative control, it repeats the interpolation with $U = 1$ (`E_sin_gauge`).
4. If a wannier90 `*_hr.dat` exists (other than Olla-DFT's own `WANNIER_hr.dat`), `leer_hr` uses it directly and skips the localization.
5. `wannier.report` prints mesh, neighbours and residual, windows, $\Omega_I$, the decomposed $\Omega$, centres with assignment to atom or bond (`asignar`, bond window 0.5–3.2 Å), decay of $H(\mathbf R)$, exactness on the mesh and error against DFT; `export` writes `WANNIER_hr.dat` (wannier90 format), `WANNIER_centros.dat`, `WANNIER_bandas.dat`, `WANNIER.txt` and optionally `WANNIER_dos.dat`; `plot` draws Wannier over DFT and the $\Omega$ trace.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Energies $\varepsilon_n(\mathbf k)$ | `seedname.eig` from pw2wannier90.x | absolute eV; `wannier.leer_eig` |
| Projections $A_{mn}(\mathbf k)$ | `seedname.amn` | `wannier.leer_amn` |
| Overlaps $M^{\mathbf k,\mathbf b}_{mn}$ | `seedname.mmn` | `wannier.leer_mmn` |
| Cell, mesh, excluded bands | `seedname.nnkp` (written by Olla-DFT) | `wannier._leer_nnkp` |
| External $H(\mathbf R)$ | `seedname_hr.dat` from wannier90 | `wannier.leer_hr` |
| DFT validation bands | `out_bandas/*.xml` from step 4 | `qeout.read_xml`, channel 0 |
| Trial orbitals $(l, m_r)$ | table `wannier.ORBITALES` | Table 3.1/3.2 of the wannier90 manual |
| Tolerances | `TOL_COMPLETITUD 1e-5`, `TOL_PESO 1e-8`, `TOL_EXACTA 1e-6` eV | module constants |
| High-symmetry path | seekpath (`wannier.camino_denso`) | 30 points per segment (`--points`) |

**Limits and pitfalls.**

- The `--window` and `--frozen` windows are compared with the **absolute** energies of the `.eig` (not relative to $E_F$), as in wannier90.
- With disentanglement and no frozen window, nothing has to be reproduced exactly; the report warns: "Sin ventana congelada no hay ninguna banda que la interpolación tenga que reproducir exactamente… Si quieres que la valencia salga exacta, pásala en --frozen".
- The `auto` projections ($s$ and $p$ per atom) fail for transition metals (the $d$ are missing) and for strongly covalent bonds; the report always says so.
- If $H(\mathbf R)$ at the edge of the superlattice exceeds 5 % of $H(0)$: "H(R) apenas ha decaído al borde de la superred: la base no está localizada".
- Without `--dft-bands` (or `out_bandas`) the report warns: "No has comparado con bandas de DFT. Que la interpolación reproduzca la malla es trivial".
- Only spin channel 0 of the DFT bands is read; the workflow is not designed for `nspin = 2` or SOC (pw2wannier90 supports them, but `prepare` does not write `nspin`).
- The minimization is gradient descent with line search, not wannier90's conjugate gradient: it may need more steps and may stop in a local minimum.
- Very anisotropic meshes may admit no shells satisfying completeness: "no encuentro un conjunto de capas de vecinos que cumpla la condición de completitud con esta malla".
- `--collect` without the structure as first argument fails: "para analizar hace falta la estructura".
- The interpolated DOS carries no spin factor 2 (the file header says so and asks to multiply by 2 to compare with spin-unpolarized dos.x) and is only valid within the energy range covered by the Wannier functions.

**References.**

- N. Marzari, D. Vanderbilt, *Phys. Rev. B* **56**, 12847 (1997) — maximally localized Wannier functions. DOI 10.1103/PhysRevB.56.12847.
- I. Souza, N. Marzari, D. Vanderbilt, *Phys. Rev. B* **65**, 035109 (2001) — disentanglement. DOI 10.1103/PhysRevB.65.035109.
- N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, D. Vanderbilt, *Rev. Mod. Phys.* **84**, 1419 (2012) — review. DOI 10.1103/RevModPhys.84.1419.
- G. Pizzi et al., *J. Phys.: Condens. Matter* **32**, 165902 (2020) — Wannier90 v3 (`.nnkp`, `.amn`, `.mmn`, `_hr.dat` formats). DOI 10.1088/1361-648X/ab51ff.
- P.-O. Löwdin, *J. Chem. Phys.* **18**, 365 (1950) — symmetric orthonormalization.

---

### `olla-dft topology` — Chern number and Wilson loops

**What it answers.** Whether the occupied subspace of a Wannier model, on a two-dimensional section of the Brillouin zone, has a non-zero Chern number (an integer topological invariant), and how the hybrid Wannier centres (Wilson loops) evolve across that section.

**Background for non-experts.** Besides their energies, bands have a "geometry": when traversing a closed loop in k-space, the occupied states accumulate a phase (Berry phase) that does not depend on how the phases of each state are chosen. Summing that phase over a whole 2D section of the Brillouin zone yields an integer, the Chern number, which does not change under smooth deformations of the system: it is *topological*. A non-zero Chern number implies dissipationless edge currents (quantum anomalous Hall effect). The *Wilson loop* is the "slice by slice" version: for every $k_2$ one computes the product of overlaps along $k_1$; the phases of its eigenvalues are the positions (modulo 1) of the hybrid Wannier functions, and their "pumping" as $k_2$ varies is another way of seeing the Chern number.

**Formulas.**

Mesh and states (`qekit/modules/topology.py: kmesh`, `analyze`): $\mathbf k_{ij}$ with $k_a = i/n_1$, $k_b = j/n_2$ and the third coordinate fixed at `--fixed` (mod 1), in the `--plane` (`xy`, `xz`, `yz`); the eigenvectors $|u_n(\mathbf k)\rangle$ come from `wannier.interpolar(..., vectores=True)`.

Unitary links and discrete Berry curvature (`topology._unitary_overlap`, `invariants_from_vectors`), with $V(\mathbf k)$ the $N_w\times N_{\text{occ}}$ matrix of occupied eigenvectors:

$$
O_\mu(\mathbf k) = V^\dagger(\mathbf k)\,V(\mathbf k+\hat\mu), \qquad
Q_\mu = u\,v^\dagger \ \text{(unitary part of } O_\mu = u\,s\,v^\dagger), \qquad
U_\mu(\mathbf k) = \frac{\det Q_\mu(\mathbf k)}{|\det Q_\mu(\mathbf k)|}
$$

$$
F_{12}(\mathbf k) = \arg\!\left[U_1(\mathbf k)\,U_2(\mathbf k+\hat 1)\,U_1^*(\mathbf k+\hat 2)\,U_2^*(\mathbf k)\right], \qquad
C = \frac{1}{2\pi}\sum_{\mathbf k} F_{12}(\mathbf k)
$$

- $\hat\mu$: mesh step in direction $\mu$ (periodic). $F_{12} \in (-\pi, \pi]$ per plaquette (rad); $C$ is rounded to the nearest integer and the residual $|C - \mathrm{round}(C)|$ is reported.
- The smallest singular value of all $O_\mu$ (`min_overlap`) is also reported; if $< 10^{-6}$, a warning about a too-coarse mesh.

Wilson loops (`invariants_from_vectors`):

$$
W(k_2) = \prod_{i=0}^{n_1-1} Q_1(k_1^{(i)}, k_2), \qquad
x_n(k_2) = \frac{\arg\lambda_n\!\left[W(k_2)\right]}{2\pi} \bmod 1
$$

- $x_n$: sorted hybrid Wannier centres, in fractions of the lattice vector along direction 1.

Section gaps: $E_g^{\text{dir}} = \min_{\mathbf k}[\varepsilon_{N_{\text{occ}}+1} - \varepsilon_{N_{\text{occ}}}]$, $E_g^{\text{ind}} = \min_{\mathbf k}\varepsilon_{N_{\text{occ}}+1} - \max_{\mathbf k}\varepsilon_{N_{\text{occ}}}$. $E_g^{\text{dir}} > $ `--gap-tol` (1e-8 eV) is required.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_topology` requires exactly one of `--occupied N` or `--fermi EV`, and a `--grid` of at least 3×3 (40×40 by default).
2. `topology.resolve_model` accepts a `*_hr.dat` or a folder containing `WANNIER_hr.dat` (or a single `*_hr.dat`; with several, error "indica el archivo exacto").
3. `wannier.leer_hr` reads $H(\mathbf R)$, $\mathbf R$ and degeneracies; `wannier.interpolar` diagonalizes on the section mesh.
4. With `--fermi`, it counts the states with $\varepsilon < E_F$ at every k; if the count varies, error: "el nivel de Fermi corta bandas… El sistema es metálico en esta sección y el Chern de 'las ocupadas' no está definido".
5. `invariants_from_vectors` computes curvature, Chern and Wilson loops.
6. `topology.report` prints gaps, discrete and integer Chern, residual and minimum overlap; `export` writes `TOPOLOGY_curvature.dat` (flux per plaquette), `TOPOLOGY_wilson.dat` (centres vs. $k_2$) and `TOPOLOGY.txt`; `plot` draws the flux map and the centres.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $H(\mathbf R)$, $\mathbf R$, $\deg(\mathbf R)$ | `WANNIER_hr.dat` (Olla-DFT) or `seedname_hr.dat` (wannier90) | `wannier.leer_hr` |
| Eigenvectors on the mesh | `wannier.interpolar` | fractional coordinates, phase $e^{2\pi i\mathbf k\cdot\mathbf R}$ |
| Occupation | `--occupied` or `--fermi` (user) | never guessed |
| Gap tolerance | `--gap-tol` | 1e-8 eV |

**Limits and pitfalls.**

- 2D sections only: for a 3D material `--fixed` must be scanned by hand; the Chern number of a section is the invariant of a 2D Chern insulator or of a slice.
- "La señal cambia al invertir la orientación del plano" (report text): the sign of $C$ depends on the `(a, b)` order of the chosen plane.
- $\mathbb Z_2$ is not computed: "no se asigna un Z2 automático sin comprobar simetría de reversión temporal". With time-reversal symmetry the Chern number is always 0; the exported Wilson loops let one read $\mathbb Z_2$ by eye, but the code does not do it.
- If the direct gap closes on the mesh: "el subespacio ocupado no está aislado… El número de Chern no está definido".
- If the discrete Chern does not close to an integer within $10^{-6}$: "refina la malla y revisa la localización del modelo Wannier".
- The result inherits every defect of the Wannier model (bad projections, undecayed $H(\mathbf R)$).

**References.**

- T. Fukui, Y. Hatsugai, H. Suzuki, *J. Phys. Soc. Jpn.* **74**, 1674 (2005) — discrete Chern number on a mesh. DOI 10.1143/JPSJ.74.1674.
- R. Yu, X. L. Qi, A. Bernevig, Z. Fang, X. Dai, *Phys. Rev. B* **84**, 075119 (2011) — Wilson loops and hybrid centres.
- A. A. Soluyanov, D. Vanderbilt, *Phys. Rev. B* **83**, 235401 (2011) — hybrid Wannier centres and invariants.
- D. Vanderbilt, *Berry Phases in Electronic Structure Theory* (Cambridge, 2018).
- X.-L. Qi, Y.-S. Wu, S.-C. Zhang, *Phys. Rev. B* **74**, 085308 (2006) — test model used in `tests/test_topology.py`.

---
### `olla-dft berry` — Berry-phase polarization, Born charges

**What it answers.** How much the electric polarization of an insulating crystal changes when going from a reference structure (usually the centrosymmetric one) to the polar one — the spontaneous polarization of a ferroelectric — and how much effective charge "moves" when an atom is displaced (Born effective charge $Z^*$).

**Background for non-experts.** The polarization of a periodic solid **cannot** be computed as the dipole moment of the cell: that number depends on where the cell boundaries are cut. King-Smith and Vanderbilt showed that what is well defined is a geometric phase (Berry phase) accumulated by the occupied states when traversing the Brillouin zone along one direction. That phase is defined modulo $2\pi$, so the polarization is defined modulo a "quantum" $e\mathbf R/\Omega$: only **differences** between two structures connected by a path are measurable, exactly as in experiment (one measures the charge that flows while the structure changes, not $P$). pw.x computes that phase with `lberry = .true.` on "strings" of k-points parallel to a reciprocal vector; Olla-DFT prepares the strings correctly, follows the branch of the phase along the path, and checks the ionic part against its exact formula.

**Formulas.**

k strings (`qekit/modules/berry.py: cuerdas`): for every point $(i/n_\perp^{(1)}, j/n_\perp^{(2)})$ of the perpendicular mesh (`--kperp` 6×6), `nppstr` points (9 by default) along $\mathbf b_{\text{gdir}}$ with coordinate $l/(n_{\text{pp}}-1)$, $l = 0,\dots,n_{\text{pp}}-1$: the last point is the first one plus $\mathbf G$.

Ionic phase (`berry.fase_ionica`), in QE's units (the quantum is `MOD_TOT` = 2 if all valences are even, 1 if any is odd; `berry.modulo_de`):

$$
\varphi_{\text{ion}} = \sum_a \left[Z_a f_a^{(g)}\right]_{\bmod\, m_a}\Big|_{\bmod\, m}, \qquad m_a = \begin{cases}1 & Z_a \text{ odd}\\ 2 & Z_a \text{ even}\end{cases}
$$

- $Z_a$: valence charge of the pseudopotential of atom $a$ (electrons); $f_a^{(g)}$: fractional coordinate along `gdir`. The per-ion folding and the final folding reproduce what pw.x does; folding to $[-m/2, m/2)$ uses Fortran's `NINT` (`berry._nint`, half rounds away from zero), so that half a quantum comes out as $-1$ as in QE.

Electronic phase from Wannier centres (`berry.desde_wannier`), as an independent check:

$$
\varphi_{\text{el}} = -f_s\sum_n \bar r_n^{(g)}, \qquad f_s = 2
$$

- $\bar r_n^{(g)}$: fractional coordinate of Wannier centre $n$ along `gdir`; $f_s$ is the spin factor. The total phase is $\varphi_{\text{el}} + \varphi_{\text{ion}}$ folded.

Polarization and quantum (`berry.polarizacion`, `berry.cuanto`):

$$
P_g = \varphi\,\frac{|\mathbf R_g|}{\Omega}, \qquad
\Delta P_{\text{quantum}} = m\,\frac{|\mathbf R_g|}{\Omega}, \qquad
1\ e/\text{Å}^2 = 16.02176634\ \text{C/m}^2
$$

- $\varphi$: total phase in QE units (dimensionless, quantum $m$); $\mathbf R_g$: lattice vector `gdir` (Å); $\Omega$: volume (Å³). $P_g$ is the **projection** of $P\Omega/e$ onto $\mathbf R_g$, not the modulus of $\mathbf P$.

Branch tracking (`berry.desenrollar`): $\tilde\varphi_0 = \varphi_0$, $\tilde\varphi_i = \varphi_i + m\cdot\mathrm{round}\big((\tilde\varphi_{i-1} - \varphi_i)/m\big)$; a warning is issued if any jump $|\tilde\varphi_i - \tilde\varphi_{i-1}| > 0.25\,m$ (`FRACCION_SOSPECHOSA`).

Born effective charge (`berry.analizar`), with $\mathbf u$ the total displacement (Å) and $\mathbf B_g$ the `gdir` reciprocal vector (with $2\pi$):

$$
Z^*_{g} = \frac{2\pi\,\dfrac{d\tilde\varphi}{d\lambda}}{\mathbf u\cdot\mathbf B_g}
$$

- $d\tilde\varphi/d\lambda$: slope of the linear fit of the tracked phase versus $\lambda \in [0,1]$ (`np.polyfit` degree 1 if more than 2 points; finite difference otherwise). It is the $Z^*_{g,\hat u}$ component of the tensor. If $\mathbf u\perp\mathbf B_g$, it is not computed.

Adiabatic path (`berry._interpolar_estructuras`): positions interpolated in fractional coordinates by minimum image, $f(\lambda) = f_a + \lambda\,[(f_b - f_a) - \mathrm{round}(f_b - f_a)]$, and the cell interpolated linearly.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_berry` loads the polar structure, optionally `--reference` (centrosymmetric) or `--displace ATOM:dx,dy,dz` (Å, atom 1-based), and `--kperp`.
2. `berry.prepare` builds the list of structures (`--nlambda` 5 values of $\lambda$; a single point if there is no path), resolves pseudopotentials and cutoffs (`sweep.prepare_common(insulator=True)`) and, in every `pNN/`, writes `1_scf.in` and `2_berry.in` (`calculation='nscf'`, `occupations='fixed'`, `conv_thr 1e-10`, `nosym`, `noinv`, with `lberry`, `gdir` and `nppstr` inserted into `&CONTROL`), plus `correr.sh`/`correr.py`.
3. `--run` → `berry.correr`: pw.x on scf and berry at every point, skipping those that already contain `JOB DONE` unless `--redo`.
4. `--collect` → `berry.collect`: `leer_berry` extracts from `2_berry.out` `Ionic Phase`, `Electronic Phase`, `TOTAL PHASE`, `MOD_TOT`, `P = … (mod …) (e/Omega).bohr`, `direction of vector`, `Number of k-points per string`, `Number of different strings`; `valencias_de` reads the table "atomic species valence mass pseudopotential" from `1_scf.out`.
5. `berry.analizar`: unwraps the phases, converts to C/m², computes $\Delta P$ and, if the path is a displacement, $Z^*$; `comprobar_ionica` compares pw.x's ionic phase with $\sum Z_a f_a$ (warning if they differ by more than $10^{-4}$).
6. `berry.report` prints the table $\lambda$ / ionic / electronic / total / tracked / $P$; `export` writes `BERRY.dat` and `BERRY.txt`; `plot` draws $P(\lambda)$, the folded values from pw.x and a band one quantum wide.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Ionic, electronic, total phases, `MOD_TOT` | `pNN/2_berry.out` from pw.x (`lberry`) | `berry.leer_berry`, regular expression over the text |
| Valences $Z_a$ | `atomic species / valence` table in `1_scf.out` | `berry.valencias_de` |
| Cell, volume, $\mathbf R_g$, $\mathbf B_g$ | user structure (last one of the path) | `berry.cuanto`, `berry.analizar` |
| e/Å² → C/m² conversion | constant `berry.E_A2_A_C_M2` | 16.02176634 |
| Suspicious-jump threshold | `berry.FRACCION_SOSPECHOSA` | 0.25 of the quantum |
| Wannier centres (check) | `olla-dft wannier` | `berry.desde_wannier`, API/tests only |

**Limits and pitfalls.**

- **Insulators** only: the nscf is written with `occupations='fixed'`; in a metal the phase is undefined.
- A single point is useless: "Un solo punto. P está definida módulo el cuanto, así que este número por sí solo no significa nada".
- If a step moves the phase by more than 25 % of the quantum: "El seguimiento de la rama supone que el paso es pequeño; con saltos así, elegir la imagen más cercana es una apuesta. Sube --nlambda". If $|\Delta P| > 0.9$ quanta: "Comprueba con más puntos que no es un salto de rama disfrazado".
- Only **one component** (`--gdir`) is computed; for the vector $\mathbf P$ three runs are needed.
- If pw.x stops with "Wrong k-strings", `nosym`/`noinv` were almost certainly missing; Olla-DFT forces them, but a hand-edited input may lose them.
- In the figure, the markers "lo que escribe pw.x (plegado)" come from `berry.polarizacion_plegada`: $P = \varphi_{\text{tot}}/m \cdot \Delta P_{\text{quantum}}$ with the same `MOD_TOT` that `analizar` uses, so they coincide with the tracked branch at $\lambda = 0$ and differ from it only by integer multiples of the quantum.
- `desde_wannier` (check against Wannier centres) is not wired to the CLI; it is only used from Python or in the tests.
- No correction is applied for spin polarization or SOC (the spin factor is a fixed 2 in `desde_wannier`; pw.x handles it internally in `lberry`).

**References.**

- R. D. King-Smith, D. Vanderbilt, *Phys. Rev. B* **47**, 1651 (1993) — modern theory of polarization. DOI 10.1103/PhysRevB.47.1651.
- R. Resta, *Rev. Mod. Phys.* **66**, 899 (1994). DOI 10.1103/RevModPhys.66.899.
- N. A. Spaldin, "A beginner's guide to the modern theory of polarization", *J. Solid State Chem.* **195**, 2 (2012). DOI 10.1016/j.jssc.2012.05.010.
- D. Vanderbilt, *Berry Phases in Electronic Structure Theory* (Cambridge, 2018).

---

### `olla-dft hubbard` — Hubbard U by linear response (hp.x)

**What it answers.** How large the DFT+U parameter $U$ is for the localized ($d$ or $f$) orbitals of your system, computed by linear response with `hp.x` instead of copied from a paper, and — with `--cycle` — its self-consistent value.

**Background for non-experts.** Semilocal functionals (LDA, GGA) let an electron "see itself" (self-interaction), which over-delocalizes $d$ and $f$ orbitals and turns insulating oxides such as NiO into metals. DFT+U adds a penalty $U$ to fractional occupation of those orbitals. The value of $U$ is not a property of the element but of the system and of the *projection scheme* with which the occupations are counted. Cococcioni and de Gironcoli obtain it by measuring how the orbital occupation responds to a small perturbation of the potential: the "bare" response $\chi_0$ (without letting the rest of the system readjust) and the full one $\chi$. Their difference is the spurious curvature that $U$ must cancel. `hp.x` does that calculation with perturbation theory (DFPT) on a mesh of $\mathbf q$ vectors equivalent to a supercell. Since the $U$ obtained depends on the $U$ used in the starting scf, one must iterate until it stabilizes.

**Formulas.**

Linear response (computed by `hp.x`, not by Olla-DFT; `qekit/modules/hubbard.py`, docstring):

$$
U_I = \left(\chi_0^{-1} - \chi^{-1}\right)_{II}
$$

- $\chi_0$, $\chi$: response matrices of the occupations $n_I$ of Hubbard site $I$ to the perturbation $\alpha_J$ of the potential on site $J$, without and with self-consistent readjustment (eV⁻¹). $U_I$ in eV.

Self-consistency cycle (`hubbard.ciclo`):

$$
U^{(k+1)}_s = (1 - \mu)\,U^{(k)}_s + \mu\,U^{\text{hp}}_s\!\left[U^{(k)}\right], \qquad
\text{converged if } k \ge 1 \ \wedge\ \max_s\left|U^{\text{hp}}_s - U^{(k)}_s\right| < \text{tol}
$$

- $\mu$ = `--mixing` (1.0 by default), tol = `--tol` (0.05 eV), at most `--max-iter` = 6 rounds; $U^{(0)}_s$ = `U_SEMILLA` = $10^{-8}$ eV. The $U$ reported per species is the mean over its sites (`HubbardRun.U`).

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_hubbard` loads the structure; `--species` or, by default, `hubbard.elementos_hubbard` (those in the `ORBITAL_HUBBARD` table: 3d Sc–Zn, 4d Y–Cd, 5d Hf–Hg, 4f La–Lu); `--qgrid` 2×2×2; `--hubbard-style legacy|card` (the same selector as `gen`).
2. `hubbard.prepare` writes `scf.in` with `inputgen.build_pw_input(hubbard={s: U_seed}, hubbard_style=…, conv_thr=1e-15)`. With `legacy` (default, QE ≤ 7.0): `lda_plus_u = .true.`, `Hubbard_U(i) = 1e-8` and `U_projection_type = 'ortho-atomic'` inserted into `&SYSTEM` (`_fijar_proyeccion`). With `card` (QE ≥ 7.1): a `HUBBARD (<projection>)` card with `U El-orb 1e-8` and no `U_projection_type`, which is an error in those versions. `--projection` accepts `atomic`, `ortho-atomic`, `norm-atomic`, `wannier`, `pseudo`. And `hp.in` (`build_hp_input`: `nq1..3`, `conv_thr_chi = 1e-8`, `iverbosity = 2`). Fixed occupations unless `--metal`; `--nspin 2` and `--mag` are passed through.
3. `--cycle` → `hubbard.ciclo`: per iteration it creates `iterNN/`, runs pw.x (`runner.run_all`) and hp.x (`run_hp`, searched next to pw.x), reads `*.Hubbard_parameters.dat` (`collect` → `leer_parametros`, section "Hubbard U parameters", columns site/type/label/spin/new type/new label/U) and mixes.
4. `--collect` → `hubbard.collect` reads the first `*.Hubbard_parameters.dat` in the folder; `--intersite` adds `leer_v` (section "Hubbard V parameters", table atom 1 / atom 2 / distance in bohr / V) and writes `HUBBARD.card` with `tarjeta_hubbard` (`U El-orb value` and `V El-orb El-orb i j value`, with hp.x supercell indices and threshold `--v-threshold` 0.01 eV).
5. `hubbard.report` prints the table of $U$ per site, the cycle history and the warnings; `export` writes `HUBBARD_U.dat` and `HUBBARD_U.txt`, and suggests the line `olla-dft gen … --hubbard El=U`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $U$ per site | `<prefix>.Hubbard_parameters.dat` from hp.x | `hubbard.leer_parametros` |
| Intersite $V$, neighbour supercell | same output, section "Hubbard V parameters" | `hubbard.leer_v` |
| Corrected orbital | table `hubbard.ORBITAL_HUBBARD` | per element; `3d` if not in the table (`2p` for the second atom of a V) |
| Seed $U$ | `hubbard.U_SEMILLA` | 1e-8 eV |
| scf `conv_thr`, `conv_thr_chi` | constants in `prepare` / `build_hp_input` | 1e-15 Ry, 1e-8 |
| $\mathbf q$ mesh | `--qgrid` | 2×2×2 (8 cells) |

**Limits and pitfalls.**

- Olla-DFT **does not compute** $\chi$ or $U$: it reads them from hp.x. Without hp.x compiled (`make hp`) the command fails: "no se encontró hp.x junto a pw.x".
- By default the scf uses the `lda_plus_u`/`Hubbard_U(i)` syntax (QE ≤ 7.0); with QE ≥ 7.1 you must request `--hubbard-style card`. The `tarjeta_hubbard` docstring warns that the card "está probado contra la sintaxis documentada, no contra una corrida de QE 7.1, porque el QE de esta máquina es 6.6".
- A single round gives "U de PRIMERA ITERACIÓN. Depende del U que llevaba el scf de partida".
- With `nq = 1×1×1`: "la perturbación ve sus propias imágenes periódicas y el U sale mal. Usa al menos 2x2x2".
- The $U$ "solo vale con la MISMA proyección"; the report repeats it in every output.
- If the cycle does not converge in `--max-iter`: "Se hicieron N vueltas sin bajar de tol eV… si el número oscila arriba y abajo, baja --mixing a 0.5; si baja despacio pero siempre en el mismo sentido, sube --max-iter"; the command returns exit code 1.
- The HUBBARD-card orbital for an element outside the table is `3d` (or `2p` as the second atom of a $V$), which may be wrong (e.g. `O-2p` is fine, `S` would get `2p`).
- The indices of the $V$ pairs are in hp.x's **supercell** numbering; the card copies them verbatim, as QE requires.

**References.**

- M. Cococcioni, S. de Gironcoli, *Phys. Rev. B* **71**, 035105 (2005) — U by linear response. DOI 10.1103/PhysRevB.71.035105.
- I. Timrov, N. Marzari, M. Cococcioni, *Phys. Rev. B* **98**, 085127 (2018) — hp.x, DFPT for U. DOI 10.1103/PhysRevB.98.085127.
- I. Timrov, N. Marzari, M. Cococcioni, *Phys. Rev. B* **103**, 045141 (2021) — self-consistent U and V, ortho-atomic. DOI 10.1103/PhysRevB.103.045141.
- V. L. Campo Jr., M. Cococcioni, *J. Phys.: Condens. Matter* **22**, 055602 (2010) — DFT+U+V.
- S. L. Dudarev et al., *Phys. Rev. B* **57**, 1505 (1998) — simplified DFT+U formulation used by QE.

---

### `olla-dft align` — band alignment between two materials

**What it answers.** Where the valence band (and the conduction band) of one material sits relative to that of the other when they are brought into contact: the *offsets* $\Delta E_v$ and $\Delta E_c$ and the heterojunction type (I nested, II staggered, III broken).

**Background for non-experts.** Every periodic calculation fixes the zero of its potential arbitrarily (the $G = 0$ term of the Hartree potential), so directly subtracting the VBMs of two different calculations gives a meaningless number. There are two ways to put them on a common scale. In **vacuum mode**, each material is computed as a slab with vacuum and its VBM is measured relative to the vacuum level of its own calculation (its ionization potential); the offset is the difference. It ignores the charge transferred on forming the contact. In **interface mode** (Van de Walle and Martin) both bulks and also the interface are computed, and the macroscopically averaged electrostatic potential on each side of the interface serves as a bridge between the two scales: it is the only term that knows about the contact dipole.

**Formulas.**

Offsets (`qekit/modules/align.py: alinear`), with $E_v^{A}$ the VBM and $V^{A}_{\text{ref}}$ the reference of calculation $A$:

$$
\Delta E_v = \left(E_v^{A} - V_{\text{ref}}^{A}\right) - \left(E_v^{B} - V_{\text{ref}}^{B}\right) + \Delta\bar V, \qquad
\Delta E_c = \left(E_c^{A} - V_{\text{ref}}^{A}\right) - \left(E_c^{B} - V_{\text{ref}}^{B}\right) + \Delta\bar V
$$

- Vacuum mode: $V_{\text{ref}}$ = vacuum level (maximum of the planar potential, mean over a 20 % window around it; `fields.work_function`), $\Delta\bar V = 0$.
- Interface mode: $V_{\text{ref}}$ = mean electrostatic potential of the bulk cell ($\langle V\rangle$ of the planar average), and $\Delta\bar V = \bar V_A - \bar V_B$ measured at the interface (`align.puente_interfaz`).
- Everything in eV; the pp.x potential (`plot_num = 11`, $V_{\text{bare}} + V_H$) comes in Ry and is multiplied by `RY_EV = 13.605693122994`.

Interface bridge (`align.puente_interfaz`): planar average $\bar V(z)$ of the cube, periodic moving macroscopic average with window $w$ (`fields.macroscopic_average`; $w$ = `--window` or $L/8$), and

$$
\bar V_A = \langle \bar{\bar V}\rangle_{z \in [L/8,\, L/4]}, \qquad \bar V_B = \langle \bar{\bar V}\rangle_{z \in [5L/8,\, 3L/4]}, \qquad \Delta\bar V = \bar V_A - \bar V_B
$$

- Material $A$ is assumed to occupy the first half of the interface cell and $B$ the second.

Alignment type (`align.alinear`), on $B$'s scale (VBM of $B$ at 0): $v_A = \Delta E_v$, $c_A = E_g^{B} + \Delta E_c$, $v_B = 0$, $c_B = E_g^{B}$:

- `=` if $|\Delta E_v| < 0.05$ and $|\Delta E_c| < 0.05$ eV (`TOL_ALINEADOS`);
- I if one gap contains the other ($v_A \le v_B \wedge c_A \ge c_B$, or the reverse);
- III if $c_A \le v_B$ or $c_B \le v_A$;
- II in every other case.

**How Olla-DFT computes it.**

1. `qekit/cli.py: _cmd_align` receives folders `a` and `b`, `--interface FOLDER` (switches on interface mode), `--axis` (c by default), `--window`, `--names`.
2. `align.leer_lado` reads the XML (`qeout.read_xml`): VBM = `<highestOccupiedLevel>`, CBM = `<lowestUnoccupiedLevel>`, $E_F$; without HOMO it fails: "no da un VBM. En un metal no hay banda de valencia que alinear; y si es un aislante, al cálculo le faltan bandas vacías (nbnd) o no usó occupations='fixed'". Without LUMO the side is flagged `es_metal` and only $\Delta E_v$ is given.
3. `align._potencial` reuses `potencial.cube` or runs `pp.x` (`fields.run_pp` with `plot_num = 11`, `output_format = 6`) and reads it with `fields.read_cube`.
4. Vacuum mode: `fields.work_function` gives `v_vacuum` and the plateau flatness; interface mode: `fields.planar_average` and its mean.
5. With `--interface`, `align.puente_interfaz` computes $\Delta\bar V$.
6. `align.alinear`, `report` (table VBM/CBM/gap relative to the reference, offsets, type and which material each carrier goes to in type II), `export` (`ALINEAMIENTO.dat`, `.txt`) and `plot` (box diagram). The box positions come from `align.posiciones_en_escala_de_b` — $v_A = \Delta E_v$, $c_A = E_g^{B} + \Delta E_c$, $v_B = 0$, $c_B = E_g^{B}$ — the same convention with which `alinear` classifies the type, so that report, export and figure cannot disagree.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| VBM, CBM, $E_F$ | `<highestOccupiedLevel>`, `<lowestUnoccupiedLevel>`, `<fermi_energy>` in `prefix.xml` | `align.leer_lado`; requires fixed occupations |
| Electrostatic potential | `potencial.cube` from pp.x (`plot_num = 11`) | `fields.read_cube`; Ry → eV with 13.605693122994 |
| Vacuum level and flatness | maximum of the planar average, 20 % window | `fields.work_function` |
| Macroscopic window | `--window` or $L/8$ | `align.puente_interfaz` |
| "Aligned" threshold | `align.TOL_ALINEADOS` | 0.05 eV |
| Flatness threshold | constant in `alinear` | 0.05 eV |

**Limits and pitfalls.**

- Vacuum mode: the report always warns: "son las dos superficies AISLADAS. Al ponerlas en contacto se transfiere carga y aparece un dipolo de interfaz que desplaza el offset, típicamente entre 0.1 y 0.5 eV".
- If the vacuum plateau varies by more than 0.05 eV: "O falta vacío, o la losa tiene dipolo neto: usa --dipole al generarla. El nivel de vacío es la referencia de todo esto, así que ese error entra entero en el offset".
- Interface mode assumes $A$ lies in the first half of the cell and $B$ in the second, and uses two fixed windows ($[L/8, L/4]$ and $[5L/8, 3L/4]$); an asymmetric interface or layers of different thickness give a wrong bridge without warning.
- VBM/CBM are read from `highestOccupiedLevel`/`lowestUnoccupiedLevel`, which depend on the scf k-mesh; no band analysis is performed.
- If $A$ has no CBM (a metal or no empty bands), the figure draws its gap as a box of height $E_g^{A}$ (or 1 eV if there is no gap either) above $v_A$: a visual filler, not a datum.
- The offsets carry the systematic error of the functional; `TIPOS["="]` reminds that "con funcionales semilocales el error frente al experimento es de varias décimas".

**References.**

- C. G. Van de Walle, R. M. Martin, *Phys. Rev. B* **35**, 8154 (1987) — alignment via macroscopic potential. DOI 10.1103/PhysRevB.35.8154.
- A. Baldereschi, S. Baroni, R. Resta, *Phys. Rev. Lett.* **61**, 734 (1988) — macroscopic average. DOI 10.1103/PhysRevLett.61.734.
- L. Kleinman, *Phys. Rev. B* **24**, 7412 (1981) — the arbitrary zero of the potential in periodic calculations.
- J. Tersoff, *Phys. Rev. B* **30**, 4874 (1984) — alignment and interface dipoles.

## Mechanics, vibrations, temperature and transport

This part documents the physics Olla-DFT implements in the commands that go from the total energy to mechanical, vibrational, thermal and transport properties: from convergence tests (`converge`, `tune`) and the equation of state (`eos`), through elastic constants (`elastic`, `derived`), strain sweeps (`strain`), surfaces and layered materials (`gamma`, `layers`, `xrd`, `exfoliate`), phonons and everything derived from them (`phonons`, `qha`, `thermochem`, `kappa`, `elph`), molecular dynamics (`md`), diffusive and ballistic electronic transport (`transport`, `ballistic`) and the cost estimator (`cost`). Every section was written by reading the code in `qekit/modules/*.py` and `qekit/cli.py`, and it only lists the formulas the code really executes, with the constants and defaults exactly as written. Whenever a docstring promises something the code does not do, it is said in "Limits and pitfalls". A note on file names: the module `qekit/modules/thermo.py` does NOT contain the harmonic thermodynamics (that lives in `phonons.thermodynamics`) but the convex hull of formation energies used by the `hull` command, which is documented elsewhere.

---

### `olla-dft converge` — Convergence of cutoffs and k-mesh

**What it answers.** From which `ecutwfc`, `ecutrho` or k-point mesh does the total energy stop changing by more than a threshold (1 meV/atom by default)? It is the first question for any new system.

**Background for non-experts.** A plane-wave calculation describes the electrons with a sum of waves; `ecutwfc` sets how many waves are included (the "resolution" of the wavefunctions), `ecutrho` the resolution of the charge density, and the k-mesh how many points of the Brillouin zone are sampled. With too few waves or points the result is coarse; with too many, the calculation costs more without gaining anything. A convergence test repeats the same calculation while raising the parameter and looks at when the energy "flattens", like adjusting the zoom of a microscope until the image stops changing.

The criterion Olla-DFT uses has an important subtlety: it compares each point against the DENSEST one in the series (the last), not against its previous neighbour. Two adjacent points may look alike by chance in the middle of a curve that has not flattened yet; comparing them with each other is the usual mistake.

**Formulas.** Per-atom difference with respect to the densest point (`converge.ConvergenceRun.per_atom_diffs`):

$$
\Delta E_i = \frac{|E_i - E_{\mathrm{ref}}|}{N_{\mathrm{at}}} \times 1000
$$

- $E_i$: total energy of point $i$, in eV per cell (read from the XML and converted from Hartree with $27.211386245988$ eV).
- $E_{\mathrm{ref}}$: energy of the last point that finished (the densest).
- $N_{\mathrm{at}}$: atoms in the cell.
- $\Delta E_i$: in meV/atom.

Convergence index (`converge.ConvergenceRun.converged_index`): the first $i$ such that every $\Delta E_j$ with $j \ge i$ satisfies $\Delta E_j \le$ threshold (failed points are ignored). k-mesh from a spacing (`kpoints.kgrid_from_spacing`):

$$
n_i = \left\lceil \frac{|\mathbf{b}_i|}{k_{\mathrm{spacing}}} \right\rceil, \qquad |\mathbf{b}_i| \text{ including the factor } 2\pi
$$

- $\mathbf{b}_i$: reciprocal vectors in Å⁻¹; $k_{\mathrm{spacing}}$ in Å⁻¹ (configuration `kspacing`, default 0.20). Directions with ≥ 8 Å of vacuum get a single point.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_converge` loads the structure and calls `qekit/modules/converge.py: prepare`.
2. `sweep.prepare_common` resolves pseudopotentials and cutoffs (`pseudo.recommend_cutoffs`: the maximum declared by the UPF files; if none declares any, `ecutwfc` from the configuration (60 Ry) and `dual` (8); `ecutrho` never below $4\,\mathrm{ecutwfc}$).
3. Default series: `ecutwfc` = 30, 40, …, 100 Ry with `ecutrho = dual × ecutwfc`; `ecutrho` = 4, 6, 8, 10, 12 × ecutwfc; `kmesh` = meshes for the spacings 0.40, 0.30, 0.25, 0.20, 0.15, 0.12 Å⁻¹ (without repeats). `--values` replaces the series (for `kmesh` it accepts `8x8x8` or spacings).
4. One `pw.in` (`calculation='scf'`, `conv_thr = 1e-8`, `tstress`/`tprnfor` on) per point via `sweep.write_scf_job`, plus `run.sh` and `run.py`.
5. With `--run`, `runner.run_all` executes `pw.x`; with `--collect`, `converge.collect` reads `out/*.xml` (`qeout.read_xml`, tag `<total_energy><etot>`).
6. `converge.report` prints the table, the convergence point and the recommendation (`--ecutwfc N` or the mesh); `converge.export` writes `CONVERGENCIA.dat` and `.txt`; `converge.plot` draws $|\Delta E|$ on a log scale with the threshold band.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Total energy | pw.x XML (`output/total_energy/etot`, Hartree) | `qeout.read_xml` → eV |
| scf convergence | XML (`convergence_info/scf_conv/convergence_achieved`) | an unconverged point counts as failed under `--run` |
| Threshold | `--threshold` parameter | 1.0 meV/atom by default |
| Base cutoffs | UPF headers or `olla-dft config` | `pseudo.recommend_cutoffs` |
| Fixed k-mesh | `sweep.default_grid` | configuration `kspacing` (0.20 Å⁻¹) |
| Ry ↔ eV | `qeout.RY_EV` | 13.605693122994 eV |

**Limits and pitfalls.** It only looks at the total energy; the report warns: "convergence depends on the property: the total energy converges before stresses or phonons". If only the last point passes, it says: "Only the last point is below … there is no margin to be sure it has already flattened there". If none passes: "does NOT converge within … Extend the series towards denser values". The `energies` field of the dataclass is commented as "eV per cell", but the table is printed in Ry (divided by `RY_EV`): not a bug, just a display conversion. With `--collect` the inputs are not rewritten (`sweep.set_write_inputs(False)`), so the report describes what actually ran.

**References.** Quantum ESPRESSO manual (`pw.x`, variables `ecutwfc`, `ecutrho`, `K_POINTS`). Monkhorst and Pack, *Phys. Rev. B* 13, 5188 (1976), DOI 10.1103/PhysRevB.13.5188.

---

### `olla-dft tune` — Adaptive convergence recommendation

**What it answers.** Given an already generated `CONVERGENCIA.dat`, is the series converged, and if not, which value should be tried next?

**Background for non-experts.** It is pure post-processing of the `converge` table: it applies the same criterion ("from this point on, the whole tail stays within the threshold") and, when it is not met, proposes the next value with a sensible step instead of leaving the user to guess.

**Formulas.** Criterion (`tuning.analyze`): minimum index $i$ with $|\Delta E_j| \le$ threshold for all $j \ge i$. States: `ready` (such $i$ exists and is not the last), `confirm` (only the last passes), `extend` (none). Next value (`tuning._next_value`):

$$
v_{\mathrm{next}} = v_{\mathrm{last}} + \max\!\left(\mathrm{median}\{v_{k+1}-v_k > 0\},\; 0.10\,|v_{\mathrm{last}}|\right)
$$

- With fewer than two values, or no positive steps: $v_{\mathrm{next}} = 1.25\,v_{\mathrm{last}}$ (or $v_{\mathrm{last}}+1$ if it is not positive).

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_tune` → `qekit/modules/tuning.py: read` reads the numeric rows (column 1 value, 2 energy in Ry, 3 $\Delta E$ in meV/atom; comments and NaN are skipped).
2. `tuning.analyze` applies the criterion and picks the state and the recommended value.
3. `tuning.report` prints it; with `-o`, `tuning.export` writes a JSON (`CONVERGENCIA_RECOMENDACION.json` by default).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Value, E, ΔE | `CONVERGENCIA.dat` (from `olla-dft converge`) | `tuning.read`; ΔE taken in absolute value |
| Threshold | `--threshold` | 1.0 meV/atom if omitted; must be > 0 |

**Limits and pitfalls.** It runs nothing and reads no QE outputs: only the table. It uses the ΔE column as written, which `converge` computed against the densest point of THAT series; if points are added later, the table must be regenerated. The report reminds: "The energy property may converge before forces, phonons or tensors".

**References.** None specific; it is the same logic as `converge`.

---

### `olla-dft eos` — E–V equation of state and bulk modulus

**What it answers.** What are the equilibrium volume $V_0$, minimum energy $E_0$, bulk modulus $B_0$ and its derivative $B_0'$ of the crystal? And, if cubic, the lattice parameter $a_0$.

**Background for non-experts.** The cell is compressed and stretched a little around the starting size, the energy is computed at each volume, and a valley-shaped curve is fitted. The bottom of the valley is the equilibrium volume; the "stiffness" of the valley (its curvature) is the bulk modulus, which measures how hard you must squeeze to reduce the volume. Olla-DFT fits three different equations; if all three agree, the fit is reliable, and if they disagree, usually the range is too narrow or there are noisy points.

**Formulas.** Third-order Birch–Murnaghan (`eos.birch_murnaghan`), with $\eta = (V_0/V)^{2/3}$:

$$
E(V) = E_0 + \frac{9 V_0 B_0}{16}\left[(\eta-1)^3 B_0' + (\eta-1)^2 (6 - 4\eta)\right]
$$

Murnaghan (`eos.murnaghan`):

$$
E(V) = E_0 + \frac{B_0 V}{B_0'}\left[\frac{(V_0/V)^{B_0'}}{B_0'-1} + 1\right] - \frac{B_0 V_0}{B_0'-1}
$$

Vinet (`eos.vinet`), with $x = (V/V_0)^{1/3}$ and $\xi = \tfrac{3}{2}(B_0'-1)$:

$$
E(V) = E_0 + \frac{9 B_0 V_0}{\xi^2}\left[1 + \left(\xi(1-x) - 1\right) e^{\xi(1-x)}\right]
$$

- $V$, $V_0$: volumes in Å³; $E$, $E_0$: eV; $B_0$: eV/Å³ inside the fit, converted to GPa with `EV_A3_GPA = 160.21766208`; $B_0'$: dimensionless.
- Fit seed (`eos.fit`): parabola $E = aV^2+bV+c$ → $V_0 = -b/2a$, $B_0 = 2aV_0$, $B_0' = 4$.
- RMSE: $\sqrt{\langle (E - E_{\mathrm{fit}})^2 \rangle}/N_{\mathrm{at}}$, in eV/atom (printed in meV/atom).
- Cubic lattice parameter (`eos.fit`, field `EOSFit.a0`, with $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ measured in `prepare`): $a_0 = (V_0 \cdot V_{\mathrm{conv}}/V_{\mathrm{prim}})^{1/3}$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_eos` → `qekit/modules/eos.py: prepare`. Requires `--npoints` ≥ 5 (default 9); `--span` 0.10 (±10 % in VOLUME); `--scale` (linear centring factor) 1.0.
2. Equally spaced volume factors $f \in [c^3(1-s),\, c^3(1+s)]$; linear factor $f^{1/3}$ applied to the cell with `set_cell(..., scale_atoms=True)`.
3. Cubic or not is decided by asking spglib (`structure.symmetry_dataset`, space group ≥ 195) and $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ is stored via `structure.conventional`.
4. One `scf` (or `relax` with `--relax-ions`) per volume, all with the SAME k-mesh (`sweep.default_grid`), written by `sweep.write_scf_job` into `V_<factor>/pw.in`.
5. `--run` executes `pw.x`; `--collect`/`eos.collect` reads `etot` from the XML.
6. `eos.fit_all` fits the three equations with `scipy.optimize.curve_fit` (`maxfev=20000`); the fit is rejected if $V_0 \notin (0.6 V_{\min}, 1.4 V_{\max})$ or $B_0 \le 0$.
7. `eos.report` prints the table, the three fits, the Birch–Murnaghan result and the spread between equations; `eos.export` writes `EOS.dat` and `EOS.txt`; `eos.plot` draws $E - E_0$ with residuals.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Total energy per volume | pw.x XML (`etot`) | `qeout.read_xml` |
| Volume of each point | scaled cell (ASE) | $|\det(\mathbf{a})|$ in Å³ |
| Cubic symmetry and conventional cell | spglib via `structure` | `symmetry_dataset`, `conventional` |
| eV/Å³ → GPa | `eos.EV_A3_GPA` | 160.21766208 |
| Non-linear fit | `scipy.optimize.curve_fit` | library |

**Limits and pitfalls.** It does not relax the cell shape: only isotropic scaling (for non-cubic crystals $c/a$ stays fixed unless `--relax-ions` is used, which relaxes positions, not the cell). It warns: "V₀ falls OUTSIDE the computed range. Re-run the sweep centred on that volume" and, if the three equations differ by more than 5 %: "usually indicates missing points or a very narrow volume range". `--relax-ions` uses `calculation='relax'`, so those inputs carry `forc_conv_thr = 1e-4`.

**References.** F. Birch, *Phys. Rev.* 71, 809 (1947), DOI 10.1103/PhysRev.71.809. F. D. Murnaghan, *Proc. Natl. Acad. Sci. USA* 30, 244 (1944). P. Vinet, J. Ferrante, J. R. Smith and J. H. Rose, *J. Phys. C* 19, L467 (1986).

---

### `olla-dft elastic` — Elastic constants by stress–strain

**What it answers.** What are the elastic constants $C_{ij}$ of the crystal (or of the sheet, with `--2d`), the bulk, shear and Young moduli, the Poisson ratio, and is the structure mechanically stable?

**Background for non-experts.** A slightly deformed solid responds with a stress (force per unit area) proportional to the strain: this is the generalised Hooke's law, and the proportionality constants are the $C_{ij}$. Olla-DFT deforms the cell by ±1 % (and ±0.5 %) along each of the six independent directions (three stretches and three shears), asks `pw.x` for the stress tensor in each, and fits a straight line. Since `pw.x` gives all six stresses at once, every strain provides six equations, far fewer runs than the energy method.

In a sheet (graphene, MoS₂) stretching the vacuum direction makes no sense: only the two in-plane directions and the in-plane shear are strained, and the constants are given in N/m by multiplying by the cell height, so the vacuum cancels.

**Formulas.** Applied strain (`elastic.strain_matrix`): $\mathbf{a}' = \mathbf{a}(\mathbf{I}+\boldsymbol{\varepsilon})$ with $\varepsilon_{ii}=\delta$ for normal components and $\varepsilon_{ij}=\varepsilon_{ji}=\delta/2$ for shears (Voigt convention, $\varepsilon_4 = 2\varepsilon_{23}$). Fit (`elastic.fit`), with the sign inverted because the tensor `pw.x` writes is the opposite of the elasticity one:

$$
C_{ij} = -\frac{\partial\,\sigma^{\mathrm{pw}}_i}{\partial \varepsilon_j}\Big|_{\text{least squares}}, \qquad \sigma^{\mathrm{pw}}_i \to \sigma^{\mathrm{pw}}_i - \sigma^{\mathrm{pw}}_i(\text{ref})
$$

Voigt, Reuss and Hill averages (`elastic.moduli`), with $S = C^{-1}$:

$$
B_V = \frac{(C_{11}+C_{22}+C_{33}) + 2(C_{12}+C_{23}+C_{13})}{9}, \quad
G_V = \frac{(C_{11}+C_{22}+C_{33}) - (C_{12}+C_{23}+C_{13}) + 3(C_{44}+C_{55}+C_{66})}{15}
$$

$$
B_R = \frac{1}{(S_{11}+S_{22}+S_{33}) + 2(S_{12}+S_{23}+S_{13})}, \quad
G_R = \frac{15}{4(S_{11}+S_{22}+S_{33}) - 4(S_{12}+S_{23}+S_{13}) + 3(S_{44}+S_{55}+S_{66})}
$$

$$
B_H = \tfrac{1}{2}(B_V+B_R),\quad G_H = \tfrac{1}{2}(G_V+G_R),\quad
E = \frac{9 B_H G_H}{3B_H + G_H},\quad \nu = \frac{3B_H - 2G_H}{2(3B_H+G_H)},\quad
A^U = 5\frac{G_V}{G_R} + \frac{B_V}{B_R} - 6
$$

- $C_{ij}$, $B$, $G$, $E$ in GPa; $\nu$ and $A^U$ dimensionless; Pugh ratio $B_H/G_H$ (ductility threshold 1.75).
- Stability (generalised Born): all eigenvalues of $\tfrac{1}{2}(C+C^T)$ positive.

Sheet (`elastic.constantes_2d`, `modulos_2d`, `born_2d`): $C^{2D}_{ij} = C_{ij}\,c\times 0.1$ (GPa·Å → N/m), with $c$ the cell height;

$$
Y_x = \frac{C_{11}C_{22}-C_{12}^2}{C_{22}},\quad \nu_x = \frac{C_{12}}{C_{22}},\quad
K = \frac{C_{11}+C_{22}+2C_{12}}{4},\quad G = C_{66};\qquad
C_{11}>0,\; C_{66}>0,\; C_{11}C_{22}-C_{12}^2>0
$$

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_elastic` → `qekit/modules/elastic.py: prepare`. In 3D the structure is ALWAYS taken to spglib's standardised primitive cell (`structure.primitive`) so that the Cartesian axes coincide with the crystal-physical ones; in `--2d` it is not (it requires vacuum along $c$ via `kpoints.direcciones_con_vacio`).
2. Crystal family from the space-group number (`elastic.crystal_family`: ≥195 cubic, ≥168 hexagonal, ≥143 trigonal, ≥75 tetragonal, ≥16 orthorhombic).
3. Strains: `--delta` 0.010, `--npoints` 4 (even) → ±δ/2, ±δ. Components: the 6 Voigt ones, or (1, 2, 6) in 2D. Plus an undeformed reference cell.
4. `--ion-mode auto` (default): `scf` (fixed ions) for ε1–ε3 and `relax` for ε4–ε6; `relax`: all relaxed; `fixed`: all fixed. Deformed cells use `conv_thr = 1e-9`.
5. `pw.x` with `tstress = .true.`; `elastic.collect` reads `<stress>` from the XML (Ha/bohr³ → GPa with `qeout.HA_BOHR3_GPA = 29421.026`).
6. `elastic.fit` fits column by column with `np.polyfit(..., 1)`; `elastic.symmetrize` averages the equivalents of the family (cubic, hexagonal with $C_{66}=(C_{11}-C_{12})/2$, partial tetragonal); `elastic.moduli` computes VRH and Born.
7. `elastic.report`/`_report_2d`, `elastic.export` (`ELASTIC_C.dat`, `ELASTIC.txt`), `elastic.plot` (σ–ε lines).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Stress tensor | pw.x XML (`output/stress`, Ha/bohr³, Fortran order) | `qeout.read_xml`; requires `tstress=.true.` (always set) |
| Space group and family | spglib (`structure.symmetry_dataset`) | `elastic.crystal_family` |
| Cell height (2D) | `|a_3|` of the input cell | `ElasticRun.altura` |
| GPa·Å → N/m | `elastic.GPA_A_NM` | 0.1 |
| Assumed thickness (2D) | `--thickness` | only for the GPa equivalent, "a convention, not a measurement" |

**Limits and pitfalls.** Symmetrisation only covers cubic, hexagonal and (partially) tetragonal; trigonal, orthorhombic and monoclinic/triclinic are left with the symmetrised matrix $\tfrac{1}{2}(C+C^T)$ and nothing else. The Born criterion is the general one (eigenvalues), not the family-specific inequalities. It warns if the residual stress of the reference cell exceeds 0.5 GPa: "it is high. Relax the cell with vc-relax before computing the elastic constants". In 2D it warns that with `--ion-mode auto` the identity $C_{66}=(C_{11}-C_{12})/2$ stops holding even if the sheet is isotropic. The GPa equivalent of a sheet depends on the chosen thickness: "This thickness is a CONVENTION, not a measurement". With fewer than 3 stresses read, no fit is done.

**References.** R. Hill, *Proc. Phys. Soc. A* 65, 349 (1952), DOI 10.1088/0370-1298/65/5/307. S. I. Ranganathan and M. Ostoja-Starzewski, *Phys. Rev. Lett.* 101, 055504 (2008), DOI 10.1103/PhysRevLett.101.055504 ($A^U$ index). F. Mouhat and F.-X. Coudert, *Phys. Rev. B* 90, 224104 (2014), DOI 10.1103/PhysRevB.90.224104 (Born criteria). S. F. Pugh, *Philos. Mag.* 45, 823 (1954).

---

### `olla-dft strain` — Strain sweep: gap, energy and magnetic moment

**What it answers.** How do the energy, band gap, pressure and magnetic moment change when the cell is strained (biaxial, uniaxial, hydrostatic or shear)? What is the deformation potential $dE_{\mathrm{gap}}/d\varepsilon$ and at which strain does the gap close?

**Background for non-experts.** Stretching or compressing a crystal changes the interatomic distances and with them the electronic structure: the gap can open, close or change type, and a magnetic material can lose its moment. The deformation potential is the slope of that response, and it is what gets compared with pressure experiments or with sheets on substrates that stretch them. Olla-DFT ALWAYS applies each strain to the original cell (not to the previous point's cell, which would accumulate error) and relaxes the internal positions at each point.

**Formulas.** Strain (`strain.matriz`): $\mathbf{a}' = \mathbf{a}_0(\mathbf{I}+\boldsymbol{\varepsilon})$ with the Voigt components of each mode (`strain.MODOS`: biaxial (xx, yy), uniaxial-a/b/c, hydrostatic (xx, yy, zz), xy shear with $\varepsilon_{xy}=\varepsilon_{yx}=\varepsilon/2$). Energy minimum by a local parabola (`strain.minimo`, up to 3 points on each side of the sampled minimum): $\varepsilon^* = -b/2a$. Deformation potential (`strain.potencial_deformacion`):

$$
E_{\mathrm{gap}}(\varepsilon) \approx m\,\varepsilon + b, \qquad R^2 = 1 - \frac{\sum (y - \hat y)^2}{\sum (y-\bar y)^2}
$$

- $m$: in eV per unit strain (fraction, not per cent). Gap $= E_{\mathrm{LUMO}} - E_{\mathrm{HOMO}}$ from the XML.
- Gap closing (`strain.cierre_de_gap`): linear interpolation of the strain at which the gap crosses 0.02 eV.

2D biaxial modulus (`strain.modulo_biaxial`, biaxial mode only, points with $|\varepsilon|\le 0.03$):

$$
Y_{2D} = \frac{1}{A_0}\frac{d^2E}{d\varepsilon^2} \times 16.021766 \;\; [\mathrm{N/m}], \qquad \frac{d^2E}{d\varepsilon^2} = 2a
$$

- $A_0 = |\mathbf{a}_1\times\mathbf{a}_2|$ in Å²; $E$ in eV; 1 eV/Å² = 16.021766 N/m. It is the combination $C_{11}+2C_{12}+C_{22}$, NOT the Young modulus.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_strain` → `qekit/modules/strain.py: prepare`. `--range MIN:MAX:N` in PER CENT (default `-5:5:11`; it rejects $|\varepsilon| > 30$ %, requires $N \ge 3$ and adds $\varepsilon=0$ if missing).
2. Estimates `nbnd` with empty bands (`inputgen._estimate_nbnd`: $\lceil 1.25\,N_{\mathrm{occ}} + 4\rceil$, ×1.2+2 if `nspin=2`) so that a LUMO exists.
3. Calculation type: `relax` (default), `scf` with `--fixed-ions`, `vc-relax` with `cell_dofree` (`z`, `shape` or `2Dxy` depending on the mode) with `--relax-perp`.
4. One input per strain (`sweep.write_scf_job`, accepts `--nspin/--mag`, `--hubbard`, `--vdw`).
5. `strain.collect` reads from each XML `etot`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`, `stress` (pressure = trace/3), `magnetization/total` and `convergence_achieved`.
6. `strain.report` prints the table, the minimum, the deformation potential, the gap closing, the moment and the biaxial modulus (if there is vacuum along $c$); `strain.export` (`STRAIN.dat`, `.txt`); `strain.plot` (two panels).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Energy, HOMO, LUMO | pw.x XML (`etot`, `highestOccupiedLevel`, `lowestUnoccupiedLevel`) | `qeout.read_xml` |
| Pressure | XML (`stress`, trace/3, QE sign) | `QEResult.pressure` in GPa |
| Magnetic moment | XML (`magnetization/total`) | μ_B per cell |
| Convergence | XML (`convergence_achieved`) | rows flagged `<< SIN CONVERGER` |
| Reference area and volume | input cell | `StrainRun.area0`, `volume0` |
| eV/Å² → N/m | literal constant in `modulo_biaxial` | 16.021766 |

**Limits and pitfalls.** If the HOMO exists but the LUMO does not (no empty bands) the gap column stays empty and it warns: "No gap in the table: the calculations have no empty bands". If the minimum is not at $\varepsilon=0$ (|ε*| > 0.3 %): "The starting structure was not relaxed". With $R^2 < 0.9$: "The gap does not respond linearly in this range". Biaxial without vacuum: "if it is bulk material, perhaps you wanted 'hidrostatica'". `--relax-perp` with hydrostatic strain is rejected. Unconverged points DO enter the table (they are read from the XML even if the runner marks them failed) but it warns that they "are NOT comparable with the rest". The "gap" is that of the scf k-mesh, not the fundamental gap from a band path.

**References.** J. Bardeen and W. Shockley, *Phys. Rev.* 80, 72 (1950), DOI 10.1103/PhysRev.80.72 (deformation potentials). C. G. Van de Walle, *Phys. Rev. B* 39, 1871 (1989).

---

### `olla-dft gamma` — Surface energy and the Fiorentini–Methfessel fit

**What it answers.** How much energy per unit area does it cost to create the (hkl) surface of a crystal, $\gamma$ in J/m², and how does it converge with slab thickness?

**Background for non-experts.** Cutting a crystal leaves atoms with fewer neighbours: that costs energy, and the surface energy is that cost per unit area. It is computed with a "slab" (a few atomic layers with vacuum above and below) by subtracting what the same atoms would be worth inside the crystal. The problem is that the bulk energy comes from ANOTHER calculation, with another k-mesh, and any residual error per atom is multiplied by the number of atoms: $\gamma$ does not converge, it drifts. The way out is to fit a straight line $E_{\mathrm{slab}}(N)$ over several thicknesses: the intercept gives $2\gamma A$ and the slope a bulk energy consistent with the slabs themselves.

**Formulas.** Direct (`surfen.GammaRun.gamma_directo`):

$$
\gamma_{\mathrm{dir}}(N) = \frac{E_{\mathrm{slab}}(N) - N_{\mathrm{at}}\,E_{\mathrm{bulk}}}{2A}
$$

Fit (`surfen.ajustar`), least squares over the pairs $(N_{\mathrm{at}}, E_{\mathrm{slab}})$:

$$
E_{\mathrm{slab}}(N_{\mathrm{at}}) = 2\gamma A + N_{\mathrm{at}}\,E_{\mathrm{bulk}}^{\mathrm{fit}}
$$

- $E_{\mathrm{slab}}$: eV per slab cell; $N_{\mathrm{at}}$: atoms in the slab; $E_{\mathrm{bulk}}$: eV/atom from the separate conventional-cell calculation; $A = |\mathbf{a}_1\times\mathbf{a}_2|$ in Å² (one face); the 2 stands for the two faces (`GammaRun.caras` is always 2).
- $\gamma$ in eV/Å² → J/m² with `EV_A2_A_J_M2 = 16.021766`. Cleavage energy $= 2\gamma$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_gamma` → `qekit/modules/surfen.py: prepare`. `--miller` (default `1 0 0`), `--layers` 3,4,5,6 (at least two, minimum 3 layers), `--vacuum` 20 Å.
2. `builder.surface` cuts each slab from the CONVENTIONAL cell with `ase.build.surface`, centres the vacuum, detects whether it is symmetric ($z$ profile equal to its mirror, tol 0.3 Å) and polar (composition of the top layer ≠ bottom) and emits warnings.
3. Unless `--no-reduce` or `--fix`, `surfen.reducir_losa` replaces the slab by its spglib primitive if the $c$ axis does not change (same $\gamma$, fewer atoms).
4. k-mesh fixed with the smallest slab and reused for all (`sweep.default_grid`); the bulk (conventional cell) gets its own mesh. `scf` or `relax` (`--relax`) calculations, with options `--vdw`, `--dipole` (`dipole_correction=3`), `--nspin/--mag`.
5. `surfen.collect` reads `etot` and `convergence_achieved` from each XML and fits (`surfen.ajustar`).
6. `surfen.report` prints the direct-γ table, the drift, the fit with $R^2$ and the difference $E_{\mathrm{bulk}}^{\mathrm{fit}} - E_{\mathrm{bulk}}$; `surfen.export` (`GAMMA.dat`, `GAMMA.txt`); `surfen.plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E_{\mathrm{slab}}(N)$ | XML of each `capasNN/` (`etot`) | `qeout.read_xml` |
| $E_{\mathrm{bulk}}$ | XML of `_bulto/` (`etot`) / atoms of the conventional cell | skipped with `--no-bulk` |
| Area $A$ | cell of the thinnest slab | `np.cross(a1, a2)` |
| Symmetry / polarity | `builder.surface` (geometry) | tolerance 0.3 Å |
| eV/Å² → J/m² | `surfen.EV_A2_A_J_M2` | 16.021766 |

**Limits and pitfalls.** It is the LINEAR Fiorentini–Methfessel fit; Boettger's incremental scheme (which takes $E_{\mathrm{bulk}}$ from the difference between consecutive slabs) is not implemented, and the docstring says so. Non-symmetric slab: "γ is the AVERAGE of its two faces, not that of one". Polar: "use --dipole". Without `--relax`: "Unrelaxed: γ comes out high. Surface relaxation lowers it by 5 to 20 %". If the drift between slabs exceeds 0.05 J/m²: "It does not converge … It is the residual error of E_bulk multiplied by the number of atoms, not physics. The good value is the fitted one". $R^2 < 0.999$: "either some point lacks convergence, or the thin slabs do not yet have a bulk interior". With `--fix` the cell is not reduced (constraints refer to specific atoms).

**References.** V. Fiorentini and M. Methfessel, *J. Phys.: Condens. Matter* 8, 6525 (1996), DOI 10.1088/0953-8984/8/36/005. J. C. Boettger, *Phys. Rev. B* 49, 16798 (1994), DOI 10.1103/PhysRevB.49.16798.

---

### `olla-dft layers` — Layer detection by connectivity

**What it answers.** Is the structure layered? How many layers per cell, along which axis are they stacked, what are the basal spacing $d$ and the interlayer gap, and where would the (00l) basal peak fall in a diffractogram?

**Background for non-experts.** The geometry is not eyeballed; bonds are: two atoms are bonded if their distance does not exceed the sum of covalent radii plus a tolerance. The bond network is built respecting periodicity and the connected pieces are separated. A piece that repeats in exactly two directions is a layer; in three, a 3D framework; in one, a chain; in none, a molecule. The dimensionality is read from the "closure vectors": walking the bonds while assigning each atom a cell displacement relative to a root atom, every bond that "does not fit" contributes an integer vector, and the rank of the set of those vectors is the number of periodic directions.

**Formulas.** Bond criterion (`layers.bonds`, with `ase.neighborlist.neighbor_list` and per-atom radii $r_i + \mathrm{tol}/2$):

$$
d_{ij} \le r^{\mathrm{cov}}_i + r^{\mathrm{cov}}_j + \mathrm{tol}, \qquad \mathrm{tol} = 0.45\ \text{Å (default)}
$$

Dimensionality (`layers._components_and_dim`): $\dim = \operatorname{rank}\{\mathbf{d}\}$ with $\mathbf{d} = \mathbf{o}_a + \mathbf{S}_{ab} - \mathbf{o}_b \ne 0$. Spacings (`layers.analyze`):

$$
d_{\mathrm{basal}} = \frac{P}{n_{\mathrm{layers}}}, \qquad P = |\mathbf{a}_{\mathrm{stack}}\cdot\hat{\mathbf{n}}|, \qquad
\mathrm{gap} = \min_k\left(z^{\mathrm{bot}}_{k+1} - z^{\mathrm{top}}_k\right)
$$

Basal reflection (`layers.report`), Bragg: $2\theta = 2\arcsin\!\left(\lambda/(2 d_{\mathrm{basal}}/l)\right)$ for $l = 1, 2, 3$.

- $\hat{\mathbf{n}}$: unit normal to the plane of the two non-stacking cell vectors; $z^{\mathrm{top/bot}}_k$: centre ± thickness/2 of layer $k$ (no van der Waals radii); $\lambda$ in Å.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_layers` → `qekit/core/layers.py: analyze` (`--tol` 0.45 Å).
2. Bonds with ASE; connected components and rank with `np.linalg.matrix_rank` over the closure vectors.
3. Stacking axis: the fractional direction outside the plane spanned by the closure vectors (SVD), taking the cell vector with the largest out-of-plane component.
4. Each layer is rebuilt contiguous (BFS with Cartesian displacements) to measure centre and thickness without cell cuts.
5. `layers.report` prints layers, $d$, gap, period and the 00l reflections with the λ from `--wavelength` (default CuKα = 1.54184 Å, `xrd.wavelength_value`), labelled with the actual radiation name (`xrd.wavelength_name`: "Cu Kα", "Mo Kα1", or "λ dada" if a number was given).
6. With `--slab FILE`, `layers.make_slab` isolates the first layer: it unwraps it along the stacking axis (minimum image in fractional coordinates relative to the first atom, so a layer crossing the cell boundary is not split), replaces the stacking vector by the normal with height thickness + `--vacuum` (20 Å), centres it and `structure.convert` writes it.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Covalent radii | `ase.data.covalent_radii` | library |
| Periodic neighbour list | `ase.neighborlist.neighbor_list("ijS")` | library |
| Wavelength | `xrd.WAVELENGTHS` or a value in Å | CuKa 1.54184, MoKa 0.71073, CoKa 1.79026, … |
| Tolerance | `--tol` | 0.45 Å (`layers.DEFAULT_TOL`) |

**Limits and pitfalls.** No QE calculation: it is pure geometry. If there are no 2D components: "No layers detected … you can try a smaller --tol". The stacking axis and the normal are computed from the FIRST layer; layers with different orientations will not be noticed. `make_slab` only unwraps along the stacking axis (not in-plane), which is all that affects centring and thickness.

**References.** M. Ashton, J. Paul, S. B. Sinnott and R. G. Hennig, *Phys. Rev. Lett.* 118, 106101 (2017), DOI 10.1103/PhysRevLett.118.106101 (topological dimensionality criterion). W. H. Bragg and W. L. Bragg, *Proc. R. Soc. Lond. A* 88, 428 (1913).

---

### `olla-dft xrd` — Simulated powder diffractogram

**What it answers.** Where do the powder X-ray diffraction peaks of this structure appear, with which relative intensity and which hkl indices? Does it look like the measured diffractogram?

**Background for non-experts.** A crystal diffracts X-rays in directions fixed by Bragg's law: each family of planes with spacing $d$ produces a peak at an angle $2\theta$. The intensity depends on how the waves scattered by each atom in the cell interfere (structure factor), on how much each element scatters (atomic scattering factor, which decays with angle), and on geometric factors of the experiment (Lorentz–polarisation). Small crystallites broaden the peaks (Scherrer). Olla-DFT computes all of that and overlays, if given, an experimental diffractogram, to see at a glance whether the structural model is the right one.

**Formulas.** (`xrd.compute`) For each $hkl$ with $|\mathbf{g}| = 1/d$ inside the accessible sphere:

$$
\sin\theta = \frac{\lambda |\mathbf{g}|}{2}, \qquad s^2 = \left(\frac{\sin\theta}{\lambda}\right)^2 = \frac{|\mathbf{g}|^2}{4}
$$

$$
f(s) = Z - 41.78214\, s^2 \sum_{i=1}^{4} a_i\, e^{-b_i s^2}, \qquad f \to f\, e^{-B_{\mathrm{iso}} s^2}
$$

$$
F(hkl) = \sum_j f_j(s)\, e^{2\pi i\,(hkl)\cdot\mathbf{r}_j}, \qquad
I \propto |F|^2 \cdot \frac{1+\cos^2 2\theta}{\sin^2\theta\,\cos\theta}
$$

- $\lambda$: Å; $\mathbf{g} = (hkl)\,\mathbf{B}$ with $\mathbf{B} = (\mathbf{A}^{-1})^T$ without $2\pi$; $Z$: atomic number; $a_i, b_i$: analytical coefficients (data file taken from pymatgen, values from the *International Tables*); $B_{\mathrm{iso}}$: global temperature factor in Å² (`--biso`, 0 by default); $\mathbf{r}_j$: fractional positions.
- Multiplicity: emerges by enumerating ALL hkl and merging those coinciding in $2\theta$ (tolerance 0.02°). Intensities normalised to 100; peaks < 0.1 are dropped.

Profile (`xrd.broaden`), pseudo-Voigt with $\eta = 0.5$ and width $w$ (FWHM in ° 2θ):

$$
y(x) = \sum_p I_p\left[(1-\eta)\, e^{-\frac{(x-x_p)^2}{2\sigma^2}} + \eta\,\frac{1}{1+\left(\frac{x-x_p}{w/2}\right)^2}\right], \quad \sigma = \frac{w}{2\sqrt{2\ln 2}}, \quad
w_{\mathrm{Scherrer}} = \frac{K\lambda}{L\cos\theta},\; K = 0.9
$$

- $L$: crystallite size (`--size` in nm, converted to Å); without `--size`, $w$ = `--fwhm` (0.15°).

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_xrd` → `qekit/modules/xrd.py: compute`. With `--basis conventional` (default) the cell is standardised to the conventional one (`structure.conventional`) so the hkl match the PDF cards; `input` uses the cell as given.
2. Enumerates hkl in the box $|h_i| \le \lceil g_{\max}/|\mathbf{b}_i|\rceil + 1$, filters $g_{\min} \le |\mathbf{g}| \le g_{\max}$ (range `--tt-min` 5°, `--tt-max` 70°).
3. Factors $f_j(s)$ from `qekit/data/atomic_scattering_params.json` (`xrd.scattering_params`), vectorised structure factor, LP, removal of extinct reflections ($I < 10^{-8} I_{\max}$), merging by 2θ and the "most readable" hkl label (`Peak.label`, Friedel-oriented).
4. `xrd.broaden` generates the continuous profile (step 0.02°).
5. `xrd.read_experimental` reads `--exp` (two columns, ≥ 10 rows; subtracts the minimum and normalises to 100).
6. `xrd.report` (12 strongest peaks), `xrd.export` (`XRD.dat`, `XRD_HKL.dat`), `xrd.plot` (experimental offset +105 above the simulated), `--suite` (exchange JSON).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Coefficients $a_i, b_i$ | `qekit/data/atomic_scattering_params.json` | from pymatgen (MIT), *International Tables* |
| $Z$ | `ase.data.atomic_numbers` | library |
| Wavelength | `xrd.WAVELENGTHS` | CuKa 1.54184 Å, CuKa1 1.54056, MoKa 0.71073, CoKa 1.79026, FeKa 1.93735, CrKa 2.29100, AgKa 0.56087 |
| Scherrer constant | `xrd.SCHERRER_K` | 0.9 |
| Conventional cell | spglib (`structure.conventional`) | if it fails, the input cell is used |

**Limits and pitfalls.** There is no Rietveld refinement nor R factor: the comparison with `--exp` is purely visual (overlay). No absorption factor, preferred orientation, anomalous-dispersion correction or Kα1/Kα2 doublet (a single λ). The temperature factor is a single isotropic $B$ for all atoms. Scherrer broadening ignores strain broadening. The $f(s)$ formula is pymatgen's parametrisation (the same as the *International Tables* in the form $Z - 41.78214 s^2\sum a_i e^{-b_i s^2}$), valid for X-rays. With `--basis input` on a primitive cell, "the hkl are NOT those of the PDF card".

**References.** *International Tables for Crystallography*, Vol. C (scattering factors). P. Scherrer, *Nachr. Ges. Wiss. Göttingen* 2, 98 (1918). S. P. Ong et al., *Comput. Mater. Sci.* 68, 314 (2013), DOI 10.1016/j.commatsci.2012.10.028 (pymatgen, origin of the coefficients). B. E. Warren, *X-ray Diffraction*, Dover (1990).

---

### `olla-dft exfoliate` — Exfoliation energy

**What it answers.** How much does it cost to separate one layer from the layered crystal, in J/m² (and meV/Å², meV/atom)? Is it exfoliable?

**Background for non-experts.** The energy of the crystal (per layer) is compared with that of an isolated monolayer in vacuum. The difference per unit area is the exfoliation energy; typical layered materials lie between 0.2 and 0.6 J/m² (graphite ≈ 0.35 J/m² experimental). Interlayer cohesion is mostly van der Waals dispersion, which LDA and PBE describe poorly: without a dispersion correction the number is not comparable with experiment, and the module says so.

**Formulas.** (`exfoliate.report_result`)

$$
E_{\mathrm{exf}} = \frac{E_{\mathrm{mono}} - E_{\mathrm{bulk}}/N_{\mathrm{layers}}}{A}
$$

- $E_{\mathrm{mono}}$, $E_{\mathrm{bulk}}$: total energies in eV; $N_{\mathrm{layers}}$: layers per bulk cell detected by `layers.analyze`; $A = |\mathbf{a}_i\times\mathbf{a}_j|$ of the two non-stacking vectors, in Å². Conversion with `EV_A2_TO_J_M2 = 16.02176634`.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_exfoliate` → `qekit/modules/exfoliate.py: prepare`: `layers.analyze(atoms, tol)` (`--tol` 0.45 Å); without layers, usage error.
2. `layers.make_slab` builds the monolayer (first layer) with `--vacuum` 20 Å.
3. Bulk k-mesh from `sweep.default_grid`; the monolayer's is the same in-plane and 1 along the stacking axis.
4. Two `scf` runs (`bulk/pw.in`, `monocapa/pw.in`; `relax` for the monolayer with `--relax-slab`), both with the same `--vdw` (`grimme-d2`, `grimme-d3`, `DFT-D`, `ts-vdw`, `xdm`, `mbd`).
5. `exfoliate.collect` reads `etot` from both XML files; `exfoliate.report_result` prints the result.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E_{\mathrm{bulk}}$, $E_{\mathrm{mono}}$ | pw.x XML (`etot`) | `qeout.read_xml` |
| $N_{\mathrm{layers}}$, stacking axis | `layers.analyze` | bond connectivity |
| Area $A$ | bulk cell | two non-stacking vectors |
| eV/Å² → J/m² | `exfoliate.EV_A2_TO_J_M2` | 16.02176634 |

**Limits and pitfalls.** It assumes all layers in the cell are equivalent (divides $E_{\mathrm{bulk}}$ by $N_{\mathrm{layers}}$ and uses only the first one). Without `--vdw`: "WITHOUT van der Waals correction: PBE barely binds the layers and LDA binds by error cancellation". With pseudopotentials that look like LDA (name containing `pz`, `lda`, `pw92`) plus Grimme: "combining them counts dispersion twice". If negative: "Almost always means the vdW correction is missing or some calculation is not well converged". There is no bulk relaxation. It writes no `.dat` files (only the on-screen report).

**References.** J. H. Jung, C.-H. Park and J. Ihm, *Nano Lett.* 18, 2759 (2018), DOI 10.1021/acs.nanolett.7b04201 (exfoliation vs. interlayer binding energy). S. Grimme, *J. Comput. Chem.* 27, 1787 (2006), DOI 10.1002/jcc.20495 (D2); S. Grimme et al., *J. Chem. Phys.* 132, 154104 (2010), DOI 10.1063/1.3382344 (D3).

---

### `olla-dft phonons` — DFPT phonons: dispersion, DOS, thermodynamics, IR, Raman and electronic temperature

**What it answers.** What are the vibrational frequencies of the crystal (at Γ or across the Brillouin zone), is it dynamically stable (no imaginary frequencies), what are its harmonic zero-point energy, free energy, entropy and heat capacity, which modes are IR and Raman active (`--raman`) and, with `--tscan`, does a soft mode stabilise as the electronic temperature rises?

**Background for non-experts.** The atoms of a crystal vibrate around their equilibrium positions as if joined by springs. Density-functional perturbation theory (DFPT, what `ph.x` does) computes the stiffness of those springs (the force constants) from the electronic density, without displacing atoms by hand. From there come the frequencies of all vibrational waves (phonons). An "imaginary" frequency (negative in the output) means the structure is not at a minimum: either it was not relaxed well or it is unstable. With the phonon density of states the harmonic thermodynamics is computed: even at 0 K the atoms vibrate (zero-point energy) and on heating more modes get populated (entropy, heat capacity). The Γ modes are the ones an infrared or Raman spectrometer sees.

**Formulas.** Harmonic thermodynamics per cell from the DOS $g(\omega)$ (`phonons.thermodynamics`), with $\epsilon = \hbar\omega$ in eV, $x = \epsilon/k_BT$ (capped at 500), $n = 1/(e^x - 1)$, and $g$ renormalised to $\int g\,d\omega = 3N_{\mathrm{at}}$ (only $\omega > 1$ cm⁻¹):

$$
E_{\mathrm{ZPE}} = \int \tfrac{1}{2}\epsilon\, g\, d\omega, \qquad
F(T) = E_{\mathrm{ZPE}} + k_B T \int \ln\!\left(1 - e^{-x}\right) g\, d\omega
$$

$$
U(T) = \int \left(\tfrac{1}{2} + n\right)\epsilon\, g\, d\omega, \qquad
C_v(T) = k_B \int x^2 e^{x} n^2\, g\, d\omega, \qquad S(T) = \frac{U - F}{T}
$$

- $k_B$ = `KB_EV` = 8.617333262e-5 eV/K; cm⁻¹ → eV with `CM1_TO_EV` = 1.239841984e-4; cm⁻¹ → THz with `CM1_TO_THZ` = 0.0299792458. Trapezoidal integrals. $T$ = 0…1000 K in steps of 10.

Stokes Raman spectrum (`phonons.raman_spectrum`) from the activity $A$ (Å⁴/amu) of `dynmat.x`:

$$
I(\omega) \propto \frac{(\omega_L - \omega)^4}{\omega}\,[n(\omega,T)+1]\,A(\omega), \qquad \omega_L = \frac{10^7}{\lambda_{\mathrm{laser}}[\mathrm{nm}]}\ \mathrm{cm^{-1}}
$$

convolved with Lorentzians of FWHM 5 cm⁻¹ at $T$ = 300 K; $\omega \le 1$ cm⁻¹ is excluded. Electronic temperature (`tphonons.degauss_de_T`):

$$
\mathrm{degauss} = k_B T, \qquad k_B = 6.333621\times10^{-6}\ \mathrm{Ry/K}, \qquad \text{smearing = fermi-dirac}
$$

Stabilisation temperature (`tphonons.temperatura_de_estabilizacion`): linear interpolation of the $T$ at which the softest mode (minimum of the frequencies with $|\omega| > 10$ cm⁻¹) crosses from negative to non-negative.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_phonons` → `qekit/modules/phonons.py: prepare`. The structure is taken to the standardised primitive cell (`structure.primitive`). It writes `scf.in` (`conv_thr = 1e-12`, k-mesh from `kspacing`), `ph.in` (`tr2_ph = 1e-14`, `fildyn='dyn'`; `epsil=.true.` if `--insulator` or `--raman`; `lraman=.true.` and `trans=.true.` with `--raman`; `ldisp` with `nq` mesh = `--qgrid` or `kgrid_from_spacing(atoms, 0.6)`).
2. Mesh mode: `q2r.in` (`zasr='simple'`, `flfrc='fuerzas.fc'`), `matdyn_band.in` (seekpath path via `kpoints.get_kpath`, 30 points per segment, `q_in_band_form`, `q_in_cryst_coord`), `matdyn_dos.in` (`dos=.true.`, 12×12×12 mesh, `fldos='fonones.dos'`). Γ mode (`--gamma` or `--raman`): `dynmat.in` (`asr='simple'`).
3. `--raman` requires norm-conserving pseudopotentials (`p["type"] == "NC"`), otherwise a usage error.
4. `--run`: `runner.run_all` runs `pw.x`; `phonons.run_chain` executes `ph.x` → (`dynmat.x` | `q2r.x` → `matdyn.x` ×2), skipping steps whose `.out` already says `JOB DONE`.
5. `phonons.collect`: at Γ it reads the table `# mode [cm-1] [THz] IR [Raman depol]` from `dynmat.out` (`read_dynmat_table`); on a mesh it reads `bandas.freq` (`_read_flfrq`, `&plot` format, q in Cartesian 2π/alat) and `fonones.dos`.
6. `phonons.report` / `report_gamma_activities` (mutual exclusion rule, depolarisation 0.75), `phonons.thermodynamics`, `phonons.export` (`FONONES_GAMMA.dat` or `FONONES_BANDAS.dat`, `FONONES_DOS.dat`, `FONONES_TERMO.dat`), and `phonons.plot` only if `phonons.has_dispersion(run)` (there are `band_freqs` and `qdist`, i.e. it is not a Γ-only run).
7. `--tscan T1,T2,...` → `qekit/cli.py: _cmd_phonons_tscan` → `qekit/modules/tphonons.py: prepare`: one full chain per temperature in `T00300/` etc., with `insulator=False`, `smearing='fermi-dirac'` and `degauss = k_B T`; `tphonons.collect`, `report` (table of imaginary modes per T, monotonicity, $T_{\mathrm{stab}}$), `export` (`FONONES_T.dat`), `plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Γ frequencies, IR, Raman, depol | `dynmat.out` (`# mode` table) | `phonons.read_dynmat_table`; IR in (D/Å)²/amu, Raman in Å⁴/amu |
| Dispersion | `bandas.freq` from `matdyn.x` | `phonons._read_flfrq` |
| Phonon DOS | `fonones.dos` from `matdyn.x` | `np.loadtxt`, states/cm⁻¹ |
| High-symmetry path | seekpath (`kpoints.get_kpath`) | labels and discontinuities |
| $k_B$, cm⁻¹→eV, cm⁻¹→THz | `phonons.KB_EV`, `CM1_TO_EV`, `CM1_TO_THZ` | CODATA |
| $k_B$ in Ry/K | `tphonons.KB_RY` | 6.333621e-6 |
| Imaginary threshold | literal −5 cm⁻¹ (`phonons.report`), `tphonons.UMBRAL_IMAGINARIO` = 10 | numerical noise below that |

**Limits and pitfalls.** It is harmonic: no thermal expansion or anharmonicity (for that, `qha`). The `prepare` docstring defaults to `insulator=True`, but the CLI passes `args.insulator`, which is `False` unless `--insulator`: by default the scf uses smearing and `epsil` is NOT enabled (no LO–TO splitting). It warns: "there are imaginary (negative) frequencies. Either the structure is not relaxed, or it is unstable at Γ" and "the structure must be relaxed (vc-relax) with these same cutoffs". The absolute scale of the `matdyn` DOS depends on the mesh; it is renormalised to $3N$ in the thermodynamics. With `--raman` (which forces Γ mode even without `--gamma`) no dispersion is drawn: the CLI decides with `has_dispersion(run)`, not with the flag. The thermodynamics ignores $\omega \le 1$ cm⁻¹. In `--tscan`, "this is ELECTRONIC temperature. The ions stay still"; if the number of imaginary modes does not decrease monotonically: "Usually means the k-mesh is not converged". Only `smearing='fermi-dirac'` corresponds to a real temperature.

**References.** S. Baroni, S. de Gironcoli, A. Dal Corso and P. Giannozzi, *Rev. Mod. Phys.* 73, 515 (2001), DOI 10.1103/RevModPhys.73.515 (DFPT). M. Lazzeri and F. Mauri, *Phys. Rev. Lett.* 90, 036401 (2003), DOI 10.1103/PhysRevLett.90.036401 (Raman via 2n+1). D. Porezag and M. R. Pederson, *Phys. Rev. B* 54, 7830 (1996) (Raman intensities). A. A. Maradudin et al., *Theory of Lattice Dynamics in the Harmonic Approximation*, Academic Press (1971).

---

### `olla-dft qha` — Quasi-harmonic approximation

**What it answers.** How does the crystal expand with temperature ($V(T)$, $\alpha(T)$, $a(T)$), what is the Grüneisen parameter, $C_p$ versus $C_v$ and $B(T)$?

**Background for non-experts.** The harmonic approximation gives no expansion: if the frequencies do not depend on volume, the minimum of the free energy does not move. The QHA keeps the harmonic modes but lets their frequencies change with VOLUME. At each temperature the static energy $E(V)$ and the vibrational free energy $F_{\mathrm{vib}}(V,T)$ are added, and the volume minimising the sum is found. On heating, the vibrational free energy favours larger volumes (modes soften) and the minimum shifts: that is thermal expansion.

**Formulas.** (`qha.f_vib`, `qha.cv_modos`, `qha.run`) For each volume $V_i$ with its modes $\omega_k$ (cm⁻¹, only $\omega > 1$), $\epsilon_k = \hbar\omega_k$:

$$
F(V_i,T) = E(V_i) + \frac{1}{N_{\mathrm{cells}}}\left[\sum_k \tfrac{1}{2}\epsilon_k + k_B T \sum_k \ln\!\left(1 - e^{-\epsilon_k/k_BT}\right)\right]
$$

Minimum by a local parabola (up to 2 points on each side of the sampled minimum): $V(T) = -b/2a$ (clipped to the range), $B(T) = 2a\,V(T)\times 160.21766208$ GPa.

$$
\alpha(T) = \frac{1}{V}\frac{dV}{dT}\ (\texttt{np.gradient}), \qquad
C_v = k_B\sum_k x_k^2 \frac{e^{x_k}}{(e^{x_k}-1)^2}\ \text{(interpolated at } V(T)), \qquad
C_p = C_v + \alpha^2 B V T
$$

$$
\gamma = -\frac{d\ln\langle\omega\rangle}{d\ln V}\ \text{(straight line over } \ln V\text{, } \ln\bar\omega\text{)}, \qquad
a(T) = \begin{cases} (V \cdot V_{\mathrm{conv}}/V_{\mathrm{prim}})^{1/3} & \text{with } \texttt{--structure} \\ V_{\mathrm{prim}}^{1/3} & \text{without it} \end{cases}\ (\texttt{--cubic})
$$

- $E$ in eV, $V$ in Å³, $C_v$, $C_p$ in meV/K per cell, $\alpha$ in K⁻¹, $B$ in GPa; `KB_EV` = 8.617333262e-5, `CM1_EV` = 1.239841984e-4.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_qha` reads a TABLE (`data`): columns $V$ (Å³), $E$ (eV), $\omega_1, \omega_2, \ldots$ (cm⁻¹) per volume; values ≤ −1000 are discarded as padding.
2. With `--structure CIF`, `qha.factor_convencional` counts how many primitive cells fit in the conventional one ($N_{\mathrm{conv}}/N_{\mathrm{prim}}$ via spglib: 4 for fcc/diamond, 2 for bcc, 1 for simple cubic) and `qha.es_cubico` (space group ≥ 195) turns on `--cubic` automatically.
3. `qekit/modules/qha.py: run` with $T$ = 0 … `--tmax` (1000) in steps of `--dt` (5), `--natoms` (1), `--cells` (primitive cells per supercell of the modes, 1), `--cubic`, `factor_conv`.
4. Warnings if there are < 4 volumes or frequencies < −5 cm⁻¹ at some volume.
5. `qha.report` (at `--temp` 300 K; $a(T)$ labelled "lattice parameter (conventional cell)" or "V_prim^(1/3) (NOT the conventional lattice parameter)" according to `QHAResult.a_convencional`), `qha.export` (`QHA.dat`), `qha.plot` (V, α, $C_v$/$C_p$).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E(V)$ and $\omega_k(V)$ | user table (`data`) | from `eos` + `phonons` (or `mlip phonons`) per volume; Olla-DFT does not generate it |
| $V_{\mathrm{conv}}/V_{\mathrm{prim}}$ and cubicity | `--structure` via spglib | `qha.factor_convencional`, `qha.es_cubico` |
| Constants | `qha.KB_EV`, `CM1_EV`, `EV_A3_GPA` | CODATA; 160.21766208 |

**Limits and pitfalls.** It launches no QE calculation: it receives the table. It uses discrete frequencies (one set of modes per volume), not a DOS: the thermodynamics is done on the list it is given, so the quality depends on the supercell/mesh of those modes. The Grüneisen parameter is an average over the mean frequency, not mode by mode. The QHA "holds up to ~half the melting temperature". Without `--structure`, $a(T)$ is only $V_{\mathrm{prim}}^{1/3}$ and the report warns: "In fcc, bcc or diamond that is NOT the conventional lattice parameter (they differ by 4^(1/3) or 2^(1/3)). Pass the structure with --structure". With a single temperature, $\alpha$ is NaN and a warning is issued.

**References.** A. Togo, L. Chaput, I. Tanaka and G. Hug, *Phys. Rev. B* 81, 174301 (2010), DOI 10.1103/PhysRevB.81.174301. G. Grimvall, *Thermophysical Properties of Materials*, North-Holland (1999). P. Pavone et al., *Phys. Rev. B* 48, 3156 (1993) (Si, negative expansion).

---

### `olla-dft derived` — Debye, sound velocities and Slack from the $C_{ij}$

**What it answers.** From an already computed elastic matrix: what are the density, the sound velocities, the elastic Debye temperature, the Poisson ratio, an approximate Grüneisen parameter and an estimate of the lattice thermal conductivity?

**Background for non-experts.** The stiffness of a solid sets how fast sound waves travel through it, and that speed sets the fastest possible vibration; the Debye temperature is that maximum frequency expressed in kelvin. With it and an anharmonicity parameter (Grüneisen), the Slack model estimates how much heat the lattice conducts. Everything is post-processing: it costs no new calculation.

**Formulas.** (`derived.density`, `sound_velocities`, `debye_from_velocity`, `poisson_ratio`, `gruneisen_from_poisson`, `slack`, `cubic_directional`)

$$
\rho = \frac{\sum_i m_i}{V}, \qquad
v_l = \sqrt{\frac{B + \tfrac{4}{3}G}{\rho}}, \qquad v_t = \sqrt{\frac{G}{\rho}}, \qquad
v_m = \left[\frac{1}{3}\left(\frac{2}{v_t^3} + \frac{1}{v_l^3}\right)\right]^{-1/3}
$$

$$
\Theta_D = \frac{\hbar}{k_B}\left(6\pi^2 n\right)^{1/3} v_m, \qquad
\nu = \frac{3B - 2G}{2(3B+G)}, \qquad
\gamma = \frac{3(1+\nu)}{2(2-3\nu)}
$$

$$
\kappa_L = A\,\frac{\bar M\,\Theta_D^3\,\delta}{\gamma^2\, n_{\mathrm{at}}^{2/3}\, T}, \qquad
A = \frac{3.1\times10^{-6}}{1 - 0.514/\gamma + 0.228/\gamma^2}, \qquad
v_L^{[100]} = \sqrt{C_{11}/\rho},\ v_T^{[100]} = \sqrt{C_{44}/\rho}
$$

- $B$, $G$: Hill averages in GPa (×10⁹ to Pa); $\rho$ in kg/m³ (masses in amu × 1.66053906660e-27, $V$ in Å³ × 1e-30); $n$: atoms per m³; $\hbar$ = 1.054571817e-34 J·s, $k_B$ = 1.380649e-23 J/K; $\bar M$: mean mass in amu; $\delta = (V/n_{\mathrm{at}})^{1/3}$ in Å; $T$ in K (`--temp`, 300); $\kappa_L$ in W/(m·K).

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_derived` loads the structure (masses and volume) and `--cij` (`ELASTIC_C.dat` from `elastic`, 6×6 matrix).
2. `elastic.moduli` → $B_H$, $G_H$; `qekit/modules/derived.py: analyze` computes everything above.
3. `derived.cubic_directional` prints $v_L$, $v_T$ along [100] only if the structure is cubic according to spglib (`elastic.crystal_family`) or the tensor has cubic form (`derived.is_cubic_tensor`: $C_{11}=C_{22}=C_{33}$, $C_{12}=C_{13}=C_{23}$, $C_{44}=C_{55}=C_{66}$ and zeros elsewhere, with a 5 % or 2 GPa tolerance).
4. `derived.report`; `derived.export` writes `DERIVED.dat`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $C_{ij}$ | `ELASTIC_C.dat` (`--cij`) | `np.loadtxt`, GPa |
| Masses and volume | structure (ASE `get_masses`, `get_volume`) | amu, Å³ |
| $\hbar$, $k_B$, amu | `derived.HBAR`, `KB`, `AMU` | CODATA 2018 |
| Slack prefactor | literal 3.1e-6 and correction $(1 - 0.514/\gamma + 0.228/\gamma^2)$ | Slack / Julian |

**Limits and pitfalls.** The $\Theta_D$ is the ELASTIC one (low-temperature acoustic limit): "The one from the phonon DOS uses the whole spectrum and gives another number; they are not the same quantity" (`derived.debye_from_dos`, $\Theta_D = (\hbar/k_B)\sqrt{5\langle\omega^2\rangle/3}$, exists and is used by `crosscheck`, not by this command). The Grüneisen parameter "comes from an empirical correlation with the Poisson ratio" (Belomestnykh) and Slack "is an order-of-magnitude estimate". Negative Poisson: auxetic-material warning. The Slack κ is labelled with the temperature actually used (`Termoelastico.T`, key `kappa_Slack_<T>K` in `DERIVED.dat`). If $G \le 0$ there are no velocities.

**References.** O. L. Anderson, *J. Phys. Chem. Solids* 24, 909 (1963), DOI 10.1016/0022-3697(63)90067-2 (elastic $\Theta_D$). G. A. Slack, *Solid State Phys.* 34, 1 (1979); D. T. Morelli and G. A. Slack, in *High Thermal Conductivity Materials*, Springer (2006) (prefactor with Julian's correction). V. N. Belomestnykh and E. P. Tesleva, *Tech. Phys.* 49, 1098 (2004) (Grüneisen–Poisson).

---

### `olla-dft thermochem` — ZPE, entropy and free energy

**What it answers.** How much must be added to a DFT energy (at 0 K, without vibrations) to obtain a free energy $G(T,p)$ comparable with experiment, for a solid, an adsorbate, an ideal gas or a transition state?

**Background for non-experts.** A DFT energy is electronic and at 0 K. What is measured is a free energy at the laboratory temperature and pressure. Between the two there are three terms: the zero-point energy (modes vibrate even at 0 K), the enthalpy correction (modes get populated on heating) and the entropic term $-TS$, which for a gas molecule includes the translational (Sackur–Tetrode) and rotational (rigid rotor) entropies and can amount to about 1 eV at 500 K. Forgetting it can flip the sign of an adsorption energy.

**Formulas.** (`thermochem.zpe`, `H_vib`, `S_vib`, `Cv_vib`, `S_traslacional`, `S_rotacional`, `corregir`) With $\epsilon_k = h c\,\tilde\nu_k$, $x_k = \epsilon_k/k_BT$ (capped at 500):

$$
E_{\mathrm{ZPE}} = \tfrac{1}{2}\sum_k \epsilon_k, \quad
H_{\mathrm{vib}} = \sum_k \frac{\epsilon_k}{e^{x_k}-1}, \quad
S_{\mathrm{vib}} = k_B\sum_k\left[\frac{x_k}{e^{x_k}-1} - \ln\!\left(1-e^{-x_k}\right)\right], \quad
C_v = k_B\sum_k \frac{x_k^2 e^{x_k}}{(e^{x_k}-1)^2}
$$

$$
S_{\mathrm{trans}} = k_B\left[\ln\!\left(\frac{V}{\Lambda^3}\right) + \tfrac{5}{2}\right], \quad V = \frac{k_BT}{p}, \quad \Lambda = \frac{h}{\sqrt{2\pi m k_B T}}
$$

$$
S_{\mathrm{rot}}^{\mathrm{linear}} = k_B\left[\ln\frac{T}{\sigma\Theta_r} + 1\right], \qquad
S_{\mathrm{rot}}^{\mathrm{non\,linear}} = k_B\left[\tfrac{1}{2}\ln\frac{\pi T^3}{\sigma^2\Theta_A\Theta_B\Theta_C} + \tfrac{3}{2}\right], \qquad \Theta_i = \frac{\hbar^2}{2 I_i k_B}
$$

$$
G - E_{\mathrm{DFT}} = E_{\mathrm{ZPE}} + H_{\mathrm{corr}} - TS, \qquad
H_{\mathrm{corr}}^{\mathrm{gas}} = H_{\mathrm{vib}} + \left(\tfrac{3}{2} + n_{\mathrm{rot}} + 1\right)k_BT, \qquad S_{\mathrm{elec}} = k_B\ln(\text{multiplicity})
$$

- $\tilde\nu$ in cm⁻¹ (`C_CM` = 2.99792458e10 cm/s, `H_EVS` = 4.135667696e-15 eV·s); $m$: molecular mass (amu → kg); $p$ in Pa (`--pressure` in bar); $\sigma$: symmetry number (`--symmetry`, 1); $I_i$: principal moments of inertia (amu·Å² → kg·m²); linear if $I_1 < 10^{-3} I_3$; $n_{\mathrm{rot}}$ = 1 (linear), 1.5 (non-linear), 0 (atom). Adsorption energy (`thermochem.adsorcion`): $E_{\mathrm{ads}} = E_{\mathrm{slab+ads}} - E_{\mathrm{slab}} - nE_{\mathrm{gas}}$; $G_{\mathrm{ads}} = E_{\mathrm{ads}} + G^{\mathrm{corr}}_{\mathrm{ads}} - n\,G^{\mathrm{corr}}_{\mathrm{gas}}$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_thermochem` reads the frequencies (`_leer_frecuencias`: a one- or multi-column file —last column—, or an inline list) and, for gas, the structure with `ase.io.read`.
2. `qekit/modules/thermochem.py: limpiar_frecuencias`: separates imaginary ones ($\tilde\nu < -1$), discards $|\tilde\nu| \le 1$ (residual translations/rotations), raises soft modes to the `--floor` (e.g. 100 cm⁻¹) and emits warnings according to `--phase` (`solido`, `adsorbato`, `gas`, `transicion`).
3. `thermochem.corregir` sums the terms at `--temp` (298.15 K) and `--pressure` (1 bar), with `--multiplicity`.
4. `thermochem.report` (with `G(T)` if `--energy` is given); with `-o`, it writes `TERMOQUIMICA.txt`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Frequencies | file (`FONONES_GAMMA.dat`, last column) or list | `_leer_frecuencias` |
| Masses and geometry (gas) | `--structure` via ASE | moments of inertia |
| $h$, $k_B$, $c$, amu, $\hbar$ | `thermochem.H_EVS`, `KB_EV`, `C_CM`, `AMU_KG`, `HBAR_JS`, `KB_J` | CODATA |
| Soft-mode floor | `--floor` (`PISO_BLANDO` = 100 only as a reference) | without `--floor` nothing is raised |

**Limits and pitfalls.** Verified in the tests against NIST (H₂O, N₂, CH₄ to 0.5 %). Transition state without an imaginary mode or with more than one: explicit warning. Imaginary modes at a minimum: "the structure is NOT a minimum … They are excluded from the sums". Raised modes: "it is a CORRECTION, not a calculation: say so if you publish". Gas with a number of modes ≠ $3N-6$ (or $3N-5$): "they are counted twice with the translational and rotational terms". In the gas phase no anharmonicity or internal rotors are included; the solid carries no $pV$ term. `adsorcion` applies no vibrational corrections to the clean slab (assumes they do not change). The CLI `thermochem` command does not expose `adsorcion` (used from `adsorb`/API).

**References.** D. A. McQuarrie, *Statistical Mechanics*, University Science Books (2000). C. J. Cramer, *Essentials of Computational Chemistry*, Wiley (2004), ch. 10. O. Sackur, *Ann. Phys.* 36, 958 (1911); H. Tetrode, *Ann. Phys.* 38, 434 (1912).

---

### `olla-dft md` — Analysis of a molecular-dynamics trajectory

**What it answers.** From a `pw.x` output with `calculation='md'`: what structure does the system have (g(r), coordination numbers), do the atoms diffuse (MSD and diffusion coefficient $D$) and what is its vibrational spectrum (VDOS) including temperature and anharmonicity?

**Background for non-experts.** A molecular dynamics run is a "movie" of the atoms moving. Three functions summarise the movie: the radial distribution function $g(r)$ says how many neighbours there are at each distance (its first peak is the bond length, its area the coordination number); the mean-square displacement (MSD) says whether atoms move away from where they were (if it grows linearly, they diffuse; if it flattens, they only vibrate); and the Fourier transform of the velocity autocorrelation gives the frequencies at which they vibrate. Before trusting any of them, the initial equilibration stretch must be discarded and the trajectory must be long enough.

**Formulas.** (`dynamics.rdf`, `coordinacion`, `msd`, `difusion`, `vdos`) With minimum-image distances in fractional coordinates:

$$
g(r) = \frac{h(r)}{N_{\mathrm{steps}}\,\frac{N(N-1)}{2}\,\frac{4\pi r^2\,\Delta r}{V}}, \qquad
g_{AB}(r) = \frac{h_{AB}(r)}{N_{\mathrm{steps}}\,N_{\mathrm{pairs}}\,\frac{4\pi r^2\Delta r}{V}}, \qquad
n_{\mathrm{coord}} = \int_0^{r_{\min}} 4\pi r^2 \rho\, g(r)\, dr
$$

$$
\mathrm{MSD}(\tau) = \left\langle |\mathbf{r}_i(t+\tau) - \mathbf{r}_i(t)|^2 \right\rangle_{i,t}, \qquad
\mathrm{MSD} = 6 D \tau + b \;\Rightarrow\; D = \frac{m}{6}\times 10^{-1}\ [\mathrm{cm^2/s}]
$$

$$
C(t) = \frac{\sum_{i,\alpha}\langle v_{i\alpha}(0)v_{i\alpha}(t)\rangle}{C(0)}\ (\text{via FFT}), \qquad
\mathrm{VDOS}(\tilde\nu) = \left|\mathcal{F}\{C(t)\,w_{\mathrm{Hann}}(t)\}\right|, \quad \tilde\nu = \frac{f[\mathrm{fs^{-1}}]\times 10^{15}}{c[\mathrm{cm/s}]}
$$

- $r_{\max}$ = half the shortest cell edge (or `--rmax` if smaller); `--bins` 200; $N_{\mathrm{pairs}} = N_AN_B$ or $N_A(N_A-1)/2$; $\rho = N/V$; $r_{\min}$: first local minimum with $g<1$ after the first maximum; the $D$ fit uses only the 20–80 % stretch of the lags (lags up to $n/2$), slope $m$ in Å²/fs; velocities by `np.gradient` of the UNWRAPPED positions (no periodic jumps), not mass-weighted.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_md` → `qekit/modules/dynamics.py: leer_md`: reads from `pw.out` (text) the `ATOMIC_POSITIONS` blocks (alat, bohr, angstrom or crystal → Å), the cell from `a(i) = (...)`·alat or `CELL_PARAMETERS`, `temperature = … K`, `!    total energy` and `Time step = … femto-seconds`; discards `--skip` steps.
2. `dynamics.analizar`: `rdf`, `coordinacion` per pair, `desdoblar` + `msd` (total and per species), `difusion`, `vdos` (≥ 8 steps), temperature drift between halves.
3. `dynamics.report`, `dynamics.export` (`MD_RDF.dat`, `MD_MSD.dat`, `MD_VDOS.dat`, `MD.txt`), `dynamics.plot` (three panels).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Positions per step | `pw.out` (`ATOMIC_POSITIONS`) | `dynamics._leer_marcos`; units detected |
| Cell | `pw.out` (`a(1..3) = (...)` × alat, or `CELL_PARAMETERS`) | assumed constant |
| Time step | `pw.out` (`Time step = … a.u., X femto-seconds`) | 1 fs if absent |
| Temperature and energy | `pw.out` (`temperature =`, `!    total energy`) | K; Ry → eV with 13.605693122994 |
| bohr → Å | `dynamics.BOHR_A` | 0.529177210903 |

**Limits and pitfalls.** Constant cell: not usable for `vc-md`. If $g(r)$ is empty up to the cutoff: "the cell is too small to extract structure from it: build a supercell". Less than 2 ps: "Good for looking at the structure, not for a diffusion coefficient". $R^2 < 0.95$: "the MSD is NOT linear: no diffusion, or not enough time". Temperature drift > 15 %: "it is still equilibrating: discard more steps with --skip". The coordination number is integrated up to the FIRST MINIMUM, "a convention, not a measurement". The VDOS is not mass-weighted (it is not the phonon DOS) and takes the modulus of the spectrum, not the real part; its resolution is $1/(N_{\mathrm{steps}}\,\Delta t)$. `KB_RY` in `dynamics.py` is unused.

**References.** M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*, Oxford (2017). A. Einstein, *Ann. Phys.* 17, 549 (1905). J.-P. Hansen and I. R. McDonald, *Theory of Simple Liquids*, Academic Press (2013).

---

### `olla-dft kappa` — Lattice thermal conductivity (fc3 + BTE with phono3py)

**What it answers.** How much heat does the crystal lattice conduct, $\kappa_L(T)$ in W/(m·K), which exponent does it follow with temperature and which phonon mean free paths carry the heat (to know whether nanostructuring helps)?

**Background for non-experts.** In a perfectly harmonic crystal a phonon would travel forever and the conductivity would be infinite. What makes it finite is that a phonon can split into two or two can merge into one: that is allowed by the cubic term of the energy (the third-order force constants, fc3). They are obtained by displacing two atoms at a time in a supercell and computing the forces; with them, the phonon Boltzmann equation (in the relaxation-time approximation, RTA) gives $\kappa$. It is expensive because the number of configurations grows fast with the supercell. Olla-DFT allows the forces to be computed with `pw.x` (the real calculation) or with a learned potential (MACE, etc.) for exploration.

**Formulas.** They are solved by phono3py (`kappa.resolver`); Olla-DFT post-processes:

$$
\kappa_L^{\alpha\beta} = \frac{1}{NV}\sum_\lambda C_\lambda\, v_\lambda^\alpha v_\lambda^\beta\, \tau_\lambda, \qquad \tau_\lambda = \frac{1}{2\Gamma_\lambda}, \qquad \Lambda_\lambda = |\mathbf{v}_\lambda|\,\tau_\lambda
$$

$$
\bar\kappa = \frac{\kappa_{xx}+\kappa_{yy}+\kappa_{zz}}{3}, \qquad
\kappa \propto T^{-n}\ (n \text{ by a straight line in } \ln\kappa\text{–}\ln T,\ T \ge 200\ \mathrm{K}), \qquad
\kappa_{\mathrm{cum}}(\Lambda) = \frac{\sum_{\lambda:\Lambda_\lambda<\Lambda} w_\lambda C_\lambda \tfrac{|\mathbf{v}_\lambda|^2}{3}\tau_\lambda}{\sum_\lambda w_\lambda C_\lambda \tfrac{|\mathbf{v}_\lambda|^2}{3}\tau_\lambda}
$$

- $\Gamma_\lambda$: linewidth (THz) from phono3py; $\mathbf{v}_\lambda$: group velocity (THz·Å); $C_\lambda$: modal heat capacity; $w_\lambda$: q-point weight; $\Lambda$ in Å (reported in nm). Modes with $\Gamma = 0$ (acoustic at Γ) are discarded.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_kappa` → `qekit/modules/kappa.py: preparar`: `Phono3py(..., supercell_matrix=--dim (2x2x2), phonon_supercell_matrix=--dim-fc2, primitive_matrix="auto", symprec=1e-5)` and `generate_displacements(distance=--distance 0.03 Å)`.
2. `kappa.configuraciones` converts the displaced supercells to ASE (fc3 and, if any, fc2).
3. Forces: (a) `--model mace|chgnet|m3gnet` → `kappa.fuerzas_mlip`; (b) without `--model` → `kappa.escribir_inputs` writes one `scf` per configuration in `fc3/dNNNN/pw.in` (and `fc2/`), `conv_thr = 1e-10`, mesh from `--kspacing` 0.35 Å⁻¹, `occupations='fixed'` unless `--metal` (smearing), plus `correr.sh`; it refuses above 150 configurations without `--force`; (c) `--collect` → `kappa.leer_fuerzas` reads `<forces>` from each XML (Ha/bohr → eV/Å) and requires ALL of them.
4. `kappa.resolver`: `produce_fc3`, `produce_fc2`, symmetrisation, `mesh_numbers = --mesh (13)`, `init_phph_interaction`, `run_thermal_conductivity(temperatures=--temps 100:800:8, is_isotope=--isotopes, boundary_mfp=--grain µm ×1e4 Å or 1e6)`.
5. `kappa.recoger` stores κ (Voigt 6), Γ, velocities, $C_\lambda$, weights; `kappa.report`, `export` (`KAPPA.dat`, `KAPPA_recorrido.dat`, `KAPPA.txt`), `plot` (κ(T) log-log with a $T^{-1}$ guide; cumulative vs Λ).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Forces | pw.x XML (`output/forces`) or MLIP potential | `qeout.read_xml` / `mlip.calculator` |
| fc2, fc3, Γ, v, C, κ | `phono3py` library | `Phono3py.thermal_conductivity` (RTA) |
| Mean-free-path grid | `kappa.RECORRIDOS` | `np.logspace(0, 7, 141)` Å |
| Isotopes | phono3py (natural abundances) | `--isotopes` |

**Limits and pitfalls.** "It is RTA, not the exact solution of the Boltzmann equation. RTA underestimates κ (≈10-15 % in silicon, much more in graphene or diamond)". "Only three-phonon scattering is included". Without `--isotopes`: "Natural silicon conducts ~10 % less than isotopically pure". With a learned potential: "the absolute value may be far off: with small MACE-MP silicon gives ~51 W/mK at 300 K where experiment is ~140". Supercell ≤ 8 cells: "κ must converge in the supercell size AND the q-mesh at the same time". By default the fc2/fc3 scf runs use `occupations='fixed'` ("the right choice for insulators"); for a metal `--metal` must be passed, or the scf runs will not converge. There is no option to write/read phono3py's `fc2.hdf5`/`fc3.hdf5`: every `--collect` rebuilds everything.

**References.** A. Togo, L. Chaput and I. Tanaka, *Phys. Rev. B* 91, 094306 (2015), DOI 10.1103/PhysRevB.91.094306 (phono3py). J. M. Ziman, *Electrons and Phonons*, Oxford (1960). L. Lindsay, D. A. Broido and T. L. Reinecke, *Phys. Rev. B* 87, 165201 (2013) (RTA vs. exact solution).

---

### `olla-dft elph` — Electron–phonon coupling: λ, ω_log, Tc and τ

**What it answers.** How strongly do electrons couple to phonons ($\lambda$), what is the Allen–Dynes superconducting critical temperature with its strong-coupling corrections, and what is the phonon-limited relaxation time $\tau(T)$ that the CRTA transport lacks?

**Background for non-experts.** Electrons moving through a metal collide with the lattice vibrations; how much they collide is measured by $\lambda$, a dimensionless number. In a conventional superconductor that same coupling is what pairs the electrons, and the Allen–Dynes formula turns $\lambda$ and a typical phonon frequency ($\omega_{\log}$) into a critical temperature. The same $\lambda$ gives the time between collisions at high temperature, $\tau$. `ph.x` computes the coupling for several numerical broadenings; the good value is the one at the "plateau", where it stops depending on the broadening.

**Formulas.** (`elph.lambda_de_a2F`, `omega_log_de_a2F`, `omega_2`, `factores_correccion`, `allen_dynes`, `tau_elph`)

$$
\lambda = 2\int \frac{\alpha^2F(\omega)}{\omega}\,d\omega, \qquad
\omega_{\log} = \exp\!\left[\frac{2}{\lambda}\int \ln\omega\,\frac{\alpha^2F(\omega)}{\omega}\,d\omega\right], \qquad
\bar\omega_2 = \left[\frac{2}{\lambda}\int \omega\,\alpha^2F(\omega)\,d\omega\right]^{1/2}
$$

$$
T_c = f_1 f_2\,\frac{\omega_{\log}}{1.2}\exp\!\left[\frac{-1.04(1+\lambda)}{\lambda - \mu^*(1+0.62\lambda)}\right], \qquad
f_1 = \left[1 + \left(\frac{\lambda}{\Lambda_1}\right)^{3/2}\right]^{1/3}, \quad \Lambda_1 = 2.46(1+3.8\mu^*)
$$

$$
f_2 = 1 + \frac{(r-1)\lambda^2}{\lambda^2 + \Lambda_2^2}, \quad r = \frac{\bar\omega_2}{\omega_{\log}}, \quad \Lambda_2 = 1.82(1+6.3\mu^*)\,r, \qquad
\frac{1}{\tau} = \frac{2\pi\lambda k_B T}{\hbar}
$$

- $\omega$ in THz in `a2F.dos*`; $\omega_{\log}$, $\bar\omega_2$ in K (`THZ_K` = 47.9924 K/THz); $\mu^*$ = 0.10, 0.13, 0.16 (a range, not computed); $T_c$ = 0 if the denominator ≤ 0; $\hbar$ = `HBAR_EVS` = 6.582119569e-16 eV·s, $k_B$ = 8.617333262e-5 eV/K; $\tau$ in s. Plateau (`elph.plato`): longest stretch of ≥ 3 consecutive λ that do not differ by more than 5 % from the first; its midpoint.

**How Olla-DFT computes it.**
1. Preparation (`qekit/cli.py: _cmd_elph` without `--collect` → `qekit/modules/elph.py: prepare`): writes `1_scf.in` (mesh `--kgrid` or `kspacing`), `2_nscf.in` with `la2F = .true.` and mesh `--kgrid-nscf` (default $q_i\cdot\max(2, \lceil 2k_i/q_i\rceil)$, a multiple of the q-mesh), and `3_ph.in` with `electron_phonon='interpolated'`, `el_ph_sigma = --sigma (0.005 Ry)`, `el_ph_nsigma = --nsigma (10)`, `fildvscf='dvscf'`, `tr2_ph = 1e-12`, `ldisp` with `--qgrid` (2x2x2). Smearing `methfessel-paxton`, `--degauss` 0.02 Ry, `conv_thr = 1e-10`.
2. The user runs `pw.x` ×2, `ph.x` and, optionally, `lambda.x` (Olla-DFT has `elph.build_lambda_input` but the CLI never writes it).
3. `--collect`: `elph.leer_elph_ph` reads from `ph.out` the broadenings (`Gaussian Broadening: X Ry`) and `DOS = … states/spin/Ry`; `elph.leer_lambda_out` reads `lambda.dat` (columns σ, λ, ∫α²F, ⟨log ω⟩, N(E_F)) or the text of `lambda.out`, takes $\mu^*$ from `lambda.in` (last numeric line) and fills the per-broadening $T_c$ column: from the final `lambda omega_log T_c` table of `lambda.out` if present and of matching size, and otherwise computes it row by row with `allen_dynes(λ_i, ω_log,i, μ*, correcciones=False)`; `ElPhRun.Tc_fuente` records which of the two was done and the report prints it. `elph.leer_a2F` reads `a2F.dos*` (or `A2F.dat`) and from it λ, $\omega_{\log}$ (if missing) and, always, $\bar\omega_2$ (`elph.omega_2`, computed in the CLI), which is what enables the $f_2$ factor in the summary $T_c$ values.
4. `elph.plato`, `elph.report` (table per broadening, regime, $T_c$ for three $\mu^*$, τ at 100/300/500/800 K), `elph.export` (`ELPH.dat`, `A2F.dat`, `ELPH.txt`), `elph.plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Broadenings, N(E_F) | `ph.out` (`Gaussian Broadening`, `DOS =`) | `elph.leer_elph_ph` |
| λ, ⟨log ω⟩ per broadening | `lambda.dat` / `lambda.out` from `lambda.x` | `elph.leer_lambda_out` |
| $\alpha^2F(\omega)$ | `a2F.dos*` (or `A2F.dat`) | `elph.leer_a2F`, column 1 THz, column 2 α²F |
| THz → K | `elph.THZ_K` | 47.9924 |
| $\mu^*$ | literal (0.10, 0.13, 0.16) | empirical |
| $T_{\mathrm{Debye}}$ | `--debye` | only to mark the validity regime of τ |

**Limits and pitfalls.** The λ values in `ph.out` are NOT read (only σ and N(E_F)); λ comes from `lambda.x` or from the $\alpha^2F$. The "Tc(K)" column of the per-broadening table is the one from `lambda.x` (Allen–Dynes WITHOUT $f_1 f_2$, with the $\mu^*$ of `lambda.in`) or the same thing recomputed by Olla-DFT; the corrected $T_c$ values for the three $\mu^*$ are those of the "Critical temperature" block, and without `a2F.dos*` there is no $\bar\omega_2$ and $f_2 = 1$. No plateau: "the k-mesh is insufficient … Any lambda reported from here is arbitrary". $\mu^*$ "is empirical (0.10-0.16) and is NOT computed here". τ "holds ABOVE the Debye temperature; below it overestimates the scattering". `lambda.x` with a coarse q-mesh leaves $\omega_{\log}$ as NaN: it is recomputed from the $\alpha^2F$ if available. τ is NOT injected automatically into `transport`: the module docstring says so explicitly and gives the sequence (`transport --collect` → `elph --collect` → $\sigma(T) = [\sigma/\tau](T)\cdot\tau(T)$ by hand on the columns of `TRANSPORTE.dat`).

**References.** P. B. Allen and R. C. Dynes, *Phys. Rev. B* 12, 905 (1975), DOI 10.1103/PhysRevB.12.905. W. L. McMillan, *Phys. Rev.* 167, 331 (1968), DOI 10.1103/PhysRev.167.331. P. B. Allen, *Phys. Rev. B* 3, 305 (1971) (high-T τ). G. Grimvall, *The Electron–Phonon Interaction in Metals*, North-Holland (1981).

---

### `olla-dft transport` — Electronic transport in CRTA: Seebeck, σ/τ, κ_e/τ, Lorenz and spin

**What it answers.** From the bands on a dense mesh: what are the Seebeck coefficient $S$, the conductivity per relaxation time $\sigma/\tau$, the electronic thermal conductivity $\kappa_e/\tau$, the power factor $S^2\sigma/\tau$ and the carrier concentration as functions of chemical potential and temperature? Does Wiedemann–Franz hold? How is transport shared between the two spin channels?

**Background for non-experts.** An electron in a band moves with velocity $v = (1/\hbar)\,dE/dk$. At a given temperature only the states within a few $k_BT$ of the chemical potential take part in transport (the $-\partial f/\partial E$ "window"). Summing velocity times velocity over that window gives the conductivity; weighting additionally by $(E-\mu)$ gives the Seebeck coefficient, which measures how much voltage appears per degree of temperature difference. The constant relaxation-time approximation (CRTA) assumes all electrons collide at the same rate $1/\tau$: then $\tau$ cancels in $S$ and in the Lorenz number (real predictions) but not in σ or κ_e, which are reported divided by τ.

**Formulas.** (`transport._fd_derivative`, `transport.compute`, `lorenz`, `cancelacion`, `TransporteEspin`) With $x = (E-\mu)/k_BT$, $-\partial f/\partial E = \mathrm{sech}^2(x/2)/(4k_BT)$, weights $w_k = 1/N_k$, $V$ the cell volume:

$$
\mathbf{v}_{n\mathbf{k}} = \frac{1}{\hbar}\nabla_{\mathbf{k}}E_{n\mathbf{k}}\ (\text{periodic finite differences, } \texttt{np.gradient}), \qquad
\mathbf{s}_m = \sum_{n\mathbf{k}} w_k\, \mathbf{v}\otimes\mathbf{v}\,(E-\mu)^m\left(-\frac{\partial f}{\partial E}\right)
$$

$$
\frac{\boldsymbol\sigma}{\tau} = \frac{e}{V}\mathbf{s}_0, \qquad
\mathbf{S} = -\frac{1}{T}\,\mathbf{s}_1\mathbf{s}_0^{-1}, \qquad
\frac{\boldsymbol\kappa_e}{\tau} = \frac{e}{VT}\mathbf{s}_2 - \mathbf{S}\mathbf{S}\,\frac{\boldsymbol\sigma}{\tau}\,T, \qquad
\mathrm{PF} = \bar S^2\,\bar\sigma/\tau
$$

$$
n = \frac{N_{\mathrm{elec}} - 2\sum_{n\mathbf{k}} w_k f(E_{n\mathbf{k}})}{V}, \qquad
L = \frac{\bar\kappa_e}{\bar\sigma T}, \qquad L_0 = 2.44\times10^{-8}\ \mathrm{W\,\Omega/K^2}, \qquad
c = \frac{|\bar\kappa_e|}{|\bar\kappa_e + \bar S^2\bar\sigma T|}
$$

$$
\sigma_{\mathrm{tot}} = \sigma_\uparrow + \sigma_\downarrow, \qquad
S_{\mathrm{tot}} = \frac{S_\uparrow\sigma_\uparrow + S_\downarrow\sigma_\downarrow}{\sigma_\uparrow+\sigma_\downarrow}, \qquad
P = \frac{\sigma_\uparrow-\sigma_\downarrow}{\sigma_\uparrow+\sigma_\downarrow}, \qquad S_{\mathrm{spin}} = S_\uparrow - S_\downarrow
$$

- $e$ = 1.602176634e-19 C; $\hbar$ = 6.582119569e-16 eV·s; $k_B$ = 8.617333262e-5 eV/K; σ/τ in S/(m·s); $S$ in V/K (µV/K in the report); κ_e/τ in W/(m·K·s); $n$ in cm⁻³ (positive = holes); bar = trace/3. The "cancellation" $c$ measures which fraction survives the subtraction $\kappa^0 - S^2\sigma T$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_transport` → `qekit/modules/transport.py: prepare`: standardised primitive cell; `scf.in` (mesh from `--kspacing` or the configured one) and `nscf.in` with `K_POINTS automatic` `--grid` (16x16x16), `nosym=.true.` (full mesh), `nbnd = 2 × estimated nbnd`; `--metal` turns off `occupations='fixed'`; `--nspin 2` and `--mag EL=value` (which implies `nspin=2`) write scf and nscf with spin polarisation, required for `--spin-resolved`.
2. `--run` executes scf and nscf; `--collect` → `transport.load` reads the first `out/*.xml`, rebuilds the grid from the fractional coordinates (rejects anything that is not a complete grid) and differentiates $E(\mathbf{k})$ with `np.gradient` on the periodically wrapped mesh, converting to Cartesian with $\mathbf{B}^{-T}$.
3. `transport.compute` over $T$ = `--temperatures` (300) and 201 values of µ in $E_F \pm$ `--mu-span` (1 eV).
4. `transport.report` (best p- and n-type $S$, maximum PF), `report_lorenz`, with `--spin-resolved` it loads `spin=1` and `report_espin`; `transport.export` (`TRANSPORTE.dat`), `transport.plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Eigenvalues, k, weights | nscf XML (`ks_energies`, Hartree → eV) | `qeout.read_xml`; `weights` replaced by $1/N_k$ |
| $E_F$, $N_{\mathrm{elec}}$, volume, cell | XML (`fermi_energy`, `nelec`, `cell`) | `qeout.read_xml` |
| Constants | `transport.E_CHARGE`, `HBAR_EVS`, `KB_EV`, `L0_SOMMERFELD` | CODATA; $L_0$ = 2.44e-8 |

**Limits and pitfalls.** CRTA only: "To give σ in S/m you need a τ that comes from a fit to a measurement or from an electron–phonon calculation — Olla-DFT does not invent it". It does NOT compute ZT (that would need κ_lattice and τ) nor couple the τ from `elph` automatically. It does not interpolate bands (unlike BoltzTraP): with a mesh < 24 per side or < 12000 points it warns "INSUFFICIENT … sigma comes out as isolated spikes". The Lorenz number inside the gap suffers catastrophic cancellation: "DO NOT TRUST THIS NUMBER … X % survives"; only points with $c > 0.10$ are summarised. `--spin-resolved` on an XML with `nspin = 1` is rejected with the instruction to re-prepare with `--nspin 2 --mag EL=0.7 --run`. Two-current model: "Valid as long as spin-flip scattering is slow … it stops holding near [the Curie temperature]".

**References.** G. K. H. Madsen and D. J. Singh, *Comput. Phys. Commun.* 175, 67 (2006), DOI 10.1016/j.cpc.2006.03.007 (BoltzTraP, same CRTA formulation). N. W. Ashcroft and N. D. Mermin, *Solid State Physics*, ch. 13 (Wiedemann–Franz). N. F. Mott, *Proc. R. Soc. A* 153, 699 (1936) (two-current model).

---

### `olla-dft ballistic` — Landauer conductance with `pwcond.x`

**What it answers.** How many conduction channels does an electrode have at each energy (complex bands) and what transmission probability $T(E)$ does a nanocontact or a molecule between electrodes have? What is its conductance in units of $G_0$?

**Background for non-experts.** In a macroscopic crystal the electron collides many times (diffusive transport, `transport`). In a few-atom contact it crosses in one go: there is no conductivity, there is conductance, given by the Landauer formula: $G = G_0 T(E_F)$, with $T$ the probability of passing summed over all open "lanes". Since $T$ cannot exceed the number of lanes, the conductance comes out quantised in steps of $G_0 = 2e^2/h$, and seeing those steps is the sign that the calculation is right.

**Formulas.** (`ballistic.G0`, `CondRun`)

$$
G = G_0\,T(E_F), \qquad G_0 = \frac{2e^2}{h} = 7.748091729\times10^{-5}\ \mathrm{S}, \qquad R = \frac{1}{G} = \frac{12.906\ \mathrm{k\Omega}}{T(E_F)}, \qquad T(E) \le N_{\mathrm{channels}}(E)
$$

- $T(E_F)$: transmission at the energy closest to $E - E_F = 0$ in the window. Region limits (`ballistic.longitud_z`): $\mathrm{bdl} = |\mathbf{a}_3|_{\mathrm{electrode}}/a$, $\mathrm{bds} = |\mathbf{a}_3|_{\mathrm{scatterer}}/a$ with $a = |\mathbf{a}_1|$ (alat units): the boundary of each region is the end of ITS cell, not the height of the last atom.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_ballistic` → `qekit/modules/ballistic.py: prepare`: `comprobar_geometria` requires $\mathbf{a}_3 \parallel z$, $\mathbf{a}_{1,2} \perp z$ and, with `--scatterer`, the same in-plane cell.
2. Writes `scf_electrodo.in` (and `scf_dispersor.in`) with `insulator=False` and prefixes `electr`/`disper`, and `cond.in` (`&inputcond`: `ikind` = `--ikind` (only 0 or 1) or 1 if there is a scatterer / 0 if not; `ikind=1` without `--scatterer` and `ikind=2` are rejected with an explanation; `energy0 = --emax` (3), `denergy = -(emax-emin)/(n-1)` (`--emin` −3, `--points` 61), `ewind = 1`, `epsproj = 1e-3`, `nz1 = --nz1` (3), `bdl` = `longitud_z(electrode)`, `bds` = `longitud_z(scatterer)`, one k-point (0, 0, 1)).
3. The user runs `pw.x` and `pwcond.x`; `--collect` → `ballistic.collect` reads `trans*.dat` (E, T) or the `T_tot` lines of `cond*.out`, and `Nchannels of the left tip` per energy (maximum over k); `ikind` from the `.out` or from `cond.in`.
4. `ballistic._avisar`, `report`, `export` (`BALISTICO.dat`, `.txt`), `plot` (T(E) and channel steps).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $T(E)$ | `trans.dat` from `pwcond.x` (or `T_tot` in `cond.out`) | `ballistic.collect` |
| Open channels | `cond.out` (`Nchannels of the left tip`) | maximum per energy |
| $G_0$ | `ballistic.G0` | 7.748091729e-5 S (CODATA) |
| bdl/bds limits | cell length along $z$ (`longitud_z`) | in alat units |

**Limits and pitfalls.** `ikind=0` "is NOT the conductance. It is the number of open channels, which bounds the conductance from above". If $T > N$: "That is impossible: T <= N by construction. Check that the bdl/bds limits…". Negative transmissions: "the calculation did not converge or … the geometry of the regions is badly cut". Different electrodes (`ikind=2` of `pwcond.x`) are not supported: `--ikind` only accepts 0 and 1, and asking for 2 gives "not implemented … the third scf and 'prefixr' and 'bdr' in cond.in must be written by hand". A single transverse k-point by default (0, 0, 1). There is no `--run`: it is always run by hand.

**References.** R. Landauer, *IBM J. Res. Dev.* 1, 223 (1957); M. Büttiker, *Phys. Rev. Lett.* 57, 1761 (1986), DOI 10.1103/PhysRevLett.57.1761. H. J. Choi and J. Ihm, *Phys. Rev. B* 59, 2267 (1999), DOI 10.1103/PhysRevB.59.2267; A. Smogunov, A. Dal Corso and E. Tosatti, *Phys. Rev. B* 70, 045417 (2004) (`pwcond.x`).

---

### `olla-dft cost` — Cost estimator calibrated with your history

**What it answers.** How long will this sweep take on THIS machine, and with what uncertainty? (`cost` shows the model; `--estimate` on any sweep applies it.)

**Background for non-experts.** The time of a plane-wave calculation scales in a known way with the number of k-points, plane waves, bands and iterations. What is not known in advance is the proportionality constant of each machine. Olla-DFT takes the shape from the physics and fits the scale with the calculations the user has already indexed in `olla-dft db` (with their wall times), and measures how wrong it is by leaving one system out and predicting it with the others.

**Formulas.** (`cost.n_ondas_planas`, `trabajo`, `iteraciones`, `_ajusta`, `estimar`)

$$
N_{\mathrm{PW}} = \frac{V\,E_{\mathrm{cut}}^{3/2}}{6\pi^2}\ (\text{bohr}^3,\ \mathrm{Ry}), \qquad
w_1 = n_k\, s\, N_{\mathrm{PW}}\, n_{\mathrm{b}}, \qquad
w_2 = n_k\, s\, N_{\mathrm{PW}}\, n_{\mathrm{b}}^2
$$

$$
t = t_0 + \left(C_1 w_1 + C_2 w_2\right)\, n_{\mathrm{scf}}\, n_{\mathrm{ion}}, \qquad
\text{NNLS fit with weights } t^{-1/2},\qquad
[t/\mathrm{disp},\ t\cdot\mathrm{disp}],\ \mathrm{disp} = e^{\sigma(\ln(\mathrm{pred}/\mathrm{real}))}
$$

- $V$: volume (Å³ → bohr³ with `A3_BOHR3`); $n_k$: irreducible k-points (spglib `get_ir_reciprocal_mesh`, or the real `number of k points` from a `pw.out`); $s$: 1, 2 or 4 (`nspin`, non-collinear); $n_{\mathrm{b}}$: `nbnd` from the input or $\max(4, 2N_{\mathrm{at}})$; $n_{\mathrm{scf}}$: median of the history (14 by default); $n_{\mathrm{ion}}$: median per type (`relax` 8, `vc-relax` 12 by default). The three-coefficient fit only if there are ≥ 8 calculations and $w_1^{\max}/w_1^{\min} \ge 5$; otherwise $C_1$ = geometric median of $t/w_1$.

**How Olla-DFT computes it.**
1. `qekit/cli.py: _cmd_cost` → `qekit/modules/cost.py: calibrar(--db olla-dft.db)`: `cost.historial` queries the `calculos` table (natoms, ecutwfc, kgrid, nspin, volume, n_scf, nk, nbnd, n_bfgs, wall_s, calculation).
2. `cost._prepara` builds $w_1 n_{it}$, $w_2 n_{it}$ and $t$; `cost._ajusta` (`scipy.optimize.nnls`, or `lstsq` clipped to ≥ 0).
3. Out-of-sample validation per system (`_clave_sistema`: natoms, ecutwfc, nk, calculation, nspin) if there are ≥ 4 systems and ≥ 8 remaining calculations: bias and dispersion of $\ln(\mathrm{pred}/\mathrm{real})$; otherwise the fit residual.
4. `cost.report_modelo` prints $t_0$, $C_1$, $C_2$, iterations, accuracy and warnings. In a sweep, `cli._run_or_explain` calls `cost.estimar_barrido` (reads each `pw.in` with `descriptores_de_input`; reuses the real $n_k$ of a point already run or from the history with the same formula and mesh) and `cost.report`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Timing history | `olla-dft.db` (SQLite, table `calculos`, `wall_s`) | from `olla-dft db` |
| $N_{\mathrm{PW}}$ | volume and `ecutwfc` from the `pw.in` | `cost.n_ondas_planas`; verified against QE in the tests |
| Irreducible $n_k$ | spglib or `pw.out` (`number of k points`) | `cost.k_irreducibles`, `nk_de_salida` |
| Fit weight | `cost.EXP_PESO` | 0.5 (chosen over 63 real calculations) |

**Limits and pitfalls.** "This tool tells ten minutes from six hours … It is not a stopwatch". Without history: "Nothing to calibrate with: the calculation database is empty or stores no times". Poorly varied history (`extrapola_bien` requires ≥ 8 calculations and a range ≥ 5): "predicting a system of another size can be off by a factor of two". spglib and `pw.x` do not always see the same symmetry: "that is where a factor of two or three goes". It does not model MPI parallelism (the time with `-j N` is simply total/N), nor the cost of `ph.x`, nor memory.

**References.** M. C. Payne et al., *Rev. Mod. Phys.* 64, 1045 (1992), DOI 10.1103/RevModPhys.64.1045 (scaling of plane-wave methods). C. L. Lawson and R. J. Hanson, *Solving Least Squares Problems*, SIAM (1995) (NNLS).

## Spectra, surfaces, chemistry and quality control

This part documents the physics behind the Olla-DFT commands that go beyond the total energy: optical and X-ray spectra (`optics`, `tddft`, `corehole`, `xanes`, `xps`), analysis of the electron density and of surfaces (`charges`, `charge`, `wf`, `esm`, `surface`, `adsorb`, `interface`), defect and reaction chemistry (`defect`, `eform`, `echem`, `neb`, `hull`), structure generation with machine-learned potentials (`amorphous`, `mlip`) and the quality-control tools that check that all of the above is comparable and credible (`audit`, `db`, `doctor`, `crosscheck`, `selftest`, `suggest`, `pseudos`). Each section states which question the command answers, which formulas the code actually implements (with the responsible function), which Quantum ESPRESSO file each datum comes from, and where the limits are. Whenever the code's internal documentation promises something the code does not do, it is said under "Limits and pitfalls".

---

### `olla-dft optics` — Dielectric function, absorption and Tauc gap

**What it answers.** How does the material respond to light? It gives $\varepsilon(\omega)$, the refractive index $n$, the extinction coefficient $k$, the absorption coefficient $\alpha$, the reflectivity $R$ and an optical gap extrapolated the way it is done with a UV-Vis spectrum.

**Background for non-experts.** When light crosses a solid, the electric field of the wave pushes the electrons. If the photon energy matches what an electron needs to jump from an occupied band to an empty one, the light is absorbed. The *dielectric function* $\varepsilon(\omega) = \varepsilon_1 + i\varepsilon_2$ summarises that response: the imaginary part $\varepsilon_2$ counts how many jumps exist at each energy (absorption) and the real part $\varepsilon_1$ how much the material polarises (refraction). The two are not independent: causality (the response cannot precede the cause) ties them through the Kramers–Kronig relations, so knowing one lets you rebuild the other.

`epsilon.x` computes $\varepsilon_2$ by summing all *vertical* transitions between Kohn–Sham bands, as if each electron jumped alone, without feeling the hole it leaves behind. This is the independent-particle approximation (RPA without local fields). The gap that comes out is the functional's, usually too small; the "scissor" is the simplest correction: all transitions are shifted rigidly by $\Delta$ and $\varepsilon_1$ is rebuilt by Kramers–Kronig so as not to break causality.

**Formulas.** All in `qekit/modules/optics.py`.

Isotropic average (`optics.collect`):
$$\varepsilon_{1,2}(\omega) = \tfrac{1}{3}\left[\varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz}\right]$$

Derived optical functions (`optics.derived`):
$$|\varepsilon| = \sqrt{\varepsilon_1^2 + \varepsilon_2^2},\qquad n = \sqrt{\frac{|\varepsilon| + \varepsilon_1}{2}},\qquad k = \sqrt{\frac{|\varepsilon| - \varepsilon_1}{2}}$$
$$\alpha(E) = \frac{2\,k\,E}{\hbar c},\qquad R = \frac{(n-1)^2 + k^2}{(n+1)^2 + k^2}$$

- $E = \hbar\omega$: photon energy (eV). $\hbar c$ = `HBAR_C_EV_CM` = $1.9732698\times10^{-5}$ eV·cm, so that $\alpha$ comes out in cm⁻¹. $n$, $k$, $R$ are dimensionless. Negative radicands are clipped to zero (`np.maximum`).

Kramers–Kronig (`optics.kramers_kronig`):
$$\varepsilon_1(\omega) = 1 + \frac{2}{\pi}\,\mathcal{P}\!\int_0^{\omega_{\max}} \frac{\omega'\,\varepsilon_2(\omega')}{\omega'^2 - \omega^2}\,d\omega'$$
- $\mathcal{P}$: principal value; implemented by removing the point $\omega'=\omega$ from the trapezoidal quadrature on the uniform `epsilon.x` grid. The integral is truncated at `wmax`.

Scissor (`optics.scissor`):
$$\varepsilon_2'(E) = \varepsilon_2(E-\Delta)\left(\frac{E-\Delta}{E}\right)^2,\qquad \varepsilon_1' = \mathrm{KK}[\varepsilon_2']$$
- $\Delta$: shift in eV (`--scissor`). The factor $((E-\Delta)/E)^2$ comes from $\varepsilon_2 \propto |p|^2/\omega^2$ with untouched matrix elements $|p|^2$. It is applied to each Cartesian component and then averaged.

Tauc plot (`optics.tauc_gap`):
$$y(E) = \left(\alpha E\right)^{1/r},\qquad r = \tfrac{1}{2}\ (\text{allowed direct}),\quad r = 2\ (\text{indirect})$$
$$E_g^{\mathrm{opt}} = -\frac{b}{m}\quad\text{with}\quad y \approx m E + b\ \text{fitted on the first absorption edge}$$

**How Olla-DFT computes it.**
1. `optics.prepare` resolves pseudopotentials and cutoffs (`sweep.prepare_common`, task `optics`) and **refuses** if any is not norm-conserving (`epsilon.x` has no matrix elements for USPP/PAW).
2. Writes `scf.in` (grid from the configured `kspacing`, 0.20 Å⁻¹ by default), `nscf.in` with a dense grid (`--kspacing`, default 0.12 Å⁻¹), `nosym=.true.` and `nbnd = 3 ×` the estimate of `inputgen._estimate_nbnd` (`nbnd_factor=3.0`), and `epsilon.in` (`calculation='eps'`, `smeartype='gauss'`, `intersmear=--smear` (0.10 eV), `wmin=0`, `wmax=--wmax` (20 eV), `nw=800`).
3. With `--run`: `pw.x` scf → `pw.x` nscf (`runner.run_all`) → `optics.run_epsilon` launches `epsilon.x` (looked up next to `pw.x`).
4. `optics.collect` reads `epsr_<prefix>.dat` and `epsi_<prefix>.dat` (columns: energy, xx, yy, zz) and averages.
5. If `--scissor Δ ≠ 0`: `optics.scissor` shifts $\varepsilon_2$ and rebuilds $\varepsilon_1$ with `optics.kramers_kronig`.
6. `optics.derived` yields $n, k, \alpha, R$; `optics.tauc_gap` fits the direct and indirect gaps; `optics.report` prints $\varepsilon_1(0)$ (value at $E \approx 0.05$ eV), $n(0)$, the maximum of $\varepsilon_2$ and the gaps.
7. `optics.export` writes `OPTICS.dat` with the columns of `optics.OPTICS_COLUMNS` (`E(eV)`, `eps1`, `eps2`, `n`, `k`, `alpha(1/cm)`, `R`), named in the last comment line so that `optics.read_optics_dat` can read them by name; `optics.plot` draws the three-panel figure.

Detail of the Tauc fit (`optics.tauc_gap`): the curve is smoothed with a ~0.05 eV moving average; the noise floor is the maximum of $y$ in the first 1 % of the spectrum; the edge starts where $y$ exceeds $\max(2\cdot\text{floor}, 10^{-3}\cdot\mathrm{median}(y>0))$ and $E > 0.1$ eV; it ends at the first local maximum that triples the onset value or at most `max_span` = 1.5 eV higher; a straight line is fitted in a window of `fit_window` = 0.6 eV centred on the steepest point of that stretch. Returns `None` if there is no absorption, if the slope is not positive or if the intercept falls outside the range.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $\varepsilon_1(\omega)$ per direction | `epsr_<prefix>.dat` from `epsilon.x` | `optics.collect`, columns 1–3 |
| $\varepsilon_2(\omega)$ per direction | `epsi_<prefix>.dat` from `epsilon.x` | `optics.collect` |
| $\hbar c$ | constant `optics.HBAR_C_EV_CM` | $1.9732698\times10^{-5}$ eV·cm |
| Type of each pseudo (NC/US/PAW) | UPF header (`pseudo_type`) | via `sweep.prepare_common` |
| $\Delta$ (scissor) | parameter `--scissor` | eV; recommended exp./GW gap − DFT gap |
| Broadening | parameter `--smear` | Gaussian `intersmear`, 0.10 eV |
| Window and points | `--wmax` (20 eV), `nw=800` | fixed in `optics.prepare` |
| nscf grid | `--kspacing` (0.12 Å⁻¹) | `sweep.default_grid` |

**Limits and pitfalls.**
- It is independent-particle RPA: no local fields, no excitons. The report reminds you: *"Recuerda: RPA de partícula independiente y gap del funcional…"*.
- Without NC pseudos the command aborts: *"epsilon.x solo funciona con pseudopotenciales de NORMA CONSERVADA…"*.
- `epsilon.x` does not include phonon-assisted transitions: in an indirect semiconductor $\varepsilon_2 = 0$ below the direct gap and the "indirect" fit does **not** give the true indirect gap (`tauc_gap` docstring).
- The Kramers–Kronig integral is truncated at `wmax`: $\varepsilon_1(0)$ inherits an error if there is strong absorption above 20 eV.
- The scissor only moves the gap; it neither corrects intensities nor adds excitons.
- If the Tauc fit fails, the report prints *"no se pudo ajustar"* instead of a number.

**References.**
- J. Tauc, R. Grigorovici, A. Vancu, *Phys. Status Solidi* 15, 627 (1966) — Tauc plot.
- `epsilon.x` manual (Quantum ESPRESSO, PP package): A. Benassi, *"epsilon.x: a post-processing tool for the calculation of the dielectric properties"*.
- M. Dressel, G. Grüner, *Electrodynamics of Solids* (Cambridge, 2002) — Kramers–Kronig and optical functions.

---

### `olla-dft tddft` — Optical absorption with TDDFPT (Lanczos/Davidson)

**What it answers.** Does the absorption spectrum change when the excited electron and the hole it leaves see each other? Where are the first excitations, which are bright and which dark, and is there a bound exciton below the gap?

**Background for non-experts.** `optics` sums one-electron transitions one at a time. In reality the excited electron (charge −) and the hole (charge +) attract each other; in molecules and wide-gap insulators that attraction lowers the energy of the pair and creates an absorption peak **inside** the gap: the exciton. Time-dependent density functional theory in linear response (TDDFPT) includes that interaction through the exchange-correlation kernel. Quantum ESPRESSO solves it in two ways: with the **Lanczos** algorithm (`turbo_lanczos.x` + `turbo_spectrum.x`), which gives the whole spectrum without computing empty states, or with **Davidson** (`turbo_davidson.x`), which gives the first N excitations one by one with their energy and oscillator strength $f$. An excitation with $f \approx 0$ exists but does not absorb light: it is "dark".

**Formulas.** In `qekit/modules/tddft.py`.

Unit conversion of the inputs (`tddft.build_lanczos_input`, `build_spectrum_input`, `build_davidson_input`):
$$E_{\mathrm{Ry}} = \frac{E_{\mathrm{eV}}}{\mathrm{RY\_EV}},\qquad \mathrm{RY\_EV} = 13.605693122994\ \mathrm{eV}$$

Wavelength (`tddft.report`):
$$\lambda\,(\mathrm{nm}) = \frac{1239.84}{E\,(\mathrm{eV})}$$

Absorption onset (`TddftRun.onset`): first local maximum of $dS/dE$ exceeding 20 % of the maximum of the derivative (inflection point of the first rise).

Exciton signature (`tddft._avisar`):
$$d = E_{\mathrm{onset}} - E_g^{\mathrm{IP}},\qquad d < -\max(0.10\ \mathrm{eV},\ 2\,\eta)\ \Rightarrow\ \text{bound exciton}$$
- $E_g^{\mathrm{IP}}$: independent-particle gap supplied by the user (`--gap`); $\eta$: broadening in eV, from `--broadening` at `--collect` or, if omitted, read from `spectrum.in` (`epsil`) or `davidson.in` (`broadening`) by `tddft._broadening_de_inputs` (Ry → eV). `UMBRAL_EXCITON` = 0.10 eV; `BROADENING_DEFAULT` = 0.05 eV.

Anisotropy (`tddft._anisotropia`): $\max_E[\max_i S_i(E) - \min_i S_i(E)] / \max_{i,E} S_i(E)$ over the $x,y,z$ components.

**How Olla-DFT computes it.**
1. `tddft.prepare`: if the minimum vacuum (`_vacio_minimo`) exceeds 5 Å or `--gamma` is given, it uses `K_POINTS gamma` (the only case TDDFPT implements); otherwise an automatic grid with a warning that `turbo_*.x` will stop. Writes `scf.in` with `nosym` and `noinv`.
2. Lanczos: `lanczos.in` (`itermax=--iter` 500, `ipol=--pol` 4 → `n_ipol=3`, `ltammd` with `--tamm-dancoff`, `lrpa` with `--rpa`, `scissor=--scissor/RY_EV` if a rigid shift of the empty bands is requested; `prepare` rejects a negative scissor or one combined with `--method davidson`) and `spectrum.in` (`itermax0=itermax`, `itermax=4×itermax` for the `--extrapolation` osc/constant/no, `epsil=--broadening/RY_EV`, `units=1` (eV), `start/end/increment`).
3. Davidson: `davidson.in` with `num_eign=--states` (10), `num_init=2N`, `num_basis_max=max(80, 8N)`, `residue_conv_thr=1e-4`, `p_nbnd_virt=15`, window and `broadening` in Ry, `reference` at the window centre.
4. The user runs `pw.x` → `turbo_lanczos.x` → `turbo_spectrum.x` (or `turbo_davidson.x`) by hand.
5. `tddft.collect --collect` (with `broadening` from `--broadening` or from the inputs): Lanczos reads the first `*plot*.dat` (columns: E in eV, total S, S_x, S_y, S_z) and, from `lanczos.out`, the `itermax` and the functional. Davidson (`_collect_davidson`) reads `<prefix>.eigen` (energy in Ry → eV, total strength, strengths per direction) and the `*plot*.dat` if present.
6. `_picos` lists local maxima above 5 % of the maximum; `_avisar` compares the onset with `--gap`; `report` flags as bright the excitations with $f > 0.01$ and counts the dark ones.
7. `export` writes `TDDFT.dat`, `TDDFT_EXCITACIONES.dat`, `TDDFT.txt`; `plot` optionally overlays the `optics` spectrum (`--compare OPTICS.dat`: the CLI reads the `alpha(1/cm)` column by name with `optics.read_optics_dat` and normalises it to the TDDFPT maximum).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $S(E)$ and components | `<prefix>.plot_S.dat` (or `*plot*.dat`) from `turbo_spectrum.x` | `tddft.collect`; energy in eV because of `units=1` |
| Excitations $(E, f, f_x, f_y, f_z)$ | `<prefix>.eigen` from `turbo_davidson.x` | `_collect_davidson`; E in Ry × `RY_EV` |
| `itermax`, functional | `lanczos.out` | regex in `tddft.collect` |
| IP gap | parameter `--gap` | eV |
| Broadening $\eta$ | `--broadening` or `spectrum.in`/`davidson.in` | `tddft._broadening_de_inputs`, Ry × `RY_EV` |
| Scissor | `--scissor` (Lanczos only) | eV → Ry in `lanczos.in` |
| $\alpha(E)$ from `optics` | `OPTICS.dat`, column `alpha(1/cm)` | `optics.read_optics_dat` |
| Ry → eV | `tddft.RY_EV` | 13.605693122994 |
| Minimum vacuum | cell geometry | `tddft._vacio_minimo` |

**Limits and pitfalls.**
- With LDA/GGA the adiabatic kernel does **not** bind excitons in a solid; the report says so: *"con LDA o GGA el kernel adiabático NO liga excitones en un SÓLIDO… En MOLÉCULAS sí mejora."*
- Γ point only: with a k grid the report warns *"OJO: TDDFPT solo tiene implementado el caso gamma y se plantará al leer el input"*.
- Molecule with < 6 Å of vacuum and < 30 atoms: *"AVISO: solo hay X Å de vacío…"*.
- `--scissor` only exists in `turbo_lanczos.x`: with `--method davidson` or with a negative value the command aborts (*"--scissor solo existe en turbo_lanczos.x…"*). The TDDFPT scissor rebuilds nothing by Kramers–Kronig: QE's own code applies it to the empty bands.
- `--compare` requires an `OPTICS.dat` with the `alpha(1/cm)` column; if missing: *"'…' no tiene la columna 'alpha(1/cm)'; --compare espera el OPTICS.dat de 'olla-dft optics'."*
- If neither `--broadening` nor the inputs give the broadening, the exciton threshold stays at `UMBRAL_EXCITON` = 0.10 eV.
- The command does not launch the `turbo_*.x` executables: it only writes inputs and reads outputs.

**References.**
- D. Rocca, R. Gebauer, Y. Saad, S. Baroni, *J. Chem. Phys.* 128, 154105 (2008) — TDDFPT Lanczos.
- O. B. Malcıoğlu, R. Gebauer, D. Rocca, S. Baroni, *Comput. Phys. Commun.* 182, 1744 (2011) — turboTDDFT.
- X. Ge, S. J. Binnie, D. Rocca, R. Gebauer, S. Baroni, *Comput. Phys. Commun.* 185, 2080 (2014) — turboTDDFT 2.0 (Davidson).

---

### `olla-dft corehole` — Core-hole pseudopotentials (ld1.x)

**What it answers.** How to describe an atom from which an electron of an inner shell has been removed? It generates the pair of pseudopotentials (normal + core-hole) that `xps` and `xanes` need, with the same configuration and the same radii, and extracts the core wavefunction that `xspectra.x` reads.

**Background for non-experts.** A pseudopotential replaces the nucleus and the inner ("core") electrons by an effective potential, so that the calculation only treats the valence electrons. To simulate an X-ray spectroscopy one must remove an electron from that frozen core: this requires a different pseudopotential, generated on purpose with the atomic program `ld1.x`, in which the occupation of the core level (1s for the K edge, 2p for L₂,₃, etc.) is one less. Since the core has one electron less, the declared valence charge `z_valence` rises by exactly 1: that unit **is** the hole. The two pseudos must be generated together with the same parameters, because comparing energies made with pseudos from different families means nothing.

**Formulas.** This module evaluates no physical formulas; it builds `ld1.x` inputs from explicit rules in `qekit/core/atomconf.py` and `qekit/modules/corehole.py`:

- Electronic configuration by Aufbau filling (`atomconf.aufbau`, Madelung order `ORDEN`) with the exceptions in `atomconf.EXCEPCIONES` (Cr, Cu, Nb, Mo, Ru, Rh, Pd, Ag, La, Ce, Gd, Pt, Au).
- Core/valence partition (`atomconf.particion`): valence = $s,p$ shell of $n_{\max}$ + any partially filled $d$/$f$ + filled $d$ of the previous row; with `--semicore`, also $(n-1)s,(n-1)p$.
- Hole (`atomconf.config_hueco`): occupation of the `BORDES[edge]` level reduced by 1.0; rejected if the level is not in the core.
- Pseudisation channels (`atomconf.canales_pseudo`): the valence plus one **unoccupied** channel (occupation −2) for every missing $l \le 2$, with $n = \max(n_{\max}, l+1)$; with `--projectors 2` a second projector per channel labelled $n+1$ with occupation −1.
- Cutoff radius per row (`corehole.RCUT_FILA`): {1: 1.0, 2: 1.3, 3: 1.7, 4: 2.0, 5: 2.2, 6: 2.4} bohr; `rcutus = 1.25 · rcut` only if `pseudotype=3`.
- Reference energies of unbound channels: `E_CANAL_VACIO` = 0.15 Ry; second projector `E_SEGUNDO_PROYECTOR` = 0.05 Ry.

**How Olla-DFT computes it.**
1. `corehole.generar`: validates the element (H..Rn), forces `pseudotype=3` if 2 projectors are requested, obtains partition, channels and `rcut` (or `--rcut`).
2. `corehole.input_ld1` writes `ld1_base.in` and `ld1_hueco.in` (`iswitch=3`, `rel=--rel`, `beta=0.3`, `dft=--functional` (PBE), `tm=.true.`, `lloc` = highest $l$ among the channels, `lgipaw_reconstruction=.true.`, `author='Olla-DFT'`). The empty channels are also added to the all-electron configuration (`_con_canales_vacios`) because `ld1.x` requires them to exist.
3. `corehole._correr_ld1` runs `ld1.x < ld1_X.in > ld1_X.out` (unless `--only-inputs`) and fails if `Error in routine` appears.
4. `corehole.leer_upf` reads from each UPF `element`, `z_valence`, `mesh_size`, `pseudo_type`, `functional`, `wfc_cutoff`, `rho_cutoff` and the `PP_GIPAW_CORE_ORBITAL` labels.
5. `corehole.verificar` applies the checks in the following table; `report` and `export` (`PSEUDOS_HUECO.txt`) list them. The exit code is 1 if there is any `FALLA`.
6. With `--core-wfc UPF`: `corehole.core_wfc` extracts the core wavefunctions in the `filecore` format of `xspectra.x` (one block per orbital, separated by a blank line, in UPF order) and verifies that the number of points matches `mesh_size`.

| Check (`corehole.verificar`) | Criterion | Flag |
|---|---|---|
| Difference in `z_valence` | exactly +1 (tolerance 1e-6) | FALLA otherwise |
| Radial meshes | `mesh_size` equal in both UPFs | FALLA otherwise |
| Hole orbital | present among the `PP_GIPAW_CORE_ORBITAL` of the core-hole UPF | FALLA otherwise |
| Functional | same in both UPFs | FALLA otherwise |
| Projectors | warning if only one per channel (XSpectra recommends two) | warning |
| Ghost states, logarithmic derivatives, transferability | **not checked** | explicit warning |

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Electronic configuration | `atomconf.aufbau` + `EXCEPCIONES` | rule, not experimental data |
| Edge level | `atomconf.BORDES` | K=1s, L1=2s, L23=2p, M1=3s, M23=3p, M45=3d |
| `z_valence`, `mesh_size`, type, functional | `PP_HEADER` of the generated UPF | `corehole.leer_upf` |
| Core orbitals | `PP_GIPAW_CORE_ORBITAL.n` sections of the UPF | `leer_upf`, `core_wfc` |
| Radial mesh | `PP_R` of the UPF | `core_wfc` |
| Cutoff radius | `RCUT_FILA` or `--rcut` | bohr |

**Limits and pitfalls.**
- The M edges (`M1`, `M23`, `M45`) exist in `atomconf.BORDES` and serve to generate the core-hole pseudo (XPS), but `xspectra.x` only implements K, L1, L2, L3 and L23: `olla-dft xanes` rejects them (`xanes.validar_borde`).
- The report warns: *"NO verificado automáticamente: estados fantasma, derivadas logarítmicas y transferibilidad… el cutoff del pseudo anterior NO sirve para este."* The cutoff must be reconverged with `olla-dft converge`.
- With `--projectors 2` the pseudo comes out ultrasoft and *"casi siempre hay que ajustar --rcut a mano hasta que ld1.x converja"*.
- `ld1.x` is not built by default in QE (`make ld1`).

**References.**
- A. Dal Corso, *Comput. Mater. Sci.* 95, 337 (2014) — pslibrary and `ld1.x`.
- N. Troullier, J. L. Martins, *Phys. Rev. B* 43, 1993 (1991) — TM pseudisation (`tm=.true.`).
- C. J. Pickard, F. Mauri, *Phys. Rev. B* 63, 245101 (2001) — GIPAW reconstruction.

---

### `olla-dft xanes` — X-ray absorption near the edge (xspectra.x)

**What it answers.** What is the shape of the XANES/NEXAFS spectrum of a given atom at a given edge, with a core hole and a given polarisation, and how much does it depend on the field direction?

**Background for non-experts.** An X-ray photon knocks an electron out of a deep level (1s at the K edge) and sends it to the empty states. The dipole selection rule only allows final states with angular momentum $l \pm 1$: from 1s one sees the empty $p$ states **of that atom**. The spectrum is, in essence, the density of empty states projected on the absorber, which is why it is local, element-selective and sensitive to oxidation state and coordination. The hole left by the electron attracts the empty states and shifts the edge, so the absorbing atom is described with the core-hole pseudopotential from `corehole`, and since the ejected electron is assumed to leave the system, the cell carries total charge +1 (full core hole approximation, FCH). `xspectra.x` computes the cross-section with the Lanczos method and continued fractions without building the empty states.

**Formulas.** In `qekit/modules/xanes.py`.

Powder average (`xanes.collect`):
$$\sigma(E) = \tfrac{1}{3}\left[\sigma_x(E) + \sigma_y(E) + \sigma_z(E)\right]$$

Minimum distance between absorber images (`xanes.distancia_imagen_minima`):
$$d_{\min} = \min_{(i,j,k)\neq 0,\ |i|,|j|,|k|\le 1}\left|i\,\mathbf{a} + j\,\mathbf{b} + k\,\mathbf{c}\right|$$

Operational onset (`xanes.onset`): first energy at which $\sigma \ge 0.5\,\sigma_{\max}$. Anisotropy (`_anisotropia`): $\max_E[\mathrm{ptp}_i\,\sigma_i(E)]/\max\sigma$; highlighted if $> 0.1$.

**How Olla-DFT computes it.**
1. `xanes.validar_borde` (also in `_cmd_xanes`) normalises `--edge` and only accepts `BORDES_XSPECTRA` = K, L1, L2, L3, L23; M edges are rejected with an explicit message. `BORDE_COREHOLE` says which `--edge` of `corehole` generates the hole for each edge (L2 and L3 share the 2p hole = `L23`). `xanes.prepare` locates the `--element`/`--site` atom, moves it to the **first** position of the list and declares it as a separate species with a three-letter label (`etiqueta_excitada`, e.g. `Sih`; the QE limit is `CHARACTER(LEN=3)`).
2. `sweep.prepare_common` (task `xanes`, excluding the core-hole UPF) and `inputgen.build_pw_input` write `scf.in` with `tot_charge = 1.0`; `_marcar_absorbedor` adds the excited species to `ATOMIC_SPECIES`, changes the label of the first atom and increments `ntyp` by 1.
3. `corehole.core_wfc` extracts `<El>.wfc` from the core-hole UPF (`PP_GIPAW_CORE_ORBITAL` sections).
4. `xanes.build_xspectra_input` writes `xspectra_pol.in` (or `xspectra_x/y/z.in` with `--average`): `calculation='xanes_dipole'`, `edge=--edge`, `xiabs=1`, `xepsilon=--polarization`, `xniter=2000`, `xcheck_conv=10`, `xerror=0.001`; `&plot`: `xnepoint=1000`, `xgamma=--broadening` (0.8 eV), `xemin=-10`, `xemax=30`, `terminator`, `cut_occ_states=.true.`; `&pseudos`: `filecore`, `r_paw(1)=--r-paw` (3.0); `&cut_occ`: `cut_desmooth=0.1`, `cut_stepl=0.01`; k grid at the end.
5. The report measures $d_{\min}$ and warns if it is below `DIST_MINIMA` = 8 Å.
6. The user runs `pw.x -in scf.in` and `xspectra.x -in xspectra_*.in`.
7. `xanes.collect --collect` reads all `xanes_*.dat` (columns E − E_F, σ), averages if there are several, and reads `xgamma` from the *"Broadening parameter (in eV)"* comment in the first file.
8. `report` gives the 50 % onset, main maximum, peaks (> 5 % of the maximum), anisotropy; `export` writes `XANES.dat` and `XANES.txt`; `plot` the figure.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $\sigma(E)$ per polarisation | `xanes_<dir>.dat` from `xspectra.x` | `xanes._leer_dat` |
| Broadening `xgamma` | header of `xanes_*.dat` | regex *"Broadening parameter"* |
| Core wavefunction | core-hole UPF (`PP_GIPAW_CORE_ORBITAL`) | `corehole.core_wfc` |
| Total charge +1 | fixed in `xanes.prepare` | `tot_charge=1.0` |
| Polarisation | `--polarization` (1 0 0) or the `EJES` axes with `--average` | Cartesian vector (`xcoordcrys=.false.`) |
| k grid | `--kspacing` → `sweep.default_grid` | also in `xspectra.in` |
| $d_{\min}$ | cell vectors | `distancia_imagen_minima` |

**Limits and pitfalls.**
- The energy axis is relative to the Fermi level, not photon energy: *"Para comparar con un experimento se alinea el borde y se compara la FORMA."*
- Supercell warning: *"AVISO: X Å es poco. Con condiciones periódicas el hueco de core ve sus propias imágenes…"* (threshold 8 Å).
- Single polarisation: *"UNA sola polarización. En un cristal anisótropo el espectro depende de la dirección…"*.
- The onset (`xanes.onset`) is the first point where σ reaches 50 % of the **global maximum**: a weak pre-edge before the white line does not count as the onset (the docstring now states this).
- M edges: *"xspectra.x solo calcula bordes K y L (K, L1, L2, L3, L23); los bordes M no están implementados en QE, aunque 'olla-dft corehole' pueda generar el pseudo con ese hueco."*
- `distancia_imagen_minima` only looks at the 26 neighbouring cells: for very oblique cells it can overestimate $d_{\min}$.
- Without `--core-hole` the command aborts: *"falta --core-hole con el UPF de hueco de core. Sin él se calcularía el espectro del estado fundamental…"* and suggests `olla-dft corehole <El> --edge <BORDE_COREHOLE[edge]>`.

**References.**
- M. Taillefumier, D. Cabaret, A.-M. Flank, F. Mauri, *Phys. Rev. B* 66, 195107 (2002) — XSpectra, Lanczos with continued fractions.
- C. Gougoussis, M. Calandra, A. P. Seitsonen, F. Mauri, *Phys. Rev. B* 80, 075102 (2009) — XSpectra with PAW/GIPAW.
- O. Bunău, M. Calandra, *Phys. Rev. B* 87, 205105 (2013) — L₂,₃ edges.

---

### `olla-dft xps` — Initial-state core-level shifts (initial_state.x)

**What it answers.** By how much does the core-level energy of each atom shift relative to the others of its species? It is the theoretical counterpart of the chemical shift in an XPS spectrum.

**Background for non-experts.** XPS measures the energy needed to eject a core electron. An atom surrounded by electronegative neighbours has its core more tightly bound (positive shift) than one in a metallic environment. The **initial-state** approximation computes only how the potential felt by the core electron changes *before* it is removed; it ignores the relaxation of the other electrons around the hole (the *final state*), which can amount to several tenths of an eV. That is why what comes out are **relative** shifts between sites, not absolute binding energies. `initial_state.x` needs two species of the same element in the input — the normal one and one with a core hole — because it defines the shift from `delta_zv = zv(excited) − zv(normal)`; if both are the same it returns zeros without warning.

**Formulas.** In `qekit/modules/xps.py`. The shift is computed by `initial_state.x`; Olla-DFT only reads and rearranges it:

$$\Delta_i = \text{shift}_i^{\mathrm{TOTAL}},\qquad \Delta_i^{\mathrm{rel}} = \Delta_i - \min_j \Delta_j,\qquad \text{spread} = \max_i\Delta_i - \min_i\Delta_i$$

Cancellation indicator (`xps.report`):
$$\frac{\max_{c}\,\mathrm{ptp}(\text{contribution}_c)}{\text{spread}} > 20 \Rightarrow \text{numerical-cancellation warning}$$

- $\Delta_i$: shift of atom $i$ in eV, read from the line `atom i type t shift = … Ry, = … eV` of the *TOTAL* section. The contributions $c$ (Fermi, local, non-local, ionic, core-correction, Hubbard…) are read from the *"The X contribution to shift"* sections.

**How Olla-DFT computes it.**
1. `xps.prepare` reads `--core-hole EL=file.UPF` (repeatable). For each element: `_verificar_par` requires the normal and the core-hole UPF to be different files and `z_valence` to differ by exactly +1 (`qekit.core.pseudo.z_valence`).
2. `inputgen.build_pw_input` writes `scf.in` with the extra species (`extra_species`) declared in `ATOMIC_SPECIES` **without** any atom using them; `_copiar_pseudos` copies the core-hole UPF into `pseudo_dir`.
3. `xps.build_input` writes `initial_state.in` with `excite(t_normal) = t_hole` (1-based indices in `ATOMIC_SPECIES` order); `excite(t)=t` is rejected.
4. `structure.symmetry_dataset` counts orbits of equivalent atoms; if there is only one it warns that everything will come out zero.
5. The user runs `pw.x -in scf.in` and `initial_state.x -in initial_state.in > initial_state.out`.
6. `xps.collect --collect` parses `initial_state.out` with `_RE_SECCION` and `_RE_ATOMO`, takes the eV column, and sets `equivalentes=True` if all $|\Delta_i| < 10^{-6}$ eV.
7. `report` tabulates shifts, shift relative to the minimum, spread, decomposition per contribution and the cancellation warning; `export` writes `XPS_CORE.dat`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Shift per atom and contribution | `initial_state.out` | `xps.collect`, regex `atom N type T shift = X Ry, = Y eV` |
| `z_valence` normal and core-hole | UPF headers | `pseudo.z_valence` in `_verificar_par` |
| Inequivalent sites | spglib via `structure.symmetry_dataset` | `equivalent_atoms` |
| Excited species label | `xanes.etiqueta_excitada` | 3 characters |
| Symbols per atom | input structure | `atoms.get_chemical_symbols()` |

**Limits and pitfalls.**
- **Initial state only.** There is no ΔSCF and no final state; the module docstring now states it explicitly: the core-hole UPF is used *only as the "excited species" that initial_state.x needs to define the shift, not to relax the system around the hole*. The report refers you on: *"las energías de enlace absolutas necesitan un ΔSCF con hueco de core."*
- Spread < 0.1 eV: *"Por debajo de ~0.1 eV el corrimiento no es concluyente: la relajación de estado final… es del mismo orden."*
- Without `--core-hole` only `scf.in` is written and the report explains that `initial_state.x` would return zeros.
- All atoms equivalent: *"AVISO: todos los átomos son equivalentes por simetría, así que todos los corrimientos van a salir exactamente cero."*
- Large cancellation: *"CUIDADO con la cancelacion… baja conv_thr (1e-10 o menos) y sube la malla k antes de creerte la tercera cifra."*
- The error messages point to `olla-dft corehole <El> --edge K` to generate the consistent pair.

**References.**
- E. Pehlke, M. Scheffler, *Phys. Rev. Lett.* 71, 2338 (1993) — initial vs final state in core-level shifts.
- L. Köhler, G. Kresse, *Phys. Rev. B* 70, 165405 (2004) — core-level binding energies with a core hole.
- `initial_state.x` documentation (Quantum ESPRESSO, PP package).

---

### `olla-dft charges` — Löwdin charges, on-grid Bader and density difference

**What it answers.** How much electronic charge "belongs" to each atom, and where does the density accumulate or deplete when a bond or an adsorption forms?

**Background for non-experts.** The electron density is continuous; sharing it among atoms requires a rule. **Löwdin** projects the states onto orthogonalised atomic orbitals (done by `projwfc.x`); it is cheap and depends on the orbital basis of the pseudopotential. **Bader** uses no orbitals: it divides space into "basins" by following the steepest ascent of the density from each point to a maximum, like rainwater running down slopes into each valley, but reversed. The **density difference** $\rho_{AB} - \rho_A - \rho_B$ shows, point by point, what changed when the two parts were joined.

**Formulas.** In `qekit/modules/charges.py`.

Löwdin (`charges.read_lowdin`, `report_lowdin`):
$$q_i^{\mathrm{net}} = Z_i^{\mathrm{val}} - Q_i^{\mathrm{Löwdin}}$$

On-grid Bader (`charges.bader`): for each grid point the neighbour $\nu$ (out of 26) that maximises the slope
$$s_\nu = \frac{\rho(\mathbf{r}+\mathbf{d}_\nu) - \rho(\mathbf{r})}{|\mathbf{d}_\nu|}$$
is chosen and the chain is followed up to a local maximum (path compression, max. 64 iterations). Each maximum is assigned to the nearest atom with periodic images. Then
$$Q_i = \sum_{\mathbf{r}\in\Omega_i}\rho(\mathbf{r})\,\Delta V,\qquad V_i = N_i\,\Delta V_{\mathrm{Å}^3},\qquad \Delta V_{\mathrm{Å}^3} = \frac{V_{\mathrm{cell}}}{n_1 n_2 n_3},\quad \Delta V = \frac{\Delta V_{\mathrm{Å}^3}}{a_0^3}$$
- $\rho$: density from the `.cube` in e/bohr³ (what `pp.x` writes, `density_units="e/bohr3"` by default; `"e/A3"` is also accepted); `charges._voxel_volume` returns the voxel volume in the units of the density (bohr³, with $a_0$ = `fields.BOHR` = 0.529177210903 Å) so that $\rho\,\Delta V$ is a number of electrons, and in Å³ to report the basin volumes.

Density difference (`charges.difference`, `report_difference`), with the same $\Delta V$ in bohr³:
$$\Delta\rho = \rho_{\mathrm{total}} - \sum_p \rho_p,\qquad Q_{\mathrm{net}} = \sum \Delta\rho\,\Delta V,\qquad Q_{\mathrm{acc}} = \sum_{\Delta\rho>0}\Delta\rho\,\Delta V$$

**How Olla-DFT computes it.**
1. If the structure is given, `charges.valence_from_pseudos` reads `z_valence` from the UPFs in `--pseudo-dir` (or the configured `pseudo_dir`) via `pseudo.resolve`; if any UPF cannot be read it returns `None`, the CLI warns (*"no pude leer z_valence de los UPF…"*) and the "neta" column stays `n/d`.
2. `--lowdin projwfc.out`: `charges.read_lowdin` looks for `Atom #  i: total charge = q` and `Spilling Parameter:`; with the structure it adds symbols and the net charge $Z^{\mathrm{val}} - Q$.
3. `--bader density.cube` (needs the structure): `fields.read_cube` reads the cube (`plot_num=0` from `pp.x`), `charges.bader` partitions and compares the sum of basins with the total integral; `report_bader` also compares the integral with $\sum_i Z_i^{\mathrm{val}}$ and warns if they differ by more than 5 %.
4. `--difference total.cube part1.cube …`: `charges.difference` requires identical grids and subtracts; `report_difference` gives net charge, accumulated charge and the extrema of the planar profile (`fields.planar_average`, axis `--axis`); `plot_difference` draws the profile.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Löwdin charges, spilling | `projwfc.x` output | regex `_RE_LOWDIN`, `_RE_SPILL` |
| $\rho(\mathbf{r})$ | `.cube` from `pp.x` (`plot_num=0`, `output_format=6`) | `fields.read_cube` |
| Atomic positions (Bader) | structure `file` | `atoms.positions` |
| $Z^{\mathrm{val}}$ per atom | `z_valence` from the UPFs (`--pseudo-dir` or configuration) | `charges.valence_from_pseudos` → `pseudo.resolve` |
| $a_0$ (bohr → Å) | `fields.BOHR` | 0.529177210903 |
| Profile axis | `--axis` (0/1/2) | `fields.planar_average` |

**Limits and pitfalls.**
- On-grid Bader: *"Hereda el sesgo de malla del método (centésimas de electrón); para números finos usa la variante near-grid del código `bader` de Henkelman."* Warning if the sum of basins differs from the integral by more than 1e-3 e: *"la malla del cube es demasiado gruesa."*
- If the grid integral does not match $\sum Z^{\mathrm{val}}$ (by more than 5 %): *"Revisa que el cube sea la densidad de valencia completa (plot_num=0) y que los UPF de --pseudo-dir sean los del cálculo."* A cube already in e/Å³ must be declared with `density_units="e/A3"` (Python only; the CLI assumes e/bohr³).
- Löwdin: spilling > 0.05 → *"AVISO: por encima de ~0.05 la base atómica no describe bien los estados"*. Useful to compare atoms, not as an absolute charge.
- Without a readable `--pseudo-dir` the "neta" column is `n/d` with a warning; the UPFs must be those of the calculation, because $Z^{\mathrm{val}}$ depends on the pseudo (semicore or not).
- The $\Delta\rho$ profile is plotted in e/bohr³ (cube units), not e/Å³.
- `--difference` requires the same cell, FFT grid and cutoffs: *"las rejillas no coinciden… la resta no significa nada."*

**References.**
- R. F. W. Bader, *Atoms in Molecules: A Quantum Theory* (Oxford, 1990).
- G. Henkelman, A. Arnaldsson, H. Jónsson, *Comput. Mater. Sci.* 36, 354 (2006) — on-grid Bader.
- P.-O. Löwdin, *J. Chem. Phys.* 18, 365 (1950).

---

### `olla-dft charge` — pp.x scalar fields and planar profile

**What it answers.** How are the charge density, the spin density, the ELF or the electrostatic potential of a finished calculation distributed along an axis?

**Background for non-experts.** `pp.x` extracts from the already computed wavefunctions and density a scalar field on the 3D grid. Averaging it over the planes perpendicular to an axis gives a 1D "profile" that is easy to read: where the layers of a slab are, where spin accumulates, where the vacuum is.

**Formulas.** `fields.planar_average`:
$$\bar f(z_k) = \frac{1}{n_1 n_2}\sum_{i,j} f(i,j,k),\qquad z_k = k\,|\mathbf{h}_3|$$
- $\mathbf{h}_3$: grid step along the chosen axis (Å). The other axes are obtained by permutation.

**How Olla-DFT computes it.**
1. `_cmd_charge`: if `<name>.cube` does not exist (or with `--rerun`), `fields.run_pp` writes `pp_<field>.in` with the `plot_num` from `fields.PLOTS` (density 0, vtotal 1, spin 6, elf 8, potential 11), `iflag=3`, `output_format=6`, and runs `pp.x` (looked up next to `pw.x`); requires `JOB DONE`.
2. `fields.read_cube` reads origin, axes (bohr → Å if $n>0$) and values.
3. `fields.planar_average` along `--axis` (a/b/c); `PERFIL_PLANAR.dat` and the figure `perfil_<name>` are written.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| 3D field | `<name>.cube` from `pp.x` | `pp.x` units: e/bohr³ (density), Ry (potentials) |
| `prefix` | XML of the calculation | `qeout.read_xml(...).prefix` |
| `plot_num` | table `fields.PLOTS` | 0, 1, 6, 8, 11 |
| Bohr → Å | `qeout.BOHR_ANG` | 0.529177210903 |

**Limits and pitfalls.**
- The profile is exported in the raw cube units (no Ry → eV conversion here; `wf` does it).
- Needs `pp.x` compiled (`make pp`); if missing: *"no se encontró pp.x junto a pw.x…"*.
- The command does not interpret the field: it only averages and plots it. The `.cube` opens in VESTA for isosurfaces.

**References.** `pp.x` documentation (INPUT_PP, Quantum ESPRESSO).

---

### `olla-dft wf` — Work function from the vacuum level

**What it answers.** How much energy does it take to remove an electron from a surface into vacuum? $\Phi = V_{\mathrm{vac}} - E_F$.

**Background for non-experts.** In a slab with vacuum, the electrostatic potential flattens far from the material: that plateau is the "vacuum level", the energy of an electron at rest outside the solid. The work function is the distance from the Fermi level (the last occupied level) to that plateau. If the plateau is not flat, either the vacuum is short or the slab has a net dipole that tilts the potential.

**Formulas.** `fields.work_function`:
$$\bar V(z) = \mathrm{RY\_EV}\cdot\overline{V_{\mathrm{pp}}}(z),\qquad V_{\mathrm{vac}} = \frac{1}{2h+1}\sum_{k=-h}^{h}\bar V\big(z_{i^\ast + k}\big),\qquad \Phi = V_{\mathrm{vac}} - E_F$$
$$\text{flatness} = \max_{k}\bar V - \min_{k}\bar V\ \text{in the same window}$$
- The index window $\{i^\ast + k\}$ is given by `fields.vacuum_window` when the atomic positions are known (the CLI passes them from the XML): it is the central 20 % of the widest **atom-free** gap along the axis (measured periodically in fractional coordinates), with $h = \max(2, 0.1\,f_{\mathrm{gap}} N_z)$. Without positions it falls back to the blind criterion: $i^\ast = \arg\max_z \bar V$ and $h = \max(2, N_z/10)$ (±10 % of the cell around the maximum). $E_F$ in eV from the XML; `RY_EV` = 13.605693122994.

**How Olla-DFT computes it.**
1. `_cmd_wf`: if `potencial.cube` does not exist, `fields.run_pp(path, "potential", ...)` runs `pp.x` with `plot_num=11` ($V_{\mathrm{bare}} + V_H$).
2. `fields.read_cube` and `qeout.read_xml` (for `fermi`, from the `fermi_energy` tag in Ha → eV).
3. `fields.work_function(cube, E_F, axis, positions=qe.positions)` averages in the plane, locates the vacuum plateau (`vacuum_window`) and computes $\Phi$ and the flatness; the report states over which $z$ range it was evaluated.
4. `report_wf`, `export_wf` (`WF.dat` with header `Phi_eV`, `V_vacio_eV`, `E_Fermi_eV`, `planitud_eV` and the profile) and `plot_profile`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $V(\mathbf{r})$ | `potencial.cube` (`pp.x`, `plot_num=11`, Ry) | `fields.read_cube` |
| $E_F$ | `pw.x` XML, tag `fermi_energy` | `qeout.read_xml`, Ha → eV |
| Ry → eV | `qeout.RY_EV` | 13.605693122994 |
| Atomic positions (for `vacuum_window`) | `pw.x` XML (`atomic_positions`) | `qeout.read_xml(...).positions` |
| Axis | `--axis` (c by default) | `_AXES` |

**Limits and pitfalls.**
- Warning if flatness > 0.05 eV: *"la meseta de vacío varía más de 0.05 eV. El vacío es insuficiente o hay un dipolo neto; aumenta el vacío (o usa una losa simétrica)…"*.
- Without positions (use from Python) the plateau is searched blindly around the potential maximum: *"con poco vacío la ventana puede pisar la cola del potencial atómico"* (docstring). The CLI always passes the positions from the XML.
- It does not apply a dipole correction by itself: a polar slab gives two different vacuum levels and this command takes the higher one. For polar slabs the calculation must be generated with `--dipole` (`gen`, `eform`) or `esm` must be used.
- If the XML has no `fermi_energy` (fixed occupations): *"el XML no trae energía de Fermi (¿terminó el scf?)"*.

**References.**
- N. D. Lang, W. Kohn, *Phys. Rev. B* 3, 1215 (1971) — work function in the jellium model.
- L. Bengtsson, *Phys. Rev. B* 59, 12301 (1999) — dipole correction in slabs.

---

### `olla-dft esm` — Charged surfaces with the effective screening medium

**What it answers.** What are the work function, the capacitance and the potential of zero charge of a slab (neutral or charged) without periodic images or the compensating background contaminating the result?

**Background for non-experts.** A charged slab in a periodic cell is an ill-posed problem: QE spreads a uniform background of opposite charge over the whole volume, vacuum included, and the energy depends on the cell size without converging to anything. The **ESM** (Effective Screening Medium) replaces periodicity along $z$ by an explicit boundary condition: the Poisson equation is solved inside the cell and matched to an analytic solution outside. Three variants: `bc1` (vacuum on both sides, neutral slabs; the vacuum level is zero by construction), `bc2` (two metal plates: a capacitor, admits a field) and `bc3` (vacuum/metal: an electrode that receives the counter-charge). With `bc2`/`bc3` the distance to the electrode is no longer a convergence parameter but **physics**: it fixes the capacitance.

**Formulas.** In `qekit/modules/esm.py`.

Centring (`esm.centrar`): $z_i \leftarrow z_i - \tfrac{1}{2}(z_{\min}+z_{\max})$ (ESM measures $z$ from the cell centre).

Vacuum level (`esm.nivel_vacio`): average of $V_{\mathrm{tot}}(z)$ from the `.esm1` in the region $|z| > t/2 + m$, with $t$ the slab thickness and a margin $m$ that starts at `MARGEN_VACIO` = 2 Å and grows in 0.5 Å steps (up to `margen_max` = 8 Å) until the standard deviation of the potential drops below `tol` = 1e-3 eV; with `bc3` only the $z<0$ side.

$$\Phi = V_{\mathrm{vac}} - E_F$$

Capacitance (`esm.capacitancia`), linear fit $q = C' V + b$:
$$C = \frac{dq}{dV}\,\frac{1}{A}\cdot 1.602176634\times10^{3}\quad[\mu\mathrm{F/cm^2}],\qquad R^2 = 1 - \frac{\sum(q-\hat q)^2}{\sum(q-\bar q)^2}$$
- $q$ in e per cell, $V$ in V (eV/e), $A$ = cell area in Å² (`|(\mathbf a\times\mathbf b)_z|`); `E_A2_A_UF_CM2` = $1.602176634\times10^{3}$ converts e/(Å²·V) to µF/cm².

Linearity (`esm.linealidad`): $\max|P - \hat P| / (\max P - \min P) \le$ `tol` = 0.02.

Potential of zero charge (`esm.potencial_de_carga_cero`): linear interpolation of $\Phi(q)$ at $q = 0$.

Grand canonical (`esm.gran_canonico`, library only): $\Omega = E + q\,\Phi$.

**How Olla-DFT computes it.**
1. `esm.comprobar`: rejects `bc1` with charge (*"bc1 es vacío por los dos lados… la energía diverge"*) and cells not orthogonal in $z$; warns if vacuum < `VACIO_MINIMO` = 6 Å, if the slab was not centred and, with `bc2/bc3` and charge, that the vacuum is physics.
2. `esm.prepare` centres the slab, computes thickness, vacuum and area, and writes one `scf` per charge in `q00/`, `q01/`… (`inputgen.build_pw_input`, `conv_thr=1e-8`, `mv` smearing with `degauss=0.02`, grid $n_1\times n_2\times 1$, `tot_charge=q`) and inserts into `&SYSTEM`: `assume_isolated='esm'`, `esm_bc`, `esm_nfit=--nfit` (4), `esm_w=--esm-w` if ≠ 0, `esm_efield=--field` only with `bc2`. Writes `run.sh`.
3. `--run` or by hand: `pw.x` in each folder.
4. `esm.collect` reads from each folder the XML (`total_energy`, `fermi`) and the `<prefix>.esm1` (`esm.leer_esm1`: z (Å), charge (e/Å), $V_H$, $V_{\mathrm{loc}}$, $V_{\mathrm{tot}}$ in eV); `nivel_vacio` and $\Phi$.
5. `esm.report`: table $q, E, E_F, V_{\mathrm{vac}}, \Phi$; with `bc1` checks $|V_{\mathrm{vac}}| < 10^{-3}$ eV; with several charges, capacitance from $V_{\mathrm{vac}}(q)$ (cell voltage) and, if $\Phi(q)$ is linear, also from $\Phi(q)$ with the PZC.
6. `export` (`ESM.dat`, `ESM_perfil_qNN.dat`, `ESM.txt`) and `plot` (profiles and $q$ vs $\Phi$).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $V_{\mathrm{tot}}(z)$, charge$(z)$ | `<prefix>.esm1` written by `pw.x` with ESM | `esm.leer_esm1`, columns 0–4 |
| $E$, $E_F$ | `pw.x` XML | `qeout.read_xml` (Ha → eV) |
| Area $A$ | cell vectors $\mathbf a,\mathbf b$ | `esm.prepare` |
| µF/cm² factor | `esm.E_A2_A_UF_CM2` | $e/(10^{-8}\,\mathrm{cm})^2$ |
| Charges | `--charge` (list) | e per cell |
| Field | `--field` (Ry/a.u.) | `bc2` only |

**Limits and pitfalls.**
- *"Con bc2 o bc3 la capacitancia depende de la distancia al contraelectrodo: es una capacitancia DE ESTE MONTAJE, no una propiedad del material."*
- Energies with net charge are not comparable with each other: *"la energía de ESM incluye la interacción con la carga imagen del electrodo, que crece como q²."*
- If $\Phi(q)$ is not linear: *"Φ(q) = V_vac − E_F NO es una recta… no doy un potencial de carga cero sobre ella."*
- `gran_canonico` (Ω = E + qΦ) exists in the module but **no command uses it**; the "grand canonical" of the module title is not exposed in the CLI.
- The slab is centred automatically; if the user had already centred it at $c/2$ (ASE) the warning explains why it was re-centred.
- The calculation always uses smearing (`insulator=False`): intended for metals/electrodes.

**References.**
- M. Otani, O. Sugino, *Phys. Rev. B* 73, 115407 (2006) — ESM.
- N. Bonnet, T. Morishita, O. Sugino, M. Otani, *Phys. Rev. Lett.* 109, 266101 (2012) — constant potential with ESM.

---

### `olla-dft echem` — Computational hydrogen electrode: HER and OER

**What it answers.** What potential must be applied so that all the steps of hydrogen evolution (HER) or oxygen evolution (OER) become downhill, and how far is it from the equilibrium potential (overpotential)?

**Background for non-experts.** Computing a solvated proton is a very hard problem. The trick of the computational hydrogen electrode (CHE) is to notice that, at 0 V versus the standard hydrogen electrode and pH 0, the pair $\mathrm{H^+ + e^-}$ has the same free energy as $\tfrac12\mathrm{H_2(g)}$, which can be computed. Every step that releases a $(\mathrm{H^+ + e^-})$ is evaluated that way, and the potential $U$ and the pH enter afterwards as additive terms. The step with the largest $\Delta G$ is the "limiting" one: the potential that makes it exergonic is the limiting potential, and its distance to the equilibrium one is the overpotential. This is thermodynamics of intermediates: no kinetic barriers and no solvent.

**Formulas.** In `qekit/modules/echem.py`.

Dependence on $U$ and pH (`Echem.dG`):
$$\Delta G_i(U, \mathrm{pH}) = \Delta G_i(0,0) - eU - k_B T\ln 10\cdot\mathrm{pH} = \Delta G_i(0,0) - e\,U_{\mathrm{RHE}}$$
$$U_{\mathrm{RHE}} = U_{\mathrm{SHE}} + k_B T\ln 10\cdot\mathrm{pH}\quad(\text{`echem.u_rhe`; } 0.0592\,\mathrm{pH\ V\ at\ 298\ K})$$
- $k_B$ = `KB_EV` = $8.617333262\times10^{-5}$ eV/K; $T$ = `--temperature` (298.15 K); $U$ = `-U` in V **versus SHE** (at pH 0 it coincides with RHE); the pH term is exactly the SHE → RHE conversion, so on the RHE scale the $\Delta G$ do not depend on pH. One electron per step.

HER (`echem.her`):
$$\Delta G_{\mathrm{H^*}} = E_{\mathrm{ads}}(\mathrm{H}) + c_{\mathrm{H}},\qquad \text{steps: } (+\Delta G_{\mathrm{H^*}},\ -\Delta G_{\mathrm{H^*}})$$
- $E_{\mathrm{ads}}(\mathrm{H})$: `--her`, referred to $\tfrac12\mathrm{H_2}$ (eV); $c_{\mathrm{H}}$ = ZPE − TΔS = 0.24 eV by default (`CORRECCIONES`).

OER (`echem.oer`), with $G_X = E_{\mathrm{ads}}(X) + c_X$:
$$\Delta G_1 = G_{\mathrm{OH}},\quad \Delta G_2 = G_{\mathrm{O}} - G_{\mathrm{OH}},\quad \Delta G_3 = G_{\mathrm{OOH}} - G_{\mathrm{O}},\quad \Delta G_4 = 4.92\ \mathrm{eV} - (\Delta G_1+\Delta G_2+\Delta G_3)$$
- $c_{\mathrm{OH}} = 0.35$, $c_{\mathrm{O}} = 0.05$, $c_{\mathrm{OOH}} = 0.40$ eV by default; `DG_AGUA_TOTAL` = 4.92 eV (experimental, $2\mathrm{H_2O} \to \mathrm{O_2} + 2\mathrm{H_2}$).

Limiting potential and overpotential (`Echem.U_limitante`, `Echem.sobrepotencial`):
$$U_L = \max_i \Delta G_i(0,0)/e,\qquad \eta = U_L - U_{\mathrm{eq}},\quad U_{\mathrm{eq}}^{\mathrm{OER}} = 1.229\ \mathrm{V},\ U_{\mathrm{eq}}^{\mathrm{HER}} = 0$$
- $\eta$ is returned **with sign**: positive = at $U_{\mathrm{eq}}$ the limiting step is still uphill (with the profiles built here it never comes out negative; it only could with a `dG_total` different from the experimental one).

Scaling relation (`echem.escala_ooh_oh`, OER) and its limit (`echem.sobrepotencial_minimo_escala`):
$$\Delta_{\mathrm{sc}} = G_{\mathrm{OOH}} - G_{\mathrm{OH}}\ \text{(compared with `ESCALA_OOH_OH` = 3.2 ± 0.2 eV)},\qquad \eta_{\min} = \frac{\Delta_{\mathrm{sc}}}{2} - \frac{\Delta G_{\mathrm{total}}}{4} = 0.37\ \mathrm{V}$$
- If OOH* and OH* are separated by a fixed $\Delta_{\mathrm{sc}}$, steps 2 and 3 add up to $\Delta_{\mathrm{sc}}$ and the worse one cannot drop below $\Delta_{\mathrm{sc}}/2$ = 1.6 eV; against 4.92/4 = 1.23 V that leaves ~0.37 V.

Pourbaix-like grid (`echem.pourbaix`, library only): $\Delta G_{\lim}(U,\mathrm{pH}) = \max_i\Delta G_i(0,0) - eU - k_BT\ln10\cdot\mathrm{pH}$ over $U\in[-0.5,2]$ V and pH $\in[0,14]$.

**How Olla-DFT computes it.**
1. `_cmd_echem` requires exactly one of `--her E` or `--oer OH=..,O=..,OOH=..`; `--corrections X=eV` overrides the thermal corrections.
2. `echem.her` or `echem.oer` build the list of steps $(\text{name}, \Delta G_i)$; `oer` warns if $\Delta G_4 < 0$ and if tabulated corrections were used.
3. The user's `U` (vs SHE) and `pH` are set; `echem.report` also prints $U_{\mathrm{RHE}}$ (`Echem.U_rhe`) if pH ≠ 0, tabulates $\Delta G(0)$ and $\Delta G(U,\mathrm{pH})$, the limiting step, $U_L$ (vs RHE), signed $\eta$, the descriptor $\Delta G_{\mathrm{H^*}}$ (HER) or the scaling relation and the $\eta_{\min}$ it imposes (OER).
4. `export` writes `ECHEM.dat` and `ECHEM.txt`; `plot` draws the staircase diagram at $U = 0$, $U_{\mathrm{eq}}$ and $U_L$.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E_{\mathrm{ads}}$ of H, OH, O, OOH | parameters `--her`, `--oer` | eV, referred to H₂O and ½H₂ (from `adsorb`) |
| ZPE − TΔS corrections | `echem.CORRECCIONES` or `--corrections` | H 0.24, OH 0.35, O 0.05, OOH 0.40 eV (Nørskov et al.) |
| Total $\Delta G$ of water | `echem.DG_AGUA_TOTAL` | 4.92 eV, experimental |
| $U_{\mathrm{eq}}$ | `echem.U_EQ_OER`, `U_EQ_HER` | 1.229 V, 0 V |
| $k_B$ | `echem.KB_EV` | $8.617333262\times10^{-5}$ eV/K (CODATA) |
| Universal $\Delta_{\mathrm{sc}}$ | `echem.ESCALA_OOH_OH` | 3.2 eV (Man et al. 2011) |

**Limits and pitfalls.**
- *"El CHE es termodinámica de intermedios: NO hay barreras cinéticas, ni disolvente explícito, ni doble capa."*
- `-U` is versus the **SHE** (CLI help: *"a pH 0 es el mismo que frente al RHE; el pH lo convierte"*); $U_L$ and $\eta$ are on the RHE scale. For the HER $U_L = |\Delta G_{\mathrm{H^*}}|$, so $\eta \ge 0$ always.
- Fourth step by difference: *"El cuarto paso sale NEGATIVO… o hay un error en las referencias, o tu superficie liga los intermedios muchísimo."*
- `pourbaix()` is not wired to any command: the "Pourbaix diagram" of the module title is not produced from the CLI.

**References.**
- J. K. Nørskov, J. Rossmeisl, A. Logadottir, L. Lindqvist, J. R. Kitchin, T. Bligaard, H. Jónsson, *J. Phys. Chem. B* 108, 17886 (2004) — CHE. DOI: 10.1021/jp047349j.
- J. K. Nørskov, T. Bligaard, A. Logadottir, J. R. Kitchin, J. G. Chen, S. Pandelov, U. Stimming, *J. Electrochem. Soc.* 152, J23 (2005) — HER volcano.
- I. C. Man et al., *ChemCatChem* 3, 1159 (2011) — OER scaling relation.

---

### `olla-dft adsorb` — Adsorption sites and adsorption energy

**What it answers.** On which inequivalent sites of a surface can a molecule sit, and how much does the system gain (or lose) by doing so on each?

**Background for non-experts.** A molecule on a metal sits on top of an atom (*top*), over the midpoint between two (*bridge*) or over the centre of a triangle of atoms (*hollow*; on fcc(111) there are two: with or without an atom underneath in the second layer). Many of those sites are copies by symmetry, so they are grouped by their "fingerprint": the sorted list of distances to their 24 nearest neighbours counting all layers. The adsorption energy is a subtraction of three total energies that only makes sense if the three calculations share cell, cutoffs, k grid and pseudos; that is why they are generated together.

**Formulas.** `thermochem.adsorcion` (called from `AdsorbRun.energias_ads`):
$$E_{\mathrm{ads}} = E(\text{slab}+\text{mol}) - E(\text{slab}) - n\,E(\text{mol})$$
- All in eV; $n$ = number of molecules (`n_mol`, 1). Negative = favourable.

Geometry after relaxation (`adsorb.collect`):
$$h = \min_{a\in\mathrm{ads}} z_a - \max_{s\in\mathrm{slab}} z_s,\qquad d_{\mathrm{contact}} = \min_{a,s}|\mathbf r_a - \mathbf r_s|$$

Site fingerprint (`adsorb._huella`): sorted distances to the $k$ = `N_VECINOS_HUELLA` = 24 nearest atoms (with periodic replicas); two sites are the same if $\max|\Delta d| <$ `TOL_HUELLA` = 0.05 Å.

**How Olla-DFT computes it.**
1. `adsorb.prepare` requires vacuum along $c$ (`kpoints.direcciones_con_vacio`) and loads the molecule (`cargar_molecula`: file or ASE G2 database).
2. `adsorb.sitios`: exposed layer = atoms within `TOL_CAPA` = 0.6 Å of the extreme $z$; *top* over each of them; *bridge* between pairs closer than `R_VECINO` = 3.6 Å; *hollow* at the centroids of the Delaunay triangulation (triangles with a side > 1.6·3.6 Å are discarded); they are brought into the cell and deduplicated by fingerprint; labelled `top1`, `bridge1`, `hollow1`…
3. With `--rotations N` and a polyatomic molecule, each site is repeated with rotations of $360k/N$ degrees around $z$.
4. `sweep.prepare_common` is resolved on the **union** slab + molecule (same pseudos and cutoffs for everything). `_losa/`, `_molecula/` (molecule centred in the **same** cell) and one folder per site are written (`adsorb.colocar`: atom `--anchor` at `--height` = 2.0 Å above the site), all `relax` unless `--fixed-ions`; `run.sh`. With `--dipole`, `dipole_correction=3` enters the **three** calculations (`inputgen.build_pw_input`: `tefield`, `dipfield`, `edir=3`, `emaxpos`/`eopreg` at the centre of the vacuum gap via `inputgen._region_vacio`, `eamp=0`); without ≥ 5 Å of vacuum it aborts.
5. `--run`/`--collect`: `adsorb.collect` reads the XMLs (`qeout.read_xml`), energies, convergence, height and contact.
6. `adsorb.report`: table sorted by $E_{\mathrm{ads}}$, best site, diagnosis by ranges (>0: does not bind; > −0.30 eV: weak physisorption; < −2 eV: probable reaction/dissociation or atomic chemisorption), difference with the second (< 50 meV: indistinguishable). `export`: `ADSORCION.dat/.txt`; `plot`: bars.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E$(slab), $E$(mol), $E$(slab+mol) | `pw.x` XML in `_losa/`, `_molecula/`, `<site>/` | `total_energy` (Ha → eV) |
| Relaxed positions | XML (`atomic_positions`) | height and contact |
| Molecule | file or `ase.build.molecule` | `--mol` |
| Neighbour radius, tolerances | `adsorb.R_VECINO`, `TOL_CAPA`, `N_VECINOS_HUELLA`, `TOL_HUELLA` | 3.6 Å, 0.6 Å, 24, 0.05 Å |
| Initial height | `--height` | 2.0 Å |
| vdW correction | `--vdw` | passed to `inputgen.build_pw_input` |
| Dipole correction | `--dipole` | `dipole_correction=3` in slab, molecule and slab+molecule |

**Limits and pitfalls.**
- Without `--vdw`: *"AVISO: sin corrección de van der Waals. En fisisorción… la energía sale cerca de cero y la geometría desligada."*
- Without `--dipole` on the `top` face: *"Sugerencia: una molécula adsorbida en una sola cara deja la losa polar. Con --dipole se cancela el dipolo artificial a través del vacío."* The sawtooth is put into all three calculations on purpose: *"si la referencia se calcula sin corregir, la resta arrastra el error."*
- $E_{\mathrm{ads}} > 0$ with fixed ions: *"lo más probable es que la altura inicial… no sea la de equilibrio y estés midiendo la repulsión."*
- The reference is the molecule exactly as given: with `--mol H` the reference is the **atom**, not ½H₂ (the report warns about it for $|E_{\mathrm{ads}}| > 2$ eV).
- The isolated molecule is computed with the same k grid as the slab (deliberate consistency, not a separate box).
- Site enumeration is geometric: it does not detect sites over second layers or reconstructions.

**References.**
- B. Hammer, J. K. Nørskov, *Adv. Catal.* 45, 71 (2000) — adsorption on metal surfaces.
- S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* 132, 154104 (2010) — DFT-D3.

---

### `olla-dft surface` — Cutting an (hkl) slab with vacuum

**What it answers.** Given a crystal, what does the $(hkl)$ surface slab with $N$ layers and vacuum look like, is it symmetric, is it polar, and how much real vacuum is left?

**Background for non-experts.** A surface is simulated with a "slab": a few atomic layers parallel to the $(hkl)$ plane and, above them, enough vacuum so that the slab does not see its periodic copy. If the two faces are not the same (*polar* slab), an artificial dipole appears across the vacuum and shifts the work functions; QE corrects it with `dipfield`. The vacuum that matters is the one between atoms, not between cell edges.

**Formulas.** `builder.surface`:
$$t = z_{\max} - z_{\min},\qquad v_{\mathrm{real}} = c - t$$
- Symmetric: the sorted profile $z_i - \bar z$ coincides with its mirror within `tol` = 0.3 Å. Polar: composition of the top layer ≠ that of the bottom layer (atoms within `tol` of the extreme).

**How Olla-DFT computes it.**
1. `structure.conventional` → `ase.build.surface(base, miller, layers, vacuum=vacuum/2, periodic=True)` and `slab.center(vacuum=vacuum/2, axis=2)`.
2. `builder.surface` computes thickness, real vacuum, number of atomic planes (`_planos_z`, tolerance 0.3 Å), symmetry and polarity; with `--fix N` it marks the atoms of the $N$ lowest planes in two ways (`_fijar_capas`): the array `slab.arrays['qekit_fijo']` and an ASE `FixAtoms` constraint. `inputgen.fixed_atoms` reads either and writes `0 0 0` in the third column of `ATOMIC_POSITIONS`.
3. Warnings: > 1.5 atoms per plane (cell is a multiple of the minimal one), real vacuum < 10 Å, polar slab, < 4 layers, freezing all the planes.
4. `report_slab` and, with `-o`, `structure.convert` writes CIF/POSCAR/XYZ. If there are fixed atoms and the format does not keep them (`structure.conserva_fijos`: only POSCAR/CONTCAR/`.vasp` store them as *Selective dynamics*), the CLI warns and recommends `builder.FORMATO_CON_FIJOS` (POSCAR or `.vasp`) or `olla-dft gamma --fix`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Slab | `ase.build.surface` on the conventional cell | `--miller`, `--layers` (6), `--vacuum` (15 Å) |
| Conventional cell | spglib via `structure.conventional` | reference for the hkl indices |
| Atomic planes | distinct $z$ heights (tol 0.3 Å) | `builder._planos_z` |

**Limits and pitfalls.**
- *"la losa es POLAR… Añade 'dipfield = .true.' y 'edir = 3' al input, o corta una losa simétrica."*
- `--fix` is lost when exporting to CIF or XYZ: *"el CIF no tiene dónde ponerlo, así que al volver a cargarlo se relajaría todo. Escribe la losa en POSCAR (o .vasp)…"*. Only POSCAR keeps the `FixAtoms` constraint, which `inputgen.fixed_atoms` translates into `0 0 0`.
- Polarity detection only compares compositions of the extreme layers: a slab with terminations of the same composition but different geometry is not flagged.
- Cutting on the conventional cell may give a surface cell larger than the minimal one (a warning is issued).

**References.**
- P. W. Tasker, *J. Phys. C* 12, 4977 (1979) — polar surfaces.
- ASE: A. H. Larsen et al., *J. Phys.: Condens. Matter* 29, 273002 (2017).

---

### `olla-dft defect` — Building a point defect

**What it answers.** What do the perfect supercell and the supercell with a vacancy, a substitution or an interstitial look like, and what is the formation-energy formula that will have to be filled in?

**Background for non-experts.** A point defect is modelled by repeating the primitive cell $n_1\times n_2\times n_3$ times and modifying one atom. The supercell must be large so that the defect does not interact with its periodic images. This command only builds the two structures and writes the formula with its terms; `eform` does the calculation.

**Formulas.** `builder.formation_energy_text` writes:
$$E_f = E(\text{defect}) - E(\text{perfect}) \pm \mu(\cdot)\ \ [+\,q(E_F + E_v) + E_{\mathrm{corr}}]$$
- vacancy: $+\mu(\text{species that leaves})$; substitution: $+\mu(\text{leaves}) - \mu(\text{enters})$; interstitial: $-\mu(\text{enters})$.

**How Olla-DFT computes it.**
1. `structure.primitive` → `repeat(supercell)` (default 2×2×2).
2. `builder.defect`: vacancy (`del d[site]`), substitution (`d[site].symbol = new`), interstitial (fractional position `--position` of the supercell; warns if it ends up < 1.0 Å from a neighbour, minimum-image distance).
3. Warning if the shortest side of the supercell < 10 Å.
4. `report_defect` and writing of `perfecto.cif` and `defecto.cif` in `--outdir`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Primitive cell | spglib via `structure.primitive` | basis of the supercell |
| Site, species, position | `--site`, `--new-element`, `--position` | 0-based indices in the supercell |

**Limits and pitfalls.** *"la supercelda mide X Å en su lado más corto: el defecto se ve con sus imágenes periódicas. Para energías de formación conviene ≥ 10-12 Å."* It relaxes nothing and computes no energies; the `--site` index refers to the repeated supercell, not to the input crystal.

**References.** C. Freysoldt, B. Grabowski, T. Hickel, J. Neugebauer, G. Kresse, A. Janotti, C. G. Van de Walle, *Rev. Mod. Phys.* 86, 253 (2014).

---

### `olla-dft eform` — Formation energy of charged defects

**What it answers.** How much does it cost to form the defect in each charge state, how does it vary with the Fermi level, where are the charge-transition levels and what is the finite-size correction?

**Background for non-experts.** Forming a defect costs an energy that depends on three things: where the atoms come from or go to (chemical potential $\mu$, fixed by the synthesis conditions), where the electrons come from or go to (Fermi level $\varepsilon_F$, measured from the valence-band maximum) and an artefact: a charged periodic cell interacts with its own images and with the neutralising background that QE adds. That artefact is corrected with the electrostatic energy of a point charge in a lattice of image charges (Makov–Payne) screened by the dielectric constant, or with the Lany–Zunger version that includes a shape term. The point where two lines $E_f(q)$ cross is a transition level: the Fermi level at which the defect changes charge.

**Formulas.** In `qekit/modules/defects.py`.

Formation energy (`DefectRun.E_f`):
$$E_f[D^q](\varepsilon_F) = E[D^q] - E[\mathrm{perf}] - \sum_i n_i\mu_i + q\,(\varepsilon_{\mathrm{VBM}} + \varepsilon_F) + E_{\mathrm{corr}}(q) + q\,\Delta V$$
- $n_i$: atoms **added** of species $i$ (−1 for the one that leaves); $\mu_i$: `--mu EL=eV` (for an elemental crystal, $\mu = E[\mathrm{perf}]/N$ automatically, `asignar_mu_elemental`); $\varepsilon_{\mathrm{VBM}}$ = `highestOccupiedLevel` of the perfect supercell (eV); $\varepsilon_F \in [0, E_g]$; $\Delta V$: potential alignment (`--dv` or `--align`).

Madelung constant by Ewald summation (`defects.madelung_xi`, `constante_madelung`):
$$\xi = \sum_{\mathbf R\neq 0}\frac{\mathrm{erfc}(\eta R)}{R} + \frac{4\pi}{V}\sum_{\mathbf G\neq 0}\frac{e^{-G^2/4\eta^2}}{G^2} - \frac{2\eta}{\sqrt\pi} - \frac{\pi}{\eta^2 V},\qquad \alpha_M = -\xi\,L,\quad L = V^{1/3}$$
- $\eta = \sqrt\pi / V^{1/3}$; real- and reciprocal-space cutoffs set by `tol` = 1e-10. Gives $\alpha_M = 2.8372974$ for the simple cubic lattice.

Image correction (`defects.correccion_imagen`):
$$E_{\mathrm{MP}} = \frac{k_e\,q^2\,\alpha_M}{2\,\varepsilon\,L},\qquad E_{\mathrm{LZ}} = E_{\mathrm{MP}}\left[1 + c_{\mathrm{sh}}\left(1 - \frac{1}{\varepsilon}\right)\right]$$
- $k_e$ = `KE` = 14.399645 eV·Å; $\varepsilon$ = `--epsilon`; $c_{\mathrm{sh}}$ = `C_SHAPE` = −0.35 (single value; LZ give −0.369 sc, −0.343 fcc, −0.342 bcc). `--correction` ∈ {`ninguna`, `makov-payne`, `lany-zunger`}.

Alignment (`defects.alineamiento`): $\Delta V = f\,\langle \bar V_{\mathrm{def}}(z) - \bar V_{\mathrm{perf}}(z)\rangle$ averaged over the 25 % of the cell opposite to the point of largest $|\Delta V - \mathrm{median}|$, with its standard deviation; $f$ = `UNIDADES_POTENCIAL[unidades_cube]` converts the `pp.x` cubes (`plot_num=11`, Ry, `unidades_cube="Ry"` by default, $f$ = `RY_EV`) to eV; with `"eV"`, $f = 1$. The result (`dV`, `sigma`, `perfil`) is always in eV.

Transition levels (`defects.niveles_transicion`), one entry per pair of consecutive charges $a<b$ (sorted by $q$), with the flag `dentro` = $0 \le \varepsilon \le E_g$; it is **not** filtered by the lower envelope (for the observable levels cross it with `envolvente`):
$$\varepsilon(a/b) = \frac{E_f(a, 0) - E_f(b, 0)}{b - a}$$

**How Olla-DFT computes it.**
1. `defects.prepare`: requires `--epsilon` if there are charges ≠ 0 and correction ≠ `ninguna`; builds the cells with `builder.defect`; resolves pseudos on the union of species.
2. Parity: if `--insulator` and some charge state leaves an odd number of electrons (`defects.electrones` with the `z_valence` of the UPFs), it activates `nspin=2` in **all** states with `tot_magnetization` 1 (odd) or 0 (even).
3. Writes `_perfecto/` (scf) and `qm1/`, `qp0/`, `qp1/`… (`relax` unless `--fixed-ions`, `tot_charge=q`, same estimated `nbnd`) and `run.sh`.
4. `--run`/`--collect`: `defects.collect` reads energies, convergence, `homo` (VBM) and `lumo` (gap) of the perfect cell; `--mu`; `--align POT_DEF POT_PERF` or `--dv`.
5. `report`: table $q$, $E$, $E_{\mathrm{corr}}$, $E_f(\varepsilon_F=0)$, $E_f(\varepsilon_F=E_g)$; transition levels, flagging those outside the gap; lower envelope (`envolvente`) and stable charges across the gap. `export`: `FORMACION.dat` (table and $E_f(\varepsilon_F)$ at 51 points); `plot`: $E_f$ vs $\varepsilon_F$.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E[D^q]$, $E[\mathrm{perf}]$ | `pw.x` XML | `total_energy` |
| $\varepsilon_{\mathrm{VBM}}$, $E_g$ | XML of the perfect supercell | `highestOccupiedLevel`, `lowestUnoccupiedLevel` |
| $\mu_i$ | `--mu` or $E[\mathrm{perf}]/N$ (elemental) | eV/atom |
| $\varepsilon$ | `--epsilon` | e.g. $\varepsilon_1(0)$ from `optics` |
| $\alpha_M$ | Ewald sum over the real cell | `defects.madelung_xi` |
| $k_e$, $c_{\mathrm{sh}}$ | `defects.KE`, `defects.C_SHAPE` | 14.399645 eV·Å, −0.35 |
| $\Delta V$ | two potential `.cube` files (`--align`, Ry → eV) or `--dv` (already in eV) | `defects.alineamiento`, `UNIDADES_POTENCIAL` |
| Electrons per cell | `z_valence` of the UPFs | `defects.electrones` |

**Limits and pitfalls.**
- Without `--epsilon`: *"la constante dieléctrica es lo que apantalla la interacción del defecto con sus imágenes; sin ella la corrección sale ε veces de más."*
- With `--correction ninguna`: *"SIN CORREGIR: las E_f de los estados cargados están sistemáticamente bajas, y el error crece con q²."*
- `--dv` is given directly in eV (not converted); `--align` assumes `pp.x` cubes in Ry and the report says so: *"entra en E_f como q·ΔV = … eV por unidad de carga (el potencial de pp.x viene en Ry y se pasó a eV)"*. If $\sigma_{\Delta V} > 0.3\,|\Delta V|$: *"el defecto todavía se nota en la zona 'lejana', o sea que la supercelda es pequeña."*
- The listed transition levels include crossings between states that are never the most stable; the report marks *"<< fuera del gap"* those outside $[0, E_g]$, but a level inside the gap between two states that are not on the envelope is not observable either.
- Without a VBM (metal, no empty bands): *"No pude leer el VBM… E_f de los estados cargados no está definida."*
- Missing $\mu$ in a compound: *"FALTA el potencial químico… las DIFERENCIAS entre cargas y los niveles de transición sí valen, el valor absoluto de E_f no."*
- The correction only removes the leading $\propto q^2/L$ term; side < 10 Å with charge: warning.

**References.**
- G. Makov, M. C. Payne, *Phys. Rev. B* 51, 4014 (1995).
- S. Lany, A. Zunger, *Phys. Rev. B* 78, 235104 (2008); *Modelling Simul. Mater. Sci. Eng.* 17, 084002 (2009).
- C. Freysoldt, J. Neugebauer, C. G. Van de Walle, *Phys. Rev. Lett.* 102, 016402 (2009).
- C. Freysoldt et al., *Rev. Mod. Phys.* 86, 253 (2014). DOI: 10.1103/RevModPhys.86.253.

---

### `olla-dft interface` — Heterostructures and lattice mismatch

**What it answers.** Which common supercell allows stacking two 2D materials (or two slabs) with the least possible strain, how large is that strain and what does the initial structure look like?

**Background for non-experts.** Two crystal lattices almost never fit. To put them in the same periodic cell one must look for integer multiples of the vectors of each that resemble each other and stretch one of the two. That strain is the number that decides whether the calculation describes the material or a stretched version of it: 1 % is tolerable, 8 % is already another material.

**Formulas.** In `qekit/modules/interface.py`.

Candidate supercells (`_celdas_candidatas`): $\mathbf A' = M\mathbf a$, $\mathbf B' = N\mathbf b$ with $M, N \in \mathbb Z^{2\times2}$, $|M_{ij}|,|N_{ij}| \le$ `--max-index` (4), $\det > 0$, grouped by determinant (the areas must match within $2\cdot$`tol`).

Strain (`_deformacion`):
$$\boldsymbol\epsilon = B'^{-1}A' - I,\qquad \epsilon_{\max} = \max_{ij}|\epsilon_{ij}| \le \texttt{--tol}\ (0.05)$$

Lagrange–Gauss reduction (`reducir_2d`) so as not to repeat the same lattice with different bases; tie-break by "simplicity" of $M, N$ (`_simplicidad`: sum of |entries|, maximum, negatives, non-zeros).

Initial separation (`separacion_vdw`): $d_0 = 0.85\,(r_1 + r_2)$ with van der Waals radii from `R_VDW` (Bondi; 2.0 Å if missing).

With `--strain both`: target cell $= (w A' + v B')/(w+v)$ with $w = n_1\,|\det \mathbf a|$, $v = n_2\,|\det\mathbf b|$.

**How Olla-DFT computes it.**
1. `interface.buscar`: enumerates, filters by atoms (`--max-atoms` 200) and strain, deduplicates by $(n_1, n_2, \text{reduced shape}, \epsilon_{\max})$, sorts by $(\epsilon_{\max}, N_{\mathrm{at}}, \text{simplicity})$ and returns the `--top` (10) best. `--list` only prints them.
2. `interface.emparejar` chooses `--index` and `construir`: `ase.build.make_supercell` for each material, the in-plane cell is taken to the target dragging fractional positions (`_supercelda_deformada`), material 2 is stacked at `--separation` (or $d_0$) above material 1, `--shift` is applied (fractions of the common cell), `--vacuum` (20 Å) is added and the cell is centred.
3. Warnings: $\epsilon_{\max} > 3\,\%$, vdW separation as a starting point, registry not optimised.
4. `export`: `<name>.cif` and `<name>.txt`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| In-plane vectors | cells of `file1`, `file2` (rows 0–1, columns 0–1) | `interface._plano` |
| vdW radii | table `interface.R_VDW` | Å; `R_VDW_DEFECTO` = 2.0 |
| Search limits | `--max-index`, `--tol`, `--max-atoms` | 4, 0.05, 200 |

**Limits and pitfalls.**
- The **largest component** $\max|\epsilon_{ij}|$ of the matrix is reported, not a norm or an average: *"una deformación de 0 % en una dirección y 6 % en la otra no es '3 %'."*
- *"La deformación es del X %. Por encima de ~3 % no se está modelando el material sino una versión estirada de él."*
- The separation is a starting point: *"con un funcional sin corrección de dispersión la distancia de equilibrio saldrá demasiado grande."*
- *"El REGISTRO… no está optimizado. Dos apilamientos distintos pueden diferir en decenas de meV por átomo."*
- $c$ is assumed to be the normal and the cell a slab; the actual strain with `--strain both` does not coincide with the reported $\boldsymbol\epsilon$ (which is that of taking B to A).

**References.**
- A. Bondi, *J. Phys. Chem.* 68, 441 (1964) — van der Waals radii.
- P. Lazić, *Comput. Phys. Commun.* 197, 324 (2015) — CellMatch, lattice matching.

---

### `olla-dft neb` — Reaction barriers with neb.x

**What it answers.** What is the minimum-energy path between reactant and product and how high is the activation barrier (forward and backward)?

**Background for non-experts.** Between two energy minima there is a "mountain pass": the transition state. The nudged elastic band (NEB) stretches a chain of images between reactant and product, joined by springs, and relaxes each image perpendicular to the path until the chain rests in the valley. The climbing image (CI) pushes the highest image up to the exact pass; without it the barrier is underestimated.

**Formulas.** In `qekit/modules/neb.py`, `neb.collect`:
$$E_a^{\rightarrow} = E_{\max} - E_1,\qquad E_a^{\leftarrow} = E_{\max} - E_N,\qquad \Delta E = E_N - E_1$$
- Energies in eV relative to the first image (column 2 of `<prefix>.dat`); if `neb.out` contains `activation energy (->)`/`(<-)`, those are used. Conversion to kJ/mol: × 96.485.

**How Olla-DFT computes it.**
1. `neb.comprobar_extremos`: same number and **order** of atoms, same cell (tol 1e-4), non-identical structures; if it fails, it aborts.
2. `neb.build_neb_input` writes `neb.in`: `&PATH` with `string_method='neb'`, `nstep_path=--nstep` (50), `ds=1`, `opt_scheme='broyden'`, `num_of_images=--images` (7), `k_max=0.3`, `k_min=0.2`, `CI_scheme='auto'` (or `'no-CI'` with `--no-ci`), `path_thr=--path-thr` (0.05 eV/Å); `pw.x` engine trimmed from `inputgen.build_pw_input` (no positions or cell); `FIRST_IMAGE`/`LAST_IMAGE` in Å with `0 0 0` on the `--fix` atoms; `CELL_PARAMETERS`.
3. The user runs `neb.x -inp neb.in > neb.out`.
4. `neb.collect --collect`: reads `<prefix>.dat` (s, E, F), `<prefix>.int` (interpolation), and from `*.out` the barriers, convergence (`convergence achieved`), iterations, `CI_scheme` and the images with *"scf convergence NOT achieved on image"*.
5. `report`: barriers, table per image, warning if the interpolated maximum falls more than 0.4 steps from any image; `export` (`NEB.dat`, `NEB.txt`); `plot`.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $s$, $E$, $F$ per image | `<prefix>.dat` from `neb.x` | `neb.collect` |
| Interpolated curve | `<prefix>.int` from `neb.x` | optional |
| Barriers, convergence, iterations, CI | `neb.out` (regex) | take precedence over the own computation |
| eV → kJ/mol | 96.485 in `neb.report` | — |

**Limits and pitfalls.**
- *"Esta barrera es ELECTRÓNICA, a 0 K y sin energía de punto cero."* Thermal corrections in `thermochem`.
- Without CI: *"esta barrera es una COTA INFERIOR."* Few images (< 5): warning.
- Images with unconverged scf: *"El scf NO convergió en la(s) imagen(es)…: por eso el perfil sale dentado."*
- The endpoints must be relaxed with the same parameters; the module does not check this.

**References.**
- G. Henkelman, B. P. Uberuaga, H. Jónsson, *J. Chem. Phys.* 113, 9901 (2000) — climbing-image NEB. DOI: 10.1063/1.1329672.
- G. Henkelman, H. Jónsson, *J. Chem. Phys.* 113, 9978 (2000) — improved tangent.

---

### `olla-dft amorphous` — Amorphous solid by melt-quench with an MLIP

**What it answers.** How to generate an amorphous structure of given composition and density, and what coordination and first-neighbour distances does it have?

**Background for non-experts.** A glass is not drawn: it is manufactured by heating the material until it melts and cooling it so fast that it has no time to crystallise. On the computer the quench is millions of times faster than in the laboratory, so the result is more disordered and somewhat less dense than the real one. Here the dynamics is done with a machine-learned interatomic potential (MACE by default), not with DFT, because thousands of steps are needed; the resulting structure is a starting point that must then be relaxed with `pw.x`.

**Formulas.** In `qekit/modules/amorphous.py`.

Edge of the cubic cell (`celda_para_densidad`) and density (`densidad_de`):
$$L = \left(\frac{\sum_i m_i\,u}{\rho}\right)^{1/3}\times 10^{8},\qquad \rho = \frac{\sum_i m_i\,u}{V}$$
- $m_i$ in amu; $u = 1.66053906660\times10^{-24}$ g; $\rho$ in g/cm³; $V$ in Å³ (× $10^{-24}$ cm³).

Quench rate (`Protocolo.velocidad_temple`):
$$\dot T = \frac{T_{\mathrm{melt}} - T_{\mathrm{final}}}{N_{\mathrm{quench}}\,\Delta t}$$
- By default $(3000 - 300)\,\mathrm{K}/(1000\times 1\ \mathrm{fs}) = 2.7\times10^{15}$ K/s.

Coordination (`coordinaciones`): $Z_{ab} = \frac{1}{N_a}\sum_{i\in a}\#\{j\in b: d_{ij} < 1.25\,(r_a^{\mathrm{cov}} + r_b^{\mathrm{cov}})\}$ with minimum image; mean first-neighbour distance with the same cutoff (`distancia_media`).

**How Olla-DFT computes it.**
1. `formula_a_simbolos` expands `SiO2` × `--units` (8).
2. `empaquetar` places atoms at random (seed `--seed`) rejecting distances < `--min-dist` × (sum of covalent radii), `FACTOR_MINIMO` = 0.75; up to 20000 attempts per atom; error if they do not fit.
3. `fundir_y_templar` (unless `--pack-only`): calculator `mlip.calculator(--model)`; Maxwell–Boltzmann velocities at `--melt` (3000 K); ASE `Langevin` with `friction=0.02` and `--dt` (1 fs); `--melt-steps` (500) at $T_{\mathrm{melt}}$; quench in 20 segments of $N_{\mathrm{quench}}/20$ steps lowering the thermostat temperature linearly to `--final` (300 K); `--anneal-steps` (200) at $T_{\mathrm{final}}$. $E$ and $T$ are recorded every 10 steps (`traza.dat`).
4. Warnings: final temperature $> 2.5\,T_{\mathrm{final}} + 200$ K (the thermostat did not follow the ramp) and $\dot T > 10^{13}$ K/s.
5. `report` (density, protocol, coordinations, distances, final $T$) and `export` (`amorfo.cif`, `AMORFO.dat`, `AMORFO.txt`).

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| Masses and covalent radii | `ase.data.atomic_masses`, `covalent_radii` | — |
| amu → g | local constant $1.66053906660\times10^{-24}$ | CODATA 2018 |
| Energies and forces | MLIP potential (`mlip.calculator`) | MACE-MP-0 small, CHGNet or M3GNet |
| Target density | `--density` | g/cm³ |
| Protocol | `--melt`, `--final`, `--melt-steps`, `--quench-steps`, `--anneal-steps`, `--dt` | K, steps, fs |

**Limits and pitfalls.**
- *"Esta estructura viene de un potencial aprendido, NO de DFT… relájala con 'olla-dft gen -p relax'… y compara varias realizaciones (--seed distintas)."*
- The default protocol is an **exploration** one: 2.7×10¹⁵ K/s, and the report warns about it (*"Velocidad de temple X K/s. Un vidrio de verdad se enfría a 1-100 K/s"*). The docstring and the `--quench-steps` help say so: 27 000 steps bring it down to 10¹⁴ K/s, ten times more to 10¹³ K/s, where the warning disappears.
- NVT dynamics at fixed volume: the final density is the imposed one, it is not relaxed.
- With `friction=0.02` and fast ramps the system may end up liquid: *"El sistema acabó a X K, no a los Y K pedidos."*
- Requires `torch` + the model package (not dependencies of Olla-DFT).

**References.**
- I. Batatia et al., *MACE-MP-0* (arXiv:2401.00096, 2023).
- ASE Langevin: A. H. Larsen et al., *J. Phys.: Condens. Matter* 29, 273002 (2017).

---

### `olla-dft mlip` — Pre-relaxation, volume scan and phonon screening with a machine-learned potential

**What it answers.** Before spending DFT: what is a nearly relaxed geometry, where approximately is the $E(V)$ minimum, and does the structure have imaginary frequencies?

**Background for non-experts.** A machine-learned interatomic potential (MLIP) gives energies and forces thousands of times cheaper than DFT. It does not replace `pw.x` — it is trained on PBE data from Materials Project and describes *another* energy surface — but it helps reach the DFT calculation with the geometry almost ready, bound the range of an equation of state, and detect before DFPT that a structure is not at a minimum.

**Formulas.** In `qekit/modules/mlip.py`.

Relaxation (`mlip.relax`): ASE BFGS until $f_{\max} <$ `--fmax` (0.01 eV/Å) or `--steps` (300), with `FrechetCellFilter` if the cell is relaxed. Pressure:
$$P = -\tfrac{1}{3}\,\mathrm{tr}\,\boldsymbol\sigma\times 160.21766208\ \ [\mathrm{GPa}]$$

Volume scan (`mlip.volume_scan`): 15 scales in $[1-s, 1+s]$, $s$ = `--span` (0.10); parabola $E = aV^2 + bV + c$:
$$V_0 = -\frac{b}{2a},\qquad B_0 \approx 2aV_0\times160.21766208\ \mathrm{GPa},\qquad \text{scale} = (V_0/V)^{1/3}$$

Finite-difference phonons (`mlip.phonon_check`, `frequencies`): Hessian $H_{i\alpha,j\beta} = -\partial F_{j\beta}/\partial u_{i\alpha}$, central differences with $\delta$ = 0.01 Å in a `--supercell` (2×2×2), symmetrised; dynamical matrix $D = H/\sqrt{m_im_j}$;
$$\omega = \mathrm{sign}(\lambda)\sqrt{|\lambda|}\times 521.4708\ \mathrm{cm^{-1}}$$
- $\lambda$: eigenvalues of $D$ in eV/(Å²·amu); imaginary if $\omega < -5$ cm⁻¹.

**How Olla-DFT computes it.**
1. `mlip.calculator` loads MACE (`mace_mp(model=--size, default_dtype='float64')`), CHGNet or M3GNet; if the package is missing it explains what to install.
2. `relax`: initial/final forces and pressure, maximum displacement, volume change; warnings if it does not converge or if some atom moved > 0.5 Å. Writes the structure (`relajado_mlip.cif`) and `MLIP_PROCEDENCIA.json` (`write_provenance`) so that `audit` knows it is not DFT.
3. `scan`: `report_scan` suggests `olla-dft eos --scale X --span 0.04`; warns if the minimum falls outside the range.
4. `phonons`: `report_phonon`; exit code 1 if there are imaginary modes.

**Where each datum comes from.**

| Datum | Source | Detail |
|---|---|---|
| $E$, $F$, $\sigma$ | MLIP calculator | `mace_mp`, `CHGNetCalculator`, `PESCalculator` |
| eV/Å³ → GPa | 160.21766208 | local constant |
| $\sqrt{\mathrm{eV/(Å^2\,amu)}}$ → cm⁻¹ | 521.4708 | `CONV` |
| Masses | `atoms.get_masses()` (ASE) | amu |

**Limits and pitfalls.**
- *"ESTO NO ES EL RESULTADO FINAL. El modelo está entrenado con datos PBE… no mezcles sus energías con las de QE."* Example from the report: Si, MACE 5.464 Å vs LDA 5.402 Å.
- `phonon_check` diagonalises the **full** dynamical matrix of the supercell: the Γ modes of the primitive cell come out plus those of the q points the supercell folds onto Γ (the docstring states this). It is not a dispersion.
- The $B_0$ of the scan comes from a parabola: *"sirve para saber el orden de magnitud, no para reportarlo."*
- Without `torch`/`mace-torch`: *"para usar 'mace' hace falta instalar 'mace-torch'… Ocupa algo más de 1 GB."*

**References.**
- I. Batatia, D. P. Kovács, G. N. C. Simm, C. Ortner, G. Csányi, *NeurIPS* 35 (2022) — MACE.
- B. Deng et al., *Nat. Mach. Intell.* 5, 1031 (2023) — CHGNet.
- C. Chen, S. P. Ong, *Nat. Comput. Sci.* 2, 718 (2022) — M3GNet.

---

### `olla-dft audit` and `olla-dft db` — Comparability between calculations and local index

**What it answers.** Can the total energies of this set of calculations be subtracted? And `db`: which calculations do I have, with which parameters and what came out?

**Background for non-experts.** Two QE total energies can only be subtracted if they come from the same "recipe": same functional, same pseudopotentials, same cutoffs and same treatment of occupations. Otherwise the difference is a perfectly well-formed number without meaning, and QE does not warn. The audit computes a fingerprint with those parameters and groups: more than one group = not comparable. The k grid is treated separately as a warning, comparing the **density** of k points, which is what is comparable between different cells.

**Formulas.** `audit.kdensity`:
$$\rho_k = \frac{n_1 n_2 n_3}{(2\pi)^3 / V}\quad[\text{points}/\text{Å}^{-3}]$$

Fingerprint (`qeout.QEResult.fingerprint` + `origen`): (origin, functional, {element: UPF}, `ecutwfc`, `ecutrho`, `smearing`, `degauss`, `occupations`, `nspin`).

**Implemented rules.**

| Rule | Where | Effect |
|---|---|---|
| DFT vs MLIP origin enters the fingerprint | `audit.audit` (reads `MLIP_PROCEDENCIA.json` via `mlip.read_provenance`) | different groups: NOT COMPARABLE |
| Functional, pseudos, ecutwfc, ecutrho, smearing, degauss, occupations, nspin | `_campos`/`ETIQUETAS` | the differing ones are listed |
| Unconverged SCF | only `scf/relax/vc-relax/md/vc-md` with `converged=False` | "NO CONVERGIERON — sus energías no sirven" |
| `nscf`/`bands` | by calculation type | "Sin energía utilizable" |
| Disparate k density | $\max\rho_k/\min\rho_k > 2$ | WARNING, not an incompatibility |
| Folder without its own XML but with children | `audit.collect` | the children are audited (a sweep) |

**How Olla-DFT computes it.**
1. `audit.collect(paths)`: for each folder reads the MLIP mark and the XML (`qeout.read_xml`).
2. `audit.audit`: groups by fingerprint, lists differences, unconverged ones and those without energy.
3. `audit.report`; exit code 1 if not comparable. `--index` registers them in `olla-dft.db`.
4. `db folder/…` indexes (`audit.index`, `INSERT OR REPLACE` by absolute path); `db --query "SELECT …"` (SELECT only); `db --formula/--calculation/--gap-min/--gap-max` (`audit.search`); `db --export` (JSON); with no arguments, `audit.summary`.

**Where each datum comes from.** Everything from the `pw.x` XML (`qeout.read_xml`): functional, `pseudo_files`, cutoffs, smearing, occupations, `nspin`, energy (Ha → eV), volume, pressure, maximum force, `homo/lumo` → gap, magnetisation, convergence, SCF steps, `nk` (k points used), `nbnd`, BFGS steps, wall time; plus `MLIP_PROCEDENCIA.json` if present. Columns of the `calculos` table in `audit.ESQUEMA`.

**Limits and pitfalls.**
- The fingerprint includes neither the k grid nor the cell: *"un bulk y una losa necesitan mallas distintas por construcción."*
- It compares the **names** of the UPFs, not their content: two different files with the same name pass.
- `hull` and `thermo.from_runs` rely on this audit and refuse to mix origins.
- `db --query` only accepts `SELECT`; old databases are migrated by adding `nk`, `nbnd`, `n_bfgs` (`_migrar`).

**References.** Quantum ESPRESSO manual (`qes` XML schema); K. Lejaeghere et al., *Science* 351, aad3000 (2016) — why pseudos and cutoffs fix the energy reference.

---

### `olla-dft hull` — Formation energies and convex hull

**What it answers.** Is each phase stable against decomposing into the others, and how much energy per atom is it above the convex hull?

**Background for non-experts.** The formation energy per atom is plotted against composition. The lowest curve that envelops all points from below (convex hull) joins the stable phases; any phase above it gains energy by decomposing into the two (or three) hull phases surrounding it, and that vertical distance is $E_{\mathrm{hull}}$. It is energy at 0 K without entropy: a phase 25 meV/atom above is sometimes synthesised anyway.

**Formulas.** In `qekit/modules/thermo.py`.
$$E_f = \frac{E(\text{compound}) - \sum_i n_i\,\mu_i}{N},\qquad \mu_i = \min_{\text{pure phases of } i}\frac{E}{N}$$
$$E_{\mathrm{hull}} = E_f - E_{\mathrm{hull\ line}}(\mathbf x)$$
- Binary (`_casco`): lower envelope by monotone chain over $x$ and linear interpolation. Ternary or higher: `scipy.spatial.ConvexHull` in $(x_1,\dots,x_{n-1}, E_f)$, keeping facets whose normal points downward in energy (`eq[-2] < 0`), and $E_{\mathrm{hull\ line}}$ is obtained by barycentric coordinates inside the facet (`Delaunay.find_simplex`).

**How Olla-DFT computes it.**
1. `audit.collect` + `audit.audit`; if not comparable, it prints the audit and refuses unless `--force`.
2. `thermo.from_runs`: discards `nscf/bands`, no-energy or unconverged runs; refuses to mix DFT and MLIP; formula with `ase.Atoms`; elemental references = lowest energy per atom of the pure phases (warning if any is missing).
3. `_casco`; `report` with metastability threshold `--threshold` (0.025 eV/atom): ESTABLE / metaestable / inestable / fuera del dominio.
4. `export` (`CASCO_CONVEXO.dat`); `plot` only for binaries.

**Where each datum comes from.** Total energies and symbols from the XML of each folder (`qeout.read_xml`); element order from `--elements` or alphabetical.

**Limits and pitfalls.**
- *"Esto es energía a 0 K, sin punto cero ni entropía."*
- Without elemental references: *"hay que calcular cada elemento puro en su fase estable, con los mismos parámetros."*
- `--force` builds the hull with non-comparable calculations at the user's own risk.
- A pure element with several phases: the lowest is the reference; the others come out with $E_f > 0$.
- The plot is only for binaries.

**References.** S. P. Ong, L. Wang, B. Kang, G. Ceder, *Chem. Mater.* 20, 1798 (2008); W. Sun et al., *Sci. Adv.* 2, e1600225 (2016) — metastability scale.

---

### `olla-dft doctor` — pw.x convergence diagnostics

**What it answers.** Is this calculation usable and, if the SCF did not converge, is it because of charge sloshing (mix less) or slowness (mix more or more steps)?

**Background for non-experts.** The self-consistent cycle mixes the new density with the old one. If it mixes too much, the charge "sloshes" from one side of the cell to the other (oscillation, typical in slabs and metals) and the error goes up and down; if it mixes too little, the error always goes down but slowly. The two remedies are opposite, so the module looks at the **shape** of the `estimated scf accuracy` curve.

**Implemented rules** (`diagnose._clasificar`, only if not converged):

| Condition | Diagnosis | Advice |
|---|---|---|
| < 8 iterations | `pocos_datos` | raise `electron_maxstep` to ≥ 100 |
| (≥ 6 points and > 25 % of rises after the first 2 iterations) **or** the error grows > 5× in one iteration | `oscilacion` | `mixing_beta = max(0.05, β/3)`, `mixing_mode='local-TF'`, `mixing_ndim=12` |
| dropped < 3 orders of magnitude in total | `estancada` | check `starting_magnetization`, smearing, distances |
| otherwise, with β ≥ 0.6 | `lenta` | `electron_maxstep = 300` (do not raise β) |
| otherwise, with β < 0.6 | `lenta` | `mixing_beta = min(0.7, max(1.75β, 0.3))`, `electron_maxstep = 300` |

Problems from the XML (`diagnose.diagnose`): unconverged SCF; residual force > 0.05 eV/Å; $|P| > 1$ GPa in `scf/relax/vc-relax`; `Error in routine`. Relaxation: warning if the energy rose in more than $N/3$ steps.

**How Olla-DFT computes it.**
1. `qeout.find_xml` + `read_xml` (convergence, steps, error, forces, pressure, magnetisation, timings).
2. `diagnose.find_stdout` looks for the file containing `Program PWSCF`; `read_scf_history` splits the stdout into SCF cycles with `_ciclos_scf` (each `iteration #  1` opens one; in a `relax` there is one per ionic step), stores `n_ciclos` and extracts **from the last cycle only** `estimated scf accuracy`, `total energy`, `beta`, `convergence has been achieved` / `convergence NOT achieved`; `read_trajectory` reads the `!    total energy`, `Total force`, `P=` lines from the whole file.
3. `report` and `plot` (SCF accuracy on a log scale and energy per ionic step). Exit code 1 if there are problems. `--system` delegates to `health.check` (installation).

**Where each datum comes from.** XML (`converged`, `n_scf_steps`, `scf_error`, `max_force` in eV/Å, `pressure` in GPa, `wall_time`) and `pw.x` stdout (regex `_RE_ACC`, `_RE_ETOT`, `_RE_ITER`, `_RE_FORCE`, `_RE_PRESS`, `_RE_WARN`, `_RE_MAXSTEP`). $\beta$ defaults to 0.4 if not found.

**Limits and pitfalls.** In a `relax` only the last SCF cycle is diagnosed (the report says so: *"en el último de N ciclos SCF (uno por paso iónico; se diagnostica solo el último)"*); an intermediate cycle that oscillated is not seen. The thresholds (0.05 eV/Å, 1 GPa) are fixed. It does not detect symmetry or pseudopotential problems.

**References.** D. D. Johnson, *Phys. Rev. B* 38, 12807 (1988) — Broyden mixing; G. Kresse, J. Furthmüller, *Phys. Rev. B* 54, 11169 (1996) — charge sloshing and `local-TF`.

---

### `olla-dft crosscheck` — The same quantity by two independent routes

**What it answers.** Do two physically independent routes to the same quantity agree? If not, something is wrong in one of them.

**Background for non-experts.** Comparing against the literature detects errors in one module, but not a shared systematic bias. Computing $B_0$ from the equation of state and from the elastic constants, or the band gap and the Tauc gap, are routes that share no code: if they agree, it is hard for both to be wrong in the same way.

**Implemented checks** (`crosscheck.run`; relative deviation $|b-a|/|a|$, or absolute if $a = 0$):

| # | Quantity | Route A | Route B | Tolerance | Data |
|---|---|---|---|---|---|
| 1 | $B_0$ | `EOS.txt` (line with `B0` and `GPa`) | $B_{\mathrm{Hill}}$ from `ELASTIC_C.dat` (`elastic.moduli`) | 5 % | both files |
| 2 | $v_L[100]$, $v_T[100]$ | $\sqrt{C_{11}/\rho}$, $\sqrt{C_{44}/\rho}$ (`derived.cubic_directional`) | LA/TA slope at Γ from `FONONES_BANDAS.dat` | 10 % | Cij, bands, masses, volume |
| 3 | $\Theta_D$ | sound velocities (`derived.debye_from_velocity`) | second moment of `FONONES_DOS.dat` | 30 % (different definitions) | Cij, DOS, N |
| 4 | optical gap | `--gap-bandas` | `--gap-tauc` | 6 % | parameters |
| 5 | $C_v$ at 1500 K | $3Nk_B$ (Dulong–Petit) | $k_B\int x^2 e^x/(e^x-1)^2\,g(\omega)\,d\omega$ with $g$ normalised to $3N$ (`_cv_alta_T`) | 3 % | DOS |
| 6 | number of modes | $3N$ | $\int g(\omega)\,d\omega$ | 5 % | DOS |
| 7 | $\kappa_L$ | `KAPPA.dat` at ~300 K | Slack model from Cij (`derived.slack`) | 60 % | KAPPA, Cij |
| 8 | Berry phase | `BERRY.dat` (column 3 at charge 0) | $-2\sum_n (\bar r_n\cdot b)/2\pi$ from `WANNIER_centros.dat`, same branch mod 2 | 0.05 | both, cell |
| 9 | work function | `ESM.dat` (Φ at $q = 0$) | `WF.dat` (`Phi_eV`) | 5 % | both |
| 10 | $B_0$ (third route) | `EOS.txt` | $-\tfrac{1}{3}\,dP/d\epsilon$ from `STRAIN.dat` (kbar → GPa × 0.1) | 10 % | STRAIN (hydrostatic) |

Constants: `KB_EV` = $8.617333262\times10^{-5}$ eV/K; cm⁻¹ → eV: $1.239841984\times10^{-4}$.

**How Olla-DFT computes it.** `crosscheck._cargar` searches recursively for the result files in the project folder; with `-f structure` it takes masses, volume, N and cell; `run` executes every check for which data exist; `report` marks OK/FALLA with the diagnosis of what to look at first. Exit code 1 if any fails.

**Limits and pitfalls.** *"Un cruce que falla NO dice cuál de los dos caminos está mal."* Check 3 compares different definitions of $\Theta_D$ (*"coincidir al 1 % sería sospechoso"*); check 10 is only valid if the sweep was hydrostatic; check 2 blames the q grid before the Cij. Checks 8–10 swallow any exception silently (`except Exception: pass`).

**References.** R. Hill, *Proc. Phys. Soc. A* 65, 349 (1952); G. A. Slack, *Solid State Phys.* 34, 1 (1979); R. D. King-Smith, D. Vanderbilt, *Phys. Rev. B* 47, 1651 (1993).

---

### `olla-dft selftest` — Validation against known physics

**What it answers.** Does Olla-DFT reproduce measured, published or exact values, and not just what it says about itself?

**Background for non-experts.** Unit tests compare the code with itself. Here each test computes a quantity with a known answer (Ewald constants, Sackur–Tetrode entropy, $T_c$ of aluminium, topological invariants…) and checks it against that reference and its source. `--quick` (default) runs those that do not need `pw.x`; `--full` adds those that do; `--mlip` the one that needs MACE.

**Tests and references** (`selftest.PRUEBAS`; relative deviation, or absolute if the reference is 0):

| Key | Quantity | Reference | Tol. | Source (as stated in the code) | Function tested |
|---|---|---|---|---|---|
| `madelung` | $\alpha_M$ simple cubic | 2.8372974 | 1e-5 | classical Ewald value | `defects.constante_madelung` |
| `lorenz` | $L/L_0$ of the free-electron gas | 1.0 | 12 % | Sommerfeld limit | `transport.compute`, `lorenz` |
| `npw` | plane waves of Si at 30 Ry | 725 | 6 % | what `pw.x` reports (V = 39.5 Å³) | `cost.n_ondas_planas` |
| `sackur` | $S_{\mathrm{trans}}$ of N₂ at 298 K | 150.4 J/(mol·K) | 1 % | Sackur–Tetrode, NIST-JANAF | `thermochem.S_traslacional` |
| `allen_dynes` | $T_c$ of Al (λ=0.44, ω_log=270 K, µ*=0.12) | 1.18 K | 12 % | Allen–Dynes 1975, exp. | `elph.allen_dynes` |
| `allen_dynes_mu` | $T_c(0.10)/T_c(0.12)$ | 1.56 | 5 % | exponential dependence on µ* | `elph.allen_dynes` |
| `born2d` | $Y_{2D}$ with C11=352, C12=60 N/m | 341.8 N/m | 1 % | $Y = C_{11} - C_{12}^2/C_{11}$ (graphene DFT) | `elastic.modulos_2d` |
| `gap_invariante` | ΔE_v of a material with itself | 0 eV | 1e-9 | exact identity | `align.alinear` |
| `ewald_escala` | $\lvert\alpha(3) - \alpha(30)\rvert$ | 0 | 1e-6 | scale invariance | `defects.constante_madelung` |
| `chern_qwz` | $C$ of the Qi–Wu–Zhang model (m=−1) | −1 | 1e-10 | PRB 74, 085308 (2006) | `topology.invariants_from_vectors` |
| `umklapp` (`--mlip`) | exponent $n$ in $\kappa\propto T^{-n}$ of Si | 1.0 | 25 % | Umklapp law above $\Theta_D$ | `kappa.*` with MACE |
| `her_pt` | $\Delta G_{\mathrm{H^*}}$ with $E_{\mathrm{ads}} = -0.33$ eV | −0.09 eV | 5 % | Nørskov 2005, Pt(111) | `echem.her` |
| `oer_ruo2` | η with ΔG(OH,O,OOH)=(0.77, 2.16, 3.87) | 0.48 V | 10 % | Man et al. 2011 | `echem.oer` |
| `escala_oer` | ΔG(OOH) − ΔG(OH) of the RuO₂ profile | 3.2 eV | 10 % | universal scaling relation | `echem.oer` + `echem.escala_ooh_oh` |
| `escala_eta_min` | $\eta_{\min}$ = Δ/2 − ΔG_total/4 | 0.37 V | 2 % | Man et al. 2011 | `echem.sobrepotencial_minimo_escala` |
| `fonon_si` (`--full`) | optical ω(Γ) of Si | 520 cm⁻¹ | 10 % | exp. Raman 520.7 cm⁻¹ | `phonons.*` with `ph.x` |
| `wannier_si` (`--full`) | Si–Si Wannier centre | 1.17563 Å | 2 % | $\sqrt3\,a/8$ with a = 5.43 Å | `wannier.*` |
| `condensador` (`--full`) | slope of $1/C$ vs $d$ / $(1/\varepsilon_0)$ | 1.0 | 6 % | parallel-plate capacitor electrostatics | `esm.*` bc3 Al(111) |
| `born_si` (`--full`) | $Z^*$ of Si | 0 e | 0.05 | acoustic sum rule | `berry.*` |
| `gamma_al` (`--full`) | γ of Al(111) | 1.10 J/m² | 25 % | Vitos 1998 (1.20), exp. 1.14 | `surfen.*` |
| `bulk_si` (`--full`) | $B$ of Si by strain | 95 GPa | 15 % | LDA 93–97 (Nielsen & Martin 1985), exp. 98 | `strain.*` |
| `sitio_h_al` (`--full`) | $E_{\mathrm{ads}}$(top) − $E_{\mathrm{ads}}$(hollow), H/Al(111) | 5.6 eV | 60 % | hollow < bridge < top ordering | `adsorb.*` |

**How Olla-DFT computes it.** `selftest.ejecutar` filters by `--only`, `--full`, `--mlip`; creates a temporary folder (`--keep` to preserve it); runs each `fn(ctx)` and times it; `report` lists value, reference, deviation, tolerance and source. Exit code 1 if any fails or errors. `--list` prints the table without running anything.

**Limits and pitfalls.** *"Las que salen MAL no siempre son un fallo del código: una tolerancia ajustada, un pseudopotencial distinto o un cutoff bajo también las mueven."* The `--full` tests depend on the pseudos in `--pseudo-dir` and on `pw.x`/`ph.x` working.

**Note on `qekit/modules/uncertainty.py`.** It has no command of its own. It offers `propagate(f, values, sigmas)` — propagation in quadrature with central derivatives, $\sigma_f^2 = \sum_i (\partial f/\partial x_i)^2\sigma_i^2$, relative step $10^{-6}$, independent inputs — and `weighted_mean` — weighted mean with $w_i = 1/\sigma_i^2$ and $\sigma = (\sum w_i)^{-1/2}$. No module in this part calls it; only `validation`/`results` check that declared uncertainties are finite and non-negative.

**References.** P. B. Allen, R. C. Dynes, *Phys. Rev. B* 12, 905 (1975); X.-L. Qi, Y.-S. Wu, S.-C. Zhang, *Phys. Rev. B* 74, 085308 (2006); L. Vitos et al., *Surf. Sci.* 411, 186 (1998); O. H. Nielsen, R. M. Martin, *Phys. Rev. B* 32, 3792 (1985).

---

### `olla-dft suggest` — Parameters from your own history

**What it answers.** Based on the calculations that already converged with these elements, which `ecutwfc`, dual, k density and `electron_maxstep` should be used?

**Background for non-experts.** With a few dozen calculations there is no point in training anything: similar calculations are looked up (share elements, similar size) and what worked for them is examined, always stating how many cases back each number.

**Implemented rules** (`recommend.similares`, `recommend.sugerir`):

| Rule | Detail |
|---|---|
| Similarity | only calculations with `convergido`; score = Jaccard of elements $\lvert A\cap B\rvert/\lvert A\cup B\rvert$; × 0.5 if $N_{\mathrm{at}}$ differs by more than a factor of 2 |
| `ecutwfc` | **maximum** among the similar ones (not the mean), with range |
| dual | maximum of `ecutrho/ecutwfc` |
| k density | median of `kdensity` (points/Å⁻³) |
| `electron_maxstep = 300` | if the median of `n_scf` > 40 |
| `mixing_beta = 0.3` + `local-TF` | if the structure is a slab (vacuum along $c$ > 8 Å), general rule, 0 cases |
| Confidence | high ≥ 8 cases, medium ≥ 3, low < 3 |
| No history | refers to the cutoffs of the UPF itself / SSSP |

**How Olla-DFT computes it.** `_cmd_suggest` loads the structure, reads `SELECT * FROM calculos` from `--db` (`olla-dft.db`), detects whether it is a slab and calls `recommend.sugerir`; `report` prints value, number of cases and reason.

**Where each datum comes from.** `calculos` table of `olla-dft.db` (`audit.index`): `formula`, `natoms`, `ecutwfc`, `ecutrho`, `kdensity`, `n_scf`, `convergido`.

**Limits and pitfalls.** *"No sustituyen a una prueba de convergencia: 'olla-dft converge' sigue siendo la forma de saberlo de verdad."* With "low" confidence the report marks *"UN SOLO CASO: tómalo como indicio"* for 1 case and *"SOLO n CASOS"* for 2. It does not invent cutoffs without history.

**References.** G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, N. Marzari, *npj Comput. Mater.* 4, 72 (2018) — SSSP.

---

### `olla-dft pseudos` — Choosing pseudopotentials with criteria

**What it answers.** Of the UPFs available for each element, which are usable for the task (optics, spin-orbit, XANES, DFT+U, phonons) and which one is advisable?

**Background for non-experts.** A folder usually holds several pseudopotentials per element, from different families and functionals. Picking the first one alphabetically fails silently: a scalar-relativistic pseudo with `lspinorb` gives a zero splitting, an ultrasoft one with `epsilon.x` gives a whole wrong spectrum, and mixing functionals between elements invalidates the total energy. The selector applies hard requirements (which discard) and preferences (which rank) and explains every decision.

**Implemented rules** (`pseudos.TAREAS`, `pseudos.evaluar`):

| Task | Hard requirement | Preference |
|---|---|---|
| `optics` | type ∈ {NC} | — |
| `soc` | relativistic = `full` (except elements with Z < 19: note and −0.5 points) | — |
| `xanes` | UPF with `PP_GIPAW` sections | — |
| `hubbard` | — | +0.15 × `z_valence` (semicore) |
| `fonones` | — | +2.0 if type ∈ {NC, US} |
| `general` | — | — |
| all | functional equal to `--functional` (aliases PBE/`SLA PW PBX PBC`, PZ/LDA, PBEsol, BLYP, revPBE) | +max(0, (90 − ecutwfc)/30); −0.5 without a declared cutoff; +1.0 US/PAW with `--cheap`; +0.3 with GIPAW; +0.2 if `full` |

Final order: non-discarded first, then descending points, then name. Coherence across elements (`pseudos.coherencia`): warning if functionals are mixed, if NC is mixed with US/PAW (the ultrasoft dual rules) and if the suggested cutoffs differ by more than 2.5×.

**How Olla-DFT computes it.**
1. `pseudos.candidatos`: `pseudo.find_for_element` (`.UPF` files whose name starts with the symbol) and `pseudos.leer` (type, functional normalised by `_funcional`/`NOMBRE_CORTO`, relativistic, `z_valence`, suggested cutoffs, GIPAW, size).
2. `pseudos.evaluar` and `elegir`; `report` with table and discarded ones; `report_coherencia` if there is more than one element; prints the `--pseudo EL=file` line for reuse.
3. The same selector is used by `sweep.prepare_common` in every command (`pseudo.resolve` → `_elegir` with the task) and `_coherencia_de_funcional` re-selects to unify the functional (preference PBE > PBEsol > revPBE > PZ > BLYP).

**Where each datum comes from.** UPF header (first 20–30 kB): `pseudo_type`, `functional`, `relativistic`, `z_valence`, `wfc_cutoff`/`rho_cutoff` (or their v1 equivalents), presence of `PP_GIPAW`. `Z_SOC` = 19 in `pseudos.py`.

**Limits and pitfalls.** *"Esto es una recomendación, no una verdad… hay que converger el cutoff con 'olla-dft converge'."* Type/functional are inferred by regex from the header: a UPF without those fields shows as `?` and is not discarded. The suggested cutoffs declared by the UPF are a starting point, not a convergence.

**References.** M. J. van Setten et al., *Comput. Phys. Commun.* 226, 39 (2018) — PseudoDojo; A. Dal Corso, *Comput. Mater. Sci.* 95, 337 (2014) — pslibrary; G. Prandini et al., *npj Comput. Mater.* 4, 72 (2018) — SSSP.
