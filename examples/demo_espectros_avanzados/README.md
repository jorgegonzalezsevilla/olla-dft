# Spectra and advanced modules: XANES, Hubbard U, electron–phonon, band unfolding and VDOS

REAL results from five modules (`xanes`, `hubbard`, `elph`, `unfold`,
`md`), with what to look at in each one.

### XANES (`xanes_Si.png`, `XANES.dat`)

Silicon K edge. All three polarizations are plotted and cannot be told
apart because they coincide: silicon is cubic, and that they agree to
0.04 % is the check that polarization is handled correctly.

- Edge at +1.1 eV above E_F
- Maximum at +3.7 eV
- Features at +10.4, +12.5 and +17.2 eV

Reproduce with:

    olla-dft corehole Si --edge K -o pseudos --functional PZ --rcut 1.6
    olla-dft xanes Si8.cif --element Si --core-hole pseudos/Si.hueco1s.UPF --average -o xanes --ecutwfc 40
    # run pw.x and then xspectra.x on the three xspectra_*.in
    olla-dft xanes Si8.cif --collect -o xanes --element Si --edge K

### Hubbard U (`HUBBARD_U.dat`)

Hubbard U of NiO by linear response (hp.x), ortho-atomic projection. The
self-consistency cycle gives:

    iter 0:  5.4429 eV      <- what a ONE-shot calculation reports
    iter 1:  3.9429
    iter 2:  4.1323
    iter 3:  4.1087 eV      <- the self-consistent one

A 1.33 eV difference. That is the whole argument for the module.

    olla-dft hubbard NiO.cif -o hub --qgrid 2x2x2          # prepares scf + hp.x
    # run pw.x and hp.x
    olla-dft hubbard NiO.cif --collect -o hub --qgrid 2x2x2 # reads the ONE-shot U
    olla-dft hubbard NiO.cif --cycle -o hub --qgrid 2x2x2   # scf -> hp.x -> scf loop until converged

### Electron-phonon coupling (`ELPH_Al_lambda.dat`)

Electron-phonon coupling of aluminium, 2×2×2 q mesh. The column that
matters is not a single one: it is the series against the broadening, and
one has to read the PLATEAU (here around 0.35). That λ rises from 0.018 to
0.35 and then settles is normal; if it did NOT settle, the k mesh would be
insufficient and any number would be arbitrary.

From it comes τ(300 K) = 11.4 fs, which replaces the constant τ of the CRTA
approximation in the transport module.

    olla-dft elph Al.cif --qgrid 2x2x2 -o elph        # prepares scf, nscf and ph.x
    # run pw.x and ph.x
    olla-dft elph Al.cif --collect -o elph

### Band unfolding (`UNFOLD.dat`)

Unfolding of a 2× silicon supercell WITHOUT a defect. The weights come out
exactly 1.0 or 0.0 (0.5/0.5 at the zone-boundary degeneracies), and the
sum over bands is exactly 4.000 at all 21 k points. That is the theorem:
with N = 2, half the weight survives.

On a supercell WITH a defect the weights get spread out, and that smearing
is the physical result.

    olla-dft unfold carpeta_bandas_supercelda/ prim.cif -o . --bands 8 --format png

### VDOS from molecular dynamics (`MD_VDOS.dat`)

Vibrational density of states of silicon at 900 K from a molecular-dynamics
trajectory. Unlike harmonic phonons, it includes anharmonicity and
temperature. The frequency resolution is set by the trajectory length, and
the module reports it: with 0.24 ps it is 138 cm⁻¹, which is coarse. A
20 ps trajectory would give 1.7 cm⁻¹.

    olla-dft md md.out -o . --skip 50 --no-plot

### Files

| File | What it is |
|---|---|
| `XANES.dat` | σ(E) of the Si K edge: average and the three polarizations |
| `xanes_Si.png` | XANES spectrum figure |
| `HUBBARD_U.dat` | per-site U of NiO as written by `hubbard --collect` (one shot, 5.44 eV) |
| `ELPH_Al_lambda.dat` | λ, ∫α²F, ⟨log ω⟩ and N(E_F) versus broadening for Al |
| `UNFOLD.dat` | unfolded spectral weight: path distance, E−E_F and weight |
| `MD_VDOS.dat` | VDOS of Si at 900 K from the trajectory (250 steps, dt = 0.97 fs) |
