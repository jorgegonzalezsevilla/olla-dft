# Ethylene with TDDFPT and ballistic conductance of an Al wire

Optical excitations of a molecule with `tddft` (turbo_davidson.x) and
complex bands / conduction channels of a wire with `ballistic` (pwcond.x).

### TDDFPT (`tddft_C2H4.png`, `TDDFT_etileno.eigen`)

Ethylene with TDDFPT (turbo_davidson.x, PZ functional). Six excitations:

    n    E (eV)    f (osc.)   pol.
    1    6.4955     0.02893    x     <- the π->π*, the only truly bright one
    2    7.1554     0.00000    z     <- dark
    3    7.1726     0.00011    z
    4    7.2074     0.01237    x
    5    7.4331     0.00005    z
    6    7.4947     0.00001    z

What to look at: four of the six are DARK. They exist as excited states
and do not show up in an absorption spectrum. A calculation that only
reported "the first excited state is at 6.50 eV" would be hiding that the
next four contribute nothing to the spectrum.

The experimental π->π* of ethylene is at 7.66 eV. 6.50 eV is low, and that
is expected: LDA underestimates this transition. The `.eigen` file carries
the energies in RYDBERG; reading them as eV would give 0.48 eV, which is
absurd, and it is an easy mistake to make.

Reproduce with:

    olla-dft corehole C --plain -o ps --functional PZ --rcut 1.3
    olla-dft corehole H --plain -o ps --functional PZ --rcut 1.0
    olla-dft tddft c2h4.cif -o run --method davidson --states 6 --pseudo-dir ps --ecutwfc 40 --emax 12
    # run pw.x and turbo_davidson.x
    olla-dft tddft --collect -o run --method davidson --gap 6.6

### Ballistic conductance (`BALISTICO.dat`)

Complex band structure of a monatomic aluminium wire (pwcond.x, ikind=0).
The "canales" column is the number of open channels at each energy, which
bounds the conductance from above: G ≤ channels × G0.

    below -0.3 eV:          0 channels  (does not conduct)
    from -0.3 to 0.2 eV:    1 channel   (the s band)
    near +0.3 eV:           3 channels  (the degenerate p bands come in)

That jump from 1 to 3, and not to 2, is the signature of a double
degeneracy: the wire's p_x and p_y come in together.

With G0 = 2e²/h = 7.748e-5 S, a single perfectly transmitting channel is
12.906 kΩ. For the actual conductance you need a scattering region in the
middle (`--scatterer`) and ikind=1.

    olla-dft ballistic Al_hilo.cif -o . --ikind 0      # prepares scf + pwcond.x
    # run pw.x and pwcond.x
    olla-dft ballistic --collect -o . --no-plot

### Files

| File | What it is |
|---|---|
| `TDDFT_etileno.eigen` | raw turbo_davidson.x output: energy (Ry) and total and per-component oscillator strengths |
| `tddft_C2H4.png` | absorption spectrum with the six excitations marked |
| `BALISTICO.dat` | open channels versus E−E_F of the Al wire, as written by `ballistic --collect` |
