# Graphite: layers, diffractogram and exfoliation energy

Layered-materials module (`layers`, `xrd`, `exfoliate`), demonstrated on
graphite. The structure is `../grafito.cif`.

    olla-dft layers ../grafito.cif
      -> 2 layers, basal d = 3.356 Å, (002) expected at 26.54° (Cu Kα)

    olla-dft xrd ../grafito.cif --size 18 --exp demo_experimental.xy
      -> simulated diffractogram + comparison with an "experimental" one
         (the "experimental" one here is synthetic, generated only for the
         demo: the same graphite with a 2 % larger basal spacing, noise
         and background)

    olla-dft exfoliate ../grafito.cif --run
      -> E_exf(LDA) = 25.8 meV/atom = 0.157 J/m² (LDA literature value)

The XRD is verified numerically against pymatgen (Si and graphite, every
peak within Δ2θ < 0.05° and ΔI < 1.5).

### Files

| File | What it is |
|---|---|
| `XRD_HKL.dat` | list of graphite reflections: 2θ, d, intensity and (hkl), λ = 1.54184 Å |
| `demo_experimental.xy` | synthetic "experimental" diffractogram (2θ, I) to exercise `--exp` |
| `xrd.pdf`, `xrd.png` | simulated diffractogram against the synthetic one |
