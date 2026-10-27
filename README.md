# Diatomic MO Lab

An interactive molecular orbital calculator for homonuclear and heteronuclear diatomic molecules, built for general chemistry students. Pick two atoms (H through Ne, plus Cl) to see the MO energy-level diagram, bond order, magnetism, and a picture of every orbital.

The whole tool is one self-contained web page, `index.html`. It runs in any modern browser with no installation, and works on phones.

## What students can do

- **Build MO diagrams** for any pair of atoms and charges from −2 to +2, with electrons filled by the Aufbau principle and Hund's rule.
- **Look at orbitals** as a phase-colored slice through the bond axis (nodes marked) or as a rotatable 3D dot cloud.
- **Compare two molecules** side by side on the same energy scale.
- **Follow guided activities** that set up each step and pose a question with a hidden answer:
  1. Bond order across period 2 (extended Hückel, 7 steps)
  2. Oxygen and its ions (extended Hückel, 5 steps)
  3. Polar bonds: CO and HF (Hartree–Fock, 5 steps)
  4. Bond length and energy (Hartree–Fock, 5 steps)

## Two calculation methods

**Extended Hückel** (Hoffmann, 1963) uses Slater-type valence orbitals. The overlap integrals are computed numerically for each bond length, and the Hamiltonian follows the Wolfsberg–Helmholz formula, H<sub>ij</sub> = K·S<sub>ij</sub>·(H<sub>ii</sub> + H<sub>jj</sub>)/2 with K = 1.75. The method is fast and qualitative. Sliders control s–p mixing and K.

**Hartree–Fock** solves the Roothaan equations self-consistently in the STO-3G basis, with core orbitals included:

- Closed-shell molecules use RHF; open-shell molecules use ROHF.
- Integrals use the McMurchie–Davidson scheme.
- The SCF starts from a superposition of atomic densities and is accelerated with DIIS.

Because Hartree–Fock gives total energies, it can optimize the bond length and draw a potential-energy curve.

The calculations run in a background Web Worker, so the page stays responsive.

## Validation

`tests/test_engines.js` checks both engines:

- **Hartree–Fock energies.** Total energies for 17 molecules and ions agree with the independent [PyQuante2](https://github.com/rpmuller/pyquante2) program to better than 10<sup>−6</sup> hartree.
- **Rotation invariance.** Rotating the molecule off the z axis leaves the energy unchanged to 10<sup>−9</sup> hartree. This exercises every p-type integral.
- **Bond lengths.** Optimized lengths reproduce published STO-3G values: H<sub>2</sub> 0.712 Å, N<sub>2</sub> 1.134 Å, CO 1.146 Å, HF 0.956 Å.
- **Extended Hückel.** Bond orders and magnetism for the period-2 series and the O<sub>2</sub> ions.

Known model behavior, flagged on the page when it occurs:

- Extended Hückel puts σ<sub>2p</sub> slightly above π<sub>2p</sub> in O<sub>2</sub>.
- Hartree–Fock/STO-3G puts σ<sub>2p</sub> above π<sub>2p</sub> in F<sub>2</sub>.
- Restricted Hartree–Fock cannot dissociate bonds correctly.

## Repository layout

| Path | Contents |
| --- | --- |
| `index.html` | The built, self-contained page (served by GitHub Pages) |
| `src/engine.js` | Extended Hückel engine |
| `src/hf.js` | Hartree–Fock engine (integrals, RHF/ROHF, optimizer, energy curves) |
| `src/sto3g.json` | STO-3G basis set for H–Ne and Cl |
| `src/app.js` | Interface, diagrams, orbital rendering, compare view, activities |
| `src/template.html` | Page markup and styles |
| `src/build.py` | Assembles `index.html` from the sources |
| `tests/test_engines.js` | Regression tests |

## Building and testing

```sh
python3 src/build.py          # writes index.html (and src/hf.built.js for the tests)
node tests/test_engines.js    # runs the test suite
```

## Credits

Extended Hückel parameters are from the Hoffmann and Alvarez tables. The STO-3G basis is from Hehre, Stewart and Pople (1969), as distributed with PyQuante2. Fonts are IBM Plex, loaded from Google Fonts.
