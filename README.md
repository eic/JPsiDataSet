# Exclusive diffractive J/psi datasets (lAger)

Running instructions for producing the exclusive J/psi samples. Run all steps inside eic-shell. Dataset details, validation and caveats are in [validation/README.md](validation/README.md).

## Beam energies

| Beam | Configuration | lAger version | Afterburner preset |
|---|---|---|---|
| 5x41 | jpsi-5on41.disp-jpsi-00-electron.json | 3.6.1 | ip6_hiacc_41x5 |
| 5x100 | jpsi-5on100.disp-jpsi-00-electron.json | 3.6.1 | ip6_hiacc_100x5 |
| 10x100 | jpsi-10on100.disp-jpsi-00-electron.json | 3.6.1 | ip6_hiacc_100x10 |
| 10x130 | jpsi-10on130.disp-jpsi-00-electron.json | 3.6.1 | ip6_ep_130x10 |
| 10x250 | jpsi-10on250.disp-jpsi-00-electron.json | 3.6.1 | ip6_ep_250x10 |
| 18x275 | jpsi-18on275.disp-jpsi-00-electron.json | 3.6.1 | ip6_hiacc_275x18 |
| 9x130 | jpsi-9on130.disp-jpsi-00-electron.json | 3.7.0 | ip6_ep_130x9 |
| 9x275 | jpsi-9on275.disp-jpsi-00-electron.json | 3.7.0 | ip6_hiacc_275x9 |

The high-acceptance (hiacc) preset is used where the afterburner provides one; 130 and 250 GeV protons only have the ip6_ep preset.

## lAger versions

The samples for the first six beam energies were generated with lAger 3.6.1 (configurations tagged lAger3.6.1-1.0); the 9x130 and 9x275 samples with lAger 3.7.0.

3.7.0 adds a CLAS12 phi electroproduction process and compiler fixes. The lAger code for the J/psi model (oleksii_2vmp), the decay and the call to PHOTOS for the radiative decay is unchanged, so the same configurations work with both versions.

PHOTOS caveat: the radiative J/psi decay depends on the PHOTOS library lAger is built against. PHOTOS uses its own random numbers, which the lAger run number does not set, so runs with different -r do not have statistically independent radiative decays. With the current eic-shell (PHOTOS 3.64) the 9x130 and 9x275 samples are reproduced exactly except for the radiative decay, which radiates slightly more often (events without an FSR photon go from 55.3% to 53.6%). See [validation/README.md](validation/README.md).

eic-shell 26.09 and later ships lAger 3.7.0, which crashes at startup with ROOT 6.40. To run it, rebuild lAger with the "Unknown" particle in src/lager/core/pdg.cc registered with stable=false instead of true. lAger 3.6.1 is not part of eic-shell and contains the same line. See [validation/README.md](validation/README.md) for details.

## 1. Generate with lAger

Get lAger from https://eicweb.phy.anl.gov/monte_carlo/lager and run it with the configuration for the beam:

./bin/lager -c jpsi-9on130.disp-jpsi-00-electron.json -r 1 -o output_directory

The number of events is set by "lumi" (in fb^-1) in the configuration, and -r sets the run number, which is also the random seed.

## 2. Apply beam effects with the afterburner

Use abconv from https://github.com/eic/afterburner with the preset for the beam:

abconv -p "ip6_ep_130x9" generator_output.hepmc

## 3. Convert to hepmc3.tree.root

Build https://github.com/eic/hepmc3ascii2root and convert the afterburner output:

cmake -S . -B build && cmake --build build

./build/hepmc3ascii2root input_file output_file
