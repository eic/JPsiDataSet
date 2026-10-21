# Dataset details and validation

Exclusive electroproduction of J/psi off the proton (e p -> e' J/psi p', J/psi -> e+ e- with radiative decay) at 9x130 and 9x275 GeV. The datasets are part of Exclusive_Diffractive_JPSI_ABCONV_9x130 and Exclusive_Diffractive_JPSI_ABCONV_9x275. The running instructions are in the [top-level README](../README.md).

## Datasets

All samples use run number (random seed) 1, 1 < Q^2 < 10 GeV^2 and 0.001 < y < 1.

Original dataset:

| Beam | Afterburner preset | Luminosity | Events | Total accepted cross section [nb] | Partial accepted cross section with BR [nb] | Acceptance [%] |
|---|---|---|---|---|---|---|
| 9x130 | ip6_ep_130x9 | 1 fb^-1 | 19,034 | 3.187455e-01 | 1.903229e-02 | 99.994747 |
| 9x275 | ip6_hiacc_275x9 | 2.5 fb^-1 | 68,977 | 4.620712e-01 | 2.759027e-02 | 99.998550 |

Location: root://dtn-eic.jlab.org//work/eic3/users/Olaiya/9x130_GeV/ and root://dtn-eic.jlab.org//work/eic3/users/Olaiya/9x275_GeV/

New dataset (Sept 2026):

| Beam | Afterburner preset | Luminosity | Events | Total accepted cross section [nb] |
|---|---|---|---|---|
| 9x130 | ip6_ep_130x9 | 52 fb^-1 | 982,976 | 3.165866e-01 |
| 9x275 | ip6_hiacc_275x9 | 35 fb^-1 | 961,584 | 4.601208e-01 |

Location: root://dtn-eic.jlab.org//work/eic3/users/Olaiya/new_9x130_GeV/ and root://dtn-eic.jlab.org//work/eic3/users/Olaiya/new_9x275_GeV/

The new and original datasets are not statistically independent, see [Overlap between datasets](#overlap-between-datasets).

Generator level plots: c1.pdf c2.pdf c4.pdf c5.pdf c6.pdf (not included in this repository).

## Configuration files and lAger versions

jpsi-9on130.disp-jpsi-00-electron.json and jpsi-9on275.disp-jpsi-00-electron.json are the lAger configurations of the new dataset. The original dataset used the same files with "lumi" set to 1 and 2.5; these versions are in the git history (git log -p jpsi-9on130.disp-jpsi-00-electron.json).

The configurations for the other beam energies (5x41, 5x100, 10x100, 10x130, 10x250, 18x275) were made from the eic examples in https://eicweb.phy.anl.gov/monte_carlo/lager/-/tree/master/examples/eic for the earlier samples generated with lAger 3.6.1, and are tagged lAger3.6.1-1.0.

Both the original and the new 9x130 and 9x275 datasets were generated with lAger 3.7.0, built locally with modifications.

The new 9x130 and 9x275 datasets were processed with afterburner (abconv) 0.2.1.

Afterburner presets: the presets for 9x130 and 9x275 match the beam parameters stored in the new dataset files, and ip6_ep_130x10 matches the beam parameters of an existing 10x130 sample. The other presets in the top-level README follow the original instruction to use high acceptance where the afterburner provides it; they were not checked against existing samples.

## Validation and caveats

Checked in Oct 2026 against lAger 3.7.0 built in eic_xl:26.10-stable (ROOT 6.40.04, PHOTOS 3.64), using the same configurations and seed.

lAger 3.7.0 fix for ROOT 6.40: the lAger installed in eic-shell 26.09 and later crashes at startup, because ROOT 6.40 no longer accepts a negative PDG code with stable=true. The fix is to register the "Unknown" particle in src/lager/core/pdg.cc with stable=false instead of true; this does not change the generated events.

Production level: reproduced exactly. For the original 9x130 sample the number of events, the cross section and Q2, y, W, t and the J/psi, scattered electron and recoil proton momenta are identical event by event. For the first 20,000 events of the new 9x130 and 9x275 samples, the Lorentz-invariant quantities Q2, y, W, t and M(J/psi) are identical event by event after the afterburner. An independent seed (run 2, 9x130) agrees statistically.

Radiative J/psi decay (PHOTOS): not reproduced. The current build radiates slightly more often: events with no FSR photon go from 55.3% to 53.6% and the mean number of FSR photons from 0.596 to 0.615, with a correspondingly larger low-mass tail in M(e+e-). This is seen at both beam energies, and with two seeds at 9x130. The original and new datasets were produced with the same locally modified lAger build and a different PHOTOS than the current eic-shell, so their radiative decay cannot be reproduced exactly with the current software.

### Statistical independence of the radiative decay

lAger does not seed PHOTOS, which uses its own random numbers. Runs with different run numbers therefore do not have independent radiative decays. With the current build, runs 1 and 2 at 9x130 have the same number of FSR photons in 73.5% of events (correlation 0.55), compared with about 41% expected for independent samples (41.5% between the original sample and run 1, which use different PHOTOS builds). The production kinematics of the two runs are independent. Samples from several runs combined for more statistics have correlated FSR.

### Overlap between datasets

lAger produces the same sequence of events for the same seed, so a larger luminosity only adds events at the end. The original 9x130 sample (19,034 events) is identical, including the radiative decay, to the first 19,034 events of the new 9x130 sample. The original 9x275 sample is most likely the first 68,977 events of the new 9x275 sample; this could not be checked because its generator output is not available. The original and new samples must not be combined or treated as independent.

## Validation scripts

Run inside eic-shell.

compare.py - compares generator-level HepMC3 ASCII samples (kinematics, decay leptons, FSR) and writes overlay plots. The first sample is the reference.

python3 compare.py OUTDIR reference=old.hepmc patched=new.hepmc

compare_invariants.py - event-by-event comparison using only Lorentz-invariant quantities, so afterburned events can be compared with generator output.

python3 compare_invariants.py newset=afterburned.hepmc patched=generator.hepmc

fsr_corr.py - checks whether the FSR photon multiplicity is correlated between runs with different seeds.

python3 fsr_corr.py reference.hepmc seed1.hepmc seed2.hepmc

tree2ascii.cc - converts the first N events of a .hepmc3.tree.root file back to HepMC3 ASCII.

g++ -O2 tree2ascii.cc -o tree2ascii $(HepMC3-config --cflags --libs) -lHepMC3rootIO $(root-config --cflags --libs)

./tree2ascii input.hepmc3.tree.root output.hepmc 20000
