# ARC

Additive relaxed clock for BEAST 2: branch rates are Gamma distributed with a variance that grows linearly with branch duration, so the clock stays consistent when a branch is subdivided.

This repository contains the BEAST2 package implementation of the Additive Relaxed Clock (ARC) model described in:

- Additive Uncorrelated Relaxed Clock Models for the Dating of Genomic Epidemiology Phylogenies. Xavier Didelot, Igor Siveroni and Erik Volz.  *Molecular Biology and Evolution*, Volume 38, Issue 1, January 2021, Pages 307–317, https://doi.org/10.1093/molbev/msaa193

 
The example experiments (xml files, log files and statistics) used by the paper can be found in https://github.com/xavierdidelot/ARC-examples

## Install

Requires BEAST 2.7. Install through BEAUti (`File > Manage Packages`), or build with `ant build` and unpack `dist/ARC.v1.0.0.zip` into your BEAST package directory.

`examples/ARCsimulatedEpidemic.xml` is a runnable example. In XML the clock class is `arc.evolution.branchratemodel.ARClockModel`; before BEAST 2.7 it was `beast.evolution.branchratemodel.ARClockModel`.

Licensed under GPL-3, see `COPYING`.
