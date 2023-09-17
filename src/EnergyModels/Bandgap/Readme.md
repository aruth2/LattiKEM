# Bandgap Model
The bandgap module defines an energy model for evaluating the total energy of charge carriers distributed within a composition-dependent band-edge spectrum.  

The most-thorough description of the model is given here: https://pubs.acs.org/doi/10.1021/acs.jpcc.3c04708

A visual depiction of the rearrangement of iodine ions (red) and bromine ions (green) in a mixed halide perovskite lattice under illumination is shown below originally from: https://doi.org/10.1021/acsenergylett.8b01369

https://github.com/aruth2/LattiKEM/assets/8935880/9b2d6dd2-78de-435b-87dc-16d61afe3fb6
In the video, the charge carriers are not directly shown, but the rough location where they accumulated was indicated by highlighting the atoms. The bandgap module is now able to directly output the carrier density around each coordination atom using the **saveChargeDensity** setting. 


The steps of the trajectory are mapped to a photoillumination experiment by the bandgap module. The number of photocarriers within the supercell is changed up to two times. The "external conditions" of the trajectory holds the number of photocarriers at each step of the simulation. Initially, the number of photocarriers within the supercell is zero. For steps greater than or equal to **biasTurnOn** but less than **biasSwitch** the number of photocarriers is **numExcitations**. For steps greater than or equal to **biasSwitch** the number of photocarriers is **numExcitationsSecondStep**.

"mixedhalide.c" implements the bandgap model for mixing of iodine and bromine on the perovskite X-site sublattice. The number of formula units is **sizex** x **sizey** x **sizez**. The lattice can contain or not contain periodic boundary conditions in each direction using the **pbcMask**. A mask of "111" (decimal 7) gives periodic boundary conditions in all 3 directions. A mask of "100" gives periodic boundaries in z but not in y or x directions. The iodine fraction of the halide sublattice is **iodineRatio** and the vacancy fraction is **vacancyRatio**. Hopping barriers can be separately defined for the two halides. For iodine this is **iHopEnergy** and for bromine it is **brHopEnergy**. 

 

For each configuration of the lattice, the bg_energy function calculates a discrete band edge energy spectrum. First, the local compositions around each formula unit of the lattice are computed using a cluster approximation. The clusters are identified by expanding outward from a central atom by **coordinationNumber** shells. For instance, coordinationNumber=1 specifies nearest neighbors, coordinationNumber=2 specifies nearest neighbors and next-nearest neighbors, etc. The composition of each cluster is then fed to a "bandgapFunction". mixedhalide.c defines this "bandgapFunction" using a Vegard's law relationship between the bandgap and the iodine and bromine content of the cluster. The constant offset is **pureIBandgap**. The difference between pure bromine and pure iodine compositions is **bandgapDifference**, and the bowing term is **bowingParameter**. 

Once the band edge excitation spectrum has been generated, the states are occupied using either Maxwell-Boltzmann or Fermi-Dirac statistics based on the value of **thermalDistribution**. Absorption spectra, photoemission spectra, and a band-edge density of states are calculated on a continuous energy spectrum from the discrete band edge energy spectrum. The photoemission spectrum factors in the occupancy of the band edge states, while the absorption and DOS do not. The total energy of the photocarriers is computed by summing the weight and energy of each band edge state.

"bandgap.c" also allows for an interatomic interaction energy between nearest-neighbor pairs. This energy is directly analogous to a "mixing enthalpy". "mixedhalide.c" implements this functionality and exposes it through the **IBrRepulsiveEnergy**.


