# Bandgap Model
The "bandgap.c" source file defines an energy model for evaluating the total energy of charge carriers distributed within a composition-dependent band-edge spectrum.  

The most-thorough description of the model is given here: https://arxiv.org/abs/2307.06268

"mixedhalide.c" implements the bandgap model for mixing of iodine and bromine on the perovskite X-site sublattice.  

The bandgap.c library describes a photoluminescence experiment where the number of photocarriers within the supercell is changed twice. The "external conditions" functionality of lattiKEM is used to expose this change. Initially, the number of photocarriers within the supercell is zero. For steps greater than or equal to *biasTurnOn* but less than *biasSwitch* the number of photocarriers is *numExcitations*. For steps greater than or equal to *biasSwitch* the number of photocarriers is *numExcitationsSecondStep*. 

For each configuration of the lattice, the bg_energy function calculates a discrete band edge energy spectrum. First, the local compositions around each formula unit of the lattice are computed using a cluster approximation. The clusters are identified by expanding outward from a central atom by *coordinationNumber* shells. For instance, coordinationNumber=1 specifies nearest neighbors, coordinationNumber=2 specifies nearest neighbors and next-nearest neighbors, etc. The composition of each cluster is then fed to a "bandgapFunction". mixedhalide.c defines this "bandgapFunction" using a Vegard's law relationship between the bandgap and the iodine and bromine content of the cluster.

Once the band edge excitation spectrum has been generated, the states are occupied using either Maxwell-Boltzmann or Fermi-Dirac statistics based on the value of *thermalDistribution*
