# Avian hypoblast study

- **Aim and contexts**: The pipelines contained in this repository aim at quantitatively analyzing avian hypoblast morphogenesis (1) and genetic patterning (2). 
- **Languages**: Fiji macro and Matlab (2023).
- **Writer**: Aurélien Villedieu (aurelien.villedieu@pasteur.fr)

## (1) Analysis of hypoblast tissue flows and comparison with epiblast tissue flows
### Average-hypoblast-tissue-flows
- **Input**: Timelapse movies of hypoblast morphogenesis (imaged using transgenic quail lines expressing fluorescent reporters) analyzed by Particle Image Velocimetry (PIV) as described in *Saadaoui & al., Science (2020)*.
- **Aim**: Generate archetypal maps of hypoblast tissue flows by averaging PIV data of different animals.

### Compare-hypoblast-and-epiblast-flows
- **Input**: Timelapse movies of both epiblast morphogenesis (memGFP signal, analyzed by PIV) and hypoblast dynamics (reported by grafted hypoblast cells expressing tdTomato-myosin, analyzed by manual tracking).
- **Aim**: Generate archetypal maps of epiblast and hypoblast flows and compare them.


## (2) Analysis of hypoblast patterning
### Staining-signal-projection
- **Input**: Hyperstacks(x,y,z) images of staining of fixed samples (immunostaining and/or HCR-RNA-FISH) on hypoblast and epiblast sides.
- **Aim**: Project the most superficial signal of each hyperstack to specifically extract epiblast (resp. hypoblast) signals.
- **Output**: Projected signals(x,y)

### Staining-signal-quantification
- **Input**: Projected signals(x,y) for *NODAL* HCR-RNA-FISH on the epiblast and hypoblast sides, for many embryos of different timings.
- **Aim**: Generate archetypal maps of *NODAL* mRNA localization in the epiblast and hypoblast at different timings, by averaging signals between animals of the same timing.
