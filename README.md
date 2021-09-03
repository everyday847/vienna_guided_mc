# What's in this repository?
This repository contains scripts to take a ribosomal RNA and alter its sequence in ways that should lower its 'delta': the difference in energy between that sequence's minimum free energy (MFE) structure and the experimentally determined secondary structure.

# Why?
An initial hypothesis raised by players of the [Eterna](https://eternagame.org) video game was that they would be able to improve ribosome function by minimizing this difference, i.e., by getting as close as possible to "solving" the secondary structure design problem for the ribosome. (In simpler puzzles, when the target structure _is_ the MFE structure, the puzzle is solved. This is not always possible, and the experimental structure of the ribosomal RNAs may not admit a perfect solution.)

Unsurprisingly, while lowering this "delta" helped in some cases, one really significant effect was that more mutations would tend to have a deleterious effect. That's not too surprising: more mutations means more opportunities to mess something up. So how well were players doing? We needed a reasonably smart control that would attempt to make mutations in a principled way with the same strategy of "delta minimization."

The scripts in this repository either run Metropolis Criterion Monte Carlo simulations on the 23S and 16S rRNA in a secondary structure aware fashion or analyze the results of those simulations.
