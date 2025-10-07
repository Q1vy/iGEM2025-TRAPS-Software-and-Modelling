# iGEM2025-TRAPS-Software
Repository for Software and Modelling done for the TU-Dresden team for iGEM2025 competition.
General topic: RNA-binding Proteins 
- predicting binding affinities
- sequence comparison for target selection



### Content description:
#### Main target search
 given a target sequence (of an RNA to bind), can find subsequences of a given length (binding sites) that are as unique compared to reference sequences (transcriptom).
- New: 
  - weighted analysis: each transcript in the transcriptome can have its own statistical weight e.g. concentration
  - combinatorial optimization for tetramers: does not only find the best pumbys, but can optimize a selection of 4 pumbys.
  - Total runtime < 5min 
- Next:
  - validate and fix unweighted workflow (mostly focused on weighted recently)
  - combinatorial optimization for k-mers
  - secondary RNA structure -> each potential binding site gets its own weight
  - convert from current functional style to a more object oriented use experience
  - reintroduce restriction enzyme avoidence
  - switch away from csv for sequence comparison result storage: array representation is a bad
Related files
 - target search test: demonstration and development notebook for new algorithms and functions on random sequence data. Finished features get automatically written into main_search.py to be later imported into target search main.
 - target search alt: quick evaluatoin of our first selection of binding sites based on unweighted transcriptome data -> seems to hold up quiet well.
 - sequence reader: demonstration and development notebook for reading, decoding and preprocessing sequence data. Finished features get automatically written into sequence_reader.py to be later imported into target search main.

#### Modeling
simple Boltzman statistics, to predict the binding interactions of RNA binding proteins with their specific target, and off-target interactions with the transcriptome. 
- Next: 
  - include high valency binding dynamics -> increased affinity due to spatial correlation of binding partners
    - either model probabilistic, stepwise binding dynamics
    - or as effective on bond with the energy of the sum of all the high valency interactions 
  - include aspects of condensates or polymer physics

### Science communication / teaching
- Todo: simpler and shorter notebooks explaining aspects of the algorithms and datastructures used here

### Git-orga
- diffs of notebooks are messy, look for alternative

