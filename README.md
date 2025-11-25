#  SEA-STAR - Team TU-Dresden 2025 Software Tool
side effect aware target site ranking

> Additionally: the code for the image analysis we used for our experimental result of our project can be found in the image-analysis branch.

## Abstract
Our software compares sequence data of a target RNA and a (quantitative) transcriptome to evaluate which binding sites on the target (= queries) have the lowest probability for undesirable binding interactions with the transcriptome (= reference library). It reads transcriptome sequences from FASTA and can add expression atlas data (by EMBL-EBI) to use quantitative transcriptome dataset for binding affinity estimations. It is designed with sequence based probes of any length and adjustable binding strength (e.g. sgRNA of Cas13 or Pumilio RBPs) in mind, allowing users to rank the specificity of single queries, or query groups (when using multiple probes simultaneously to minimize overlapping off-target effects). We also show how to use secondary structure based binding site availability scores from vienna RNA within out tool.

## How does it work
Have a look at our [team website](https://2025.igem.wiki/tu-dresden/software) and in the 'statmech_binding_model.ipynb' for more details.


## Installation
 - download / clone this repository
 - set up python dependencies using uv, conda or pip
 - Have fun!

 >For new versions and contributions please look [here](https://github.com/Q1vy/iGEM2025-TRAPS-Software-and-Modelling)

## Getting started
> Most core utilities are defined in main_search.py and sequence_reader.py. <br>
> The jupyter notebooks contain implementations on how to use them.

Implementation notebooks by simplicity
- target_search_mini: using random sequence data to show the simplest way of ranking single binding site of a target RNA against a transcriptome without expression levels
- target_search_mCherry: using quantitative sequence data AND [Vienna RNA](http://rna.tbi.univie.ac.at/) guide scores to show a simple way of ranking single binding site of a target RNA against a transcriptome dataaset containing expression levels
- target_search_dev has definitions, explanations and usage (on random data) of all functions in full detail, also showing an implementation of evaluating of groups of queries simultaneously to avoid overlapping off-target effects, which is not yet part of main.py
- Pumby and Cas13 are (older) notebooks we used to choose our sequences for our project. Both include an evaluation of query groups, and CAS13 also uses ViennaRNA guidescores.
- eval shows how to estimate tbe affinity of a single (or a few) binding sites, not searching through all
## Authors and acknowledgment
Seastar: [Paul Stark](https://github.com/Q1vy)

## Recent
 - new V2 in progress, goal: instead of a fixed k=[0, 1, 2], make it dynamic and implement group evaluation into .py library file
