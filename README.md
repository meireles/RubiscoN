# Linking genes to ecosystems one carboxylation at a time

### Points of contact

- Dudu, jose.meireles@maine.edu
- Ya, yangya@umn.edu
- Ethan, eebutler@umn.edu
- Artur, stefa066@umn.edu

---

## Rationale

TODO

## Approaches

1. **Genetic / biochemical approach**. Use rbcL's DNA -> AA sequence to infer biochemical properties of RUBISCO -- such as the number of binding sites -- that are relevant for global models (see the [Farquhar model](https://pubmed.ncbi.nlm.nih.gov/24306196/)).

2. **Community ecology approach**. Relate N/RUBISCO to solid estimates of plant performance (N use efficiency, growth rates, etc.). N/RUBISCO could help explain local community dynamics and the productivity of a plot. This will benefit from good time series and carefully collected data (not from a big messy database). We're looking at you [nutnet](https://nutnet.org/)!

3. **Functional ecology / Macroevolutionary approach**. I.e. good old pattern finding. Here we would describe the distribution of N/RUBISCO across vascular plants and estimate how it correlates w/ other functional traits after taking the phylogeny into account.


## History of the idea and why it isn't bogus

Ethan stated that the N content in RUBISCO is a key parameter in global vegetation models. Ya said that we have a crap ton of rbcL data from which we can estimate how much N is in RUBISCO.

These quantities may not be the same. That's because **rbcL gives us N per RUBISCO molecule (N/RUBISCO)** while global models may care about the **total amount of N in RUBISCO**, which is equals **N/RUBISCO * number of RUBISCO molecules**.

Regardless, N/RUBISCO can be a cool trait on its own because it gets at the how much N a molecule of RUBISCO costs. The plant with higher N/RUBISCO will always pay a heavier price to make the stuff than one with lower N/RUBISCO.

**N/RUBISCO is like LMA while total N in RUBISCO is like leaf dry weight.**

---

## Keeping Our Sanity

We are using the R package [`targets`](https://github.com/ropensci/targets) to describe and run the project pipeline. Please check out the [`targets` book](https://books.ropensci.org/targets) to learn more.

This doesn't mean that the project is reproducible across platforms. For example, the MAFFT binaries were compiled for macOS, so the MSA step won't work if you're on another platform.


### Structure

* `_targets.R` runs the pipeline

* `bin` is where the binaries of our dependencies live

* `R` holds the R files with function definitions, **not scripts**

* `raw_data` contains data that is either inconvenient or impossible to access using an API (e.g. TRY)

* `data` contains processed data that need to be read/written by a binary, e.g. the AA fasta files touched by MAFFT.

* `RubiscoN.Rproj` use RStudio and double click this guy to open the project

---

## Data


### Phylogeny

We currently using the megatree (v. 1.1) from [Smith & Brown 2018](https://bsapubs.onlinelibrary.wiley.com/doi/full/10.1002/ajb2.1019), which can be downloaded from [github](https://github.com/FePhyFoFum/big_seed_plant_trees)

### rbcL

Gathered every rubisco aa sequence from [UniPROT(https://www.uniprot.org/) using its API.

## Trait data

### Important

** The zipped TRY data ("13403.zip") is too large for github even with git-lfs**
Therefore, file "13403.zip" was gitignored (not tracked by git or pushed to github) and saved to our google drive folder instead.
The `try_raw_data_file` target to fail if you need to run the pipeline from the very beginning (which you shouldn't). Simply copy the "13403.zip" to the `raw_data/TRY/` directory to solve that issue.

### Trait list
Here is a list requested traits and their TRY IDs

1. Leaf Area: 1, 3108, 3109, 3110, 3111, 3112, 3113, 3114
2. Specific Leaf Area: 11, 3115, 3116, 3117, 3085, 3086
3. Leaf lifespan: 12
4. Leaf nitrogen (N) content per leaf area: 50
5. Leaf nitrogen (N) content per leaf dry mass: 14
6. Leaf photosynthesis pathway: 22
7. Leaf photosynthesis carboxylation capacity (Vcmax) per leaf area (Farquhar model): 186
8. Leaf photosynthesis carboxylation capacity (Vcmax) per leaf dry mass (Farquhar model): 185
9. Leaf phosphorus (P) content per leaf area: 51
10. Leaf phosphorus (P) content per leaf dry mass: 15

---
