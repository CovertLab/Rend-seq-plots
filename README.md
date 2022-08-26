# Rend-seq

This repository contains the data and code that is necessary to
visualize the end-enriched RNA sequencing (Rend-seq) data reported
by [Lalanne et al. (2018)](https://www.sciencedirect.com/science/article/pii/S0092867418302873). The plots were used as supporting figures
for Sun et al., "Cross-evaluation of E. coli’s operon structures
via a whole-cell model suggests alternative cellular benefits
for low- versus high-expressing operons" (2022), manuscript in preparation.

You can reach out to [Gwanggyu Sun](ggsun@stanford.edu) for any questions
about the code or data in this repository.

This repository was developed with Python 3.8.7. After installing Python, you can install
all of the required packages by running

```shell
pip install -r requirements.txt
```

To regenerate all Rend-seq figures used in Sun et al., run

```shell
bash runscripts/paper_figures.sh
```
