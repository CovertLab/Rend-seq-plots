#! /bin/bash
# Figure 2E (rplKAJL-rpoBC operon)
python plot_read_counts.py --range 4176370 4187660 --ymax 5 --filename 2E.pdf

# Figure 3I (oppABCDF operon)
python plot_read_counts.py --range 1299016 1304842 --ymax 4 --filename 3I.pdf

# Figure 3J (cmk-rpsA-ihfB operon)
python plot_read_counts.py --range 960380 963395 --ymax 5 --filename 3J.pdf
