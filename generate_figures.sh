#! /bin/bash
# Figure 2E (rplKAJL-rpoBC operon)
python plot_read_counts.py --range 4176370 4187660 --ymax 5 --filename 2E.pdf

# Figure 3H (oppABCDF operon)
python plot_read_counts.py --range 1299016 1304842 --ymax 4 --filename 3H.pdf

# Figure 3I (potABCD operon)
python plot_read_counts.py --range -1180960 -1184910 --rev --wig 0 --ymax 4 --filename 3I.pdf

# Figure S2B (potABCD operon, treated with monophosphate-specific 5'-exonuclease)
python plot_read_counts.py --range -1180960 -1184910 --rev --wig 1 --ymax 4 --filename S2B.pdf
