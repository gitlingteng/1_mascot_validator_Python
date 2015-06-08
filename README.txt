Our approach is to mine the Mascot .DAT file to extract information potentially useful for automated validation.

The raw data are separated into two distinct zones: the high Mascot score peptides, most with low precursor mass error, and the low Mascot score peptides, most with high precursor
mass error. As the Mascot score increases from 0 to 35, the variance of the precursor mass errors of all peptide matches above this score falls dramatically (see also supplemental Fig. 1). We determined cutoffs for precursor mass error that would encompass 95% of all peptides (dashed lines) and 95% of peptides with Mascot peptide scores over 35 (solid lines).

The program has 2 parts. First, mascotparer is a small package to parse the mascot .DAT file,remove the useless info, and save the 
important one such as peptide, mass charge, protein info, peptide intensity as dictionary object for further input in the 2nd part Validator 1.
 
Then Validator 1 sought to find all16O/18O pairs in the Mascot summary file (“.DAT file”). The program iterates through all queries looking for identical top scoring peptides found in both 16O and 18O forms (a “16O/18O pair”) of the peptides with low Mascot peptide scores. The output will be exported to .csv file and .xls file.
