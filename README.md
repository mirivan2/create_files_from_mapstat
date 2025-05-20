# create_files_from_mapstat

This repository contains two Python scripts designed to process ".mapstat" files from KMA output for ResFinder database and generates normalized depth for downstream analysis.

---

1. "create_files_from_mapstat.py"

1.1. Purpose: 
Processes ".mapstat" files to compute log10-transformed normalized sequence depth for ARG-level data, from ResFinder output.

1.2. Usage:
python3 create_files_from_mapstat.py <input_folder> <refdata_file> <output_prefix>

1.3. Arguments:
"<input_folder>": Directory containing ".mapstat" files.
"<refdata_file>": Reference data file (e.g., "ResFinder_*.refdata.txt") with gene lengths.
"<output_prefix>": Prefix for output files (e.g., "P2_four").

1.4. Outputs:
"pivot_<output_prefix>_log10_norm_seq_dep.tsv"
"long_<output_prefix>_log10_norm_seq_dep.tsv"

1.5. Notes:
The script:
a. Filters out ARGs below a threshold coverage (e.g., < 0.9).
b. Adds missing gene lengths manually for known entries.
c. Merges gene variants by averaging normalized depth.

---

2. "create_files_from_mapstat_genomic.py"

2.1. Purpose:  
Processes genomic ".mapstat" files to compute log10-transformed normalized depth at species level, focusing on bacterial taxa.

2.2. Usage:
python3 create_files_from_mapstat_genomic.py <input_folder> <refdata_file> <output_prefix>

2.3. Arguments:
"<input_folder>": Directory containing ".mapstat" files.
"<refdata_file>": Reference data file with taxonomic annotations.
"<output_prefix>": Prefix for output files (e.g., "P2_eight").

2.4. Outputs:
"pivot_<output_prefix>_log10_nom_seq_dep.tsv"
"long_<output_prefix>_log10_nom_seq_dep.tsv"

2.5. Notes:
The script:
a. Filters out bacterial species below a threshold coverage (e.g., < 0.5).
b. Only includes entries classified as "Bacteria".
c. Requires a separate file "sequence_lengths.tsv" for reference sequence lengths.

Developer: Emilie Egholm Bruun Jensen, National Food Institute, Technical University of Denmark, emegbr@food.dtu.dk
