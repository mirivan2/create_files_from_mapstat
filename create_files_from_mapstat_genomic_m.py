import sys
import os
import pandas as pd
import numpy as np

## Usage:
# python3 create_files_from_mapstat_genomic.py /home/projects/cge/data/projects/2022/7000/KMA/output_triplicates_P2/genomic_20220524_eight/data_ind10/ /home/projects/cge/people/emijen/tools/genomic_20220524.master.refdata P2_eight_ind10


# List to store dataframes
dataframes = []

# Load refdata
refdata = pd.read_csv(sys.argv[2], sep='\t', header=0)

# Load seq lengths
seq_len_df = pd.read_csv("/home/projects/cge/people/emijen/tools/sequence_lengths.tsv",sep='\t',names=["id","id_len"])

# Iterate through all files in the folder
for filename in os.listdir(sys.argv[1]):
    
    file_path = os.path.join(sys.argv[1], filename)
    
    # Read the file skipping the first 6 lines and using the first row as column headers
    df = pd.read_csv(file_path, sep='\t', skiprows=6, header=0)
    df = df[["# refSequence","refCoveredPositions","bpTotal"]]
    df["# refSequence"] = df["# refSequence"].str.split(" ").str[0] 

    # Merge df with refdata using the column '# refSequence' from df and the column 'id' from refdata
    merged_df = pd.merge(df, refdata[["id","superkingdom_name","species_name"]], left_on='# refSequence', right_on='id').drop("id",axis=1)

    # Only include bacteria
    merged_df = merged_df[merged_df["superkingdom_name"] == "Bacteria"].drop("superkingdom_name",axis=1)

    # Add length data
    merged_df = pd.merge(merged_df, seq_len_df[["id", "id_len"]], left_on='# refSequence', right_on='id').drop("id",axis=1)

    # Print intermediate file
    merged_df.to_csv("int_1_"+filename, sep="\t", index=False)

    # Calculate coverage and remove everything below 0.1
    merged_df["coverage"] = merged_df["refCoveredPositions"] / merged_df["id_len"]
    removed = len(merged_df) - len(merged_df[merged_df["coverage"] >= 0.1])
    print("Removed {0} taxa from file {1}, because coverage < 0.1".format(removed,filename))
    merged_df = merged_df[merged_df["coverage"] > 0.1]

    # Calculate normalised sequence depth
    merged_df["nom_seq_dep"] = merged_df["bpTotal"] / merged_df["id_len"].sum()

    # Merge to species level
    merged_df = merged_df.groupby('species_name')['nom_seq_dep'].mean().reset_index()
    print(merged_df)

    # Log10 norm seq depth
    # merged_df["log10_nom_seq_dep"] = np.log10(merged_df["nom_seq_dep"])
    # merged_df = merged_df.drop("nom_seq_dep",axis=1)

    # Add sample name
    merged_df["sample"] = filename

    # Print intermediate file
    merged_df.to_csv("int_2_"+filename, sep="\t", index=False)

    # Append the dataframe to the list
    dataframes.append(merged_df)


# Make final dataframes
master_df = pd.concat(dataframes)

pivoted_df = master_df.pivot(index='sample', columns='species_name', values='nom_seq_dep')
pivoted_df.to_csv("pivot_"+sys.argv[3]+"_nom_seq_dep.tsv", sep="\t")

long_df = master_df.pivot(index='species_name', columns='sample', values='nom_seq_dep').reset_index()
long_df.to_csv("long_"+sys.argv[3]+"_nom_seq_dep.tsv", sep="\t",index=False)


