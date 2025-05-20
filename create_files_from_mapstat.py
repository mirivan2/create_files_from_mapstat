import sys
import os
import pandas as pd
import numpy as np

## Usage:
# python3 create_files_from_mapstat.py data/ data/ResFinder_20220825.refdata.txt P2_four


# List to store dataframes
dataframes = []

# Load refdata
refdata = pd.read_csv(sys.argv[2], sep='\t', skiprows=1, header=0)

# Insert the data for the genes we did not have data on in the refdata
refdata.loc[refdata['# id'] == 'Cfr(E)_1_NG_070225', 'id_len'] = 1035
refdata.loc[refdata['# id'] == 'tet(O/32/O)_6_NG_048124', 'id_len'] = 2120
refdata.loc[refdata['# id'] == 'bleO_1_AF051917', 'id_len'] = 399
refdata.loc[refdata['# id'] == 'tet(X6)_1_MN507533', 'id_len'] = 1137
refdata.loc[refdata['# id'] == 'erm(50)_1_LC473083', 'id_len'] = 792
refdata.loc[refdata['# id'] == 'dfrA36_1_CP038791', 'id_len'] = 588
refdata.loc[refdata['# id'] == 'tet(X5)_1_CP040912', 'id_len'] = 1167


# Iterate through all files in the folder
for filename in os.listdir(sys.argv[1]):
    if filename.endswith(".mapstat"):
        file_path = os.path.join(sys.argv[1], filename)
        
        # Read the file skipping the first 6 lines and using the first row as column headers
        df = pd.read_csv(file_path, sep='\t', skiprows=6, header=0)
        df = df[["# refSequence","refCoveredPositions","bpTotal"]]

        # Merge df with refdata using the column '# refSequence' from df and the column '# id' from refdata
        merged_df = pd.merge(df, refdata[["# id", "id_len"]], left_on='# refSequence', right_on='# id').drop("# id",axis=1)

        # Check if some genes do not have length provided and print if NA:
        if len(merged_df[merged_df["id_len"].isna()]) > 0:
            print("The following genes do not have any length information")
            print(merged_df[merged_df["id_len"].isna()],"\n")

        # Calculate coverage and remove everything below 0.85
        merged_df["coverage"] = merged_df["refCoveredPositions"] / merged_df["id_len"]
        removed = len(merged_df) - len(merged_df[merged_df["coverage"] >= 0.85])
        print("Removed {0} rows from file {1}, because coverage < 0.85".format(removed,filename))
        merged_df = merged_df[merged_df["coverage"] > 0.85]

        # Calculate normalised sequence depth
        merged_df["nom_seq_dep"] = merged_df["bpTotal"] / merged_df["id_len"].sum()

        # Merge together gene variants and use the mean of the nom_seq_dep as the new value
        merged_df["gene"] = merged_df["# refSequence"].apply(lambda x: x.split("_")[0])
        merged_df = merged_df.drop(["# refSequence","refCoveredPositions","bpTotal","id_len", "coverage"], axis=1)
        merged_df = merged_df.groupby('gene')['nom_seq_dep'].mean().reset_index()
        
        # Log10 norm seq depth
        merged_df["log10_nom_seq_dep"] = np.log10(merged_df["nom_seq_dep"])
        merged_df = merged_df.drop("nom_seq_dep",axis=1)

        # Add sample name
        merged_df["sample"] = filename

        # Append the dataframe to the list
        dataframes.append(merged_df)


# Make final dataframes
master_df = pd.concat(dataframes)

pivoted_df = master_df.pivot(index='sample', columns='gene', values='log10_nom_seq_dep')
pivoted_df.to_csv("pivot_"+sys.argv[3]+"_log10_nom_seq_dep.tsv", sep="\t")

long_df = master_df.pivot(index='gene', columns='sample', values='log10_nom_seq_dep').reset_index()
long_df.to_csv("long_"+sys.argv[3]+"_log10_nom_seq_dep.tsv", sep="\t",index=False)

