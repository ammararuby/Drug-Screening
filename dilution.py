import sys
import pandas as pd
import numpy as np


def usage():
    print("""
    Wrong input parameters.

    Usage example:

    python Dilution_matrix.py 3 200 mixing_matrix.csv dose_matrix.csv stocks.csv

    - 3 is the number of replicates
    - 200 is the assay volume
    - mixing_matrix.csv is the file containing the experimental drug combinations coefficient (coded concentrations)
    - dose_matrix.csv is the file containing the dosing information for each drug (real concentations)
    - stock.csv is the file containing the initial stock concentration for each drug

    """)


def parse_args(args):
    if len(args) != 6:
        usage()
        sys.exit()

    num_replicates = int(args[1])
    assay_volume = int(args[2])
    mixing_matrix = pd.read_csv(args[3], delimiter=",", header=0, index_col=0)
    dose_matrix = pd.read_csv(args[4], delimiter=",", header=0, index_col=0)
    stock = pd.read_csv(args[5], delimiter=",", header=0, index_col=0)
    return num_replicates, assay_volume, mixing_matrix, dose_matrix, stock


def basic_metrics(mixing_matrix, num_replicates, assay_volume):
    num_drugs, num_mixtures = mixing_matrix.shape[0], mixing_matrix.shape[1]
    # If it starts at 1 or 0
    min_num_doses = mixing_matrix.min().min()
    max_num_doses = mixing_matrix.max().max()
    # Add 10% to the base assay volume
    base_volume = assay_volume * num_replicates * 1.1 / num_drugs
    return num_drugs, num_mixtures, min_num_doses, max_num_doses, base_volume


def count_doses(mixing_matrix, min_num_doses, max_num_doses):
    # Initialize
    num_doses = max_num_doses + 1 - min_num_doses
    num_drugs = len(mixing_matrix.index)
    dose_count_matrix = pd.DataFrame(
        index=mixing_matrix.index, #rows, i.e. drugs
        columns=range(min_num_doses, max_num_doses + 1), #columns, i.e. number of doses
        data=np.zeros((num_drugs, num_doses), dtype=int))

    for drug_name in dose_count_matrix.index:
        for dose_name in range(min_num_doses, max_num_doses + 1):
            drug = mixing_matrix.loc[drug_name]
            # drug.values is of type numpy.ndarray 
            # drug.values == dose_name returns a list of boolean
            dose_count_matrix.at[drug_name, dose_name] = (drug.values == dose_name).sum()
    return dose_count_matrix


### Constructing Volume Matrix ###
def volume_matrix(dose_count_matrix, ):
    
    
    #volumes = count_matrix.drop(['Drugs'], axis = 1) * base_volume
    #volume_matrix = pd.concat([data[['Drugs']], volumes], axis = 1)




def main():
    num_replicates, assay_volume, mixing_matrix, dose_matrix, stock = parse_args(sys.argv)
    num_drugs, num_mixtures, min_num_doses, max_num_doses, base_volume = basic_metrics(mixing_matrix, num_replicates, assay_volume)
    dose_count_matrix = count_doses(mixing_matrix, min_num_doses, max_num_doses)
    print(dose_count_matrix)


if __name__ == '__main__':
    main()