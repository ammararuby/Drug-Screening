import sys
import pandas as pd
import numpy as np


def usage():
    print("""
    Wrong input parameters.

    Usage example:

    python Dilution_matrix.py 3 200 mixing_matrix.csv dose_matrix_stock.csv

    - 3 is the number of replicates
    - 200 is the assay volume
    - mixing_matrix.csv is the file containing the experimental drug combinations where each row represents a specific drug and each 
        column represent a different drug combination. Drug concentrations are represented as coefficient or coded concentrations. 
    - dose_matrix_stock.csv is the file containing the dosing information for each drug (real concentations corresponding to coded concentrations) and 
        the initial stock concentration for each drug


        NOTE: if some drugs have fewer than the maximum number of concentrations, always use the highest coded concentrations possible
        and leave empty lower concentrations. In example, drug 1 has only 1 dose and drug 2 has 3 doses
        i.e.
               Dose 1   Dose2   Dose3
        Drug 1  0       0       5
        Drug 2  1       2       5
    """)


# Parse arguments
def parse_args(args):
    if len(args) != 5:
        usage()
        sys.exit()

    num_replicates = int(args[1])
    assay_volume = int(args[2])
    mixing_matrix = pd.read_csv(args[3], delimiter=",", header=0, index_col=0)
    dose_matrix_stock = pd.read_csv(args[4], delimiter=",", header=0, index_col=0)
    return num_replicates, assay_volume, mixing_matrix, dose_matrix_stock


# Define basic variables used in downstream code
def basic_metrics(mixing_matrix, num_replicates, assay_volume, dose_matrix_stock):
    num_drugs, num_mixtures = mixing_matrix.shape[0], mixing_matrix.shape[1]
    # If it starts at 1 or 0
    min_num_doses = mixing_matrix.min().min()
    max_num_doses = mixing_matrix.max().max()
    dose_matrix = dose_matrix_stock.drop(['Stock'], axis = 1) 
    #dose_matrix = dose_matrix_stock[0:(max_num_doses-1)]
    # Add 10% to the base assay volume
    base_volume = assay_volume * num_replicates * 1.1 / num_drugs
    return num_drugs, num_mixtures, min_num_doses, max_num_doses, base_volume, dose_matrix


# Construct matrix that counts the number of occurences for each drug at each dose present in all mixtures
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


# Construct volume matrix. Calculate the volume of each drug dose needed.
def gen_volume_matrix(dose_count_matrix, base_volume):
    volume_matrix = dose_count_matrix*int(base_volume)
    return volume_matrix
    

# Construct table of dilution factors. Table containing the dilution factor needed for serial dilutions, i.e. the ratio
# between each drug dose and the next higest dose. Rhe highest drug dose should be diluted from the stock.
def dilution_factors(dose_matrix_stock, mixing_matrix, min_num_doses, max_num_doses):
    # Initialize
    num_doses = max_num_doses + 1 - min_num_doses
    num_drugs = len(mixing_matrix.index)
    dilution_factor_matrix = pd.DataFrame(
        index=mixing_matrix.index)

    cols = dose_matrix_stock.columns
    for i in range(0, len(cols) - 1):
        first = dose_matrix_stock.iloc[:, i]
        second = dose_matrix_stock.iloc[:, i + 1]
        
        dose_name = cols[i]
        dilution_factor_matrix[dose_name] = second / first
    return dilution_factor_matrix
        



def main():
    num_replicates, assay_volume, mixing_matrix, dose_matrix_stock = parse_args(sys.argv)
    num_drugs, num_mixtures, min_num_doses, max_num_doses, base_volume, dose_matrix = basic_metrics(mixing_matrix, num_replicates, assay_volume, dose_matrix_stock)
    dose_count_matrix = count_doses(mixing_matrix, min_num_doses, max_num_doses)
    volume_matrix = gen_volume_matrix(dose_count_matrix, base_volume)
    dilution_matrix = dilution_factors(dose_matrix_stock, mixing_matrix, min_num_doses, max_num_doses)
    print (dilution_matrix)


if __name__ == '__main__':
    main()