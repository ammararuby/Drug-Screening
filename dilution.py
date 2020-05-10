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
    min_num_doses = 1
    max_num_doses = mixing_matrix.max().max()
    # define dosage matrix by removing stock column and multiplying drug doses by number of drug to account 
    # for diluting drugs with other drugs in final mixing step
    dose_matrix = dose_matrix_stock.drop(['Stock'], axis = 1) * num_drugs
    dose_matrix['Stock'] = dose_matrix_stock['Stock']
    # Add 10% to the base assay volume
    base_volume = assay_volume * num_replicates * 1.1 / num_drugs
    return num_mixtures, min_num_doses, max_num_doses, base_volume, dose_matrix


# Construct matrix that counts the number of occurences for each drug at each dose present in all mixtures
def count_doses(mixing_matrix, min_num_doses, max_num_doses):
    # Initialize 
    num_drugs = len(mixing_matrix.index)
    num_doses =  max_num_doses + 1 - min_num_doses
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
    # add a buffer volume of 50 ul for experimental purposes 
    volume_matrix = dose_count_matrix * int(base_volume) + 50
    return volume_matrix
    

# Construct table of dilution factors. Table containing the dilution factor needed for serial dilutions, i.e. the ratio
# between each drug dose and the next higest dose. Rhe highest drug dose should be diluted from the stock.
def calc_dilution_factors(dose_matrix_stock, mixing_matrix, min_num_doses, max_num_doses):
    dilution_factor = pd.DataFrame(
        index=mixing_matrix.index)

    cols = dose_matrix_stock.columns
    for i in range(0, len(cols) - 1):
        first = dose_matrix_stock.iloc[:, i]
        second = dose_matrix_stock.iloc[:, i + 1]
        
        dose_name = cols[i]
        dilution_factor[dose_name] = first/second
    return dilution_factor


# Calculate serial dilution volumes
def calc_dilution_volumes(dilution_factor, volume_matrix, min_num_doses, max_num_doses):
    drug_names = volume_matrix.index
    dilution_factors = volume_matrix.columns
    dilution_volumes_stock = pd.DataFrame(index=drug_names, columns=dilution_factors)
    total_volume = pd.DataFrame(index=drug_names, columns=dilution_factors)
    for i in range(0, len(drug_names)):
        vol_needed = 0
        for j in range(0, len(dilution_factors)):
            current_dilution = dilution_factor.iloc[i, j]
            current_volume = volume_matrix.iloc[i, j]
            total_volume.at[drug_names[i], dilution_factors[j]] = current_volume + vol_needed
            vol_needed = current_dilution * (current_volume + vol_needed)
            dilution_volumes_stock.at[drug_names[i], dilution_factors[j]] = vol_needed
    dilution_volumes_medium = total_volume - dilution_volumes_stock
    return dilution_volumes_stock, dilution_volumes_medium


def calc_mixture_medium(mixing_matrix, base_volume, num_mixtures):
    cols = mixing_matrix.columns
    mixture_medium = pd.Series(index=mixing_matrix.columns, dtype=float)
    for col in cols:
        mixture_medium.at[col] = (mixing_matrix.loc[:, col] == 0).sum() * base_volume
    return mixture_medium


def main():
    num_replicates, assay_volume, mixing_matrix, dose_matrix_stock = parse_args(sys.argv)
    num_mixtures, min_num_doses, max_num_doses, base_volume, dose_matrix = basic_metrics(mixing_matrix, num_replicates, assay_volume, dose_matrix_stock)
    dose_count_matrix = count_doses(mixing_matrix, min_num_doses, max_num_doses)
    volume_matrix = gen_volume_matrix(dose_count_matrix, base_volume)
    dilution_factors = calc_dilution_factors(dose_matrix, mixing_matrix, min_num_doses, max_num_doses)
    dilution_volumes, dilution_volumes_medium = calc_dilution_volumes(dilution_factors, volume_matrix, min_num_doses, max_num_doses)
    mixture_medium = calc_mixture_medium(mixing_matrix, base_volume, num_mixtures)

    dilution_volumes.to_csv('dilution_volumes.csv')
    dilution_volumes_medium.to_csv('dilution_volumes_medium.csv')
    mixture_medium.to_csv('mixture_medium.csv')


if __name__ == '__main__':
    main()
