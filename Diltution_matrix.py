### Loading Libraries ###
import pandas as pd
import numpy as np

### This script replicates the excel analysis ###

### Toy Data Construction ###
np.random.seed(123)
data = pd.DataFrame({"Drugs":np.arange(1,6)})
data = pd.concat([data, pd.DataFrame(np.random.randint(0, 7, size = (5,6)))], axis = 1)
data.columns = ['Drugs'] + list(np.arange(1,7))


### Initializations ###
num_replicates = 3
stock_concentration = pd.DataFrame({"Drugs":np.arange(1,6),
                                    "Stock":np.random.randint(100, 500, size = 5)})
assay_volume = 200


### Creating tuple of important values ###
def basic_metrics(data, num_replicates, assay_volume):
    num_drugs, num_mixtures = data.shape[0], data.shape[1] - 1
    num_doses = data.drop(["Drugs"], axis = 1).max().max()
    base_volume = assay_volume * num_replicates * (1 + 0.1) / num_drugs
    return(num_drugs, num_mixtures, num_doses, base_volume)


### Constructing Count Matrix ###
num_drugs, num_mixtures, num_doses, base_volume = basic_metrics(data, num_replicates, assay_volume)
def row_count(data):
    df = data[["Drugs"]]
    for i in range(0, num_mixtures + 1):
        df[str(i)] = np.sum(data == i)
    return(df)

count_matrix = data.apply(row_count, axis = 1)
count_matrix['total'] = count_matrix.drop(["Drugs"], axis = 1).sum(axis = 1)


### Constructing Volume Matrix ###
volumes = count_matrix.drop(['Drugs'], axis = 1) * base_volume
volume_matrix = pd.concat([data[['Drugs']], volumes], axis = 1)


### Constructing Concentration Matrix ###
np.random.seed(123)
doses = pd.DataFrame(np.zeros((num_drugs, num_doses)))
doses.columns = np.arange(1, num_doses + 1)
concentration_matrix = pd.concat([stock_concentration, doses], axis = 1)
concentration_matrix[num_doses] = concentration_matrix['Stock'] - np.random.randint(10, 40, size = num_drugs)
for i in reversed(range(1, num_doses)):
    concentration_matrix[i] = concentration_matrix[i + 1] / np.round(np.random.uniform(1.5, 2.5), 2)


### Constructing Dilution Matrix ###
dilution_matrix = pd.concat([concentration_matrix.drop(['Stock'], axis = 1),
                             concentration_matrix[['Stock']]], axis = 1)
dilution_matrix = pd.melt(dilution_matrix,
                          id_vars = ['Drugs'],
                          value_vars = ['Stock'] + list(np.arange(1, num_doses + 1)),
                          var_name = 'Dose', value_name = "Concentration").sort_values(['Drugs', 'Dose'])
dilution_matrix.reset_index(inplace = True, drop = True)
dilution_matrix['lag'] = dilution_matrix.groupby(["Drugs"]).shift(-1)['Concentration']
dilution_matrix['dilution_factor'] = dilution_matrix['lag'] / dilution_matrix['Concentration']
dilution_matrix = dilution_matrix.drop(['lag'], axis = 1).fillna(0)

### Constructing Volume Matrix ###
buffer_volume = 50
volume_matrix.loc[:, '0':str(num_doses)] = volume_matrix.loc[:, '0':str(num_doses)] + buffer_volume
volume_matrix = volume_matrix.drop(['total'], axis = 1)
volume_matrix['Stock'] = stock_concentration['Stock']

# Converting to long format:
volumes = pd.melt(volume_matrix.drop(['0'], axis = 1),
                  id_vars = ['Drugs'],
                  value_vars = ['Stock'] + list(np.arange(1, num_doses + 1).astype('str')),
                  var_name = 'Dose', value_name = "Volume").sort_values(['Drugs', 'Dose'])[["Volume"]]

volumes.reset_index(inplace = True, drop = True)

### Constructing Final Matrix ###
final_matrix = pd.concat([dilution_matrix, volumes], axis = 1)
final_matrix.loc[final_matrix['Dose'] == 'Stock', 'Volume'] = pd.NA
final_matrix = final_matrix.sort_values(['Drugs', 'Dose'], ascending = [1, 0])
final_matrix.reset_index(inplace = True, drop = True)
final_matrix['volume_needed_backend'] = final_matrix['Volume'] / final_matrix['dilution_factor']
final_matrix['volume_needed_from'] = final_matrix['volume_needed_backend'].shift(-1)
final_matrix.loc[final_matrix['Dose'] == 1, 'volume_needed_from'] = 0
final_matrix['total_volume'] = final_matrix['Volume'] + final_matrix['volume_needed_from']
final_matrix['volume_stock'] = final_matrix['total_volume'] / final_matrix['dilution_factor']
final_matrix['volume_medium'] = final_matrix['total_volume'] - final_matrix['volume_stock']
final_matrix = final_matrix.drop(['volume_needed_backend'], axis = 1)

#Output is in the long format

if __name__ == '__main__':
    print(final_matrix)
