import pandas as pd
import numpy as np


# Set a random seed for reproducibility
seed = 42
np.random.seed(seed)

# Load the data
df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt', delim_whitespace=True)

# Calculate the original sum of weights
original_sum = df['weight'].sum()

# Shuffle the 'weight' column
shuffled_weights = np.random.permutation(df['weight'].values)

# Replace the original 'weight' column with shuffled weights
df['weight'] = shuffled_weights

# Calculate the new sum of weights
new_sum = df['weight'].sum()

# Save the new DataFrame to a new text file
df.to_csv('cust_weights_PACA_2037_shuffle.txt', sep=' ', index=False)

# Check if the sums are equal
if original_sum == new_sum:
    print("The sum of the weights remains the same.")
else:
    print("The sum of the weights has changed!")