import pandas as pd

# Load the data into DataFrames
cust1 = pd.read_csv("outputs/PACA_2km/cust_weights_2km.txt", delim_whitespace=True)
cust2 = pd.read_csv("outputs/PACA_2km/cust_weights_2km_fake.txt", delim_whitespace=True)

# sum the weighjt of the two files and print
print(cust1["weight"].sum())
print(cust2["weight"].sum())


# Merge the data on coordinates x_LAMB93 and y_LAMB93 with an outer join
merged = pd.merge(
    cust1,
    cust2,
    on=["x_LAMB93", "y_LAMB93"],
    how="outer",
    indicator=True,
    suffixes=('_cust1', '_cust2')
)

# Filter rows that are only in one of the files
diff_points = merged[merged["_merge"] != "both"]

# Select relevant columns: coordinates and weight
relevant_data = diff_points[[
    "x_LAMB93", "y_LAMB93", 
    "weight_cust1", "weight_cust2", 
    "_merge"
]]


# Replace NaN weights with 0 for clarity
relevant_data.fillna({'weight_cust1': 0, 'weight_cust2': 0}, inplace=True)

# Print the results
print(relevant_data)

# Save the results to a file for convenience
relevant_data.to_csv("diff_points_coordinates_weights.txt", sep="\t", index=False)

