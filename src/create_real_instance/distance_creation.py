import pandas as pd


GRID_TYPE = "2km"
dist_type = "metres"


def MergeAndFilter_dist_matrix():

    # Load distance matrix
    dist_matrix = pd.read_csv(
        f'outputs/PACA_{GRID_TYPE}/dist_files/dist_matrix_{dist_type}_{GRID_TYPE}_5435.txt', sep=' ', header=0
    )
    print("Distance matrix (original):")
    print(dist_matrix.tail())

    # Load extension file
    extension = pd.read_csv(
        '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Reseau/2km_grid/mat_2km_add_15pts_mesure/mat_15pts_mesure.csv',
        sep=';', header=None
    )
    extension.columns = ['customer', 'location', 'w1', 'w2', 'distance_metre', 'distance_secondes']
    print("Extension (original):")
    print(extension.tail())

    # Rename distance column to match
    extension = extension.rename(columns={'distance_metre': 'distance'})

    # Ensure both datasets have the same column types
    dist_matrix[['customer', 'location']] = dist_matrix[['customer', 'location']].astype(int)
    extension[['customer', 'location']] = extension[['customer', 'location']].astype(int)

    # Select only necessary columns from extension
    extension = extension[['customer', 'location', 'distance']]
    print("Extension (filtered and renamed):")
    print(extension.tail())

    # Union both datasets
    combined_data = pd.concat([dist_matrix[['customer', 'location', 'distance']], extension], ignore_index=True)
    print("Combined dataset:")
    print(combined_data.tail())

    # Create expected pairs
    unique_customers = combined_data['customer'].nunique()
    expected_pairs = pd.MultiIndex.from_product([range(1, unique_customers + 1), range(1, unique_customers + 1)], names=['customer', 'location'])

    # Create combined pairs
    combined_pairs = combined_data[['customer', 'location']].drop_duplicates().set_index(['customer', 'location'])

    # Find missing pairs
    missing_pairs = expected_pairs.difference(combined_pairs.index)
    print(f"Missing customer-location pairs: {len(missing_pairs)}")

    # Convert missing pairs to DataFrame
    missing_pairs_df = pd.DataFrame(missing_pairs.tolist(), columns=['customer', 'location'])

    # Merge missing pairs with inverse pairs in the original data
    inverse_pairs_df = combined_data[['customer', 'location', 'distance']].rename(columns={'customer': 'location', 'location': 'customer'})
    missing_pairs_with_inverse = missing_pairs_df.merge(inverse_pairs_df, on=['customer', 'location'], how='left')

    # Check if inverse distances exist
    inverse_pairs_to_add = missing_pairs_with_inverse[missing_pairs_with_inverse['distance'].notnull()]

    # Add inverse pairs to the combined data
    if not inverse_pairs_to_add.empty:
        combined_data = pd.concat([combined_data, inverse_pairs_to_add[['customer', 'location', 'distance']]], ignore_index=True)

    # check if we have number of customers * number of customers
    print(f"Number of customers: {unique_customers}")
    print(f"Number of customer-location pairs: {len(combined_data)}")
    if len(combined_data) == unique_customers * unique_customers:
        print("Number of customer-location pairs is correct.")
    else:
        print("Number of customer-location pairs is incorrect.")

    # Save the updated dataset
    file_save = f'outputs/PACA_{GRID_TYPE}/dist_files/combined_dist_matrix_{dist_type}_{GRID_TYPE}_{unique_customers}.txt'
    combined_data.to_csv(file_save, sep=" ", index=False)
    print(f"Updated dataset saved as {file_save}.")

def reindex_distance_matrix(dist_matrix_file, customers_file):
    print('Reindexing Distance Matrix')

    # Load the original distance matrix (5450 x 5450)
    dist_matrix = pd.read_csv(dist_matrix_file, sep=' ', header=0)
    max_dist_matrix_index = dist_matrix['customer'].max()  # Get max index of dist_matrix
    print(f"Original Distance Matrix has a maximum customer index of {max_dist_matrix_index}")
    print("Original Distance Matrix:")
    # print(dist_matrix.head())
    # print(dist_matrix.tail())


    ### print some specific pair of customers and their distances
    # print(dist_matrix[(dist_matrix['customer'] == 5439) & (dist_matrix['location'] == 5440)])
    # print(dist_matrix[(dist_matrix['customer'] == 5441) & (dist_matrix['location'] == 5450)])
    # print(dist_matrix[(dist_matrix['location'] == 5439) & (dist_matrix['customer'] == 5440)])
    # print(dist_matrix[(dist_matrix['location'] == 5441) & (dist_matrix['customer'] == 5450)])


    # Load the filtered customer file with 5282 points
    customers = pd.read_csv(customers_file, sep=' ', header=0)
    max_customer_index = customers['customer'].max()  # Get max index of filtered customer file
    print(f"Filtered Customers file has a maximum customer index of {max_customer_index}")
    print("Filtered Customers:")
    print(customers.tail())

    # Create a mapping from `customer_old` (from dist_matrix) to `customer` (from the filtered customers)
    customer_mapping = dict(zip(customers['customer_old'], customers['customer']))
    print("Customer Mapping (old to new):")
    # print(list(customer_mapping.items())[:5])  # Print first 5 mappings
    print(list(customer_mapping.items())[-5:])  # Print last 5 mappings
    # last five
    # print(list(customer_mapping.items())[2500:2505])
    print(list(customer_mapping.items())[4200:4205])

    # Get the list of customers that appear in both the filtered list and the original dist_matrix
    valid_customer_old = customers['customer_old'].values
    dist_matrix_filtered = dist_matrix[dist_matrix['customer'].isin(valid_customer_old)]
    dist_matrix_filtered = dist_matrix_filtered[dist_matrix_filtered['location'].isin(valid_customer_old)]

    # Reindex the filtered distance matrix based on the new `customer` indices
    dist_matrix_filtered['customer'] = dist_matrix_filtered['customer'].map(customer_mapping)
    dist_matrix_filtered['location'] = dist_matrix_filtered['location'].map(customer_mapping)

    # Now we have the long format: location, customer, distance
    dist_matrix_long = dist_matrix_filtered[['location', 'customer', 'distance']]

    # Ensure the matrix has only valid customer-location pairs (asymmetric: customer to location and vice versa)
    dist_matrix_long = dist_matrix_long.dropna(subset=['location', 'customer'])

    # check if we have number of customers * number of customers
    unique_customers = dist_matrix_long['customer'].nunique()
    print(f"Number of customers: {unique_customers}")
    print(f"Number of customer-location pairs: {len(dist_matrix_long)}")
    if len(dist_matrix_long) == unique_customers * unique_customers:
        print("Number of customer-location pairs is correct.")
    else:
        print("Number of customer-location pairs is incorrect.")


    # Save the long format dataset to a file if desired
    file_save = f'outputs/PACA_{GRID_TYPE}/dist_files/dist_matrix_{dist_type}_{GRID_TYPE}_{max_customer_index}.txt'
    dist_matrix_long.to_csv(file_save, sep=" ", index=False)
    print(f"Reindexed dist matrix saved as {file_save}.")

    # Print the head of the long format distance matrix to check
    print("Reindexed Long Format Distance Matrix:")
    # print(dist_matrix_long.head())
    # print(dist_matrix_long.tail())

    # print(dist_matrix_long[(dist_matrix_long['customer'] == 5279) & (dist_matrix_long['location'] == 5280)])
    # print(dist_matrix_long[(dist_matrix_long['customer'] == 5281) & (dist_matrix_long['location'] == 5282)])
    # print(dist_matrix_long[(dist_matrix_long['location'] == 5279) & (dist_matrix_long['customer'] == 5280)])
    # print(dist_matrix_long[(dist_matrix_long['location'] == 5281) & (dist_matrix_long['customer'] == 5282)])





# create main function
if __name__ == "__main__":
    MergeAndFilter_dist_matrix()


    dist_matrix_file = f'outputs/PACA_{GRID_TYPE}/dist_files/combined_dist_matrix_{dist_type}_{GRID_TYPE}_5450.txt'
    customers_file = f'outputs/PACA_{GRID_TYPE}/cust_weights_{GRID_TYPE}.txt'
    reindex_distance_matrix(dist_matrix_file, customers_file)