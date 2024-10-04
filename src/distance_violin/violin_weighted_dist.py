import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

# Step 1: Parse CUSTOMER ASSIGNMENTS from sol.txt
# Assuming sol.txt is structured as you posted, parse the data
def parse_customer_assignments(file_path):
    customer_assignments = []
    with open(file_path, 'r') as f:
        in_customer_section = False
        for line in f:
            if "CUSTOMER ASSIGNMENTS" in line:
                in_customer_section = True
                continue
            if in_customer_section:
                if '->' not in line:
                    break
                if line.strip().startswith("customer"):  # Skip the header line
                    continue
                
                parts = line.split('->')
                customer_part = parts[0].strip()
                locations_part = parts[1].strip()

                customer = int(customer_part.split()[0])  # Extract customer ID
                demand = float(customer_part.split('(')[1].split(')')[0])  # Extract total customer demand

                # Handle cases where customer is assigned to multiple locations
                locations_assignments = locations_part.split()
                for i in range(0, len(locations_assignments), 2):
                    location = int(locations_assignments[i])
                    assigned_demand = float(locations_assignments[i + 1].split('(')[1].split(')')[0])
                    customer_assignments.append((customer, demand, location, assigned_demand))

    return pd.DataFrame(customer_assignments, columns=['Customer', 'TotalDemand', 'Location', 'AssignedDemand'])


# Step 2: Parse distances from dist.txt
def parse_distances(file_path):
    return pd.read_csv(file_path, sep=' ', header=0)

# Step 3: Merge dataframes to calculate wi * dij
def calculate_weights(assignments_df, distances_df):
    merged_df = pd.merge(assignments_df, distances_df, left_on=['Customer', 'Location'], right_on=['customer', 'location'])
    merged_df['wi_dij'] = merged_df['AssignedDemand'] * merged_df['distance']
    return merged_df

def calculate_distances_unweighted(assignments_df, distances_df):
    merged_df = pd.merge(assignments_df, distances_df, left_on=['Customer', 'Location'], right_on=['customer', 'location'])
    merged_df['wi_dij'] = merged_df['distance']
    return merged_df


# Step 4: Create the violin plot
def plot_violin(data):
    plt.figure(figsize=(12, 8))  # Increase size for better visibility
    sns.violinplot(x='Location', y='wi_dij', data=data, inner="quartile")  # Use inner="quartile" for better visualization
    plt.title('Distribution of wi * dij (Assigned Demand * Distance) by Location')
    plt.ylabel('wi * dij')
    plt.xlabel('Location')
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.tight_layout()  # Adjust layout to prevent clipping of labels
    plt.show()


# Step 4: Create a function to plot a single violin for all wi * dij values
def plot_single_violin(data):
    # Extract all wi * dij values into a single series
    all_wi_dij = data['wi_dij'].values  # Collect all the wi * dij values

    # Create a DataFrame for plotting
    plot_data = pd.DataFrame({'wi_dij': all_wi_dij})

    plt.figure(figsize=(10, 6))  # Set figure size
    sns.violinplot(x=[1] * len(plot_data), y='wi_dij', data=plot_data)  # Use a single x value for the violin
    plt.title('Distribution of wi * dij (Assigned Demand * Distance) Across All Customers and Locations')
    plt.ylabel('wi * dij')
    plt.xlabel('All Locations')  # Label for x-axis
    plt.xticks([])  # Remove x-ticks since we have only one value
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.show()


def plot_distribution(data):
    plt.figure(figsize=(10, 6))
    sns.histplot(data['wi_dij'], bins=30, kde=True, color='blue', edgecolor='black')
    plt.title('Distribution of wi * dij Values')
    plt.xlabel('wi * dij')
    plt.ylabel('Frequency')
    plt.show()

# Function to plot combined distributions
def plot_combined_distribution(data1, data2, label1, label2):
    plt.figure(figsize=(12, 6))
    
    # Plot histograms
    sns.histplot(data1, bins=30, color='blue', label=label1, kde=True, stat="density", alpha=0.5)
    sns.histplot(data2, bins=30, color='orange', label=label2, kde=True, stat="density", alpha=0.5)
    
    plt.title('Comparison of Distribution of wi * dij Values')
    plt.xlabel('wi * dij')
    plt.ylabel('Density')
    plt.legend()
    plt.show()

def compare_distributions(wi_dij_1, wi_dij_2):
    # Check for empty datasets
    if wi_dij_1.empty:
        print("Warning: wi_dij_1 is empty!")
    if wi_dij_2.empty:
        print("Warning: wi_dij_2 is empty!")

    # Box plot
    plt.figure(figsize=(12, 6))
    
    # Create a combined list of the two datasets
    data = [wi_dij_1, wi_dij_2]
    
    # Ensure data is flattened
    flat_data = [d.flatten() for d in data]
    
    # Create the boxplot
    sns.boxplot(data=flat_data, palette=["blue", "orange"])
    plt.xticks([0, 1], ['Solution 1', 'Solution 2'])
    plt.title('Box Plot Comparison of wi * dij Values')
    plt.ylabel('wi * dij')
    plt.show()

    # KDE plot with subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # KDE for Solution 1
    sns.kdeplot(wi_dij_1, ax=axes[0], color='blue', label='Solution 1', fill=True)
    axes[0].set_title('KDE of Solution 1')
    axes[0].set_xlabel('wi * dij')
    
    # KDE for Solution 2
    if not wi_dij_2.empty:  # Only plot if not empty
        sns.kdeplot(wi_dij_2, ax=axes[1], color='orange', label='Solution 2', fill=True)
    else:
        axes[1].text(0.5, 0.5, 'No data available for Solution 2', 
                      horizontalalignment='center', verticalalignment='center', 
                      transform=axes[1].transAxes)

    axes[1].set_title('KDE of Solution 2')
    axes[1].set_xlabel('wi * dij')
    
    plt.tight_layout()
    plt.show()

# def plot_distributions_step(wi_dij_1, wi_dij_2):
#     # Combine the two datasets into one DataFrame with a 'Solution' label
#     data_combined = pd.concat([pd.DataFrame({'wi_dij': wi_dij_1, 'Solution': 'Solution 1'}),
#                                pd.DataFrame({'wi_dij': wi_dij_2, 'Solution': 'Solution 2'})])

#     # Step plot with displot
#     sns.displot(data_combined, x="wi_dij", hue="Solution", element="step", kind="kde", fill=False)

#     plt.title('Distribution of wi * dij for Solution 1 and Solution 2')
#     plt.xlabel('wi * dij')
#     plt.ylabel('Density')
#     plt.show()


def plot_distributions_kde(wi_dij_1, wi_dij_2):
    # Combine the two datasets into one DataFrame with a 'Solution' label
    data_combined = pd.concat([pd.DataFrame({'wi_dij': wi_dij_1, 'Solution': 'Solution 1'}),
                               pd.DataFrame({'wi_dij': wi_dij_2, 'Solution': 'Solution 2'})])

    # Use kdeplot instead of displot
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=data_combined, x="wi_dij", hue="Solution", fill=False)

    plt.title('Distribution of wi * dij for Solution 1 and Solution 2 (KDE)')
    plt.xlabel('wi * dij')
    plt.ylabel('Density')
    plt.show()
def plot_distributions_hist(wi_dij_1, wi_dij_2):
    # Combine the two datasets into one DataFrame with a 'Solution' label
    data_combined = pd.concat([pd.DataFrame({'wi_dij': wi_dij_1, 'Solution': 'Solution 1'}),
                               pd.DataFrame({'wi_dij': wi_dij_2, 'Solution': 'Solution 2'})])

    # Use histplot with step elements
    plt.figure(figsize=(10, 6))
    sns.histplot(data=data_combined, x="wi_dij", hue="Solution", element="step", stat="density", common_norm=False)

    plt.title('Distribution of wi * dij for Solution 1 and Solution 2 (Histogram)')
    plt.xlabel('wi * dij')
    plt.ylabel('Density')
    plt.show()


def plot_distributions_kde_intervals(wi_dij_1, wi_dij_2):
    # Combine the two datasets into one DataFrame with a 'Solution' label
    data_combined = pd.concat([pd.DataFrame({'wi_dij': wi_dij_1, 'Solution': 'Solution 1'}),
                               pd.DataFrame({'wi_dij': wi_dij_2, 'Solution': 'Solution 2'})])
    
    # Define the intervals (you can adjust these according to your data distribution)
    intervals = [(0, 5000), (5000, 15000), (15000, data_combined['wi_dij'].max())]

    # Create subplots with 1 row and 3 columns
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for i, (lower, upper) in enumerate(intervals):
        # Subset data for the current interval
        interval_data = data_combined[(data_combined['wi_dij'] >= lower) & (data_combined['wi_dij'] < upper)]
        
        # KDE plot for the current interval
        sns.kdeplot(data=interval_data, x="wi_dij", hue="Solution", ax=axes[i], fill=False)
        
        # Set titles and labels for each plot
        axes[i].set_title(f'KDE Plot for Interval [{lower}, {upper})')
        axes[i].set_xlabel('wi_dij')
        axes[i].set_ylabel('Density')

    # Show the plot
    plt.tight_layout()
    plt.show()


def plot_distributions_log(wi_dij_1, wi_dij_2):
    # Apply log transformation to the data
    log_wi_dij_1 = np.log1p(wi_dij_1)  # log(1+x) to handle zero values
    log_wi_dij_2 = np.log1p(wi_dij_2)

    # Combine the datasets into one DataFrame with a 'Solution' label
    data_combined = pd.concat([pd.DataFrame({'log_wi_dij': log_wi_dij_1, 'Solution': 'Solution 1'}),
                               pd.DataFrame({'log_wi_dij': log_wi_dij_2, 'Solution': 'Solution 2'})])

    # KDE plot with displot
    sns.displot(data_combined, x="log_wi_dij", hue="Solution", kind="kde", fill=False)

    plt.title('Log-Transformed KDE of wi * dij for Solution 1 and Solution 2')
    plt.xlabel('log(1 + wi * dij)')
    plt.ylabel('Density')
    plt.show()


def plot_outliers_boxplot(wi_dij_1, wi_dij_2):
    # Combine datasets and create labels
    data_combined = pd.concat([pd.DataFrame({'wi_dij': wi_dij_1, 'Solution': 'Solution 1'}),
                               pd.DataFrame({'wi_dij': wi_dij_2, 'Solution': 'Solution 2'})])

    # Boxplot to show distribution and outliers
    sns.boxplot(data=data_combined, x="Solution", y="wi_dij", showfliers=True, whis=1.5)

    plt.title('Comparison of Outliers between Solutions 1 and 2')
    plt.ylabel('wi * dij')
    plt.ylim(0, None)  # Adjust y-limits to focus on outliers if needed
    plt.show()

def plot_outliers_violin(wi_dij_1, wi_dij_2):
    # Combine datasets and create labels
    data_combined = pd.concat([pd.DataFrame({'wi_dij': wi_dij_1, 'Solution': 'Solution 1'}),
                               pd.DataFrame({'wi_dij': wi_dij_2, 'Solution': 'Solution 2'})])

    # Violin plot to show full distribution with focus on outliers
    sns.violinplot(data=data_combined, x="Solution", y="wi_dij", cut=0)  # cut=0 limits tails to observed data

    plt.title('Violin Plot Comparison of Outliers')
    plt.ylabel('wi * dij')
    plt.ylim(0, None)  # Adjust y-limits to focus on outliers if needed
    plt.show()

def plot_outliers_scatter(wi_dij_1, wi_dij_2, threshold=2e3):
    # Filter for extreme values
    extreme_wi_dij_1 = wi_dij_1[wi_dij_1 > threshold]
    extreme_wi_dij_2 = wi_dij_2[wi_dij_2 > threshold]

    # Create DataFrames with labels
    data_extreme = pd.concat([pd.DataFrame({'wi_dij': extreme_wi_dij_1, 'Solution': 'Solution 1'}),
                              pd.DataFrame({'wi_dij': extreme_wi_dij_2, 'Solution': 'Solution 2'})])

    # Scatter plot to compare extreme values
    sns.scatterplot(data=data_extreme, x="Solution", y="wi_dij", hue="Solution", marker="o")

    plt.title(f'Scatter Plot of Outliers (wi * dij > {threshold})')
    plt.ylabel('wi * dij')

    # add count of extreme values
    plt.text(0, threshold, f'Count: {len(extreme_wi_dij_1)}', ha='center', va='center', color='red')
    plt.text(1, threshold, f'Count: {len(extreme_wi_dij_2)}', ha='center', va='center', color='red')

    # # print id of extreme values
    # for i, val in enumerate(extreme_wi_dij_1):
    #     plt.text(0, val, f'{i}', ha='center', va='center', color='red')
    # for i, val in enumerate(extreme_wi_dij_2):
    #     plt.text(1, val, f'{i}', ha='center', va='center', color='red') 


    plt.show()


def plot_qq(wi_dij_1, wi_dij_2):
    plt.figure(figsize=(8, 6))

    # Q-Q plot comparing two distributions
    stats.probplot(wi_dij_1, dist="norm", plot=plt, label='Solution 1')
    stats.probplot(wi_dij_2, dist="norm", plot=plt, label='Solution 2')

    plt.title('Q-Q Plot: Solution 1 vs Solution 2')
    plt.legend()
    plt.show()


def plot_outliers_ecdf(wi_dij_1, wi_dij_2):
    sns.ecdfplot(wi_dij_1, label='Solution 1')
    sns.ecdfplot(wi_dij_2, label='Solution 2')

    plt.title('ECDF of wi * dij for Solutions 1 and 2')
    plt.xlabel('wi * dij')
    plt.ylabel('ECDF')
    plt.legend()
    plt.show()


def plot_distance_outliers_scatter(dist_1, dist_2, dist_3, dist_4, threshold=3600):
    # Step 1: Combine distance datasets and filter for distances greater than the threshold
    data_combined_distances = pd.concat([pd.DataFrame({'distance': dist_1, 'Cover type': 'No Cover'}),
                                         pd.DataFrame({'distance': dist_2, 'Cover type': 'Cover EPCI'}),
                                         pd.DataFrame({'distance': dist_3, 'Cover type': 'Cover Canton'}),
                                         pd.DataFrame({'distance': dist_4, 'Cover type': 'Cover Commune'})])
    
    # Step 2: Filter distances above the threshold (e.g., 3600)
    outliers_data = data_combined_distances[data_combined_distances['distance'] > threshold]
    
    # Step 3: Count the number of outliers for each Cover type in the correct order
    cover_order = ['No Cover', 'Cover EPCI', 'Cover Canton', 'Cover Commune']
    outlier_counts = outliers_data['Cover type'].value_counts().reindex(cover_order).fillna(0).to_dict()

    # Step 4: Create a scatter plot for the outliers without the legend
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=outliers_data, x="Cover type", y="distance", hue="Cover type", s=100, legend=False)

    # Step 5: Add titles and labels
    plt.title(f'Outliers in Distance (> {threshold}) between Cover Models (Cinema p = 192)')
    plt.ylabel('Distance (dij)')
    plt.xlabel('Cover type')

    # Step 6: Annotate the plot with the number of outliers in the correct x positions
    for i, cover_type in enumerate(cover_order):
        count = int(outlier_counts[cover_type])
        plt.text(x=i, 
                 y=outliers_data['distance'].max() + 100,  # Adjust for label placement
                 s=f'Count: {count}',
                 horizontalalignment='center',
                 color='red',
                 fontweight='bold')

    # Step 7: Show the plot
    plt.tight_layout()
    plt.show()


def plot_outliers_boxplot_and_line(wi_dij_1, wi_dij_2, wi_dij_3, wi_dij_4, dist_1, dist_2, dist_3, dist_4):
    # Step 1: Combine wi_dij datasets for the boxplot
    data_combined_wi_dij = pd.concat([pd.DataFrame({'value': wi_dij_1, 'Cover type': 'No Cover'}),
                                      pd.DataFrame({'value': wi_dij_2, 'Cover type': 'Cover EPCI'}),
                                      pd.DataFrame({'value': wi_dij_3, 'Cover type': 'Cover Canton'}),
                                      pd.DataFrame({'value': wi_dij_4, 'Cover type': 'Cover Commune'})])

    # Step 2: Combine distance datasets for the second plot (same structure as wi_dij)
    data_combined_distances = pd.concat([pd.DataFrame({'value': dist_1, 'Cover type': 'No Cover'}),
                                         pd.DataFrame({'value': dist_2, 'Cover type': 'Cover EPCI'}),
                                         pd.DataFrame({'value': dist_3, 'Cover type': 'Cover Canton'}),
                                         pd.DataFrame({'value': dist_4, 'Cover type': 'Cover Commune'})])

    # Step 3: Create subplots with two rows: one for boxplots (wi_dij), one for line plot (distance)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), gridspec_kw={'height_ratios': [2, 1]})

    # Step 4: Plot boxplot for wi_dij values
    sns.boxplot(data=data_combined_wi_dij, x="Cover type", y="value", showfliers=True, whis=1.5, ax=ax1)
    ax1.set_title('Comparison of wi * dij between Solutions of Covers Models (Cinema p = 192)')
    ax1.set_ylabel('wi * dij')
    ax1.set_ylim(0, None)  # Adjust y-limits to focus on outliers if needed

    # Step 5: Add a title between the plots
    # fig.suptitle('Comparison of Solutions: wi * dij and Distance', fontsize=16, y=0.92)

    # Step 6: Plot boxplot for distances (same x-axis, just different y-values)
    sns.boxplot(data=data_combined_distances, x="Cover type", y="value", showfliers=True, whis=1.5, ax=ax2)
    ax2.set_title('Comparison of Distances between Solutions of Covers Models (Cinema p = 192)')
    ax2.set_ylabel('Distance (dij)')
    ax2.set_ylim(0, None)  # Adjust y-limits if needed

    # Step 7: Adjust the layout and show the plot
    plt.tight_layout()
    plt.show()



# File paths
# sol_file = 'data/solutions_service_cinema/test_paca_cinema_canton_p_192_EXACT_CPMP_cover_canton.txt'
sol_file = 'data/solutions_service_cinema/test_paca_cinema_p_192_EXACT_CPMP.txt'
dist_file = '/home/falbuquerque/Documents/projects/GeoAvignon/PMPSolver/data/PACA_jul24/dist_matrix_minutes_2037.txt'

# # Execute the steps
# assignments_df = parse_customer_assignments(sol_file)
# distances_df = parse_distances(dist_file)
# merged_df = calculate_weights(assignments_df, distances_df)

# print(distances_df['distance'].describe())
# print(distances_df[distances_df['distance'] < 0])

# print(merged_df['wi_dij'].describe())
# print(merged_df[merged_df['wi_dij'] < 0])# Plot the violin plot

# print(merged_df.head())

# # plot_violin(merged_df)
# # plot_single_violin(merged_df)
# plot_distribution(merged_df)

# # calulate the sum of wi * dij for each location
# sum_wi_dij = merged_df.groupby('Location')['wi_dij'].sum()
# # summ all
# sum_wi_dij = sum_wi_dij.sum()
# print(sum_wi_dij)



# File paths for both solutions
sol_file_1 = 'data/solutions_service_cinema/test_paca_cinema_p_192_EXACT_CPMP.txt'
sol_file_2 = 'data/solutions_service_cinema/test_paca_cinema_EPCI_p_192_EXACT_CPMP_cover_EPCI.txt'
sol_file_3 = 'data/solutions_service_cinema/test_paca_cinema_canton_p_192_EXACT_CPMP_cover_canton.txt'
sol_file_4 = 'data/solutions_service_cinema/test_paca_cinema_commune_p_192_EXACT_CPMP_cover_commune.txt'
dist_file = '/home/falbuquerque/Documents/projects/GeoAvignon/PMPSolver/data/PACA_jul24/dist_matrix_minutes_2037.txt'  # Ensure to use the correct distance file

# Step 1: Parse both solution files
assignments_df_1 = parse_customer_assignments(sol_file_1)
assignments_df_2 = parse_customer_assignments(sol_file_2)
assignments_df_3 = parse_customer_assignments(sol_file_3)
assignments_df_4 = parse_customer_assignments(sol_file_4)
distances_df = parse_distances(dist_file)

# Step 2: Calculate weighted distances for both datasets
merged_df_1 = calculate_weights(assignments_df_1, distances_df)
merged_df_2 = calculate_weights(assignments_df_2, distances_df)
merged_df_3 = calculate_weights(assignments_df_3, distances_df)
merged_df_4 = calculate_weights(assignments_df_4, distances_df)


# Example usage (assuming wi_dij and distance values are from your merged DataFrames):
# plot_outliers_boxplot_and_line(
#     merged_df_1['wi_dij'], merged_df_2['wi_dij'], merged_df_3['wi_dij'], merged_df_4['wi_dij'],
#     merged_df_1['distance'], merged_df_2['distance'], merged_df_3['distance'], merged_df_4['distance']
# )


plot_distance_outliers_scatter(merged_df_1['distance'], merged_df_2['distance'], merged_df_3['distance'], merged_df_4['distance'], threshold=3600)

exit()


print(merged_df_1.head())


# Extract the wi * dij values for comparison
wi_dij_1 = merged_df_1['wi_dij']
wi_dij_2 = merged_df_2['wi_dij']
wi_dij_3 = merged_df_3['wi_dij']
wi_dij_4 = merged_df_4['wi_dij']


print(wi_dij_1.head())
exit()



# Step 3: Plot the combined distributions
# plot_combined_distribution(wi_dij_1, wi_dij_2, 'Solution 1', 'Solution 2')
# compare_distributions(wi_dij_1, wi_dij_2)
# plot_distributions_step(wi_dij_1, wi_dij_2)
# plot_distributions_hist(wi_dij_1, wi_dij_2)
# plot_distributions_kde(wi_dij_1, wi_dij_2)
# plot_distributions_kde_intervals(wi_dij_1, wi_dij_2)    
# plot_distributions_log(wi_dij_1, wi_dij_2)  

plot_outliers_boxplot(wi_dij_1, wi_dij_2)

# plot_outliers_violin(wi_dij_1, wi_dij_2)

plot_outliers_scatter(wi_dij_1, wi_dij_2, threshold=3600)

# compare average value and stdev 
print(wi_dij_1.mean(), wi_dij_1.std())
print(wi_dij_2.mean(), wi_dij_2.std())


# plot_qq(wi_dij_1, wi_dij_2)

# plot_outliers_ecdf(wi_dij_1, wi_dij_2)