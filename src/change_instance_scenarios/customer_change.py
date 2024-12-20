import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_data(file_path):
    """Load customer data from a file."""
    return pd.read_csv(file_path, sep='\s+')

def shuffle_weights(df, seed=42):
    """Shuffle the weights of the customers."""
    np.random.seed(seed)
    shuffled_weights = np.random.permutation(df['weight'].values)
    df['weight'] = shuffled_weights
    return df

def create_constant_weights(df):
    """Create constant weights for customers while maintaining the sum of weights."""
    total_weight = df['weight'].sum()
    num_customers = len(df)
    constant_weight = total_weight / num_customers
    df['weight'] = constant_weight
    return df

def create_random_weights(df, seed=42):
    """Create random weights for customers while maintaining the same statistical features and the sum of weights."""
    np.random.seed(seed)
    total_weight = df['weight'].sum()
    mean_weight = df['weight'].mean()
    std_weight = df['weight'].std()
    min_weight = df['weight'].min()
    max_weight = df['weight'].max()
    num_customers = len(df)
    
    # Generate random weights with the same mean and standard deviation
    random_weights = np.random.normal(loc=mean_weight, scale=std_weight, size=num_customers)

    # Clip the weights to maintain the same amplitude (range)
    random_weights = np.clip(random_weights, min_weight, max_weight)
    
    # Adjust the weights to maintain the same sum
    random_weights = random_weights / random_weights.sum() * total_weight
    
    df['weight'] = random_weights
    return df

def save_data(df, file_path):
    """Save the DataFrame to a file."""
    df.to_csv(file_path, sep=' ', index=False)

def print_stats(df, label):
    """Print the statistical features of the weights."""
    print(f"{label} - Mean: {df['weight'].mean()}, Std: {df['weight'].std()}, Min: {df['weight'].min()}, Max: {df['weight'].max()}, Sum: {df['weight'].sum()}")

def plot_distributions(df_original, df_shuffled, df_constant, df_random):
    """Plot the distribution of weights side by side."""
    fig, axs = plt.subplots(1, 4, figsize=(20, 5), sharey=True)
    
    axs[0].hist(df_original['weight'], bins=30, color='blue', alpha=0.7)
    axs[0].set_title('Original Weights')
    
    axs[1].hist(df_shuffled['weight'], bins=30, color='green', alpha=0.7)
    axs[1].set_title('Shuffled Weights')
    
    axs[2].hist(df_constant['weight'], bins=30, color='orange', alpha=0.7)
    axs[2].set_title('constant Weights')
    
    axs[3].hist(df_random['weight'], bins=30, color='red', alpha=0.7)
    axs[3].set_title('Random Weights')
    
    for ax in axs:
        ax.set_xlabel('Weight')
        ax.set_ylabel('Frequency')
    
    plt.tight_layout()
    plt.show()

def main():
    # Set parameters
    seed = 42
    GRID_TYPE = '2km'
    path_customer = f'outputs/PACA_{GRID_TYPE}/cust_weights_{GRID_TYPE}.txt'
    output_shuffle_path = f'outputs/PACA_{GRID_TYPE}/cust_weights_shuffle_{GRID_TYPE}.txt'
    output_constant_path = f'outputs/PACA_{GRID_TYPE}/cust_weights_constant_{GRID_TYPE}.txt'
    output_random_path = f'outputs/PACA_{GRID_TYPE}/cust_weights_random_{GRID_TYPE}.txt'

    # Load the data
    df = load_data(path_customer)

    # Calculate the original sum of weights
    original_sum = df['weight'].sum()

    # Print stats for the original weights
    print_stats(df, "Original Weights")

    # Shuffle the weights
    df_shuffled = shuffle_weights(df.copy(), seed)
    save_data(df_shuffled, output_shuffle_path)
    print_stats(df_shuffled, "Shuffled Weights")

    # Create constant weights
    df_constant = create_constant_weights(df.copy())
    save_data(df_constant, output_constant_path)
    print_stats(df_constant, "constant Weights")

    # Create random weights
    df_random = create_random_weights(df.copy(), seed)
    save_data(df_random, output_random_path)
    print_stats(df_random, "Random Weights")

    # Calculate the new sums of weights
    new_sum_shuffled = df_shuffled['weight'].sum()
    new_sum_constant = df_constant['weight'].sum()
    new_sum_random = df_random['weight'].sum()

    # Check if the sums are equal
    if original_sum == new_sum_shuffled:
        print("The sum of the weights remains the same after shuffling.")
    else:
        print("The sum of the weights has changed after shuffling!")

    if original_sum == new_sum_constant:
        print("The sum of the weights remains the same after creating constant weights.")
    else:
        print("The sum of the weights has changed after creating constant weights!")

    if original_sum == new_sum_random:
        print("The sum of the weights remains the same after creating random weights.")
    else:
        print("The sum of the weights has changed after creating random weights!")

    # Plot the distributions
    plot_distributions(df, df_shuffled, df_constant, df_random)

if __name__ == "__main__":
    main()