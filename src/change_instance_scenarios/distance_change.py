import pandas as pd
import os


GRID_TYPE = "2km"
dist_type = "seconds"


def modify_and_save_distances(input_file, output_dir, percentage_change):
    """
    Modify distances in the input file by a specified percentage and save the result.
    
    Args:
    - input_file (str): Path to the input CSV file containing the distance matrix.
    - output_dir (str): Directory to save the modified file.
    - percentage_change (float): Percentage by which to increase or decrease the distances.
    """
    # Load the data
    dist_matrix = pd.read_csv(input_file, sep=' ', header=0)

    print("Distance matrix (original):")    
    print(dist_matrix.head())


    # Ensure the required columns exist
    if 'distance' not in dist_matrix.columns:
        raise ValueError("The input file must contain a 'distance' column.")
    
    # Modify the distance values
    dist_matrix['distance'] = dist_matrix['distance'] * (1 + percentage_change / 100)
    
    # Prepare the output file name
    modifier = "increase" if percentage_change > 0 else "decrease"
    percentage_str = f"{abs(percentage_change):.0f}"
    base_name = os.path.basename(input_file).replace('.txt', '')
    output_file = os.path.join(output_dir, f"{base_name}_{modifier}_{percentage_str}percent.txt")
    
    # Save the modified data to a new file
    dist_matrix.to_csv(output_file, index=False)
    

    print("Distance matrix (modified):")
    print(dist_matrix.head())

    print(f"File saved: {output_file}")

# Example usage
percentage_change = 20  # Increase or decrease percentage
input_file = f'outputs/PACA_{GRID_TYPE}/dist_matrix_{dist_type}_{GRID_TYPE}.txt'  # Path to the input file
output_dir = f'outputs/PACA_{GRID_TYPE}/'  # Directory to save modified files

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Run the function
modify_and_save_distances(input_file, output_dir, percentage_change)
