import numpy as np
import glob
import os

# Function to read parameters from parameters.dat
def read_parameters_from_dat(file_path):
    params = {}

    with open(file_path, "r") as f:
        for line in f:
            key_value = line.strip().split(maxsplit=1)
            if len(key_value) == 2:
                key, value = key_value
                params[key] = value.strip()

    # Convert to correct types
    n_var = int(params["nVar"])
    nLevels = int(params["nLevels"])
    nx = 2**(nLevels - 1)  # Compute number of parcels

    # Handle space-separated variable names from the 'varName' line
    variable_names = params["varName"].split()

    return n_var, nx, variable_names

#------------------------------------------------------------------------------

def main():
    # Define the folder for saving processed data
    output_dir = "processed_data"
    os.makedirs(output_dir, exist_ok=True)

    # Read parameters from `parameters.dat`
    params_file_path = "parameters.dat"
    n_var, nx, variable_names = read_parameters_from_dat(params_file_path)

    # Find all data files
    data_files = sorted(glob.glob("../data/rlz_00001/Data_*.dat"))
    n_files = len(data_files)

    if n_files == 0:
        raise FileNotFoundError("No data files found in ../data/")

    print(f"Found {n_files} data files for processing.")

    # Initialize arrays for mean and variance calculations
    sum_values = np.zeros((nx, n_var))      # Sum of values for each parcel
    sum_sq_values = np.zeros((nx, n_var))   # Sum of squared values for variance calculation
    count = np.zeros((nx, n_var))           # Count occurrences for each parcel

    # Process each data file
    for file_path in data_files:
        print(f"Processing file: {file_path}")

        try:
            # Load data (each row corresponds to a parcel, each column is a variable)
            data = np.loadtxt(file_path)

            # Check if data has the correct shape
            if data.shape[1] != n_var:
                print(f"Warning: {file_path} has {data.shape[1]} columns, expected {n_var}. Skipping.")
                continue

            num_rows = data.shape[0]  # Number of valid parcels in this file
            
            # Accumulate sum and sum of squares for each parcel
            sum_values[:num_rows, :] += data
            sum_sq_values[:num_rows, :] += data ** 2
            count[:num_rows, :] += 1  # Increment count only for present rows

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")
            continue

    # Safely compute mean and variance per parcel across all data files
    mean_values = np.divide(sum_values, count, out=np.zeros_like(sum_values), where=(count != 0))
    variance_values = np.divide(sum_sq_values, count, out=np.zeros_like(sum_sq_values), where=(count != 0)) - (mean_values**2)

    # Ensure non-negative variances (avoid small numerical errors)
    variance_values = np.maximum(variance_values, 0)

    # Save results to processed_data/ folder
    mean_file = os.path.join(output_dir, "means_per_parcel.dat")
    var_file = os.path.join(output_dir, "variances_per_parcel.dat")

    # Prepare header with space-separated variable names
    header = " ".join(variable_names)

    np.savetxt(mean_file, mean_values, header=header, fmt="%15.8e")
    np.savetxt(var_file, variance_values, header=header, fmt="%15.8e")

    print(f"Mean and variance per parcel saved in: {output_dir}/")

if __name__ == "__main__":
    main()

