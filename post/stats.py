import numpy as np
import glob
import os

# Function to read parameters from parameters.dat
def read_parameters_from_dat(file_path):
    params = {}

    with open(file_path, "r") as f:
        for line in f:
            key, value = line.split(maxsplit=1)
            params[key] = value.strip()

    # Convert to correct types
    n_var = int(params["nVar"])
    nLevels = int(params["nLevels"])
    nx = 2**(nLevels - 1)  # Compute number of parcels

    variable_names = params["varName"].split()  # Assuming space-separated variable names
    # variable_names = params["varName"].split(",")  # If comma-separated, uncomment this

    return n_var, nx, variable_names

#------------------------------------------------------------------------------

def main():
    # Define the folder for saving processed data
    output_dir = "processed_data"
    os.makedirs(output_dir, exist_ok=True)  #  Create the folder if it doesn't exist

    # Read parameters from `parameters.dat`
    params_file_path = "parameters.dat"
    n_var, nx, variable_names = read_parameters_from_dat(params_file_path)

    # Find all data files
    data_files = sorted(glob.glob("../data/rlz_00001/Data_*.dat"))  # Ensure correct step ordering
    n_files = len(data_files)  # Number of data files (equivalent to steps)

    if n_files == 0:
        raise FileNotFoundError("No data files found in ../data/")

    print(f" Found {n_files} data files for processing.")

    # Initialize arrays for mean and variance calculations
    sum_values = np.zeros((nx, n_var))  # Sum of values for each parcel
    sum_sq_values = np.zeros((nx, n_var))  # Sum of squared values for variance calculation
    count = np.zeros((nx, n_var))  # Count occurrences for each parcel

    # Process each data file
    for file_path in data_files:
        print(f"Processing file: {file_path}")

        # Load data (each row corresponds to a parcel, each column is a variable)
        data = np.loadtxt(file_path)

        # Accumulate sum and sum of squares for each parcel
        sum_values += data
        sum_sq_values += data**2
        count += 1  # Each parcel is counted once per file

    # Compute mean and variance per parcel across all data files
    mean_values = sum_values / count
    variance_values = (sum_sq_values / count) - (mean_values**2)

    # Ensure non-negative variances (avoid small numerical errors)
    variance_values = np.maximum(variance_values, 0)

    # Save results to processed_data/ folder
    mean_file = os.path.join(output_dir, "means_per_parcel.dat")
    var_file = os.path.join(output_dir, "variances_per_parcel.dat")

    header = " ".join(variable_names)

    np.savetxt(mean_file, mean_values, header=header, fmt="%15.8e")
    np.savetxt(var_file, variance_values, header=header, fmt="%15.8e")

    print(f" Mean and variance per parcel saved in: {output_dir}/")

if __name__ == "__main__":
    main()

