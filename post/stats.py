import numpy as np
import glob
from utils import read_parameters  # Import the utility function for reading parameters

#------------------------------------------------------------------------------

def main():
    # Reading input parameters
    yaml_file_path = '../data/InputParameters.yaml'
    params = read_parameters(yaml_file_path)

    n_var = params['params'].get('nVar')  # Number of variables to process
    variable_names = params['variable_names']  # Get variable names from input parameters

    # Ensure the number of variable names matches the number of variables (n_var)
    if len(variable_names) != n_var:
        raise ValueError("The number of variable names does not match the number of variables (nVar)")

    # Find all data files
    data_files = glob.glob("../data/Data_*.dat")

    # Initialize arrays to hold mean and variance data
    mean_values = np.zeros((len(data_files), n_var))
    variance_values = np.zeros((len(data_files), n_var))

    # Process each column in the data files
    for column_index in range(n_var):
        all_data = np.array([])

        # Load and append data from each file
        for file_path in data_files:
            data = np.loadtxt(file_path)
            all_data = np.append(all_data, data[:, column_index])

        # Calculate mean and variance
        mean_values[:, column_index] = np.mean(all_data)
        variance_values[:, column_index] = np.var(all_data)

    # Create aligned header strings for variable names
    header_format = "{:<15}" * n_var
    aligned_header = header_format.format(*variable_names)

    # Save mean and variance to files with aligned variable names in the header
    np.savetxt("mean_values.dat", mean_values, header=aligned_header, fmt="%15.8e")
    np.savetxt("variance_values.dat", variance_values, header=aligned_header, fmt="%15.8e")

if __name__ == "__main__":
    main()

