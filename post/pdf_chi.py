

import numpy as np
import glob
import matplotlib.pyplot as plt
from utils import read_parameters, compute_pdf  # Import the utility functions
#------------------------------------------------------------------------------

def main():
    # Reading input parameters
    yaml_file_path = '../data/InputParameters.yaml'
    params = read_parameters(yaml_file_path)

    length_scale = params['params'].get('LengthScale', 0.03125)  # Default to 0.03125 if not found
    n_var = params['params'].get('nVar')  # Default to None if not found

    # Find all data files
    data_files = glob.glob("../data/Data_*.dat")

    plt.figure()

    # Process each column in the data files
    for column_index in range(n_var):
        # Initialize arrays to hold mixed fraction data
        mixed_fraction = np.array([])

        # Load and append data from each file
        for file_path in data_files:
            data = np.loadtxt(file_path)
            mixed_fraction = np.append(mixed_fraction, data[:, 3])

        # Calculate the average values
        average_values = (mixed_fraction[::2] + mixed_fraction[1::2]) / 2

        # Calculate scalar dissipation rate (chi)
        chi_values = np.zeros(len(average_values))
        for i in range(0, len(chi_values) - 1, 2):
            chi_values[i] = 2 * 1 * ((average_values[i + 1] - average_values[i]) / length_scale) ** 2

        # Filter and log-transform chi values
        positive_chi = chi_values[chi_values > 0]
        log_chi = np.log10(positive_chi)

        # Compute PDF for log-transformed chi values
        log_chi_min = -10
        log_chi_max = np.max(log_chi)
        pdf, bin_centers = compute_pdf(log_chi, log_chi_min, log_chi_max, nbins=100)

        # Plot the PDF
        plt.plot(bin_centers, pdf, label=f'Column {column_index + 1}')

    # Add labels and title
    plt.xlabel('log10(chi)')
    plt.ylabel('Probability Density')
    plt.title('PDF of log10(chi) for All Columns')
    plt.legend()
    plt.savefig("chi.pdf")

if __name__ == "__main__":
    main()

