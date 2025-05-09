
import numpy as np
import glob
import matplotlib
matplotlib.use('TkAgg')  # Use a standard GUI backend
import matplotlib.pyplot as plt

# Utility function to read parameters from a .dat file
def read_parameters(file_path):
    params = {}
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                key, value = parts
                params[key] = value.strip()

    # Convert to correct types
    domain_length = float(params.get('domainLength', 0.03125))
    n_var = int(params.get('nVar', 1))
    variable_names = params.get('varName', '').split()
    i_batchelor = list(map(int, params.get('i_batchelor', '').split()))

    return domain_length, n_var, variable_names, i_batchelor

# Utility function to compute PDF
def compute_pdf(data, data_min, data_max, nbins=100):
    bins = np.linspace(data_min, data_max, nbins)
    hist, bin_edges = np.histogram(data, bins=bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return hist, bin_centers

#------------------------------------------------------------------------------

def main():
    # Reading input parameters from .dat file
    dat_file_path = 'parameters.dat'
    domain_length, n_var, variable_names, i_batchelor = read_parameters(dat_file_path)

    # Find all data files
    data_files = sorted(glob.glob("../data/rlz_00001/Data_*.dat"))

    if len(data_files) == 0:
        print("Error: No data files found.")
        return

    # Process each variable in the data files
    for column_index in range(n_var):
        # Calculate the length scale for the current variable
        length_scale = domain_length / (2 ** i_batchelor[column_index])

        # Initialize arrays to hold mixed fraction data
        mixed_fraction = np.array([])

        # Load and append data from each file
        for file_path in data_files:
            try:
                data = np.loadtxt(file_path)
                if data.shape[1] <= column_index:
                    print(f"Warning: Column {column_index} not found in {file_path}")
                    continue
                mixed_fraction = np.append(mixed_fraction, data[:, column_index])
            except Exception as e:
                print(f"Error reading file {file_path}: {e}")
                continue

        # Calculate the average values between adjacent pairs
        if len(mixed_fraction) % 2 != 0:
            mixed_fraction = mixed_fraction[:-1]  # Ensure an even number of elements

        average_values = (mixed_fraction[::2] + mixed_fraction[1::2]) / 2

        # Calculate scalar dissipation rate (chi)
        chi_values = np.zeros(len(average_values) - 1)
        for i in range(len(chi_values)):
            delta = (average_values[i + 1] - average_values[i]) / length_scale
            chi_values[i] = 2 * (delta ** 2)

        # Filter positive chi values and apply log transformation
        positive_chi = chi_values[chi_values > 0]
        if len(positive_chi) == 0:
            print(f"No positive chi values for variable {variable_names[column_index]} (Column {column_index + 1})")
            continue

        log_chi = np.log10(positive_chi)

        # Compute PDF for log-transformed chi values
        log_chi_min = -10
        log_chi_max = np.max(log_chi)
        pdf, bin_centers = compute_pdf(log_chi, log_chi_min, log_chi_max, nbins=100)

        # Plot the PDF for this variable
        plt.figure(figsize=(8, 5))
        plt.plot(bin_centers, pdf, label=f'{variable_names[column_index]} (Column {column_index + 1})')
        plt.xlabel('log10(chi)')
        plt.ylabel('Probability Density')
        plt.title(f'PDF of log10(chi) for {variable_names[column_index]}')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()

        # Save the plot with a specific filename
        plot_filename = f"chi_{variable_names[column_index]}.pdf"
        plt.savefig(plot_filename)
        plt.close()
        print(f"Saved PDF plot for variable {variable_names[column_index]} as {plot_filename}")

if __name__ == "__main__":
    main()


