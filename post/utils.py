# utils.py

import numpy as np
import yaml

#-------------------------Function to read the parameters ------------------------

def read_parameters(yaml_file_path):
    """Read parameters from a YAML file."""
    with open(yaml_file_path, 'r') as file:
        parameters = yaml.safe_load(file)
    return parameters

#---------------------Function to calculate probability density function-----------

def compute_pdf(x, xmin, xmax, nbins=100):
    """Compute the probability density function (PDF) of the data."""
    dxbin = (xmax - xmin) / nbins
    xbins = np.linspace(xmin + dxbin / 2, xmax - dxbin / 2, num=nbins)

    ibin = np.floor((x - xmin) / dxbin)
    np.clip(ibin, 0, nbins - 1, out=ibin)

    P = np.zeros(nbins)

    for i in range(nbins):
        P[i] = np.size(np.where(ibin == i))

    P = P / len(x) / dxbin

    return P, xbins

