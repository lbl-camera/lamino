import click
import numpy as np
import tifffile
import toml
from pathlib import Path


def read_ovf_file(file_path):
    """
    Reads an OVF file and returns the data as a numpy array.
    
    Parameters:
        file_path (str or Path): Path to the OVF file.

    Returns:
        np.ndarray: The data from the OVF file as a numpy array.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # get the header information
    header = {}
    for line in lines:
        if line.startswith("#"):
            key_value = line[1:].strip().split(":", 1)
            if len(key_value) == 2:
                key, value = key_value
                header[key.strip()] = value.strip()
        else:
            break
    xdims = int(header.get("xnodes", 0))
    ydims = int(header.get("ynodes", 0))
    zdims = int(header.get("znodes", 0))
    
    size = xdims * ydims * zdims
    if size == 0:
        raise ValueError("The OVF file does not contain valid dimensions in the header.")

    data_start_index = None
    for i, line in enumerate(lines):
        if 'Begin: Data Text' in line:
            data_start_index = i + 1
            break
    
    if data_start_index is None:
        raise ValueError("The OVF file does not contain a '# Data' section.")
    
    # allocate a numpy array to hold the data
    vector_data = np.zeros((size, 3), dtype=np.float32)

    # Read the data into a numpy array
    data_lines = lines[data_start_index:]
    data = []
    idx = 0 
    for line in data_lines:
        if line.strip() and not line.startswith("#"):
            vector_data[idx] = np.fromstring(line.strip(), sep=' ')
            idx += 1
    
    xcomp = vector_data[:, 0].reshape((xdims, ydims, zdims))
    ycomp = vector_data[:, 1].reshape((xdims, ydims, zdims))
    zcomp = vector_data[:, 2].reshape((xdims, ydims, zdims))
    return xcomp, ycomp, zcomp

def save_as_tiff(data, output_path):
    """
    Saves a numpy array as a TIFF file.
    
    Parameters:
        data (np.ndarray): The data to save.
        output_path (str or Path): Path to the output TIFF file.
    """
    if len(data) != 3:
        raise ValueError("Data must be a list or tuple of three numpy arrays (xcomp, ycomp, zcomp).")

    if len(data) != len(output_path):
        raise ValueError("The number of data arrays and output paths must be the same.")

    for i in range(3):
        vdata = np.transpose(data[i], (2, 1, 0))
        tifffile.imwrite(output_path[i], vdata.astype(np.float32))

@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output_path', type=click.Path(), default=None, help='Path to the output TIFF file.')
def main(input_file, output_path):
    """
    Command-line interface for converting OVF files to TIFF format.
    
    Parameters:
        input_file (str): Path to the input OVF file.
        output_path (str): Path to the output TIFF files. If not provided, it will be generated based on the input file name.
    """
    data = read_ovf_file(input_file)
    
    # if output_file is not provided, save in current directory with the same name as input file but with .tiff extension
    if output_path is None:
        output_path = Path.cwd()
        base_name = Path(input_file).stem
    else:
        output_path = Path(output_path)
        base_name = Path(input_file).stem
    output_file = [output_path / f"{base_name}{i}.tiff" for i in range(3)]
    
    save_as_tiff(data, output_file)
    print(f"Converted {input_file} to {output_file}")

if __name__ == "__main__":
    main()
