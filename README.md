Moment Maps Generation with CASA for ALMA Data

## Overview

This repository provides scripts and documentation for generating moment maps (Moment 0, 1, and 2) for various molecular lines (CO J=0-1, 2-1, 3-2, 4-3, HCO+, CI, C+, H2O) from ALMA data, specifically for the Circum-Nuclear Disk (CND) scale. We utilize the Common Astronomy Software Applications (CASA) package for data processing.

## Prerequisites

- **CASA**: Ensure you have CASA installed on your system. You can download it from the [CASA website](https://casa.nrao.edu/).
- **Python**: Python 3.x is recommended for running the scripts.
- **ALMA Data**: Obtain the ALMA datasets you intend to analyze.

## Installation

1. Clone this repository to your local machine:
   ```bash
   git clone https://github.com/yourusername/repo-name.git
   cd repo-name
   ```

2. Ensure all dependencies are installed. This includes any necessary Python packages.

## Usage

### Generating Moment Maps

1. **Prepare your ALMA data**: Ensure your data is calibrated and in the correct format for CASA.

2. **Run the Moment Map Scripts**: Use the provided scripts to generate moment maps. The scripts are organized by molecular line.

   - To generate moment maps for CO J=0-1, execute:
     ```bash
     casa --nologger --nogui -c co_j01_moments.py
     ```

   - Replace `co_j01_moments.py` with the appropriate script for other molecular lines (e.g., `co_j21_moments.py`, `hcop_moments.py`, etc.).

3. **Adjust Parameters**: Modify the scripts to adjust parameters such as the moment order or any specific settings related to your data.

### Output

The scripts will generate moment maps and save them in the specified output directory. You will find:

- **Moment 0**: Integrated intensity map.
- **Moment 1**: Velocity field map.
- **Moment 2**: Velocity dispersion map.

## Example

An example script for CO J=2-1 might look like this:

```python
# co_j21_moments.py
from tasks import *
from casatools import image

# Load your ALMA data
vis = 'your_alma_data.ms'
output_dir = 'output/'

# Generate Moment 0
immoments(vis=vis, outfile=output_dir + 'co_j21_mom0.fits', moments=[0])

# Generate Moment 1
immoments(vis=vis, outfile=output_dir + 'co_j21_mom1.fits', moments=[1])

# Generate Moment 2
immoments(vis=vis, outfile=output_dir + 'co_j21_mom2.fits', moments=[2])
```

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please contact Mojtaba Raouf(mailto:raouf@strw.leidenuniv.nl).

---

Feel free to adapt this README to suit your specific project details and requirements!
