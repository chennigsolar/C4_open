# Introduction
C⁴ open is a command-line interface program for the current carrying capacity calculation of buried electrical cables 
acc. to IEC 60287 and IEC 60853-1/3. 

It is possible to calculate the current carrying capacities of cable arrangements for static and cyclic loads. 
In the latter case, a factor M according to IEC 60853-1/3 is calculated and applied to the previously calculated 
current carrying capacity.

In case of groups of cables, C⁴ open provides different possibilities to calculate the mutual heating factor F. 
For a single bundle of cables of two, three and four, a simple calculation tool can be used. 
For more parallel bundles of cables or groups of single cables, the arrangements can be imported and analyzed 
from .dxf files. Another way to calculate the factor F is the use of the 'create_arrangement' function from the module
'arrangements' Please see the separate jupyter notebook for further explanations.

C⁴ open comprises a data set of electric low voltage and medium voltage cables 
according to the German VDE standards, which are collected in a 'cable_data.db'. It can be found in the 
directory './data'.

# Pre-requisites
C⁴ open is written in Python 3.12. It uses the following libraries:
- numpy
- pandas
- matplotlib
- ezdxf
- openpyxl
- xlrd
- scipy

# Installation
You can install C⁴ open using pip:

```shell
pip install git+https://github.com/chennigsolar/C4_open.git
```

# Use
The following code snippet shows a simple example for calculating the current-carrying capacity of a system with three single-core cables under static load:

```python
import pandas as pd
from c4_open.project import Project
from c4_open.cables import Cable

# Simple example for calculating the current-carrying capacity 
# of a system with three single-core cable under static load

project_parameters = {'name': 'simple example',
                      'calc_case': 'ac_sc',
                      'F': 1,               # Mutual heating factor F
                      'L': 0.7,             # Depth of laying
                      'N': 3,               # Total number of parallel cables (or pipes)
                      'U_0': 12000.0,       # Phase-ground voltage (only for ac cables)
                      'U_n': 20000.0,       # Phase-phase voltage (only for ac cables)
                      'f': 50.0,            # Mains frequency
                      'cable_type': 'NA2XS2Y_1x150_20kV',  # Cable type according to 'cable_data.xlsx'
                      'deltatheta_x': 15.0, # Critical temperature rise of soil above drying of soil occurs 
                      'rho_T4': 1.0,        # Thermal resistivity of soil
                      'theta_amb': 20.0,    # Ambient temperature          
                      }

project = Project(**project_parameters)
cable = Cable(project, external_resistance_method='three_trefoil')

print(pd.Series(cable.get_result()))
```
Please see the example jupyter notebooks on the GitHub repository for further explanations how to use C⁴ open.

# License
This project is licensed under the GNU GPLv3 License - see the LICENSE.md file for details.

# Contribution
Please feel free to contribute to this project.