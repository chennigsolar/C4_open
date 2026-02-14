from pathlib import Path
import os

GENERAL_SETTINGS = {'rho_T4_dryoutzone': 2.5}

# Set thermal resistivity of the dryout zone
rho_T4_dryoutzone = GENERAL_SETTINGS['rho_T4_dryoutzone']
rho_T4_dryoutzone = float(rho_T4_dryoutzone)

# Get the project directory
current_file_path = Path(__file__).resolve()
project_directory = current_file_path.parent.parent.parent

# Get path to cable datat base
# cable_database_path = os.path.join(current_file_path.parent, 'data', 'cable_data_.xlsx')
cable_database_path = os.path.join(current_file_path.parent, 'data', 'cable_data.db')