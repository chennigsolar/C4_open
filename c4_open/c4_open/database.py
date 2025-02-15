"""
This module contains functions to read the cable database and return it as a pandas dataframe.
It also contains a function to add new cables to the database.
"""

import pandas as pd
import sqlite3
from c4_open.config import cable_database_path


def get_cable_database():
    """
    This function reads the cable database returns it as a pandas
    dataframe.
    """
    conn = sqlite3.connect(cable_database_path)
    cable_database = pd.read_sql_query('SELECT * FROM cable_data', conn, index_col='Type')
    conn.close()
    return cable_database

def get_cable_types():
    """
    This function reads the available cable types from the cable database and returns it as a pandas
    series.
    """
    conn = sqlite3.connect(cable_database_path)
    cable_database = pd.read_sql_query('SELECT * FROM cable_data', conn, index_col='Type')
    conn.close()
    return pd.Series(cable_database.index)

def add_cable_to_database(data):
    """
    This function adds a new cable to the cable database.

    Parameters
    ----------
    data : dict
        dictionary containing the cable data
    """

    # Check if cable type already exists in the database
    existing_cable_types = get_cable_types().tolist()
    if data['Type'] not in existing_cable_types:

        # Connect to SQLite database
        conn = sqlite3.connect(cable_database_path)
        cursor = conn.cursor()

        # Create a tuple of values from data dictionary
        values = tuple(data.values())

        # Execute an INSERT statement
        placeholders = ', '.join(['?' for _ in data])
        values = tuple(data.values())
        cursor.execute(f"INSERT INTO cable_data VALUES ({placeholders})", values)
        conn.commit()
        conn.close()

    else:
        raise ValueError(f"Cable type '{data['Type']}' already exists in the database")

def get_cable_data(cable_type):
    """
    This function reads the cable data for a specific cable type from the database.
    """
    conn = sqlite3.connect(cable_database_path)
    cable_data = pd.read_sql_query(f"SELECT * FROM cable_data WHERE Type = '{cable_type}'", conn)
    conn.close()
    return cable_data

def create_database_from_xlsx(xlsx):
    """
    This function creates a new database from a xlsx file.
    """
    path = "./data/cable_data.db"
    conn = sqlite3.connect(path)
    cable_database = pd.read_excel(xlsx)
    cable_database.to_sql('cable_data', conn, if_exists='replace')
    conn.close()
    return

# create_database_from_xlsx('./data/cable_data.xlsx')