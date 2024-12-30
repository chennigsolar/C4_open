# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import os


def get_delta_theta_x(m):
    """
    Calculating the critical temperature rise of soil acc. to DIN VDE 0276-1000

    Parameters
    ----------
    m : float
        load factor of the 24 h profile (average load)

    Returns
    -------
    delta_theta_x:  float
        critical temperature rise of soil

    References
    -----------
    [5] DIN VDE 0276 Starkstromkabel –
        Teil 1000: Strombelastbarkeit, Allgemeines, Umrechnungsfaktoren
    """
    delta_theta_x = -33.3 * m + 48.31
    return delta_theta_x


class LoadProfile:
    """
    The LoadProfile class provides a simple possibility to calculate and access load profile parameters calculated from
    24 h load profile values in hourly resolution read from an .xlsx. The columns comprising the load values must be
    labeled with 'load'.

    Parameters
    ----------
    path_or_list : basestring or list
        path to .xlsx with load values or list

    Raises
    ------
    KeyError
        If load profile values cannot be found in .xlsx
    """
    def __init__(self, path_or_list):
        self.path_or_list = path_or_list
        self.load_list = []
        self.load_profile = None
        self.mu = None
        self.load_factor = None
        self.delta_theta_x = None
        self.Y_0 = None
        self.Y_1 = None
        self.Y_2 = None
        self.Y_3 = None
        self.Y_4 = None
        self.Y_5 = None
        self.Y_0_hour = None
        self.load_data()
        self.analyze()

    def load_data(self):
        if isinstance(self.path_or_list, list):
            self.load_list = self.path_or_list

        elif os.path.isfile(self.path_or_list):
            load_profile_raw = pd.read_excel(self.path_or_list)
            load_profile_raw.columns = [col.lower() for col in load_profile_raw.columns.tolist()]

            try:
                self.load_list = load_profile_raw["load"].tolist()

            except KeyError:
                raise KeyError("No load profile found")

    def analyze(self):
        # Finding maximal load
        load_max = max(self.load_list)

        # Normalize load list
        normalized_load_list = []
        for i in range(0, 24):
            normalized_load = self.load_list[i]/load_max
            normalized_load_list.append(normalized_load)

        # Calculating load factor
        load_sum = 0
        for i in range(0, 24):
            load = normalized_load_list[i]
            load_sum = load_sum + load

        self.load_factor = load_sum/len(normalized_load_list)

        # Calculating Y-Values
        Y_value_list = []
        for i in range(0, 24):
            Y_value = normalized_load_list[i] * normalized_load_list[i]
            Y_value_list.append(Y_value)

        # Calculating mu factor
        Y_i_sum = 0
        for i in range(0, 24):
            Y_i = Y_value_list[i]
            Y_i_sum = Y_i_sum + Y_i

        self.mu = Y_i_sum/len(Y_value_list)

        # Find Y_values with unity load and write this Y-Values to a candidate list
        # First: Convert Y-value list to numpy array
        Y_value_array = np.array(Y_value_list)

        # Second: Find Y_values with unity load and write this Y-Values to a candidate array
        Y_0_candidate_hour_array = np.array(np.where(Y_value_array == 1))

        # Third: Convert array to list
        Y_0_candidate_hour_list = Y_0_candidate_hour_array.tolist()
        Y_0_candidate_hour_list = Y_0_candidate_hour_list[0]

        # finding Y_value candidate with the highest sum of itself plus the preceding five Y-values
        # First: Summing up Y-values for candidates
        Y_sum_list = []
        for i in Y_0_candidate_hour_list:
            Y_sum = 0
            for j in range(6):
                Y_value = Y_value_list[i-j]
                Y_sum = Y_sum + Y_value

            Y_sum_list.append(Y_sum)

        # Second: Identifying highest Y_sum
        Y_sum_max = max(Y_sum_list)

        # Third: Finding the hour corresponding to the highest Y-sum
        Y_sum_max_hour = Y_sum_list.index(Y_sum_max)

        self.Y_0_hour = Y_0_candidate_hour_list[Y_sum_max_hour]

        # Getting the Y_i values
        self.Y_0 = Y_value_list[self.Y_0_hour]
        self.Y_1 = Y_value_list[self.Y_0_hour-1]
        self.Y_2 = Y_value_list[self.Y_0_hour-2]
        self.Y_3 = Y_value_list[self.Y_0_hour-3]
        self.Y_4 = Y_value_list[self.Y_0_hour-4]
        self.Y_5 = Y_value_list[self.Y_0_hour-5]

        self.delta_theta_x = get_delta_theta_x(self.load_factor)

        # Summarize results in DataFrame
        self.load_profile = pd.DataFrame({'normalized_load': normalized_load_list,
                                          'Y_i': Y_value_list
                                          }
                                         )

    def get_results(self):
        return self.load_profile

    def get_parameters(self):
        return {'Y_0': self.Y_0,
                'Y_1': self.Y_1,
                'Y_2': self.Y_2,
                'Y_3': self.Y_3,
                'Y_4': self.Y_4,
                'Y_5': self.Y_5,
                'mu': self.mu,
                'm': self.load_factor,
                'deltatheta_x': self.delta_theta_x,
                'Y_0_hour': self.Y_0_hour,
                }

