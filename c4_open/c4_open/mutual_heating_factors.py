import math
import pandas as pd
import ezdxf

"""
This module provides utility functions for the calculation of the mutual heating factor F. The formulas for the 
different arrangements are derived from the formula in IEC 60287-2-1, Chapter 2.2.3.2 for groups of buried equally 
loaded identical cables (not touching).
"""


def get_mutual_heating_factor_for_three_cables_trefoil(L, D_e):
    """
    Calculate the mutual heating factor for three cables in trefoil

    Parameters
    ----------
    L : float
         distance from the surface of the ground to the cable axis in the same units as D_e (e.g. m)

    D_e : float
        external diameter of the cable in the same units as L (e.g. m)

    Returns
    -------
    F:  float
        mutual heating factor
    """

    F = (2 * L / D_e) ** 2
    return F


def get_mutual_heating_factor_for_three_cables_flat(L, D_e):
    """
    Calculate the mutual heating factor for three cables in flat formation with no spacing

    Parameters
    ----------
    L : float
         distance from the surface of the ground to the cable axis in the same units as D_e (e.g. m)

    D_e : float
        external diameter of the cable in the same units as L (e.g. m)

    Returns
    -------
    F:  float
        mutual heating factor
    """
    F = (2 * L / D_e) ** 2
    return F


def get_mutual_heating_factor_for_two_cables_flat(L, D_e):
    """
    Calculate the mutual heating factor for two cables in flat formation with no spacing

    Parameters
    ----------
    L : float
         distance from the surface of the ground to the cable axis in the same units as D_e (e.g. m)

    D_e : float
        external diameter of the cable in the same units as L (e.g. m)

    Returns
    -------
    F:  float
        mutual heating factor
    """
    F = 2 * L / D_e
    return F


def get_mutual_heating_factor_for_four_cables_square(L, D_e):
    """
    Calculate the mutual heating factor for four cables in square arrangement

    Parameters
    ----------
    L : float
         distance from the surface of the ground to the cable axis in the same units as D_e (e.g. m)

    D_e : float
        external diameter of the cable in the same units as L (e.g. m)

    Returns
    -------
    F:  float
        mutual heating factor
    """
    F = (2 * L / D_e) ** 2 * 2 * L / (math.sqrt(2) * D_e)
    return F


def read_dxf(path):
    doc = ezdxf.readfile(path)
    msp = doc.modelspace()

    x_coordinates = []
    y_coordinates = []
    colors = []

    for elements in msp:
        if elements.dxftype() == 'CIRCLE':
            x_coordinates.append(elements.dxf.center.x)
            y_coordinates.append(elements.dxf.center.y)
            colors.append(elements.dxf.color)

    coordinates = pd.DataFrame({'X': x_coordinates, 'Y': y_coordinates, 'Color': colors})
    return coordinates


def get_mutual_heating_factor_from_dxf(path):
    """
    Calculate the mutual heating factor for a group of cables from .dxf

    The .dxf must meet the following requirements:
    - the drawing must be in units equal to m
    - the external diameters of the cables must be drawn as circles in m
    - the soil surface must be at Y = 0 and the cables must be positioned
      below according to the depth of laying in m
    - The observed cable must have the color 'red'

    Parameters
    ----------
    path: str
         path to.dxf

    Returns
    -------
    F:  float
        mutual heating factor
    """
    # Get coordinates and coors of circles from dxf
    coordinates = read_dxf(path)
    # Generate data frame with observed cable
    coordinates_observed_cable = pd.DataFrame(coordinates)

    parallel_cables = coordinates_observed_cable[coordinates_observed_cable['Color'] != (1)].index
    coordinates_observed_cable.drop(parallel_cables, inplace=True)  # Drop parallel cables from data frame
    coordinates_observed_cable.reset_index(inplace=True, drop=True)  # Generate new index

    # Generate data frame with parallel cables
    coordinates_parallel_cables = pd.DataFrame(coordinates)

    observed_cable = coordinates_parallel_cables[coordinates_parallel_cables['Color'] == (1)].index
    coordinates_parallel_cables.drop(observed_cable, inplace=True)  # Remove (drop) observed cable from data frame
    coordinates_parallel_cables.reset_index(inplace=True, drop=True)  # Generate new index

    m = len(coordinates_parallel_cables)  # Number of rows excl. Header

    # variable for d_p is d_p0 and variable for d_p' is d_p1
    f = 1
    n = 0
    while n < m:
        d_p0 = math.sqrt((coordinates_parallel_cables.loc[n].X
                          - coordinates_observed_cable.loc[0].X) ** 2
                         + (coordinates_parallel_cables.loc[n].Y
                            - coordinates_observed_cable.loc[0].Y) ** 2
                         )

        d_p1 = math.sqrt((coordinates_parallel_cables.loc[n].X
                          - coordinates_observed_cable.loc[0].X) ** 2
                         + (-coordinates_parallel_cables.loc[n].Y
                            - coordinates_observed_cable.loc[0].Y) ** 2
                         )

        f_n = d_p1 / d_p0

        f = f * f_n
        n += 1
    return f
