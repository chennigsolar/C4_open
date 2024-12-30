"""
This module provides a function to create a cable arrangement and calculate a corresponding
mutual heating factor F. It bases on the formula given in chapter 2.2.3.2 of IEC 60287-2-1.
"""

import math

# The following functions are helper functions to create first cable systems for different
# arrangements. All functions return a tuple with x, y coordinates of the cables per system
def _trefoil_system(d_ext, depth_of_laying):
    return[(0, - depth_of_laying - math.tan(30 * math.pi / 180) * d_ext / 2), 
           (- d_ext, - depth_of_laying - math.tan(30 * math.pi / 180) * d_ext / 2),
           (- 0.5 * d_ext, - depth_of_laying + d_ext / (2 * math.cos(30 * math.pi / 180)))
           ]

def _three_flat_system(d_ext, depth_of_laying):
     return[(0, - depth_of_laying), 
           (d_ext, - depth_of_laying),
           (-d_ext, - depth_of_laying)
           ]

def _two_flat_system(d_ext, depth_of_laying):
     return[(0, - depth_of_laying), 
           (d_ext, - depth_of_laying),
           ]

def _single_cable(d_ext, depth_of_laying):
     return[(0, - depth_of_laying), 
           ]

def _square_system(d_ext, depth_of_laying):
     return[(0, - depth_of_laying - d_ext / 2),
            (0, - depth_of_laying + d_ext / 2), 
           (d_ext, - depth_of_laying - d_ext / 2),
           (d_ext, - depth_of_laying + d_ext / 2)
           ]

def _translate_system(system, d):
    """
    Helper function for the translation of a cable system in x direction.  
    
    Parameters
    ----------
    system : tuple with float
        x, y coordinates of the original cable system in m

    d : float
        amount of translation in m

    Returns
    -------
    dist_ratio : tuple with float
        x, y coordinates of the translated cable system in m
    
    """

    translated_system = []

    if isinstance(system, list):
        for cable in system:
            translated_system.append((cable[0] + d, cable[1]))
    else:
        raise ValueError ('No system to translate!')
    
    return translated_system

def _get_distance_ratio(obs, par):
    """
    Helper function for the calculation of the ratio of the distance
    from the observed to the mirrored parallel cable to the distance
    from the observed to the parallel cable.

    Parameters
    ----------
    obs : tuple with float
        x, y coordinates of observed cable in m

    par : tuple with float
        x, y coordinates of parallel cable in m

    Returns
    -------
    dist_ratio : float
        ratio of the distance from observed cable to the mirrored cable to the
        distance of the original cable and add the result to a list

    """
    x_0, y_0 = obs
    x_i, y_i = par

    delta_x = abs(x_i - x_0)
    delta_y = abs(y_i - y_0)


    dist = (delta_x ** 2 + delta_y ** 2) ** 0.5

    y_i_ = (-1) * y_i
    delta_y_ = abs(y_i_ - y_0)

    dist_ = (delta_x ** 2 + delta_y_ ** 2) ** 0.5
    
    dist_ratio = dist_ / dist

    return dist_ratio


def create_arrangement(system_arrangement, number_of_systems, d_ext, d_clear, depth_of_laying):
    """
    Create an arrangement of cables. The mutual heating factor F is calculated using the formula
    in chapter 2.2.3.2 of IEC 60287-2-1.

    Parameters
    ----------
    system_arrangement : string
        possible arrangements: 'three_flat', 'two_flat', 'trefoil', 'single', 'square'

    number_of_systems : int
        number of systems

    d_ext : float
        external diameter of cable in m

    d_clear : float
        clear distance between cable systems in m

    depth_of_laying : Float
        depth of laying in m

    Returns
    -------
    F : float
        mutual heating factor

    systems : list with float
        list with coordinates of all cable systems

    xx : list with float
        list with x coordinates of all cables

    yy : list with float
        list with y coordinates of all cables

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    systems = []

    if system_arrangement == "trefoil":
        systems.append(_trefoil_system(d_ext, depth_of_laying))
        d = 2 * d_ext + d_clear
    
    elif system_arrangement == "three_flat":
        systems.append(_three_flat_system(d_ext, depth_of_laying))
        d = 3 * d_ext + d_clear

    elif system_arrangement == "two_flat":
        systems.append(_two_flat_system(d_ext, depth_of_laying))
        d = 2 * d_ext + d_clear

    elif system_arrangement == "single":
        systems.append(_single_cable(d_ext, depth_of_laying))
        d = d_ext + d_clear

    elif system_arrangement == "square":
        systems.append(_square_system(d_ext, depth_of_laying))
        d = 2 * d_ext + d_clear

    else:
        raise ValueError("Invalid system arrangement!")
    
    if isinstance(number_of_systems, int):
        for i in range(2, number_of_systems + 1):
            if i%2 ==0:
                systems.append(_translate_system(systems[0], d * (i / 2)))
            else:
                systems.append(_translate_system(systems[0], -d * (i - 1) / 2))
    else:
        raise ValueError('Number of systems must be integer!')
    
    # Create lists with all x and y coordinates of cables
    xx = []
    yy = []
    for system in systems:
        for cable in system:
            xx.append(cable[0])
            yy.append(cable[1])
    
    
    # Calculate d_ as the ratio of the distance from observed cable to the mirrored parallel cable to the 
    # distance of the original parallel cable and add the result to a list.

    # For the case of three_flat systems and an even number of systems, the observed cable must be the second cable of the first system (the right one). This is done by
    # switching the coordinates of the first and second cable in the lists.

    if number_of_systems%2 == 0 and system_arrangement == 'three_flat':
        xx[0], xx[1] = xx[1], xx[0]
        yy[0], yy[1] = yy[1], yy[0]

    dd_ = []
    for i in range(1, len(xx)):  # Start with the second coordinates set as the first refers to the observed cable, which must not be considered
        d_ = _get_distance_ratio((xx[0], yy[0]), (xx[i], yy[i]))
        dd_.append(d_)

    # Calculate the mutual heating factor F as the product of all distance ratios

    F = 1
    for i in range(len(dd_)):
        F = F * dd_[i]

    return F, systems, xx, yy