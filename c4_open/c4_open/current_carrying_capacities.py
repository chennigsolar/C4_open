"""
This module contains functions to calculate the current carrying capacity of cables according to IEC 60287-1-1 and
IEC 60287-2-1.
"""

import math
import c4_open.config


def get_maximum_operating_temperature_resistance(R_20, alpha, theta_max):
    """
    Calculate the conductor resistance at maximum operating temperature.

    Parameters
    ----------
    R_20 : float
        conductor resistance at 20°c  in Ohm/m

    alpha : float
        temperature coefficient of conductor material in 1/K

    theta_max : float
        maximum operating temperature of conductor in °C

    Returns
    -------
    R : float
        conductor resistance at maximum operating temperature
    """
    R = R_20 * (1 + alpha * (theta_max - 20))
    return R


def get_skin_effect_factor(f, R_):
    """
    Calculate the skin effect factor acc. to IEC 60287-1-1 (Chapter 2.1.2)

    Parameters
    ----------
    f : float
        mains frequency in Hz

    R_ : float
        dc resistance of the conductor in Ohm/m

    Returns
    -------
    Y_s : float
        skin effect factor

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
        Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
        General
    """
    # In the following the variable x_s represents (x_s)**2 as used in the IEC 60287 for simplification
    # Value for k_s is taken from table 2 of IEC 60287-1-1 and is applicable for round, stranded and sector-shaped (not
    # dried and impregnated) conductors
    k_s = 1.0
    x_s = 8.0 * math.pi * f / R_ * 10 ** (-7) * k_s
    Y_s = None
    if math.sqrt(x_s) <= 2.8:  # only valid for x_s <= 2.8
        Y_s = x_s ** 2 / (192 + 0.8 * x_s ** 2)
    return Y_s


def get_proximity_effect_factor_round(f, R_, d_c, s, cond_material):
    """
    Calculate the proximity-effect factor for three-core cable or three single core cables acc. to IEC 60287-1-1
    (Chapter 2.1.4.1).

    Parameters
    ----------
    f : float
        mains frequency in Hz

    d_c : float
        diameter of conductor in mm

    s : float
        distance between conductor axes in mm

    R_ : float
        dc resistance of the conductor in Ohm/m at maximum operating temperature

    cond_material : str
        conductor material; 'Al' or 'Cu'

    Returns
    -------
    Y_p : float
        proximity effect factor

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
        Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
        General
    """
    # In the following the variable x_p represents (x_p)**2 used in the IEC 60287 for simplification
    # Value for k_p is taken from table 2 of IEC 60287-1-1
    if cond_material == 'Al':
        k_p = 0.8
    else:
        k_p = 1.0
    x_p = 8 * math.pi * f / R_ * 10 ** (-7) * k_p
    Y_p = None
    if math.sqrt(x_p) <= 2.8:  # only valid for x_s <= 2.8
        Y_p = (x_p ** 2 / (192 + 0.8 * x_p ** 2) * (d_c / s) ** 2 *
               (0.312 * (d_c / s) ** 2 + 1.18 / (x_p ** 2 /
                                                 (192 + 0.8 * x_p ** 2) + 0.27)))
    return Y_p


def get_proximity_effect_factor_sector_shaped(f, R_, A_c, t_ins):
    """
    Calculate the proximity-effect factor for a multicore cable with sector-shaped conductors
    acc. to IEC 60287-1-1 (Chapter 2.1.4.2).

    Parameters
    ----------
    f : float
        mains frequency in Hz

    R_ : float
        dc resistance of the conductor in Ohm/m at maximum operating temperature

    A_c : float
        Nominal cross-section of conductor in mm²

    t_ins : float
        Thickness of insulation in mm

    Returns
    -------
    Y_p : float
        proximity effect factor

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
        Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
        General
    """
    # In the following the variable x_p represents (x_p)**2 used in the IEC 60287 for simplification
    # Value for k_p is taken from table 2 of IEC 60287-1-1
    # The degree of compaction (deg_comp) of stranded conductors is taken from DIN VDE 0176-1000
    k_p = 1
    deg_comp = 0.915
    d_x = math.sqrt(4 * A_c / (math.pi * deg_comp))
    t = 2 * t_ins  # thickness of insulation between conductors is 2 * t_ins
    s = d_x + t
    x_p = 8 * math.pi * f / R_ * 10 ** (-7) * k_p
    Y_p = None
    if math.sqrt(x_p) <= 2.8:  # only valid for x_s <= 2.8
        Y_p = 2 / 3 * (x_p ** 2 / (192 + 0.8 * x_p ** 2) * (d_x / s) ** 2 *
                       (0.312 * (d_x / s) ** 2 + 1.18 / (x_p ** 2 /
                                                         (192 + 0.8 * x_p ** 2) + 0.27)))
    return Y_p


def get_ac_resistance(Y_s, Y_p, R_):
    """
    Calculate the AC resistance of the conductor acc. to IEC 60287-1-1
    (Chapter 2.1).

    Parameters
    ----------
    Y_p : float
        proximity effect factor

    Y_s : float
        skin effect factor

    R_ : float
        dc resistance of the conductor in Ohm/m at maximum operating temperature

    Returns
    -------
    R : float
        ac resistance of the conductor in Ohm/m at maximum operating temperature

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
        Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
        General
    """
    R = R_ * (1 + Y_s + Y_p)
    return R


def get_dielectric_loss(C, U_0, f):
    """
    Calculate the AC resistance of the conductor acc. to IEC 60287-1-1
    (Chapter 2.1).

    Parameters
    ----------
    C : float
        capacity of the cable in F/m

    U_0 : float
        phase-to-ground voltage in V

    f : float
        mains frequency in Hz

    Returns
    -------
    W_d : float
        dielectric loss

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
        Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
        General
    """
    # Value for tan_delta is taken from table 3 of IEC 60287-1-1
    tan_delta = 0.004
    W_d = 2 * math.pi * f * C * U_0 ** 2 * tan_delta
    return W_d


def get_loss_factor_trefoil(R, R_s, s, d, f):
    """
    Calculate the loss factor for sheath and screen for three single core cables in trefoil acc. to IEC 60287-1-1
    (Chapter 2.3.1).

    This function applies to cables single-core cables with sheaths bonded at both ends. Eddy current losses are
    ignored.


    Parameters
    ----------
    R : float
        ac resistance of the conductor at max. operating temperature in Ohm/m

    R_s : float
        resistance of the screen winding at 20°C in Ohm/m

    s : float
        distance between conductor axes in mm

    d : float
        mean diameter of the sheath in mm

    f : float
        mains frequency in Hz

    Returns
    -------
    (X, lambda_1):  float, float
        reactance
        loss factor for sheath and screen

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
        Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
        General
    """
    X = 2 * 2 * math.pi * f * 10 ** (-7) * math.log(2 * s / d)
    lambda_1 = (R_s / R) * 1 / (1 + (R_s / X) ** 2)
    return lambda_1, X


def get_loss_factor_three_flat(R, R_s, s, d, f):
    """
    Calculate the loss factor for sheath and screen for three single core cables in flat formation
    acc. to IEC 60287-1-1 (Chapter 2.3.2).

    This function applies to cables single-core cables with sheaths bonded at both ends. Eddy current losses are
    ignored.


    Parameters
    ----------
    R : float
        ac resistance of the conductor at max. operating temperature in Ohm/m

    R_s : float
        resistance of the screen winding at 20°C in Ohm/m

    s : float
        distance between conductor axes in mm

    d : float
        mean diameter of the sheath in mm

    f : float
        mains frequency in Hz

    Returns
    -------
    (X, lambda_1):  float, float
        reactance
        loss factor for sheath and screen

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
        Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
        General
    """
    X_1 = 2 * 2 * math.pi * f * 10 ** (-7) * math.log(2 * 2 ** (1. / 3.) * (s / d))
    lambda_1 = (R_s / R) * 1 / (1 + (R_s / X_1) ** 2)
    return lambda_1, X_1


def get_t1_single_core(rho_T, d_c, t_1):
    """
    Calculate the thermal resistance between conductor and sheath T_1 for single-core cables acc. to IEC 60287-2-1
    (Chapter 2.1.1.1).

    This function applies to cables single-core cables with s

    Parameters
    ----------
    rho_T : float
         thermal resistivity of insulation in Km/W

    d_c : float
        diameter of the conductor in mm

    t_1 : float
        thickness of insulation between conductor and sheath in mm

    Returns
    -------
    T_1:  float
        thermal resistance between conductor and sheath

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    T_1 = rho_T / (2 * math.pi) * math.log(1 + (2 * t_1 / d_c))
    return T_1


def get_t1_belted_sector_shaped(rho_T, t_1, A_c, d_a, r_1):
    """
    Calculate the thermal resistance between conductor and sheath T_1 for three-core belted cables with
    sector-shaped conductors acc. to IEC 60287-2-1 (Chapter 2.1.1.2.5).

    This function applies to three-core belted cables with sector-shaped conductors. According to note in IEC 60287-1-1,
    Chapter 1.4.1.1 it can also be used for four-core cables where the fourth conductor is an unloaded neutral or
    protective conductor.

    Parameters
    ----------
    rho_T : float
         thermal resistivity of insulation in Km/W

    t_1 : float
        thickness of insulation between conductor and sheath in mm

    A_c : float
        Nominal cross-section of conductor in mm²

    d_a : float
        external diameter of the belt insulation in mm

    r_1 : float
        radius of the circle circumscribing the conductors in mm

    Returns
    -------
    T_1:  float
        thermal resistance between conductor and sheath

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    # In this calculation insulation and filling material are considered to be from the same material.
    # The degree of compaction (deg_comp) of stranded conductors is taken from DIN VDE 0176-1000
    t = 2 * t_1  # Insulation thickness between conductors is 2 * insulation thickness
    deg_comp = 0.915
    d_x = math.sqrt(4 * A_c / (math.pi * deg_comp))
    F_2 = 1 + 3 * t / (2 * math.pi * (d_x + t) - t)
    G = 3 * F_2 * math.log(d_a / (2 * r_1))
    T_1 = rho_T / (2 * math.pi) * G
    return T_1


def get_t3(rho_T, D_a, t_3):
    """
    Calculate the thermal resistance T_3 of the outer covering acc. to IEC 60287-2-1
    (Chapter 2.1.3).

    Parameters
    ----------
    rho_T : float
         thermal resistivity of insulation in Km/W

    D_a : float
        external diameter of the armour in mm
        Note: For unarmoured cables D_a is taken as the external diameter of the component immediately beneath it, i.e.
        sheath, screen or bedding.

    t_3 : float
        thickness of serving (e.g. sheath) in mm

    Returns
    -------
    T_3:  float
        thermal resistance of outer covering

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    T_3 = rho_T / (2 * math.pi) * math.log(1 + (2 * t_3 / D_a))
    return T_3


def get_t4_from_mutual_heating_factor(rho_T, L, D_e, F):
    """
    Calculate the external thermal resistance T_4 acc. to IEC 60287-2-1
    (Chapter 2.2.3.2).

    Calculation for groups of buried cables (not touching).

    Parameters
    ----------
    rho_T : float
         thermal resistivity of soil in Km/W

    L : float
         distance from the surface of the ground to the cable axis in mm

    D_e : float
        external diameter of the cable in mm

    F : float
        mutual heating factor

    Returns
    -------
    T_4:  float
        external thermal resistance of cable

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    u = 2 * L / D_e
    T_4 = rho_T / (2 * math.pi) * math.log((u + math.sqrt(u ** 2 - 1)) * F)
    return T_4


def get_t4_for_three_trefoil(rho_T, L, D_e):
    """
    Calculate the external thermal resistance T_4 of three single core cable with non-metallic sheathes in trefoil
    acc. to IEC 60287-2-1 (Chapter 2.2.4.3.3)

    Calculation for three cables (touching).

    Parameters
    ----------
    rho_T : float
         thermal resistivity of soil in Km/W

    L : float
         distance from the surface of the ground to the cable axis in mm

    D_e : float
        external diameter of the cable in mm

    Returns
    -------
    T_4:  float
        external thermal resistance of cable

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    u = 2 * L / D_e
    T_4 = None
    if u >= 5:
        T_4 = rho_T / (2 * math.pi) * (math.log(2 * u) + 2 * math.log(u))
    return T_4


def get_t4_for_three_flat(rho_T, L, D_e):
    """
    Calculate the external thermal resistance T_4 of three single core cable with
    non-metallic sheathes in flat formation acc. to IEC 60287-2-1 (Chapter 2.2.4.2.2)

    Calculation for three cables (touching).

    Parameters
    ----------
    rho_T : float
         thermal resistivity of soil in Km/W

    L : float
         distance from the surface of the ground to the cable axis in mm

    D_e : float
        external diameter of the cable in mm

    Returns
    -------
    T_4:  float
        external thermal resistance of cable

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    u = 2 * L / D_e
    T_4 = None
    if u >= 5:
        T_4 = rho_T * (0.475 * math.log(2 * u) - 0.142)
    return T_4


def get_t4_for_two_flat(rho_T, L, D_e):
    """
    Calculate the external thermal resistance T_4 of two single core cable with
    non-metallic sheathes in flat formation acc. to IEC 60287-2-1 (Chapter 2.2.4.1.2)

    Calculation for two cables (touching).

    Parameters
    ----------
    rho_T : float
         thermal resistivity of soil in Km/W

    L : float
         distance from the surface of the ground to the cable axis in mm

    D_e : float
        external diameter of the cable in mm

    Returns
    -------
    T_4:  float
        external thermal resistance of cable

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    u = 2 * L / D_e
    T_4 = None
    if u >= 5:
        T_4 = rho_T / math.pi * (math.log(2 * u) - 0.295)
    return T_4


def get_thermal_resistance_between_cable_and_pipe(K, D_e, theta_mean):
    """
    Calculate the thermal resistance between cable and plastic pipe acc. to IEC 60287-2-1
    (Chapter 2.2.7.1).

    Calculation for cables in plastic pipes.

    Parameters
    ----------
    K : int
         number of cables in pipe

    D_e : float
        external diameter of one cable in mm

    theta_mean : float
        mean temperature of the medium filling the space between cable and pipe (estimation)

    Returns
    -------
    T_41:  float
        thermal resistance between cable and pipe

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    # for plastic pipes following values are taken from IEC 60287-2-1, table 4
    U = 1.87
    V = 0.312
    Y = 0.0037

    if K == 1:
        D_e = D_e

    elif K == 2:
        D_e = 1.65 * D_e

    elif K == 3:
        D_e = 2.15 * D_e

    elif K == 4:
        D_e = 2.50 * D_e

    else:
        raise ValueError('Invalid number of cables!')

    T_41 = U / (1 + 0.1 * (V + Y * theta_mean) * D_e)
    return T_41


def get_thermal_resistance_of_pipe(rho_T, D_o, D_d):
    """
    Calculate the thermal resistance of the pipe acc. to IEC 60287-2-1
    (Chapter 2.2.7.2).

    Parameters
    ----------
    rho_T : float
         thermal resistivity of pipe material in Km/W

    D_o : float
        external diameter of the pipe mm

    D_d : float
        internal diameter of the pipe mm

    Returns
    -------
    T_42:  float
        thermal resistance of the pipe

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    T_42 = rho_T / (2 * math.pi) * math.log(D_o / D_d)
    return T_42


def get_t4_pipe(T_41, T_42, rho_T, L, D_o, F):
    """
    Calculate the external thermal resistance T_4 of a pipe acc. to IEC 60287-2-1
    (Chapter 2.2.3.2 and 2.2.7).

    Calculation for buried pipe(s) (not touching).

    Parameters
    ----------
    T_41 : float
         thermal resistance between cable and pipe in Km/W

    T_42 : float
          thermal resistance of the pipe in Km/W

    rho_T : float
         thermal resistivity of soil in Km/W

    L : float
         distance from the surface of the ground to the pipe axis in mm

    D_o : float
        external diameter of the pipe in mm

    F : float
        mutual heating factor

    Returns
    -------
    T_4:  float
        external thermal resistance of pipe

    References
    -----------
    [2] Electric cables – Calculation of the current rating –
        Part 2-1: Thermal resistance – Calculation of thermal resistance
    """
    T_43 = get_t4_from_mutual_heating_factor(rho_T, L, D_o, F)
    T_4 = T_41 + T_42 + T_43
    return T_4


def get_current_carrying_capacity_dc(theta_max,
                                     theta_amb,
                                     delta_theta_x,
                                     rho_T4,
                                     R_,
                                     n,
                                     T_1,
                                     T_2,
                                     T_3,
                                     T_4,
                                     dry_zone=True):
    """
    Calculate the permissible current rating of a dc cable acc. to IEC 60287-1-1
    (Chapter 1.4.1.2 and 1.4.2.2).

    Calculation for cables where partial drying-out occurs or not.

    Parameters
    ----------
    theta_max : float
        maximum operating temperature of conductor in °C

    theta_amb : float
        ambient temperature in °C

    delta_theta_x : float
        critical temperature rise of soil in K

    rho_T4 : float
         thermal resistivity of soil in Km/W

    R_ : float
        dc resistance of the conductor at max. operating temperature in Ohm/m

    n : int
        number of load-carrying conductors in the cable

    T_1 : float
        thermal resistance per unit length between one conductor and sheath in Km/W

    T_2 : float
        thermal resistance per unit length of the bedding between sheath and armour in Km/W

    T_3 : float
        thermal resistance per unit length of the external serving of the cable in K.m/W

    T_4 : float
        thermal resistance per unit length between the cable surface and the surrounding
        medium in Km/W

    dry_zone : bool
         consider dry-zone around cable; default is 'True'

    Returns
    -------
    I:  float
        permissible current rating

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
    Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
    General
    """

    if dry_zone:
        I = math.sqrt(((theta_max - theta_amb) + (c4_open.config.rho_T4_dryoutzone / rho_T4 - 1) * delta_theta_x
                       ) / (R_ * (T_1 + n * T_2 + n * (T_3 + c4_open.config.rho_T4_dryoutzone / rho_T4 * T_4))))

    else:
        I = math.sqrt((theta_max - theta_amb) / (R_ * T_1 + n * R_ * T_2 + n * R_ * (T_3 + T_4)))

    return I


def get_current_carrying_capacity_ac(theta_max,
                                     theta_amb,
                                     delta_theta_x,
                                     rho_T4,
                                     R,
                                     W_d,
                                     lambda_1,
                                     lambda_2,
                                     n,
                                     T_1,
                                     T_2,
                                     T_3,
                                     T_4,
                                     dry_zone=True):
    """
    Calculate the permissible current rating of an ac cable acc. to IEC 60287-1-1
    (Chapter 1.4.1.1 and 1.4.2.1).

    Calculation for cables where partial drying-out occurs or not.

    Parameters
    ----------
    theta_max : float
        maximum operating temperature of conductor in °C

    theta_amb : float
        ambient temperature in °C

    delta_theta_x : float
        critical temperature rise of soil in K

    rho_T4 : float
         thermal resistivity of soil in Km/W

    R : float
        ac resistance of the conductor at max. operating temperature in Ohm/m

    W_d : float
        dielectric loss in W/m

    lambda_1 : float
        loss factor for sheath and screen

    lambda_2 : float
        loss factor for armouring

    n : int
        number of load-carrying conductors in the cable

    T_1 : float
        thermal resistance per unit length between one conductor and sheath in Km/W

    T_2 : float
        thermal resistance per unit length of the bedding between sheath and armour in Km/W

    T_3 : float
        thermal resistance per unit length of the external serving of the cable in K.m/W

    T_4 : float
        thermal resistance per unit length between the cable surface and the surrounding
        medium in Km/W

    dry_zone : bool
         consider dry-zone around cable; default is 'True'

    Returns
    -------
    I:  float
        permissible current rating

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
    Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
    General
    """
    if dry_zone:
        I = math.sqrt(((theta_max - theta_amb) - W_d * (0.5 * T_1 + n * (
                T_2 + T_3 + c4_open.config.rho_T4_dryoutzone / rho_T4 * T_4)) + (
                               c4_open.config.rho_T4_dryoutzone / rho_T4 - 1) * delta_theta_x) / (
                R * (T_1 + n * (1 + lambda_1) * T_2 + n * (1 + lambda_1 + lambda_2) * (
                              T_3 + c4_open.config.rho_T4_dryoutzone / rho_T4 * T_4))))

    else:
        I = math.sqrt(((theta_max - theta_amb) - W_d * (0.5 * T_1 + n * (T_2 + T_3 + T_4))) / (
                R * T_1 + n * R * (1 + lambda_1) * T_2 + n * R * (1 + lambda_1 + lambda_2) * (T_3 + T_4)))

    return I


def get_current_carrying_capacity_dc_in_pipe(theta_max,
                                             theta_amb,
                                             delta_theta_x,
                                             rho_T4,
                                             R_,
                                             n,
                                             K,
                                             T_1,
                                             T_2,
                                             T_3,
                                             T_4,
                                             dry_zone=True):
    """
    Calculate the permissible current rating of K dc cables in pipe acc. to IEC 60287-1-1
    (Chapter 1.4.1.2 and 1.4.2.2) and adjustment of T_3,

    Calculation for cables where partial drying-out occurs or not.

    Parameters
    ----------
    theta_max : float
        maximum operating temperature of conductor in °C

    theta_amb : float
        ambient temperature in °C

    delta_theta_x : float
        critical temperature rise of soil in K

    rho_T4 : float
         thermal resistivity of soil in Km/W

    R_ : float
        dc resistance of the conductor at max. operating temperature in Ohm/m

    n : int
        number of load-carrying conductors in the cable

    K : int
        number of cables in pipe

    T_1 : float
        thermal resistance per unit length between one conductor and sheath in Km/W

    T_2 : float
        thermal resistance per unit length of the bedding between sheath and armour in Km/W

    T_3 : float
        thermal resistance per unit length of the external serving of the cable in K.m/W

    T_4 : float
        thermal resistance per unit length between the cable surface and the surrounding
        medium in Km/W

    dry_zone : bool
         consider dry-zone around cable; default is 'True'

    Returns
    -------
    I:  float
        permissible current rating

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
    Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
    General
    """
    # For K multicore cables with n cores per cable in one pipe, the total number of cores is n*K,
    # but the value for T_3 must be divided by K. For single core cables applies the same but n is equal to 1.
    if dry_zone:
        I = math.sqrt(((theta_max - theta_amb) + (c4_open.config.rho_T4_dryoutzone / rho_T4 - 1) * delta_theta_x
                       ) / (R_ * (T_1 + n * T_2 + n * K * (T_3 / K + c4_open.config.rho_T4_dryoutzone / rho_T4 * T_4))))

    else:
        I = math.sqrt((theta_max - theta_amb) / (R_ * T_1 + n * R_ * T_2 + n * K * R_ * (T_3 / K + T_4)))

    return I


def get_current_carrying_capacity_ac_single_core_in_pipe(theta_max,
                                                         theta_amb,
                                                         delta_theta_x,
                                                         rho_T4,
                                                         R,
                                                         W_d,
                                                         lambda_1,
                                                         K,
                                                         T_1,
                                                         T_3,
                                                         T_4,
                                                         dry_zone=True):
    """
    Calculate the permissible current rating of unarmoured ac cables laid in pipe acc. to IEC 60287-1-1
    (Chapter 1.4.1.1 and 1.4.2.1).

    Calculation for cables where partial drying-out occurs or not.

    Parameters
    ----------
    theta_max : float
        maximum operating temperature of conductor in °C

    theta_amb : float
        ambient temperature in °C

    delta_theta_x : float
        critical temperature rise of soil in K

    rho_T4 : float
         thermal resistivity of soil in Km/W

    R : float
        ac resistance of the conductor at max. operating temperature in Ohm/m

    W_d : float
        dielectric loss in W/m

    lambda_1 : float
        loss factor for sheath and screen

    K : int
        number of cables in pipe

    T_1 : float
        thermal resistance per unit length between one conductor and sheath in Km/W

    T_3 : float
        thermal resistance per unit length of the external serving of the cable in K.m/W

    T_4 : float
        thermal resistance per unit length between the cable surface and the surrounding
        medium in Km/W

    dry_zone : bool
         consider dry-zone around cable; default is 'True'

    Returns
    -------
    I:  float
        permissible current rating

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
    Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
    General
    """
    # For K single cables the value for T_3 must be divided by K.
    if dry_zone:
        I = math.sqrt(((theta_max - theta_amb) - W_d * (0.5 * T_1 + K * (
                1 / K * T_3 + c4_open.config.rho_T4_dryoutzone / rho_T4 * T_4)) + (
                               c4_open.config.rho_T4_dryoutzone / rho_T4 - 1) * delta_theta_x) / (
                              R * (T_1 + K * (1 + lambda_1) * (1 / K * T_3 + c4_open.config.rho_T4_dryoutzone / rho_T4 * T_4))))

    else:
        I = math.sqrt(((theta_max - theta_amb) - W_d * (0.5 * T_1 + K * (1 / K * T_3 + rho_T4 / rho_T4 * T_4))
            + (rho_T4 / rho_T4 - 1) * delta_theta_x) / (R * (T_1 + K * (
            1 + lambda_1) * (1 / K * T_3 + rho_T4 / rho_T4 * T_4))))

    return I


def get_current_carrying_capacity_ac_multi_core_in_pipe(theta_max,
                                                        theta_amb,
                                                        delta_theta_x,
                                                        rho_T4,
                                                        R,
                                                        W_d,
                                                        n,
                                                        K,
                                                        T_1,
                                                        T_3,
                                                        T_4,
                                                        dry_zone=True):
    """
    Calculate the permissible current rating of unscreened ac multicore cables laid in pipe acc. to IEC 60287-1-1
    (Chapter 1.4.1.1 and 1.4.2.1).

    Calculation for cables where partial drying-out occurs or not.

    Parameters
    ----------
    theta_max : float
        maximum operating temperature of conductor in °C

    theta_amb : float
        ambient temperature in °C

    delta_theta_x : float
        critical temperature rise of soil in K

    rho_T4 : float
         thermal resistivity of soil in Km/W

    R : float
        ac resistance of the conductor at max. operating temperature in Ohm/m

    W_d : float
        dielectric loss in W/m

    n : int
        number of load-carrying conductors in the cable

    K : int
        number of cables in pipe

    T_1 : float
        thermal resistance per unit length between one conductor and sheath in Km/W

    T_3 : float
        thermal resistance per unit length of the external serving of the cable in K.m/W

    T_4 : float
        thermal resistance per unit length between the cable surface and the surrounding
        medium in Km/W

    dry_zone : bool
         consider dry-zone around cable; default is 'True'

    Returns
    -------
    I:  float
        permissible current rating

    References
    -----------
    [1] Electric cables – Calculation of the current rating –
    Part 1-1: Current rating equations (100 % load factor) and calculation of losses –
    General
    """
    # For K multicore cables with n cores per cable in one pipe, the total number of cores is n*K,
    # but the value for T_3 must be divided by K.
    if dry_zone:
        I = math.sqrt(((theta_max - theta_amb) - W_d *
                       (0.5 * T_1 + n * K * (1 / K * T_3 + c4_open.config.rho_T4_dryoutzone
                                             / rho_T4 * T_4)) + (c4_open.config.rho_T4_dryoutzone
                                                                 / rho_T4 - 1) * delta_theta_x)
                      / (R * T_1 + n * K * R * (1 / K * T_3 + c4_open.config.rho_T4_dryoutzone / rho_T4 * T_4)))

    else:
        I = math.sqrt(((theta_max - theta_amb) - W_d * (0.5 * T_1 + n * K * (1 / K * T_3 + rho_T4 / rho_T4 * T_4))
                       + (rho_T4 / rho_T4 - 1) * delta_theta_x) / (R * T_1 + n * K * R *
                                                                   (1 / K * T_3 + rho_T4 / rho_T4 * T_4)))
    return I

