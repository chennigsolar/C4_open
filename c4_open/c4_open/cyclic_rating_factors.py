import scipy.special as sp
import math
import c4_open.tables, c4_open.config


def get_cyclic_rating_factor_for_group(W,
                                       theta_amb,
                                       theta_max,
                                       D_e,
                                       N,
                                       F,
                                       L,
                                       rho_T4,
                                       mu,
                                       Y_0,
                                       Y_1,
                                       Y_2,
                                       Y_3,
                                       Y_4,
                                       Y_5):
    """
    Calculate the cyclic rating factor M for a group of cables acc. to IEC 60853-1

    Parameters
    ----------
    W : float
        total joule (current dependent) losses per cable determined using conductor resistance at
        maximum operating temperature in W

    theta_max : float
        maximum operating temperature of conductor in °C

    theta_amb : float
        ambient temperature in °C

    D_e : float
        external diameter of one cable in mm

    N : int
        total number of cables

    F : float
        mutual heating factor

    L : float
        distance from the surface of the ground to the pipe axis in mm

    rho_T4 : float
         thermal resistivity of soil in Km/W

    mu : float
        loss load factor for load current cycle under consideration

    Y_i : float
        coefficient proportional to the current-dependent losses
        in a cable between (i) and (i + 1) hours prior to the instant
        of highest conductor temperature

    Returns
    -------
    M : float
        rating factor due to a cyclic variation of load current

    k_1 : float
        ratio of cable (or duct) external surface temperature rise in a group above ambient to conductor
        temperature rise above ambient under steady-state conditions

    gamma : float
        function used in calculating cyclic rating factors for one of a group of cables

    theta_R : float
        conductor temperature rise above ambient at time i hours due to the application of a step
        function of load current equal to the sustained (100% load factor) rated current

    References
    -----------
    [3] IEC 60853-1 - Calculation of the cyclic and emergency current rating of cables
        Part 1:
        Cyclic rating factor for cables up to and including 18/30 (36) kV
    """
    # Getting value for soil thermal diffusity (delta) from dictionary in tables py. acc. to table II of IEC 60853-1
    delta = c4_open.tables.soil_thermal_diffusivity_table.get(rho_T4)
    # Calculating maximal temperature above ambient
    theta_max_above_amb = theta_max - theta_amb

    """ Calculation of k_1 and d_f """
    T_4 = rho_T4 / (2 * math.pi) * math.log(
        4 * L / D_e)  # Calculation of T_4 acc. to IEC 60853-1

    deltaT_4 = rho_T4 * math.log(F) / (2 * math.pi)  # Calculation of deltaT_4 acc. to IEC 60853-1

    k_1 = W * (
            T_4 + deltaT_4) / theta_max_above_amb  # Calculation of k_1 acc. to IEC 60853-1

    d_f = 4 * L * 0.001 / (F ** (1 / (N - 1)))  # Calculation of d_f acc. to IEC 60853-1

    """ Calculation of gamma_i values """
    # For the sake of an easier verifiability the calculation of the gamma_i-values is done element wise and
    # not by iterating though a loop
    gamma_1 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 1 * delta)) + (N - 1)
                    * -sp.expi(- (d_f ** 2 / (16 * 3600 * 1 * delta)))) \
                   / (2 * math.log(4 * L * 0.001 * F / (D_e * 0.001)))

    gamma_2 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 2 * delta)) + (N - 1)
                    * -sp.expi(- (d_f ** 2 / (16 * 3600 * 2 * delta)))) \
                   / (2 * math.log(4 * L * 0.001 * F / (D_e * 0.001)))

    gamma_3 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 3 * delta)) + (N - 1)
                    * -sp.expi(-(d_f ** 2 / (16 * 3600 * 3 * delta)))) \
                   / (2 * math.log(4 * L * 0.001 * F / (D_e * 0.001)))

    gamma_4 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 4 * delta)) + (N - 1)
                    * -sp.expi(-(d_f ** 2 / (16 * 3600 * 4 * delta)))) \
                   / (2 * math.log(4 * L * 0.001 * F / (D_e * 0.001)))

    gamma_5 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 5 * delta)) + (N - 1)
                    * -sp.expi(-(d_f ** 2 / (16 * 3600 * 5 * delta)))) \
                   / (2 * math.log(4 * L * 0.001 * F / (D_e * 0.001)))

    gamma_6 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 6 * delta)) + (N - 1)
                    * -sp.expi(-(d_f ** 2 / (16 * 3600 * 6 * delta)))) \
                   / (2 * math.log(4 * L * 0.001 * F / (D_e * 0.001)))

    """ Calculation of theta_R (i)/theta_R (oo) values """
    # Note: In the following theta_R (i)/(theta_R (oo) is called theta_Ri for simplification
    theta_R0 = 0  # by definition
    theta_R1 = 1 - k_1 + k_1 * gamma_1
    theta_R2 = 1 - k_1 + k_1 * gamma_2
    theta_R3 = 1 - k_1 + k_1 * gamma_3
    theta_R4 = 1 - k_1 + k_1 * gamma_4
    theta_R5 = 1 - k_1 + k_1 * gamma_5
    theta_R6 = 1 - k_1 + k_1 * gamma_6

    """ Calculation of factor M """
    M = 1 / (

            Y_0 * (theta_R1 - theta_R0
                        ) +

            Y_1 * (theta_R2 - theta_R1
                        ) +

            Y_2 * (theta_R3 - theta_R2
                        ) +

            Y_3 * (theta_R4 - theta_R3
                        ) +

            Y_4 * (theta_R5 - theta_R4
                        ) +

            Y_5 * (theta_R6 - theta_R5
                        ) +

            mu * (1 - theta_R6
                       )

    ) ** 0.5

    return M, k_1, [gamma_1, gamma_2, gamma_3, gamma_4, gamma_5, gamma_6], \
        [theta_R0, theta_R1, theta_R2, theta_R3, theta_R4, theta_R5, theta_R6]


def get_cyclic_rating_factor_for_single_cable(W,
                                              theta_amb,
                                              theta_max,
                                              D_e,
                                              L,
                                              rho_T4,
                                              mu,
                                              Y_0,
                                              Y_1,
                                              Y_2,
                                              Y_3,
                                              Y_4,
                                              Y_5):
    """
    Calculate the cyclic rating factor M for a single cable acc. to IEC 60853-1

    Parameters
    ----------
    W : float
        total joule (current dependent) losses per cable determined using conductor resistance at
        maximum operating temperature in W

    theta_max : float
        maximum operating temperature of conductor in °C

    theta_amb : float
        ambient temperature in °C

    D_e : float
        external diameter of one cable in mm

    L : float
        distance from the surface of the ground to the pipe axis in mm

    rho_T4 : float
         thermal resistivity of soil in Km/W

    mu : float
        loss load factor for load current cycle under consideration

    Y_i : float
        coefficient proportional to the current-dependent losses
        in a cable between (i) and (i + 1) hours prior to the instant
        of highest conductor temperature

    Returns
    -------
    M : float
        rating factor due to a cyclic variation of load current

    k : float
        ratio of cable (or duct) external surface temperature rise above ambient to conductor
        temperature rise above ambient under steady-state conditions

    beta : float
        cable (or duct) outer surface temperature attainment factor at time (i) hours, beta(00) = 1

    phi : float
        function defined in Appendix A of IEC 60853-1

    References
    -----------
    [3] IEC 60853-1 - Calculation of the cyclic and emergency current rating of cables
        Part 1:
        Cyclic rating factor for cables up to and including 18/30 (36) kV
    """
    # Getting value for soil thermal diffusity (delta) from dictionary in tables py. acc. to table II of IEC 60853-1
    delta = c4_open.tables.soil_thermal_diffusivity_table.get(rho_T4)

    # Calculating maximal temperature above ambient
    theta_max_above_amb = theta_max - theta_amb

    """ Calculation of T_4 and k """
    T_4 = float(rho_T4 / (2 * math.pi) * math.log(4 * L * 0.001 / (D_e * 0.001)))
    # Calculation of T_4 acc. to IEC 60853-1

    k = float(W * T_4 / theta_max_above_amb)  # Calculation of k acc. to IEC 60853-1

    """Calculation of beta_i values """
    beta_1 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 1 * delta))) \
                  / (2 * math.log(4 * L * 0.001 / (D_e * 0.001)))

    beta_2 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 2 * delta))) \
                  / (2 * math.log(4 * L * 0.001 / (D_e * 0.001)))

    beta_3 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 3 * delta))) \
                  / (2 * math.log(4 * L * 0.001 / (D_e * 0.001)))

    beta_4 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 4 * delta))) \
                  / (2 * math.log(4 * L * 0.001 / (D_e * 0.001)))

    beta_5 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 5 * delta))) \
                  / (2 * math.log(4 * L * 0.001 / (D_e * 0.001)))

    beta_6 = (-sp.expi(-(D_e * 0.001) ** 2 / (16 * 3600 * 6 * delta))) \
                        / (2 * math.log(4 * L * 0.001 / (D_e * 0.001)))

    """ Calculation of phi_i values """
    phi_0 = beta_1

    phi_1 = beta_2 - beta_1

    phi_2 = beta_3 - beta_2

    phi_3 = beta_4 - beta_3

    phi_4 = beta_5 - beta_4

    phi_5 = beta_6 - beta_5

    """ Calculation of B """
    B = float(Y_0 * phi_0 + Y_1 * phi_1 + Y_2 * phi_2 + Y_3 * phi_3
                   + Y_4 * phi_4 + Y_5 * phi_5)

    """ Calculation of factor M """
    M = 1 / (((1 - k) * Y_0 + k * (B + mu * (1 - beta_6))) ** 0.5)
    return M, k, [beta_1, beta_2, beta_3, beta_4, beta_5, beta_6], [phi_0, phi_1, phi_2, phi_3, phi_4, phi_5]


def get_M_1(M, k, Y_0, theta_max, theta_amb, delta_theta_x, rho_T4):
    """
    Calculate the cyclic rating factor corrected for moisture migration M_1 acc. to IEC 60853-3

    Parameters
    ----------
    M : float
        rating factor due to a cyclic variation of load current

    k : float
        ratio of cable (or duct) external surface temperature rise above ambient to conductor
        temperature rise above ambient under steady-state conditions

    Y_0 : float
        coefficient proportional to the current-dependent losses
        in a cable between (0) and (0 + 1) hours prior to the instant
        of highest conductor temperature

    theta_max : float
        maximum operating temperature of conductor in °C

    theta_amb : float
        ambient temperature in °C

    delta_theta_x : float
        critical temperature rise of soil in K

    rho_T4 : float
         thermal resistivity of soil in Km/W

    Returns
    -------
    M_1 : float
        cyclic rating factor corrected for moisture migration
    theta_SPK : float
        peak cyclic temperature rise of the cable surface, in kelvins (K)

    References
    -----------
    [4] IEC 60853-3 - Calculation of the cyclic and emergency
        current rating of cables –
        Part 3:
        Cyclic rating factor for cables of all voltages, with partial drying
        of the soil
    """
    # The cyclic rating factor M is calculated for the assumption of a constant external thermal
    # resistance T_4 at different loads. A dry out zone around the cable was not considered.
    # In the following a method acc. to IEC 60853-3 is used to consider a potential dry out zone around the
    # cable and adjust the M accordingly. The result is an adjusted factor M_1

    # E = theta_SPK/theta_c and theta_c = theta_max_above_amb and theta_SPK is the peak cyclic temperature rise
    # of the cable surface

    E = 1 - M ** 2 * (1 - k) * Y_0  # formula (8) with simplification
    # acc. to page 19 of IEC 60853-3, because thermal capacities of the cables can be ignored.

    theta_SPK = E * (theta_max - theta_amb)

    # If (theta_SPK + W_d * T_4) is not greater than the critical temperature rise of the soil, then there will be no
    # drying and the cyclic rating factor is given by M without any correction and this function will return 'None';
    # if the peak temperature rise of the cable surface (theta_SPK + Wd * T_4 ) is greater than the critical
    # value, then the corrected value, M_1, will be calculated and returned
    # For all in this programm considered calculation cases, dielectric losses can be ignored and the check can be
    # simplified to:

    if not theta_SPK > delta_theta_x:
        M_1 = None

    else:
        M_1 = M * math.sqrt((1 + k * (c4_open.config.rho_T4_dryoutzone / rho_T4 - 1)) / (
                1 + E * (c4_open.config.rho_T4_dryoutzone / rho_T4 - 1)))  # formula (11)

    return M_1, theta_SPK


def get_joule_loss(R, n, I, W_d, lambda_1):
    """
    Calculate the total joule loss of a cable

    Parameters
    ----------
    R : float
        (ac) resistance of one conductor at max. operating temperature in Ohm/m

    n : int
        number of load-carrying conductors in the cable

    I : current per conductor in A

    W_d : float
        dielectric loss in W/m

    lambda_1 : float
        loss factor for sheath and screen

    Returns
    -------
    W : float
        total joule loss of the cable in W/m
    """
    W = I ** 2 * R * n * (1 + lambda_1) + W_d
    return W


def get_joule_loss_pipe(R, n, K, I, W_d, lambda_1):
    """
    Calculate the total joule loss of a cable

    Parameters
    ----------
    R : float
        (ac) resistance of one conductor at max. operating temperature in Ohm/m

    n : int
        number of load-carrying conductors in the cable

    K : int
        number of cable in pipe

    I : current per conductor in A

    W_d : float
        dielectric loss in W/m

    lambda_1 : float
        loss factor for sheath and screen

    Returns
    -------
    W : float
        total joule loss of the cable in W/m
    """
    W = (I ** 2 * R * n * (1 + lambda_1) + W_d) * K
    return W


class FactorM:
    """
        The FactorM class provides a convenient way to calculate the cyclic rating factor of buried cables with
        transient loads. It selects the applicable calculation functions depending on the given arrangement and number
        of parallel cables or pipes/conduits and calls them in the correct sequence.

        The method FactorM.get_result() returns a dict with the results.

        Please see the example jupyter notebooks for further explanation.

        Parameters
        ----------
        cable : Cable
            an instance of the Cable class

        load_profile : LoadProfile
            an instance of the LoadProfile class

        Raises
        ------
        ValueError
            If `calc_case` is unknown
        """

    def __init__(self, cable, load_profile):

        # Reading data from cable
        self.calc_case = str(cable.calc_case)
        self.N = int(cable.N)
        self.F = float(cable.factor_F)
        self.L = float(cable.L)
        self.rho_T4 = float(cable.rho_T4)
        self.theta_amb = float(cable.theta_amb)
        self.d_out = float(cable.d_out)
        if (self.calc_case == "dc_sc_pipe"
                or self.calc_case == "ac_sc_pipe"
                or self.calc_case == "ac_sc_pipe"
                or self.calc_case == "ac_mc_pipe"):
            self.d_out = float(cable.D_o)
        else:
            self.d_out = float(cable.d_out)
        self.n = int(cable.n)
        self.theta_max = float(cable.theta_max)

        self.I_dry_zone_yes = float(cable.I_dry_zone_yes)
        self.I_dry_zone_no = float(cable.I_dry_zone_no)
        if self.calc_case == 'dc_sc' or self.calc_case == 'dc_sc_pipe':
            self.R = float(cable.R_) # For DC cables, the relevant resistance is the DC resistance
        else:
            self.R = float(cable.R)

        if (self.calc_case == 'ac_sc' or self.calc_case == 'ac_mc'
                or self.calc_case == 'ac_sc_pipe' or self.calc_case == 'ac_mc_pipe'):
            self.W_d = float(cable.W_d)
        if self.calc_case == 'ac_sc' or self.calc_case == 'ac_sc_pipe':
            self.lambda_1 = float(cable.lambda_1)

        # Reading data from load profile
        self.mu = float(load_profile.mu)
        self.Y_0 = float(load_profile.Y_0)
        self.Y_1 = float(load_profile.Y_1)
        self.Y_2 = float(load_profile.Y_2)
        self.Y_3 = float(load_profile.Y_3)
        self.Y_4 = float(load_profile.Y_4)
        self.Y_5 = float(load_profile.Y_5)
        self.delta_theta_x = float(load_profile.delta_theta_x)

        # Reading calculation case dependent data from project parameter dictionary
        if self.calc_case == 'ac_sc_pipe' or self.calc_case == 'ac_mc_pipe' or self.calc_case == 'dc_sc_pipe':
            self.K = int(cable.K)
            self.D_o = float(cable.D_o)

        self.results = None
        self.W_dry_zone_yes = None
        self.W_dry_zone_no = None
        self.M_dry_zone_yes = None
        self.M_dry_zone_no = None
        self.k_1 = None
        self.k_1_dry_zone_no = None
        self.k = None
        self.k_dry_zone_no = None
        self.M_1 = None
        self.theta_SPK = None
        self.beta = None
        self.phi = None
        self.gamma = None
        self.gamma = None
        self.theta_R_dry_zone_yes = None
        self.theta_R_dry_zone_no = None
        self.B = None
        self.report_0 = None
        self.report_1 = None
        self.report_2 = None

        self.run()

    def run(self):
        # Calculate the total joule loss for the different calculation cases
        if self.calc_case == 'ac_sc':
            self.W_dry_zone_yes = get_joule_loss(self.R,
                                                 self.n,
                                                 self.I_dry_zone_yes,
                                                 self.W_d,
                                                 self.lambda_1)
            self.W_dry_zone_no = get_joule_loss(self.R,
                                                self.n,
                                                self.I_dry_zone_no,
                                                self.W_d,
                                                self.lambda_1)

        elif self.calc_case == 'ac_mc':
            self.W_dry_zone_yes = get_joule_loss(self.R,
                                                 self.n,
                                                 self.I_dry_zone_yes,
                                                 self.W_d,
                                                 0.0
                                                 )
            self.W_dry_zone_no = get_joule_loss(self.R,
                                                self.n,
                                                self.I_dry_zone_no,
                                                self.W_d,
                                                0.0
                                                )

        elif self.calc_case == 'dc_sc':
            self.W_dry_zone_yes = get_joule_loss(self.R,
                                                 self.n,
                                                 self.I_dry_zone_yes,
                                                 0.0,
                                                 0.0
                                                 )
            self.W_dry_zone_no = get_joule_loss(self.R,
                                                self.n,
                                                self.I_dry_zone_no,
                                                0.0,
                                                0.0
                                                )

        elif self.calc_case == 'ac_sc_pipe':
            self.W_dry_zone_yes = get_joule_loss_pipe(self.R,
                                                      self.n,
                                                      self.K,
                                                      self.I_dry_zone_yes,
                                                      self.W_d,
                                                      self.lambda_1
                                                      )
            self.W_dry_zone_no = get_joule_loss_pipe(self.R,
                                                     self.n,
                                                     self.K,
                                                     self.I_dry_zone_no,
                                                     self.W_d,
                                                     self.lambda_1
                                                     )

        elif self.calc_case == 'ac_mc_pipe':
            self.W_dry_zone_yes = get_joule_loss_pipe(self.R,
                                                      self.n,
                                                      self.K,
                                                      self.I_dry_zone_yes,
                                                      self.W_d,
                                                      0.0
                                                      )
            self.W_dry_zone_no = get_joule_loss_pipe(self.R,
                                                     self.n,
                                                     self.K,
                                                     self.I_dry_zone_no,
                                                     self.W_d,
                                                     0.0
                                                     )

        elif self.calc_case == 'dc_sc_pipe':
            self.W_dry_zone_yes = get_joule_loss_pipe(self.R,
                                                      self.n,
                                                      self.K,
                                                      self.I_dry_zone_yes,
                                                      0.0,
                                                      0.0
                                                      )
            self.W_dry_zone_no = get_joule_loss_pipe(self.R,
                                                     self.n,
                                                     self.K,
                                                     self.I_dry_zone_no,
                                                     0.0,
                                                     0.0
                                                     )

        else:
            raise ValueError('Invalid Calculation Case')

        # Calculate M for the different calculation cases and for one cable/pipe or group of cables/pipes
        if (self.N == 1 and self.calc_case == 'dc_sc'
                or self.N == 1 and self.calc_case == 'ac_sc'
                or self.N == 1 and self.calc_case == 'ac_mc'):
            (self.M_dry_zone_yes,
             self.k,
             self.beta,
             self.phi) = get_cyclic_rating_factor_for_single_cable(self.W_dry_zone_yes,
                                                                   self.theta_amb,
                                                                   self.theta_max,
                                                                   self.d_out,
                                                                   self.L * 1000,
                                                                   self.rho_T4,
                                                                   self.mu,
                                                                   self.Y_0,
                                                                   self.Y_1,
                                                                   self.Y_2,
                                                                   self.Y_3,
                                                                   self.Y_4,
                                                                   self.Y_5
                                                                   )

            (self.M_dry_zone_no,
             self.k_dry_zone_no,
             self.beta,
             self.phi) = get_cyclic_rating_factor_for_single_cable(self.W_dry_zone_no,
                                                                   self.theta_amb,
                                                                   self.theta_max,
                                                                   self.d_out,
                                                                   self.L * 1000,
                                                                   self.rho_T4,
                                                                   self.mu,
                                                                   self.Y_0,
                                                                   self.Y_1,
                                                                   self.Y_2,
                                                                   self.Y_3,
                                                                   self.Y_4,
                                                                   self.Y_5
                                                                   )

        elif self.N == 1 and (self.calc_case == 'dc_sc_pipe'
                              or self.calc_case == 'ac_sc_pipe'
                              or self.calc_case == 'ac_mc_pipe'):
            (self.M_dry_zone_yes,
             self.k,
             self.beta,
             self.phi) = get_cyclic_rating_factor_for_single_cable(self.W_dry_zone_yes,
                                                                   self.theta_amb,
                                                                   self.theta_max,
                                                                   self.D_o,
                                                                   self.L * 1000,
                                                                   self.rho_T4,
                                                                   self.mu,
                                                                   self.Y_0,
                                                                   self.Y_1,
                                                                   self.Y_2,
                                                                   self.Y_3,
                                                                   self.Y_4,
                                                                   self.Y_5
                                                                   )

            (self.M_dry_zone_no,
             self.k_dry_zone_no,
             self.beta,
             self.phi) = get_cyclic_rating_factor_for_single_cable(self.W_dry_zone_no,
                                                                   self.theta_amb,
                                                                   self.theta_max,
                                                                   self.D_o,
                                                                   self.L * 1000,
                                                                   self.rho_T4,
                                                                   self.mu,
                                                                   self.Y_0,
                                                                   self.Y_1,
                                                                   self.Y_2,
                                                                   self.Y_3,
                                                                   self.Y_4,
                                                                   self.Y_5
                                                                   )

        elif (self.N > 1 and self.calc_case == 'dc_sc'
                or self.N > 1 and self.calc_case == 'ac_sc'
                or self.N > 1 and self.calc_case == 'ac_mc'):
            (self.M_dry_zone_yes,
             self.k_1,
             self.gamma,
             self.theta_R_dry_zone_yes) = get_cyclic_rating_factor_for_group(self.W_dry_zone_yes,
                                                                             self.theta_amb,
                                                                             self.theta_max,
                                                                             self.d_out,
                                                                             self.N,
                                                                             self.F,
                                                                             self.L * 1000,
                                                                             self.rho_T4,
                                                                             self.mu,
                                                                             self.Y_0,
                                                                             self.Y_1,
                                                                             self.Y_2,
                                                                             self.Y_3,
                                                                             self.Y_4,
                                                                             self.Y_5
                                                                             )

            (self.M_dry_zone_no,
             self.k_1_dry_zone_no,
             self.gamma,
             self.theta_R_dry_zone_no) = get_cyclic_rating_factor_for_group(self.W_dry_zone_no,
                                                                            self.theta_amb,
                                                                            self.theta_max,
                                                                            self.d_out,
                                                                            self.N,
                                                                            self.F,
                                                                            self.L * 1000,
                                                                            self.rho_T4,
                                                                            self.mu,
                                                                            self.Y_0,
                                                                            self.Y_1,
                                                                            self.Y_2,
                                                                            self.Y_3,
                                                                            self.Y_4,
                                                                            self.Y_5
                                                                            )

        elif (self.N > 1 and self.calc_case == 'dc_sc_pipe'
                or self.N > 1 and self.calc_case == 'ac_sc_pipe'
                or self.N > 1 and self.calc_case == 'ac_mc_pipe'):
            (self.M_dry_zone_yes,
             self.k_1,
             self.gamma,
             self.theta_R_dry_zone_yes) = get_cyclic_rating_factor_for_group(self.W_dry_zone_yes,
                                                                             self.theta_amb,
                                                                             self.theta_max,
                                                                             self.D_o,
                                                                             self.N,
                                                                             self.F,
                                                                             self.L * 1000,  # function expects mm
                                                                             self.rho_T4,
                                                                             self.mu,
                                                                             self.Y_0,
                                                                             self.Y_1,
                                                                             self.Y_2,
                                                                             self.Y_3,
                                                                             self.Y_4,
                                                                             self.Y_5
                                                                             )

            (self.M_dry_zone_no,
             self.k_1_dry_zone_no,
             self.gamma,
             self.theta_R_dry_zone_no) = get_cyclic_rating_factor_for_group(self.W_dry_zone_no,
                                                                            self.theta_amb,
                                                                            self.theta_max,
                                                                            self.D_o,
                                                                            self.N,
                                                                            self.F,
                                                                            self.L * 1000,
                                                                            self.rho_T4,
                                                                            self.mu,
                                                                            self.Y_0,
                                                                            self.Y_1,
                                                                            self.Y_2,
                                                                            self.Y_3,
                                                                            self.Y_4,
                                                                            self.Y_5
                                                                            )

        # in some cases the calculated theta_SPK does not exceed delta_theta_x and get_M_1() correctly returns 'None'.
        # If this is the case, a test should be made to check whether a recalculation of M_1 with M_dry_zone_no would
        # lead to a theta_SPK exceeding delta_theta_x. If this test is positive the previously calculated M_1,
        # considering M_dry_zone_yes, should be taken. Its value is then calculated by bypassing the conditional
        # check in get_M_1() by setting delta_theta_x to = 0,

        if self.N == 1:
            self.M_1, self.theta_SPK = get_M_1(self.M_dry_zone_yes,
                                               self.k,
                                               self.Y_0,
                                               self.theta_max,
                                               self.theta_amb,
                                               self.delta_theta_x,
                                               self.rho_T4
                                               )
            if self.M_1 is None:
                M_1_test, _ = get_M_1(self.M_dry_zone_no,
                                      self.k_dry_zone_no,
                                      self.Y_0,
                                      self.theta_max,
                                      self.theta_amb,
                                      self.delta_theta_x,
                                      self.rho_T4
                                      )
                if M_1_test is not None:

                    self.M_1, self.theta_SPK = get_M_1(self.M_dry_zone_yes,
                                                       self.k,
                                                       self.Y_0,
                                                       self.theta_max,
                                                       self.theta_amb,
                                                       0,
                                                       self.rho_T4
                                                       )

        else:
            self.M_1, self.theta_SPK = get_M_1(self.M_dry_zone_yes,
                                               self.k_1,
                                               self.Y_0,
                                               self.theta_max,
                                               self.theta_amb,
                                               self.delta_theta_x,
                                               self.rho_T4
                                               )
            if self.M_1 is None:
                M_1_test, _ = get_M_1(self.M_dry_zone_no,
                                      self.k_1_dry_zone_no,
                                      self.Y_0,
                                      self.theta_max,
                                      self.theta_amb,
                                      self.delta_theta_x,
                                      self.rho_T4
                                      )
                if M_1_test is not None:
                    self.M_1, self.theta_SPK = get_M_1(self.M_dry_zone_yes,
                                                       self.k_1,
                                                       self.Y_0,
                                                       self.theta_max,
                                                       self.theta_amb,
                                                       0,
                                                       self.rho_T4
                                                       )

    def get_result(self):

        return {'M (with dry zone)': self.M_dry_zone_yes,
                'M_1': self.M_1,
                'M (no dry zone)': self.M_dry_zone_no,
                }

    def get_report(self):
        """
        This method returns a dictionary with the results of the calculation.
        """
        if self.N>1:
            if self.M_1 is not None:
                self.report_0 = {'W': self.W_dry_zone_yes
                                 }
            else:
                self.report_0 = {'W': self.W_dry_zone_no
                                 }

            self.report_1 = {'L': self.L,
                           'd_out': self.d_out,
                           'rho_T4': self.rho_T4,

                           'N': self.N,
                           'F': self.F,
                           'theta_amb': self.theta_amb,
                           'theta_max': self.theta_max,
                           'mu': self.mu,
                           'Y_0': self.Y_0,
                           'Y_1': self.Y_1,
                           'Y_2': self.Y_2,
                           'Y_3': self.Y_3,
                           'Y_4': self.Y_4,
                           'Y_5': self.Y_5,

                           'k_1': self.k_1,
                           'gamma_1': self.gamma[0],
                           'gamma_2': self.gamma[1],
                           'gamma_3': self.gamma[2],
                           'gamma_4': self.gamma[3],
                           'gamma_5': self.gamma[4],
                           'gamma_6': self.gamma[5]}

            if self.M_1 is not None:
                self.report_2 = {'theta_R0':  self.theta_R_dry_zone_yes[0],
                           'theta_R1':  self.theta_R_dry_zone_yes[1],
                           'theta_R2':  self.theta_R_dry_zone_yes[2],
                           'theta_R3':  self.theta_R_dry_zone_yes[3],
                           'theta_R4':  self.theta_R_dry_zone_yes[4],
                           'theta_R5':  self.theta_R_dry_zone_yes[5],
                           'theta_R6':  self.theta_R_dry_zone_yes[6],

                           'M':         self.M_dry_zone_yes,
                           'theta_SPK': self.theta_SPK,
                           'M_1':       self.M_1,
                           }
            else:
                self.report_2 = {'theta_R0': self.theta_R_dry_zone_no[0],
                               'theta_R1': self.theta_R_dry_zone_no[1],
                               'theta_R2': self.theta_R_dry_zone_no[2],
                               'theta_R3': self.theta_R_dry_zone_no[3],
                               'theta_R4': self.theta_R_dry_zone_no[4],
                               'theta_R5': self.theta_R_dry_zone_no[5],
                               'theta_R6': self.theta_R_dry_zone_no[6],

                               'M': self.M_dry_zone_no,
                               'theta_SPK': None,
                               'M_1': None,
                                 }
            #return {**self.report_1, **self.report_2}

        else:
            if self.M_1 is not None:
                self.report_0 = {'W': self.W_dry_zone_yes
                                 }
            else:
                self.report_0 = {'W': self.W_dry_zone_no
                                 }

            self.report_1 = {'L': self.L,
                           'd_out': self.d_out,
                           'rho_T4': self.rho_T4,
                           'W': self.W_dry_zone_yes,
                           'theta_amb': self.theta_amb,
                           'theta_max': self.theta_max,
                           'mu': self.mu,
                           'Y_0': self.Y_0,
                           'Y_1': self.Y_1,
                           'Y_2': self.Y_2,
                           'Y_3': self.Y_3,
                           'Y_4': self.Y_4,
                           'Y_5': self.Y_5,

                           'k': self.k,
                           'beta_1': self.beta[0],
                           'beta_2': self.beta[1],
                           'beta_3': self.beta[2],
                           'beta_4': self.beta[3],
                           'beta_5': self.beta[4],
                           'beta_6':self.beta[5],

                           'phi_0': self.phi[0],
                           'phi_1': self.phi[1],
                           'phi_2': self.phi[2],
                           'phi_3': self.phi[3],
                           'phi_4': self.phi[4],
                           'phi_5': self.phi[5],

                           'B': self.B
                             }
            if self.M_1 is not None:
                self.report_2 = {
                           'M': self.M_dry_zone_yes,
                           'theta_SPK':  self.theta_SPK,
                           'M_1':        self.M_1,
                           }
            else:
                self.report_2 = {
                           'M': self.M_dry_zone_no,
                            'theta_SPK':  None,
                            'M_1':        None,
                }

        return {**self.report_0, **self.report_1, **self.report_2}