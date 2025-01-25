from c4_open.database import get_cable_database
from c4_open import mutual_heating_factors

"""
This module contains the class Project. It is a container for the required project data.
"""

class Project:
    """
    Project objects are convenient containers for the required project data.
    You can also assign a name to a location object.

    Parameters
    ----------
    name : basestring
        Name of the project
    calc_case : basestring
        Calculation case
        Possible values are: 'ac_sc_pipe', 'dc_sc_pipe', 'ac_mc_pipe', 'ac_sc', 'dc_sc', 'ac_mc'
    cable_type : basestring
        Cable type according to 'cable_data.xlsx'
    L : float
        Depth of laying
    N : int
        Total number of parallel cables (or pipes)
    F : float
        Mutual heating factor F
    U_0 : float
        Phase-ground voltage (only for ac cables)
    U_n : float
        Phase-phase voltage (only for ac cables)
    f : float
        Mains frequency
    deltatheta_x : float
        Critical temperature rise of soil after drying of soil occurs
    K : int
        Number of cables in pipe (optional)
    dpipe_in : basestring
        Inner pipe diameter (optional)
    dpipe_out : basestring
        External pipe diameter (optional)
    theta_mean : basestring
        Mean air temperature inside pipe (optional),
        a good assumption is mean of ambient temperature and max. conductor temperature
    rho_T4 : float
        Thermal resistivity of soil
    rho_pipe : basestring
        Thermal resistivity of pipe material (optional)
    theta_amb : float
        Ambient temperature

    Raises
    ------
    ValueError
        If `calc_case` is unknown

    ValueError
        If `cable_type` is unknown
    """


    def __init__(self, name,
                 calc_case,
                 cable_type,
                 L,
                 N,
                 deltatheta_x,
                 rho_T4,
                 theta_amb,
                 F = 1,

                 # Parameters for ac cables
                 U_0 = None,
                 U_n = None,
                 f = None,

                 # Parameters for cables in pipes
                 K = None,
                 dpipe_in = None,
                 dpipe_out = None,
                 theta_mean = None,
                 rho_pipe = None,
                 ):

        self.name = name
        self.calc_case = calc_case
        self.cable_type = cable_type
        self.L = L
        self.N = N
        self.F = F
        self.U_0 = U_0
        self.U_n = U_n
        self.f = f
        self.deltatheta_x = deltatheta_x
        self.K = K
        self.dpipe_in = dpipe_in
        self.dpipe_out = dpipe_out
        self.theta_mean = theta_mean
        self.rho_T4 = rho_T4
        self.rho_pipe = rho_pipe
        self.theta_amb = theta_amb

        if cable_type not in get_cable_database().index:
            raise ValueError(f"Unknown cable type '{cable_type}'")
        
        self.D_e = get_cable_database().loc[self.cable_type].d_out

        if calc_case not in ['ac_sc_pipe', 'dc_sc_pipe', 'ac_mc_pipe', 'ac_sc', 'dc_sc', 'ac_mc']:
            raise ValueError(f"Unknown calculation case '{calc_case}'")

    def get_parameters(self):
        """
        Returns a dictionary with the project parameters
        """
        return {
            'name': self.name,
            'calc_case': self.calc_case,
            'cable_type': self.cable_type,
            'L': self.L,
            'N': self.N,
            'F': self.F,
            'U_0': self.U_0,
            'U_n': self.U_n,
            'f': self.f,
            'deltatheta_x': self.deltatheta_x,
            'K': self.K,
            'dpipe_in': self.dpipe_in,
            'dpipe_out': self.dpipe_out,
            'theta_mean': self.theta_mean,
            'rho_T4': self.rho_T4,
            'rho_pipe': self.rho_pipe,
            'theta_amb': self.theta_amb
        }
    def get_delta_theta_x(self, load_profile):
        """
        Updates the detatheta_x value based on the load profile
        
        Parameters
        ----------
        load_profile : LoadProfile
            Load profile object

        """
        self.deltatheta_x = load_profile.delta_theta_x

        return self

    def get_F(self, method, path_to_dxf=None):
        """
        Wrapper for the mutual heating factor calculation
        Returns the mutual heating factor F

        Parameters
        ----------
        method : basestring
            Method for the mutual heating factor calculation
            Possible values are: 'three_cables_flat', 'two_cables_flat', 'three_cables_trefoil', 'four_cables_square', 'from_dxf'

        path_to_dxf : basestring
            Path to the dxf file (required if method is 'from_dxf')
        """

        if method == 'three_cables_flat':
            self.F = mutual_heating_factors.get_mutual_heating_factor_for_three_cables_flat(self.L * 1000, self.D_e)

        elif method == 'two_cables_flat':
            self.F = mutual_heating_factors.get_mutual_heating_factor_for_two_cables_flat(self.L * 1000, self.D_e)

        elif method == 'three_cables_trefoil':
            self.F = mutual_heating_factors.get_mutual_heating_factor_for_three_cables_trefoil(self.L * 1000, self.D_e)

        elif method == 'four_cables_square':
            self.F = mutual_heating_factors.get_mutual_heating_factor_for_four_cables_square(self.L, self.D_e)

        elif method == 'from_dxf':
            if path_to_dxf is None:
                raise ValueError("Path to dxf file must be provided")
            self.F = mutual_heating_factors.get_mutual_heating_factor_from_dxf(path_to_dxf)

        else:
            raise ValueError(f"Unknown method '{method}'")

        return self
        