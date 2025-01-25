import c4_open.current_carrying_capacities as ccc
from c4_open.database import get_cable_database, get_cable_types


class Cable:
    """
    The Cable class provides a simple interface to the different calculation functions for the determination
    of the current carrying capacity of buried cables. It selects the applicable calculation functions depending on the
    given calculation case and calls them in the required sequence.

    The method Cable.get_result() returns a dict with the required information for the calculation of the
    cyclic rating factor M.

    The method Cable.get_report() returns a report of all calculation results in form of a dictionary for documentation
    purposes.

    Please see the example jupyter notebook 'example_ac.ipynb' for further explanation.

    The implemented six calculation cases are:
     - 'ac_sc':          ac single core cable, optionally with screen
     - 'ac_sc_pipe':     as above but laid in pipe/conduit
     - 'dc_sc':          dc single core cable, optionally with screen
     - 'dc_sc_pipe':     as above but laid in pipe/conduit
     - 'ac_mc':          ac multicore cable
     - 'ac_mc_pipe':     as above but laid in pipe/conduit

    Parameters
    ----------
    project : Project
        An instance of the Project class

    external_resistance_method : str, optional
        If not specified, the method 'from_F' will be chosen.
        Valid strings are 'from_F', 'three_trefoil', 'three_flat', 'two_flat'.

    Raises
    ------
    ValueError
        If `calc_case` is unknown
    """

    def __init__(self, project, external_resistance_method='from_F'):
        self.external_resistance_method = external_resistance_method

        # Reading data from project
        # Convert arguments into float or string respectively
        self.calc_case = str(project.calc_case)
        if (self.calc_case == "ac_sc"
                or self.calc_case == "ac_mc"
                or self.calc_case == "ac_sc_pipe"
                or self.calc_case == "ac_mc_pipe"):
            self.U_n = float(project.U_n)
            self.U_0 = float(project.U_0)
            self.f = float(project.f)
        self.factor_F = float(project.F)
        self.L = float(project.L)
        self.rho_T4 = float(project.rho_T4)
        self.theta_amb = float(project.theta_amb)
        self.deltatheta_x = float(project.deltatheta_x)
        self.cable_type = str(project.cable_type)
        self.N = int(project.N)

        # Reading cable data from cable database
        cable_data = get_cable_database()
        self.cable_parameters = cable_data.loc[self.cable_type].to_dict()

        self.A_c = float(self.cable_parameters.get('A_c'))
        self.cond_mat = str(self.cable_parameters.get('cond_mat'))
        self.theta_max = float(self.cable_parameters.get('theta_max'))
        self.R__20 = float(self.cable_parameters.get('R__20'))
        self.alpha_con = float(self.cable_parameters.get('alpha_con'))
        if (self.calc_case == "ac_sc"
                or self.calc_case == "ac_mc"
                or self.calc_case == "ac_sc_pipe"
                or self.calc_case == "ac_mc_pipe"):
            self.R__scr20 = float(self.cable_parameters.get('R__scr20'))
            self.alpha_scr = float(self.cable_parameters.get('alpha_scr'))
            self.C = float(self.cable_parameters.get('C'))
        self.d_c = float(self.cable_parameters.get('d_c'))
        self.d_ins = float(self.cable_parameters.get('d_ins'))
        self.d_out = float(self.cable_parameters.get('d_out'))
        self.t_ins = float(self.cable_parameters.get('t_ins'))
        self.t_sheath = float(self.cable_parameters.get('t_sheath'))
        self.n = int(self.cable_parameters.get('n'))
        self.rho_T1 = float(self.cable_parameters.get('rho_T1'))
        self.rho_T3 = float(self.cable_parameters.get('rho_T3'))
        if self.calc_case == "ac_mc" or self.calc_case == "ac_mc_pipe":
            self.d_a = float(self.cable_parameters.get('d_a'))  # external diameter of the belt insulation
            self.r_1 = float(self.cable_parameters.get('r_1'))  # radius of the circle circumscribing the conductors

        if (self.calc_case == "dc_sc_pipe"
                or self.calc_case == "ac_sc_pipe"
                or self.calc_case == "ac_sc_pipe"
                or self.calc_case == "ac_mc_pipe"):
            self.rho_pipe = float(project.rho_pipe)
            self.D_d = float(project.dpipe_in)
            self.D_o = float(project.dpipe_out)
            self.theta_mean = float(project.theta_mean)
            self.K = int(project.K)

        self.R_ = None
        self.Y_s = None
        self.Y_p = None
        self.R = None
        self.W_d = None
        self.lambda_1 = None
        self.X = None
        self.T_1 = None
        self.T_2 = None
        self.T_3 = None
        self.T_4 = None
        self.T_41 = None
        self.T_42 = None
        self.T_4 = None
        self.I_dry_zone_yes = None
        self.I_dry_zone_no = None

        # Initialize result and report dictionaries
        self.result = {}

        self.report = {}
        self.report_0 = {}
        self.report_1 = {}
        self.report_2 = {}
        self.report_3 = {}
        self.report_4 = {}
        self.report_5 = {}
        self.report_6 = {}
        self.report_7 = {}
        self.report_8 = {}
        self.report_9 = {}
        self.report_10 = {}
        self.report_11 = {}
        self.report_12 = {}

        self.run()

    def run(self):
        self.R_ = ccc.get_maximum_operating_temperature_resistance(self.R__20, self.alpha_con, self.theta_max)
        if self.calc_case == 'ac_mc':
            self.T_4 = ccc.get_t4_from_mutual_heating_factor(self.rho_T4, self.L * 1000, self.d_out, self.factor_F)

        elif self.calc_case == 'ac_sc' or self.calc_case == 'dc_sc':
            if self.external_resistance_method == 'from_F':
                self.T_4 = ccc.get_t4_from_mutual_heating_factor(self.rho_T4, self.L * 1000, self.d_out, self.factor_F)

            elif self.external_resistance_method == 'three_trefoil':
                self.T_4 = ccc.get_t4_for_three_trefoil(self.rho_T4, self.L * 1000, self.d_out)

            elif self.external_resistance_method == 'three_flat':
                self.T_4 = ccc.get_t4_for_three_flat(self.rho_T4, self.L * 1000, self.d_out)

            elif self.external_resistance_method == 'two_flat':
                self.T_4 = ccc.get_t4_for_two_flat(self.rho_T4, self.L * 1000, self.d_out)

        elif self.calc_case == 'ac_sc_pipe' or self.calc_case == 'dc_sc_pipe' or self.calc_case == 'ac_mc_pipe':
            self.T_41 = ccc.get_thermal_resistance_between_cable_and_pipe(self.K, self.d_out, self.theta_mean)
            self.T_42 = ccc.get_thermal_resistance_of_pipe(self.rho_pipe, self.D_o, self.D_d)
            self.T_4 = ccc.get_t4_pipe(self.T_41, self.T_42, self.rho_T4, self.L * 1000, self.D_o, self.factor_F)

        if self.calc_case == 'ac_sc':
            self.Y_s = ccc.get_skin_effect_factor(self.f, self.R_)
            self.Y_p = ccc.get_proximity_effect_factor_round(self.f, self.R_, self.d_c, self.d_out, 'Al')
            self.R = ccc.get_ac_resistance(self.Y_s, self.Y_p, self.R_)
            self.W_d = ccc.get_dielectric_loss(self.C, self.U_0, self.f)
            # Diameter of screen is approx. diameter over insulation and distance between conductor axes is the external
            # cable diameter
            if self.external_resistance_method == 'three_trefoil':
                self.lambda_1, self.X = ccc.get_loss_factor_trefoil(self.R,
                                                                    self.R__scr20,
                                                                    self.d_out,
                                                                    self.d_ins,
                                                                    self.f)
            elif self.external_resistance_method == 'three_flat':
                self.lambda_1, self.X = ccc.get_loss_factor_three_flat(self.R,
                                                                       self.R__scr20,
                                                                       self.d_out,
                                                                       self.d_ins,
                                                                       self.f)
            else:
                self.lambda_1, self.X = ccc.get_loss_factor_trefoil(self.R,
                                                                    self.R__scr20,
                                                                    self.d_out,
                                                                    self.d_ins,
                                                                    self.f)
            self.T_1 = ccc.get_t1_single_core(self.rho_T1, self.d_c, self.t_ins)
            self.T_3 = ccc.get_t3(self.rho_T3, self.d_ins, self.t_sheath)

            self.I_dry_zone_yes = ccc.get_current_carrying_capacity_ac(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R,
                self.W_d,
                self.lambda_1,
                0,
                self.n,
                self.T_1,
                0,
                self.T_3,
                self.T_4,
                dry_zone=True
                )

            self.I_dry_zone_no = ccc.get_current_carrying_capacity_ac(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R,
                self.W_d,
                self.lambda_1,
                0,
                self.n,
                self.T_1,
                0,
                self.T_3,
                self.T_4,
                dry_zone=False
            )

            self.result = {'Y_s': self.Y_s,
                           'Y_p': self.Y_p,
                           'R': self.R,
                           'W_d': self.W_d,
                           'lambda_1': self.lambda_1,
                           'T_1': self.T_1,
                           'T_3': self.T_3,
                           'T_4': self.T_4,
                           "I (no dry zone)": self.I_dry_zone_no,
                           "I (with dry zone)": self.I_dry_zone_yes
                           }

        elif self.calc_case == 'dc_sc':
            self.T_1 = ccc.get_t1_single_core(self.rho_T1, self.d_c, self.t_ins)
            self.T_3 = ccc.get_t3(self.rho_T3, self.d_ins, self.t_sheath)
            self.I_dry_zone_yes = ccc.get_current_carrying_capacity_dc(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R_,
                self.n,
                self.T_1,
                0,
                self.T_3,
                self.T_4,
                dry_zone=True
            )

            self.I_dry_zone_no = ccc.get_current_carrying_capacity_dc(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R_,
                self.n,
                self.T_1,
                0,
                self.T_3,
                self.T_4,
                dry_zone=False
            )

            self.result = {'R': self.R_,
                           'T_1': self.T_1,
                           'T_3': self.T_3,
                           'T_4': self.T_4,
                           "I (no dry zone)": self.I_dry_zone_no,
                           "I (with dry zone)": self.I_dry_zone_yes
                           }

        elif self.calc_case == 'ac_mc':
            self.Y_s = ccc.get_skin_effect_factor(self.f, self.R_)
            self.Y_p = ccc.get_proximity_effect_factor_sector_shaped(self.f, self.R_, self.A_c, self.t_ins)
            self.R = ccc.get_ac_resistance(self.Y_s, self.Y_p, self.R_)
            self.W_d = ccc.get_dielectric_loss(self.C, self.U_0, self.f)
            self.T_1 = ccc.get_t1_belted_sector_shaped(self.rho_T1, self.t_ins, self.A_c, self.d_a, self.r_1)
            self.T_3 = ccc.get_t3(self.rho_T3, self.d_a, self.t_sheath)

            self.I_dry_zone_yes = ccc.get_current_carrying_capacity_ac(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R,
                self.W_d,
                0,
                0,
                self.n,
                self.T_1,
                0,
                self.T_3,
                self.T_4,
                dry_zone=True
                )

            self.I_dry_zone_no = ccc.get_current_carrying_capacity_ac(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R,
                self.W_d,
                0,
                0,
                self.n,
                self.T_1,
                0,
                self.T_3,
                self.T_4,
                dry_zone=False
            )

            self.result = {'Y_s': self.Y_s,
                           'Y_p': self.Y_p,
                           'R': self.R,
                           'W_d': self.W_d,
                           'T_1': self.T_1,
                           'T_3': self.T_3,
                           'T_4': self.T_4,
                           "I (no dry zone)": self.I_dry_zone_no,
                           "I (with dry zone)": self.I_dry_zone_yes
                           }

        elif self.calc_case == 'ac_sc_pipe':
            self.Y_s = ccc.get_skin_effect_factor(self.f, self.R_)
            self.Y_p = ccc.get_proximity_effect_factor_round(self.f, self.R_, self.d_c, self.d_out, 'Al')
            self.R = ccc.get_ac_resistance(self.Y_s, self.Y_p, self.R_)
            self.W_d = ccc.get_dielectric_loss(self.C, self.U_0, self.f)
            # Diameter of screen is approx. diameter over insulation and distance between conductor axes is the external
            # cable diameter
            if self.external_resistance_method == 'three_trefoil':
                self.lambda_1, self.X = ccc.get_loss_factor_trefoil(self.R,
                                                                    self.R__scr20,
                                                                    self.d_out,
                                                                    self.d_ins,
                                                                    self.f)
            elif self.external_resistance_method == 'three_flat':
                self.lambda_1, self.X = ccc.get_loss_factor_three_flat(self.R,
                                                                       self.R__scr20,
                                                                       self.d_out,
                                                                       self.d_ins,
                                                                       self.f)
            else:
                self.lambda_1, self.X = ccc.get_loss_factor_trefoil(self.R,
                                                                    self.R__scr20,
                                                                    self.d_out,
                                                                    self.d_ins,
                                                                    self.f)

            self.T_1 = ccc.get_t1_single_core(self.rho_T1, self.d_c, self.t_ins)
            self.T_3 = ccc.get_t3(self.rho_T3, self.d_ins, self.t_sheath)

            self.I_dry_zone_yes = ccc.get_current_carrying_capacity_ac_single_core_in_pipe(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R,
                self.W_d,
                self.lambda_1,
                self.K,
                self.T_1,
                self.T_3,
                self.T_4,
                dry_zone=True
                )

            self.I_dry_zone_no = ccc.get_current_carrying_capacity_ac_single_core_in_pipe(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R,
                self.W_d,
                self.lambda_1,
                self.K,
                self.T_1,
                self.T_3,
                self.T_4,
                dry_zone=False
            )

            self.result = {'Y_s': self.Y_s,
                           'Y_p': self.Y_p,
                           'R': self.R,
                           'W_d': self.W_d,
                           'lambda_1': self.lambda_1,
                           'T_1': self.T_1,
                           'T_3': self.T_3,
                           'T_4': self.T_4,
                           "I (no dry zone)": self.I_dry_zone_no,
                           "I (with dry zone)": self.I_dry_zone_yes
                           }

        elif self.calc_case == 'dc_sc_pipe':
            self.T_1 = ccc.get_t1_single_core(self.rho_T1, self.d_c, self.t_ins)
            self.T_3 = ccc.get_t3(self.rho_T3, self.d_ins, self.t_sheath)

            self.I_dry_zone_yes = ccc.get_current_carrying_capacity_dc_in_pipe(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R_,
                self.n,
                self.K,
                self.T_1,
                0,
                self.T_3,
                self.T_4,
                dry_zone=True
            )

            self.I_dry_zone_no = ccc.get_current_carrying_capacity_dc_in_pipe(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R_,
                self.n,
                self.K,
                self.T_1,
                0,
                self.T_3,
                self.T_4,
                dry_zone=False
            )

            self.result = {'R': self.R_,
                           'T_1': self.T_1,
                           'T_3': self.T_3,
                           'T_4': self.T_4,
                           "I (no dry zone)": self.I_dry_zone_no,
                           "I (with dry zone)": self.I_dry_zone_yes
                           }

        elif self.calc_case == 'ac_mc_pipe':
            self.Y_s = ccc.get_skin_effect_factor(self.f, self.R_)
            self.Y_p = ccc.get_proximity_effect_factor_sector_shaped(self.f, self.R_, self.A_c, self.t_ins)
            self.R = ccc.get_ac_resistance(self.Y_s, self.Y_p, self.R_)
            self.W_d = ccc.get_dielectric_loss(self.C, self.U_0, self.f)
            self.T_1 = ccc.get_t1_belted_sector_shaped(self.rho_T1, self.t_ins, self.A_c, self.d_a, self.r_1)
            self.T_3 = ccc.get_t3(self.rho_T3, self.d_a, self.t_sheath)

            self.I_dry_zone_yes = ccc.get_current_carrying_capacity_ac_multi_core_in_pipe(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R,
                self.W_d,
                self.n,
                self.K,
                self.T_1,
                self.T_3,
                self.T_4,
                dry_zone=True
                )

            self.I_dry_zone_no = ccc.get_current_carrying_capacity_ac_multi_core_in_pipe(
                self.theta_max,
                self.theta_amb,
                self.deltatheta_x,
                self.rho_T4,
                self.R,
                self.W_d,
                self.n,
                self.K,
                self.T_1,
                self.T_3,
                self.T_4,
                dry_zone=False
            )

            self.result = {'Y_s': self.Y_s,
                           'Y_p': self.Y_p,
                           'R': self.R,
                           'W_d': self.W_d,
                           'T_1': self.T_1,
                           'T_3': self.T_3,
                           'T_4': self.T_4,
                           "I (no dry zone)": self.I_dry_zone_no,
                           "I (with dry zone)": self.I_dry_zone_yes
                       }

        else:
            raise ValueError('Invalid Calculation Case')

    def get_result(self):
        return self.result

    def get_report(self):
        """
        This method returns a report of the calculation for documentation purposes in form of a dictionary
        """
        self.report_0 = {'cable_type': self.cable_type,
                       'theta_max': self.theta_max,
                       'R__20': self.R__20,
                       'alpha_con': self.alpha_con,
                       'R_': self.R_
                       }


        if self.calc_case == 'ac_mc' or self.calc_case == 'ac_mc_pipe'\
                or self.calc_case == 'ac_sc_pipe' or self.calc_case == 'ac_sc':
            self.report_1 = {'C': self.C,
                             'R__scr20': self.R__scr20,
                             'alpha_scr': self.alpha_scr
                             }

        if self.calc_case == 'ac_sc' or self.calc_case == 'ac_sc_pipe'\
                or self.calc_case == 'dc_sc' or self.calc_case == 'dc_sc_pipe':
            self.report_2 = {'d_c': self.d_c}

        self.report_3 = {'d_out': self.d_out,
                       't_ins': self.t_ins,
                       't_sheath': self.t_sheath,
                       'n': self.n,
                       'rho_T1': self.rho_T1,
                       'rho_T3': self.rho_T3
                       }

        if self.calc_case == 'ac_mc' or self.calc_case == 'ac_mc_pipe':
            self.report_4 = {'d_a': self.d_a,
                           'r_1': self.r_1
                           }
        if self.calc_case == 'ac_mc_pipe' or self.calc_case == 'dc_sc_pipe' or self.calc_case == 'ac_sc_pipe':
            self.report_5 = {'K': self.K,
                           'Dpipe_in': self.D_d,
                           'Dpipe_out': self.D_o,
                           'rho_pipe': self.rho_pipe
                           }
        if self.calc_case == 'ac_mc' or self.calc_case == 'ac_mc_pipe' or self.calc_case == 'ac_sc_pipe' or self.calc_case == 'ac_sc':
            self.report_6 = {'U_n': self.U_n,
                           'U_0': self.U_0,
                           'f': self.f
                           }
        self.report_7 = {'L': self.L,
                       'factor_F': self.factor_F,
                       'rho_T4': self.rho_T4,
                       'theta_amb': self.theta_amb,
                       'deltatheta_x': self.deltatheta_x
                       }

        if self.calc_case == 'ac_mc' or self.calc_case == 'ac_mc_pipe'  \
            or self.calc_case == 'ac_sc_pipe' or self.calc_case == 'ac_sc':
            self.report_8 = {'Y_s': self.Y_s,
                           'Y_p': self.Y_p,
                           'R': self.R,
                           'W_d': self.W_d
                           }

        if self.calc_case == 'ac_sc' or self.calc_case == 'ac_sc_pipe':

            self.report_9 = {'lambda_1': self.lambda_1,
                            }

        self.report_10 = {'T_1': self.T_1,
                       'T_3': self.T_3,
                       'T_4': self.T_4
                       }

        if self.calc_case == 'ac_sc_pipe' or self.calc_case == 'ac_mc_pipe' or self.calc_case == 'dc_mc_pipe':
            self.report_11 = {'T_41': self.T_41,
                           'T_42': self.T_42
                           }
        self.report_12 = {'I (with dry zone)': self.I_dry_zone_yes,
                        'I (without dry zone)': self.I_dry_zone_no
                        }

        self.report = {**self.report_0, **self.report_1, **self.report_2, **self.report_3, **self.report_4,
                       **self.report_5, **self.report_6, **self.report_7, **self.report_8, **self.report_9,
                       **self.report_10, **self.report_11, **self.report_12}

        return self.report