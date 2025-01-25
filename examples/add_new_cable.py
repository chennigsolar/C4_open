from c4_open.database import add_cable_to_database

data = {'Type': 'NA2XY 1x240',
            'A_c': 240,
            'cond_mat': "Al",
            'theta_max': 90,
            'R__20': 0.125e-3,
            'alpha_con': 0.004,
            'R__scr20': 0.0,
            'alpha_scr': 0.0,
            'C': 8e-10,
            'd_c': 20,
            'd_ins': 2.2,
            'd_out': 28.0,
            't_ins': 2.2,
            't_sheath': 1.8,
            'n': 1,
            'rho_T1': 3.5,
            'rho_T3': 5,
            'd_c_max': None,
            'd_c_min': None,
            'd_a': None,
            'r_1': None,
            's': None,
            'd_x': None}

add_cable_to_database(data)