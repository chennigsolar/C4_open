import unittest
from c4_open import arrangements, current_carrying_capacities as ccc, cyclic_rating_factors, mutual_heating_factors
from c4_open.load_profile import LoadProfile
from c4_open.project import Project
from c4_open.cables import Cable
from c4_open.cyclic_rating_factors import FactorM



class FunctionTests(unittest.TestCase):
    def test_skin_effect_factor(self):
        # test data for 150 mm² Al cable NA2XS2Y
        skin_effect_factor = ccc.get_skin_effect_factor(50, 2.06e-4 * (1 + 0.004 * 70))
        self.assertEqual(skin_effect_factor, 0.0011818253923873495)

    def test_proximity_effect_factor_round(self):
        # test data for 150 mm² Al cable NA2XS2Y
        proximity_effect_factor = ccc.get_proximity_effect_factor_round(50,
                                                                        2.06e-4 * (1 + 0.004 * 70),
                                                                        14.4, 38,
                                                                        'Al')
        self.assertEqual(proximity_effect_factor, 0.0004783918387023795)

    def test_proximity_effect_factor_sector_shaped(self):
        # test data for 150 mm² Al cable NA2XS2Y
        proximity_effect_factor = ccc.get_proximity_effect_factor_sector_shaped(50,
                                                                                2.06e-4 * (1 + 0.004 * 70),
                                                                                240,
                                                                                2.2)
        self.assertEqual(proximity_effect_factor, 0.00233061803637006)

    def test_ac_resistance(self):
        ac_resistance = ccc.get_ac_resistance(0.0011818253923873495,
                                              0.0007460731272216992,
                                              2.06e-4 * (1 + 0.004 * 70))
        self.assertEqual(ac_resistance, 0.0002641883482816506)

    def test_dielectric_loss(self):
        dielectric_loss = ccc.get_dielectric_loss(2.6e-10, 12000, 50)
        self.assertEqual(dielectric_loss, 0.047048491580160744)

    def test_loss_factor_trefoil(self):
        loss_factor = ccc.get_loss_factor_trefoil(2.06e-4, 7.41e-4, 38, 27.2, 50)
        self.assertEqual(loss_factor, (0.02709988897696266, 6.456075739922802e-05))

    def test_loss_factor_three_flat(self):
        loss_factor = ccc.get_loss_factor_three_flat(2.06e-4, 7.41e-4, 38, 27.2, 50)
        self.assertEqual(loss_factor, (0.04050492176837549, 7.907799800125202e-05))

    def test_t1_single_score(self):
        T_1 = ccc.get_t1_single_core(3.5, 15.8, 1.4)
        self.assertEqual(T_1, 0.09088236531068454)

    def test_t1_belted_sector_shaped(self):
        # Test data for NAYY-J 4x150 mm²
        T_1 = ccc.get_t1_belted_sector_shaped(3.5, 1.8, 150, 41, 17.2)
        self.assertEqual(T_1, 0.32215995934589886)

    def test_t3(self):
        # Test data for NAYY-J 4x150 mm²
        T_3 = ccc.get_t3(3.5, 18.4, 1.8)
        self.assertEqual(T_3, 0.09953888513954343)

    def test_t4_from_mutual_heating_factor(self):
        T_3 = ccc.get_t4_from_mutual_heating_factor(1, 700, 22, 1)
        self.assertEqual(T_3, 0.7713079060163779)

    def test_thermal_resistance_between_cable_and_pipe(self):
        T_41 = ccc.get_thermal_resistance_between_cable_and_pipe(3, 38, 55)
        self.assertEqual(T_41, 0.35881254155365827)

    def test_thermal_resistance_between_cable_and_pipe_invalid_number(self):
        with self.assertRaises(ValueError):
            ccc.get_thermal_resistance_between_cable_and_pipe(5, 38, 55)

    def test_thermal_resistance_of_pipe(self):
        T_42 = ccc.get_thermal_resistance_of_pipe(3.5, 160, 138)
        self.assertEqual(T_42, 0.08239776959571701)

    def test_t4_pipe(self):
        T_4 = ccc.get_t4_pipe(0.35881254155365827, 0.08239776959571701, 1, 700, 160, 1)
        self.assertEqual(T_4, 0.8962214758127067)

    def test_current_carrying_capacity_dc(self):
        I = ccc.get_current_carrying_capacity_dc(90,
                                                 20,
                                                 15,
                                                 1,
                                                 2.06e-4 * (1 + 0.004 * 70),
                                                 1,
                                                 0.09088236531068454,
                                                 0,
                                                 0.09953888513954343,
                                                 0.7713079060163779,
                                                 dry_zone=True)
        self.assertEqual(I, 406.910089878596)

    def test_current_carrying_capacity_ac(self):
        I = ccc.get_current_carrying_capacity_ac(90,
                                                 20,
                                                 15,
                                                 1,
                                                 0.0002641883482816506,
                                                 0.047048491580160744,
                                                 0.02113104974373334,
                                                 0,
                                                 1,
                                                 0.3291386548277138,
                                                 0,
                                                 0.0940005981829874,
                                                 1.8323328774746943,
                                                 dry_zone=True)
        self.assertEqual(I, 261.62334379175905)

    def test_current_carrying_capacity_dc_in_pipe(self):
        I = ccc.get_current_carrying_capacity_dc_in_pipe(90,
                                                         20,
                                                         15,
                                                         1,
                                                         0.000264,
                                                         1,
                                                         2,
                                                         0.092,
                                                         0,
                                                         0.100,
                                                         1.189,
                                                         dry_zone=True)
        self.assertEqual(I, 238.94109042324462)

    def test_current_carrying_capacity_ac_single_core_in_pipe(self):
        I = ccc.get_current_carrying_capacity_ac_single_core_in_pipe(90,
                                                                     15,
                                                                     15,
                                                                     1,
                                                                     0.000264,
                                                                     0.047048,
                                                                     0.016616,
                                                                     3,
                                                                     0.329,
                                                                     0.031,
                                                                     0.896,
                                                                     dry_zone=True)
        self.assertEqual(I, 226.2269780946717)

    def test_current_carrying_capacity_ac_multi_core_in_pipe(self):
        I = ccc.get_current_carrying_capacity_ac_multi_core_in_pipe(70,
                                                                    20,
                                                                    15,
                                                                    1,
                                                                    0.000249,
                                                                    0.005675,
                                                                    3,
                                                                    1,
                                                                    0.460,
                                                                    0.092,
                                                                    1.122,
                                                                    dry_zone=True)
        self.assertEqual(I, 178.3131472389703)

    def test_maximum_operating_temperature_resistance(self):
        R = ccc.get_maximum_operating_temperature_resistance(2.06e-4, 0.004, 90)
        self.assertEqual(R, 2.6368e-4)

    def test_cyclic_rating_factor_for_group(self):
        M = cyclic_rating_factors.get_cyclic_rating_factor_for_group(27.5,
                                                                     15,
                                                                     70,
                                                                     50,
                                                                     3,
                                                                     101,
                                                                     1000,
                                                                     1.2,
                                                                     0.6,
                                                                     1,
                                                                     1,
                                                                     1,
                                                                     1,
                                                                     1,
                                                                     1,
                                                                     )

        self.assertEqual(round(M[0], 2), 1.17)

    def test_cyclic_rating_factor_for_single_cable(self):
        M = cyclic_rating_factors.get_cyclic_rating_factor_for_single_cable(48.13,
                                                                            15,
                                                                            70,
                                                                            60,
                                                                            1000,
                                                                            1.2,
                                                                            0.504,
                                                                            0.992,
                                                                            0.728,
                                                                            0.640,
                                                                            0.596,
                                                                            0.593,
                                                                            0.796,
                                                                            )

        self.assertEqual(round(M[0], 2), 1.17)

    def test_M_1(self):
        M_1, _ = cyclic_rating_factors.get_M_1(1.17,
                                               0.7,
                                               0.992,
                                               70,
                                               20,
                                               25,
                                               1)

        self.assertEqual(round(M_1, 2), 1.22)

    def test_joule_loss(self):
        W = cyclic_rating_factors.get_joule_loss(0.0002642,
                                                 1,
                                                 200,
                                                 0.0470,
                                                 0.021131)
        self.assertEqual(W, 10.838312408000002)

    def test_joule_loss_pipe(self):
        W = cyclic_rating_factors.get_joule_loss_pipe(0.0002642,
                                                      1,
                                                      3,
                                                      200,
                                                      0.0470,
                                                      0.021131)
        self.assertEqual(W, 3 * 10.838312408000002)


    def test_mutual_heating_factor_for_three_cables_trefoil(self):
        F = mutual_heating_factors.get_mutual_heating_factor_for_three_cables_trefoil(700, 30)
        self.assertEqual(F, 2177.7777777777774)

    def test_mutual_heating_factor_for_three_cables_flat(self):
        F = mutual_heating_factors.get_mutual_heating_factor_for_three_cables_flat(700, 30)
        self.assertEqual(F, 2177.7777777777774)

    def test_mutual_heating_factor_for_two_cables_flat(self):
        F = mutual_heating_factors.get_mutual_heating_factor_for_two_cables_flat(700, 30)
        self.assertEqual(F, 46.666666666666664)

    def test_mutual_heating_factor_for_four_cables_square(self):
        F = mutual_heating_factors.get_mutual_heating_factor_for_four_cables_square(700, 30)
        self.assertEqual(F, 71863.00028058837)

    def test_t4_for_three_trefoil(self):
        T_4 = ccc.get_t4_for_three_trefoil(1, 700, 30)
        self.assertEqual(T_4, 1.9452295268798752)

    def test_t4_for_three_flat(self):
        T_4 = ccc.get_t4_for_three_flat(1, 700, 30)
        self.assertEqual(T_4, 2.0126842243880416)

    def test_t4_for_two_flat(self):
        T_4 = ccc.get_t4_for_two_flat(1, 700, 30)
        self.assertEqual(T_4, 1.350008668264133)
    
    def test_create_arrangement(self):
        F, _, _, _ = arrangements.create_arrangement('trefoil', 3, 0.028, 0.07, 0.7)
        self.assertEqual(round(F, 2), 5147580090.50)
        F, _, _, _ = arrangements.create_arrangement('three_flat', 2, 0.028, 0.07, 0.7)
        self.assertEqual(round(F, 2), 1828250.58)

    def test_load_profile(self):
        lp = LoadProfile('../examples/example_load_profile.xlsx')
        self.assertEqual(lp.mu, 0.3284154551254132)
        self.assertEqual(lp.load_factor, 0.3985628301965643)
        self.assertEqual(lp.delta_theta_x, 35.037857754454414)
        self.assertEqual(lp.Y_0, 1.0)
        self.assertEqual(lp.Y_1, 1.0)
        self.assertEqual(lp.Y_2, 1.0)
        self.assertEqual(lp.Y_3, 1.0)
        self.assertEqual(lp.Y_4, 1.0)
        self.assertEqual(lp.Y_5, 0.7276871934425999)
        self.assertEqual(lp.Y_0_hour, 13)

    def test_ac_sc(self):
        project = Project('Test',
                         'ac_sc',
                         'NA2XS2Y 1x240 20kV',
                         0.7,
                         3,
                         30,
                         1.0,
                         20,
                         1,
                         12000,
                         20000,
                         50,
                         )
        project.get_F('three_cables_trefoil')
        load_profile = LoadProfile('../examples/example_load_profile.xlsx')
        project.get_delta_theta_x(load_profile)
        cable = Cable(project)
        M = FactorM(cable, load_profile)
        self.assertEqual(round(cable.get_result()['I (with dry zone)'], 2), 387.34)
        self.assertEqual(round(M.get_result()['M_1'], 2), 1.22)


    def test_dc_sc(self):
        project = Project('Test',
                         'dc_sc',
                         'A2XH 1x240',
                         0.7,
                         2,
                         30,
                         1.0,
                         20,
                         1
                         )
        project.get_F('two_cables_flat')
        load_profile = LoadProfile('../examples/example_load_profile.xlsx')
        project.get_delta_theta_x(load_profile)
        cable = Cable(project)
        M = FactorM(cable, load_profile)
        self.assertEqual(round(cable.get_result()['I (with dry zone)'], 2), 462.24)
        self.assertEqual(round(M.get_result()['M_1'], 2), 1.19)

    def test_ac_mc(self):
        project = Project('Test',
                         'ac_mc',
                         'NAYY-J 4x70',
                         0.7,
                         1,
                         30,
                         1.0,
                         20,
                         1,
                         230,
                         400,
                         50,
                         )
        load_profile = LoadProfile('../examples/example_load_profile.xlsx')
        project.get_delta_theta_x(load_profile)
        cable = Cable(project)
        M = FactorM(cable, load_profile)
        self.assertEqual(round(cable.get_result()['I (no dry zone)'], 2), 177.21)
        self.assertEqual(round(M.get_result()['M (no dry zone)'], 2), 1.15)

    def test_ac_sc_pipe(self):
        project = Project('Test',
                         'ac_sc_pipe',
                         'NA2XS2Y 1x240 20kV',
                         0.7,
                         1,
                         35.038,
                         1.0,
                         20,
                         1,
                         12000,
                         20000,
                         50,
                         3,
                         138,
                         160,
                         55,
                          rho_pipe=3.5
                         )
        load_profile = LoadProfile('../examples/example_load_profile.xlsx')
        project.get_delta_theta_x(load_profile)
        cable = Cable(project)
        M = FactorM(cable, load_profile)
        self.assertEqual(round(cable.get_result()['I (no dry zone)'], 2), 371.16)
        self.assertEqual(round(M.get_result()['M (no dry zone)'], 3), 1.137)


if __name__ == "__main__":
    unittest.main()
