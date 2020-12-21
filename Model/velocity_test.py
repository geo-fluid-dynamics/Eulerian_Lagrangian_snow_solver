import unittest
import velocity
from ConstantVariables import a_eta, b_eta, eta_0, c_eta, T_fus,g, rho_i,Z_max ,D_rate
import numpy as np

class TestVelocity(unittest.TestCase):
    def test_choose_viscosity(self):
        eta_const_n1 = velocity.choose_viscosity(263, 0.1125, 'eta_constant_n1')
        self.assertEqual(eta_const_n1,eta_0 * rho_i * 0.1125/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.1125) )
        eta_T = velocity.choose_viscosity(263, 0.1125, 'eta_T')
        self.assertEqual(eta_T,eta_0 * rho_i * 0.1125/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.1125) )
        eta_phi = velocity.choose_viscosity(263, 0.1125, 'eta_phi')
        self.assertEqual(eta_phi,eta_0 * rho_i * 0.1125/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.1125) )
        eta_phiT = velocity.choose_viscosity(263, 0.1125, 'eta_phiT')
        self.assertEqual(eta_phiT,eta_0 * rho_i * 0.1125/c_eta * np.exp(a_eta *(T_fus - 263)+ b_eta *rho_i * 0.1125) )
        eta_const_n3 = velocity.choose_viscosity(263, 0.1125, 'eta_constant_n3')
        self.assertEqual(eta_const_n3, ((0.1125 * rho_i * g * Z_max)**3)/D_rate)

    def test_sigma_cont_croc(self):
        sigma = velocity.sigma_cont_croc(np.array([0.25,0.25]), np.array([200,100,100]),3, v_opt='layer_based')
        sigma_dz_bottom = 0.25*g*rho_i*200
        sigma_dz_centre = 0.25*g*rho_i*100
        sigma_dz_top = sigma_dz_centre/2
        sigma_top = sigma_dz_top
        sigma_centre = sigma_dz_centre+sigma_dz_top
        sigma_bottom = sigma_centre+sigma_dz_bottom
        self.assertEqual(sigma[-1],sigma_top)
        self.assertEqual(sigma[1], sigma_centre)
        self.assertEqual(sigma[0], sigma_bottom)

        sigma = velocity.sigma_cont_croc(np.array([0.25,0.25]), np.array([200,100,100]),3, v_opt='continuous')
        sigma_dz_bottom = 0.25*g*rho_i*200
        sigma_dz_centre = 0.25*g*rho_i*100
        sigma_dz_top = 0
        sigma_top = sigma_dz_top
        sigma_centre = sigma_dz_centre+sigma_dz_top
        sigma_bottom = sigma_centre+sigma_dz_bottom
        self.assertEqual(sigma[-1],sigma_top)
        self.assertEqual(sigma[1], sigma_centre)
        self.assertEqual(sigma[0], sigma_bottom)        

    def test_velocity(self):
        (v, v_dz) = velocity.velocity(np.array([1,1,1]), np.array([1,1,1]),np.array([1,1]),3)
        self.assertEqual(v[1],-1)
        self.assertEqual(v_dz[1],-1)

if __name__ == '__main__':
    unittest.main()