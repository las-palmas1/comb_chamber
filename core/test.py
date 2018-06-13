import unittest
from .geom import CombustionChamberGeom
from gas_turbine_cycle.gases import Air
import numpy as np


class CombChamberGeomTest(unittest.TestCase):
    def setUp(self):
        self.comb_geom = CombustionChamberGeom(
            work_fluid_in=Air(),
            Q_n=48e6,
            sigma_d=0.998,
            sigma_front=0.998,
            n_pipe=12,
            G_in=40,
            p_stag_in=15e5,
            T_stag_in=550,
            F_in=np.pi * 0.5 ** 2 / 4,
            D_d_av=0.37,
            G_fuel_in=0,
            G_fuel=0.9,
            alpha_sum=3.5,
            alpha1=1.1,
            c_d=35,
            c1=15,
            H=4.5e6,
            l1_rel=0.6,
            eta_burn=0.99
        )
        self.comb_geom.compute()

    def test_plot_pipes(self):
        self.comb_geom.plot(figsize=(8, 8))