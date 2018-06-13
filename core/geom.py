from gas_turbine_cycle.gases import IdealGas
from gas_turbine_cycle.tools.gas_dynamics import GasDynamicFunctions as gd
import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt


class Diffuser:
    def __init__(self, work_fluid: IdealGas, sigma, G, p_stag_in, T_stag_in, F_in, c_out, G_fuel, D_out_av):
        self.work_fluid = work_fluid
        self.sigma = sigma
        self.G = G
        self.G_fuel = G_fuel
        self.p_stag_in = p_stag_in
        self.T_stag_in = T_stag_in
        self.F_in = F_in
        self.c_out = c_out
        self.D_out_av = D_out_av

        self.T_stag_out = None
        self.p_stag_out = None
        self.fuel_content = None
        self.alpha = None
        self.c_p = None
        self.k = None
        self.c_in = None
        self.static_outlet = None
        self.T_in = None
        self.p_in = None
        self.rho_in = None
        self.lam_in = None
        self.static_outlet = None
        self.T_out = None
        self.p_out = None
        self.rho_out = None
        self.lam_out = None
        self.F_out = None
        self.D_out_hub = None
        self.D_out_per = None
        self.a_cr_out = None
        self.tau_out = None
        self.pi_out = None

    @classmethod
    def get_static(cls, c, T_stag, p_stag, k, R):
        a_cr = gd.a_cr(T_stag, k, R)
        lam = c / a_cr
        tau = gd.tau_lam(lam, k)
        pi = gd.pi_lam(lam, k)
        T = T_stag * tau
        p = p_stag * pi
        rho = p / (R * T)
        return T, p, rho, lam, a_cr, tau, pi

    def compute(self):
        self.T_stag_out = self.T_stag_in
        self.p_stag_out = self.p_stag_in * self.sigma
        self.fuel_content = self.G_fuel / (self.G - self.G_fuel)
        if self.fuel_content != 0:
            self.alpha = 1 / (self.work_fluid.l0 * self.fuel_content)
        else:
            self.alpha = 1e6
        self.c_p = self.work_fluid.c_p_real_func(self.T_stag_in, alpha=self.alpha)
        self.k = self.work_fluid.k_func(self.c_p)
        c_init = 0.7 * self.G / (self.F_in * self.p_stag_in / (self.T_stag_in * self.work_fluid.R))
        self.c_in = newton(
            lambda c: self.F_in * c * self.get_static(c, self.T_stag_out, self.p_stag_out,
                                                       self.k, self.work_fluid.R)[2] - self.G,
            x0=c_init
        )
        self.static_outlet = self.get_static(self.c_in, self.T_stag_in, self.p_stag_in, self.k, self.work_fluid.R)
        self.T_in = self.static_outlet[0]
        self.p_in = self.static_outlet[1]
        self.rho_in = self.static_outlet[2]
        self.lam_in = self.static_outlet[3]
        self.static_outlet = self.get_static(self.c_out, self.T_stag_out, self.p_stag_out, self.k, self.work_fluid.R)
        self.T_out = self.static_outlet[0]
        self.p_out = self.static_outlet[1]
        self.rho_out = self.static_outlet[2]
        self.lam_out = self.static_outlet[3]
        self.a_cr_out = self.static_outlet[4]
        self.tau_out = self.static_outlet[5]
        self.pi_out = self.static_outlet[6]
        self.F_out = self.G / (self.c_out * self.rho_out)
        self.D_out_hub = self.D_out_av - self.F_out / (np.pi * self.D_out_av)
        self.D_out_per = self.D_out_av + self.F_out / (np.pi * self.D_out_av)


class CombustionChamberGeom:
    def __init__(self, work_fluid_in: IdealGas, Q_n,
                 sigma_d, sigma_front, n_pipe, G_in, p_stag_in, T_stag_in, F_in,
                 D_d_av, G_fuel_in, G_fuel, alpha_sum, alpha1,
                 c_d=60, c1=12, H=4.5e6, l1_rel=0.6, eta_burn=0.99):
        self.work_fluid_in = work_fluid_in
        self.Q_n = Q_n
        self.sigma_d = sigma_d
        self.sigma_front = sigma_front
        self.G_in = G_in
        self.p_stag_in = p_stag_in
        self.T_stag_in = T_stag_in
        self.F_in = F_in
        self.G_fuel_in = G_fuel_in
        self.G_fuel = G_fuel
        self.alpha1 = alpha1
        self.c_d = c_d
        self.c1 = c1
        self.H = H
        self.l1_rel = l1_rel
        self.eta_burn = eta_burn
        self.D_d_av = D_d_av
        self.alpha_sum = alpha_sum
        self.n_pipe = n_pipe

        self.diffuser = None
        self.p_stag1 = None
        self.static1 = None
        self.T1 = None
        self.p1 = None
        self.rho1 = None
        self.G1 = None
        self.F_pipe = None
        self.volume_pipe = None
        self.l_pipe = None
        self.l1 = None
        self.d_pipe = None
        self.lam1 = None
        self.a_cr1 = None
        self.tau1 = None
        self.pi1 = None

    def compute(self):
        self.diffuser = Diffuser(
            work_fluid=self.work_fluid_in,
            sigma=self.sigma_d,
            G=self.G_in,
            p_stag_in=self.p_stag_in,
            T_stag_in=self.T_stag_in,
            F_in=self.F_in,
            c_out=self.c_d,
            G_fuel=self.G_fuel_in,
            D_out_av=self.D_d_av
        )
        self.diffuser.compute()
        self.p_stag1 = self.diffuser.p_stag_out * self.sigma_front
        self.static1 = self.diffuser.get_static(self.c1, self.T_stag_in, self.p_stag1, self.diffuser.k,
                                                self.work_fluid_in.R)
        self.T1 = self.static1[0]
        self.p1 = self.static1[1]
        self.rho1 = self.static1[2]
        self.lam1 = self.static1[3]
        self.a_cr1 = self.static1[4]
        self.tau1 = self.static1[5]
        self.pi1 = self.static1[6]
        self.G1 = self.alpha1 / self.alpha_sum * self.G_in
        self.F_pipe = self.G1 / (self.n_pipe * self.c1 * self.rho1)
        self.d_pipe = np.sqrt(4 * self.F_pipe / np.pi)
        self.volume_pipe = self.G_fuel * self.Q_n * self.eta_burn * 3600 / (self.H * self.p_stag_in * self.n_pipe)
        self.l_pipe = self.volume_pipe / self.F_pipe
        self.l1 = self.l1_rel * self.l_pipe

    def plot(self, figsize=(7, 7), fname=None):
        plt.figure(figsize=figsize)
        phi = np.linspace(0, 2 * np.pi, 100)
        plt.plot(self.diffuser.D_out_hub * np.cos(phi), self.diffuser.D_out_hub * np.sin(phi), lw=1.5, color='red')
        plt.plot(self.diffuser.D_out_per * np.cos(phi), self.diffuser.D_out_per * np.sin(phi), lw=1.5, color='red')
        plt.plot(self.diffuser.D_out_av * np.cos(phi), self.diffuser.D_out_av * np.sin(phi), lw=0.75, color='red',
                 ls='--')
        for i in range(self.n_pipe):
            x0 = self.diffuser.D_out_av * np.cos((2 * np.pi / self.n_pipe) * i)
            y0 = self.diffuser.D_out_av * np.sin((2 * np.pi / self.n_pipe) * i)
            plt.plot(x0 + self.d_pipe * np.cos(phi), y0 + self.d_pipe * np.sin(phi), lw=0.7, c='blue')
        plt.grid()
        if fname:
            plt.savefig(fname)
        plt.show()





