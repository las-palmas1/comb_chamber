from core.geom import CombustionChamberGeom
import unittest
from jinja2 import Template, FileSystemLoader, Environment, select_autoescape
from gas_turbine_cycle.gases import Air
import numpy as np
import core.templates.tests
import core.templates


class TestCombustionChamberGeomTemplate(unittest.TestCase):
    def setUp(self):
        self.comb_geom = CombustionChamberGeom(
            work_fluid_in=Air(),
            Q_n=48e6,
            sigma_d=0.998,
            sigma_front=0.998,
            n_pipe=12,
            G_in=45,
            p_stag_in=15e5,
            T_stag_in=550,
            F_in=np.pi * 0.5**2 / 4,
            D_d_av=0.5,
            G_fuel_in=0,
            G_fuel=0.9,
            alpha_sum=3.5,
            alpha1=1.1,
            c_d=60,
            c1=12,
            H=4.5e6,
            l1_rel=0.6,
            eta_burn=0.99
        )
        self.comb_geom.compute()

    def test_template(self):
        loader = FileSystemLoader(
            [
                core.templates.__path__[0],
                core.templates.tests.__path__[0]
            ]
        )
        env = Environment(
            loader=loader,
            autoescape=select_autoescape(['tex']),
            block_start_string='</',
            block_end_string='/>',
            variable_start_string='<<',
            variable_end_string='>>',
            comment_start_string='<#',
            comment_end_string='#>'
        )

        template = env.get_template('test_comb_geom_templ.tex')
        content = template.render(
            comb_geom=self.comb_geom
        )

        with open('test_comb_geom.tex', 'w', encoding='utf-8') as f:
            f.write(content)



