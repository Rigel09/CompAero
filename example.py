"""
This module shows some general examples on using the CompAero python package
"""

# from CompAero.ConicalFlow import ConicalFlowRelations

import matplotlib.pyplot as plt
import numpy as np

from CompAero.fanno_flow_relations import FannoFlowRelations
from CompAero.isentropic_relations import IsentropicRelations
from CompAero.normal_shock_relations import NormalShockRelations
from CompAero.oblique_shock_relations import ObliqueShockRelations
from CompAero.prandtl_meyer import PrandtlMeyer
from CompAero.rayleigh_flow_relations import RayleighFlowRelations

if __name__ == "__main__":
    ##########################################################################
    #################### Isentropic Example  #################################
    ##########################################################################
    ir_flow = IsentropicRelations(gamma=1.4, mach=1.5)
    print(ir_flow)

    ##########################################################################
    #################### Normal Shock Example  ###############################
    ##########################################################################
    nsr_flow = NormalShockRelations(gamma=1.4, m2=0.5130)
    print(nsr_flow)

    ##########################################################################
    #################### Oblique Shock Example  ##############################
    ##########################################################################
    osr_flow = ObliqueShockRelations(1.4, mach=3, wedge_angle=20)
    print()
    print(osr_flow)
    osr_flow.plot_theta_beta_mach_chart()

    ##########################################################################
    #################### Fanno Flow Example  #################################
    ##########################################################################
    T1 = 300
    P1 = 1
    po1 = IsentropicRelations.calc_p0_p(3, 1.4) * P1
    flow = FannoFlowRelations(1.4, po_po_st=4.23456790)
    flow.apply_pipe_parameters(0.4, 11, 0.005)
    print()
    print(flow)
    print()
    print("T2: ", T1 * flow.t2_t1)
    print("P2: ", P1 * flow.p2_p1)
    print("Po2: ", po1 * flow.po2_po1)

    # Plot length and diameter combinations to slow a mach 3 flow down to mach 2.5 flow at the exit
    f4ld = FannoFlowRelations.calc_4flst_d(1.5, 1.4)

    def diam_eqn(pipe_len: float) -> float:
        """Calculated the corresponing pipe diameter"""
        return 1 / f4ld * 4 * flow.friction_coeff * pipe_len

    pipe_length_arr = np.linspace(1e-5, 40, 100)
    diams = np.array([diam_eqn(pl) for pl in pipe_length_arr])

    plt.plot(pipe_length_arr, diams)
    plt.xlim(0, 40)
    plt.ylim(0, diams[-1])
    plt.xlabel("Pipe Length [ft]")
    plt.ylabel("Pipe Diameter [ft]")
    plt.title("Pipe Length Vs Diameter")
    plt.grid()
    plt.show()

    ##########################################################################
    #################### Rayleigh Flow Example  ##############################
    ##########################################################################
    rff_flow = RayleighFlowRelations(1.4, mach=1.5)
    print()
    rff_flow.simulate_heat_addition(1000, 275.2, 287)
    print(rff_flow)

    ##########################################################################
    #################### Prandtl Meyer Example  ##############################
    ##########################################################################
    pm_flow = PrandtlMeyer(1.4, down_stream_nu=86.27, deflection_angle=1)
    print(pm_flow)
