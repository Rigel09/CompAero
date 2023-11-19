# from CompAero.ConicalFlow import ConicalFlowRelations
from math import degrees, pi, radians

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
    temp = IsentropicRelations(gamma=1.4, mach=1.5)
    print(temp)

    ##########################################################################
    #################### Normal Shock Example  ###############################
    ##########################################################################
    flow = NormalShockRelations(gamma=1.4, m2=0.5130)
    print(flow)

    ##########################################################################
    #################### Oblique Shock Example  ##############################
    ##########################################################################
    flow = ObliqueShockRelations(1.4, mach=3, wedge_angle=60)
    print()
    print(flow)
    flow.plot_theta_beta_mach_chart()

    ##########################################################################
    #################### Fanno Flow Example  #################################
    ##########################################################################
    t1 = 300
    p1 = 1
    po1 = IsentropicRelations.calc_p0_p(3, 1.4) * p1
    flow = FannoFlowRelations(1.4, po_po_st=4.23456790)
    flow.apply_pipe_parameters(0.4, 11, 0.005)
    print()
    print(flow)
    print()
    print("T2: ", t1 * flow.t2_t1)
    print("P2: ", p1 * flow.p2_p1)
    print("Po2: ", po1 * flow.po2_po1)

    # Plot length and diameter combinations to slow a mach 3 flow down to a mach 2.5 flow at the exit
    f4ld = FannoFlowRelations.calc_4flst_d(1.5, 1.4)

    diamEqn = lambda len: 1 / f4ld * 4 * flow.friction_coeff * len

    lengths = np.linspace(1e-5, 40, 100)
    diams = np.array([diamEqn(len) for len in lengths])

    plt.plot(lengths, diams)
    plt.xlim(0, 40)
    plt.ylim(0, diams[-1])
    plt.xlabel("Pipe Length [ft]")
    plt.ylabel("Pipe Diameter [ft]")
    plt.title("Pipe Length Vs Diameter")
    plt.grid()
    plt.show()

    flow = RayleighFlowRelations(1.4, 0.5)
    print(flow)

    ##########################################################################
    #################### Rayleigh Flow Example  ##############################
    ##########################################################################
    flow = RayleighFlowRelations(1.4, mach=1.5)
    print()
    flow.simulate_heat_addition(1000, 275.2, 287)
    print(flow)

    ##########################################################################
    #################### Prandtl Meyer Example  ##############################
    ##########################################################################
    flow = PrandtlMeyer(1.4, down_stream_nu=86.27, deflection_angle=1)
    print(flow)

    ##########################################################################
    ##################### Conical Flow Example  ##############################
    ##########################################################################
    #
    # Not Currently Supported
    #

    # flow = ConicalFlowRelations(1.4, mach=3, shock_angle=30)
    # ans = flow.calculateConeFlowParameters(288.16, 287, [21, 25, 26, 24, 24.5, 29, flow.coneAngle])

    # print(ans[flow.coneAngle].mach)
    # print(ans[flow.coneAngle].tempRay_tempInit)
    # print(ans[flow.coneAngle].pressRay_pressInit)
    # print(ans[flow.coneAngle].rhoRay_rhoInit)
    # print(ans[flow.coneAngle])
