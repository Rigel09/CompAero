from CompAero.ConicalFlow import ConicalFlowRelations
from CompAero.PrandtlMeyer import PrandtlMeyer
from CompAero.RayleighFlowRelations import RayleighFlowRelations
from numpy.lib.twodim_base import diag
from CompAero.IsentropecRelations import IsentropicRelations
from CompAero.NormalShockRelations import NormalShockRelations
from CompAero.ObliqueShockRelations import ObliqueShockRelations
from CompAero.FannoFlowRelations import FannoFlowRelations

from math import radians, degrees, pi
import numpy as np
import matplotlib.pyplot as plt

import cProfile, pstats, io
from pstats import SortKey


if __name__ == "__main__":
    ##########################################################################
    #################### Isentropic Example  #################################
    ##########################################################################
    # temp = IsentropicRelations(gamma=1.4, a_aStar=1.0235)
    # print(temp)

    ##########################################################################
    #################### Normal Shock Example  ###############################
    ##########################################################################
    # flow = NormalShockRelations(gamma=1.4, m2=0.5130)
    # print()
    # print(flow)

    ##########################################################################
    #################### Oblique Shock Example  ##############################
    ##########################################################################
    # flow = ObliqueShockRelations(1.4, mach=3, shockAngle=60)
    # print()
    # print(flow)
    # flow.plotThetaBetaMachChart()

    ##########################################################################
    #################### Fanno Flow Example  #################################
    ##########################################################################
    # # Fanno Flow Example
    # t1 = 300
    # p1 = 1
    # po1 = IsentropicRelations.calc_p0_p(3, 1.4) * p1
    # flow = FannoFlowRelations(1.4, po_poSt= 4.23456790)
    # flow.applyPipeParameters(0.4, 11, 0.005)
    # print()
    # print(flow)
    # print()
    # print("T2: ", t1 * flow.t2_t1)
    # print("P2: ", p1 * flow.p2_p1)
    # print("Po2: ", po1 * flow.po2_po1)

    # # Plot length and diameter combinations to slow a mach 3 flow down to a mach 2.5 flow at the exit
    # f4ld = FannoFlowRelations.calc4FLSt_D_FromMach(1.5, 1.4)

    # diamEqn = lambda len : 1 / f4ld * 4 * flow.frictionCoeff * len

    # lengths = np.linspace(1e-5, 40, 100)
    # diams = np.array([diamEqn(len) for len in lengths])

    # plt.plot(lengths, diams)
    # plt.xlim(0, 40)
    # plt.ylim(0, diams[-1])
    # plt.xlabel("Pipe Length [ft]")
    # plt.ylabel("Pipe Diameter [ft]")
    # plt.title("Pipe Length Vs Diameter")
    # plt.grid()
    # plt.show()

    ##########################################################################
    #################### Rayleigh Flow Example  ##############################
    ##########################################################################
    # flow = RayleighFlowRelations(1.4, rho_rhoSt=0.5875)
    # print()
    # flow.simulateHeatAddition(260786.98, 275.2, 287)
    # print(flow)
    # flow.__preciscion = 4

    ##########################################################################
    #################### Prandtl Meyer Example  ##############################
    ##########################################################################
    # flow = PrandtlMeyer(1.4, dwnstreamNu=86.27, deflectionAngle=1)
    # print(flow)

    ##########################################################################
    ##################### Conical Flow Example  ##############################
    ##########################################################################
    flow = ConicalFlowRelations(1.4, mach=3, shockAngle=30)
    ans = flow.calculateConeFlowParameters(288.16, 287, [21, 25, 26, 24, 24.5, 29, flow.coneAngle])

    print(ans[flow.coneAngle].mach)
    print(ans[flow.coneAngle].tempRay_tempInit)
    print(ans[flow.coneAngle].pressRay_pressInit)
    print(ans[flow.coneAngle].rhoRay_rhoInit)
    print(ans[flow.coneAngle])

