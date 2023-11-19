from math import radians

from CompAero.oblique_shock_relations import ObliqueShockRelations as OSR

if __name__ == "__main__":
    # Relationships across an oblique shock can be determined by passing combinations of known values to the constructor
    # Valid combinations (from the documentation at /docs/ObliqueShockRelations.html)
    """
    Valid_Combinations_of_Parameters:
        1: gamma, mach, shock angle\n
        2. gamma, mach, wedge angle\n
        2: gamma, mach normal 1, shock angle \n
        3: gamma, mach behind shock, wedge angle, shock angle\n
        4: gamma, shock angle, wedge angle\n
        5: gamma, P2/P1, shock angle\n
        6: gamma, Rho2/Rho1, shock angle\n
        7: gamma, T2/T1, shock angle\n
        8: gamma, P02/P01, shock angle\n
        9: gamma, P02/P1, shock angle\n
        10: gamma, dwn_strm_mach, shock angle\n
    """

    # Demonstrations of each of these constructors

    # Mach 3 at a flow deflection of 25
    # This solves the TBM Equation for the shock angle
    flow = OSR(1.4, mach=3, wedge_angle=25)
    # print(flow)

    # Mach 3 with a shock angle of 45
    # This solves the TBM Equation for the flow deflection angle
    flow = OSR(1.4, mach=3, shock_angle=45)
    # print(flow)

    # Angles can also be passed in in radians
    flow = OSR(1.4, mach=3, shock_angle=radians(45), use_degrees=False)
    # print(flow)

    # Determine flow from the normal component of the flow ahead of the shcok and the shock angle
    flow = OSR(1.4, mn1=1.8, shock_angle=45)
    # print(flow)

    # Flow State can be determined by knowing the shock angle and wedge angle as well
    # This solves the TBM Equation for the mach number
    flow = OSR(1.4, shock_angle=43, wedge_angle=22)
    # print(flow)

    # Flow state can be determined by knowing a condition across the shock and the shock angle
    # This is availble for any of the states across the shock such as T2/T1
    flow = OSR(1.4, po2_po1=0.76, shock_angle=42)
    # print(flow)
