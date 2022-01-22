from math import pow, log, sqrt


def thrust_coefficient(gamma: float, areaRatio: float, pe_pc: float, pa_pc: float) -> float:
    """ Calculates the thrust coefficient Cf

    Args:
        gamma (float): Ratio of specific heats
        areaRatio (float): Area ratio of the nozzle Ae/At
        pe_pc (float): Ratio of exit pressure to chamber pressure (exit / chamber)
        pa_pc (float): Ratio of ambient pressure to chamber pressure (ambient / chamber)

    Returns:
        float: thrust coefficient Cf
    """
    if areaRatio < 1.0:
        raise ValueError("Area ratio cannot be less than 1")

    gp1 = gamma + 1
    gm1 = gamma - 1
    # sqrt
    first = 2 * pow(gamma, 2) / gm1
    second = pow(2 / gp1, gp1 / gm1)
    third = 1 - pow(pe_pc, gm1 / gamma)
    return sqrt(first * second * third) + (pe_pc - pa_pc) * areaRatio


def max_thrust_coefficient(gamma: float, areaRatio: float, pe_pc: float) -> float:
    """ Calculates the maximum thrust coefficient for a nozzle flow Cf_max

    Args:
        gamma (float): Ratio of specific heats
        areaRatio (float): Area ratio of the nozzle Ae/At
        pe_pc (float): Ratio of exit pressure to chamber pressure (exit / chamber)

    Returns:
        float: Max thrust coefficient Cf
    """
    return thrust_coefficient(gamma, areaRatio, pe_pc, pe_pc)


def min_thrust_coefficient(areaRatio: float) -> float:
    """ Calculates the minimum thrust coefficient for a nozzle. Values below this have a high probability of flow seperation in the nozzle

    Args:
        areaRatio (float): Area ratio of the nozzle Ae/At

    Returns:
        float: Minimum thrust coefficient Cf
    """
    if areaRatio < 1.0:
        raise ValueError("Area ratio cannot be less than 1")

    return -0.0445 * pow(log(areaRatio), 2) + 0.5324 * log(areaRatio) + 0.1843
