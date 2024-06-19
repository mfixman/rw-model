from __future__ import annotations

def get_alpha_mack(s: Individual, g: Group, sigma: float) -> float:
    return 1/2 * (1 + (2 * s.assoc - sigma))

def get_alpha_hall(s: Individual, g: Group, sigma: float, lamda: float) -> float:
    delta_ma_hall = s.delta_ma_hall or 0

    surprise = abs(lamda - sigma)
    window_term =  1 - g.xi_hall * math.exp(-delta_ma_hall**2 / 2)
    gamma = 0.99
    kayes = gamma*surprise +  (1-gamma)*s.alpha_hall

    new_error = kayes

    return new_error

def linearStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    s.alpha *= 1 + sign * 0.05
    s.assoc += s.alpha * beta * (g.prev_lamda - sigma)

def exponentialStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    if sign == 1:
        s.alpha *= (s.alpha ** 0.05) ** sign
    s.assoc += s.alpha * beta * (g.prev_lamda - sigma)

def mackStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    s.alpha_mack = g.get_alpha_mack(s, g, sigma)
    s.alpha = s.alpha_mack
    s.assoc += s.alpha * beta * (g.prev_lamda - sigma)

def hallStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    s.alpha_hall = g.get_alpha_hall(s, g, sigma, g.prev_lamda)
    s.alpha = s.alpha_hall
    s.assoc += .5 * s.alpha * abs(g.prev_lamda)

def macknhallStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    s.alpha_mack = g.get_alpha_mack(s, g, sigma)
    s.alpha_hall = g.get_alpha_hall(s, g, sigma, g.prev_lamda)
    s.alpha = (1 - abs(g.prev_lamda - sigma)) * s.alpha_mack + s.alpha_hall
    s.assoc += s.alpha * beta * (g.prev_lamda - sigma)

def dualVStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    # Ask Esther whether this is lamda^{n + 1) or lamda^n.
    rho = lamda - (sigmaE - sigmaI)

    if rho >= 0:
        s.Ve += g.betap * s.alpha * lamda
    else:
        s.Vi += g.betan * s.alpha * abs(rho)

    s.alpha = g.gamma * abs(rho) + (1 - g.gamma) * s.alpha
    s.assoc = s.Ve - s.Vi

def lepelleyStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    rho = lamda - (sigmaE - sigmaI)

    VXe = sigmaE - s.Ve
    VXi = sigmaI - s.Vi

    DVe = 0.
    DVi = 0.
    if rho >= 0:
        DVe = s.alpha * g.betap * (1 - s.Ve + s.Vi) * abs(rho)

        if rho > 0:
            s.alpha += -g.thetaE * (abs(lamda - s.Ve + s.Vi) - abs(lamda - VXe + VXi))
    else:
        DVi = s.alpha * g.betan * (1 - s.Vi + s.Ve) * abs(rho)
        s.alpha += -g.thetaI * (abs(abs(rho) - s.Vi + s.Ve) - abs(abs(rho) - VXi + VXe))

    s.alpha = min(max(s.alpha, 0.05), 1)
    s.Ve += DVe
    s.Vi += DVi

    s.assoc = s.Ve - s.Vi

def dualmackStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    rho = lamda - (sigmaE - sigmaI)

    VXe = sigmaE - s.Ve
    VXi = sigmaI - s.Vi

    if rho >= 0:
        s.Ve += s.alpha * g.betap * (1 - s.Ve + s.Vi) * abs(rho)
    else:
        s.Vi += s.alpha * g.betan * (1 - s.Vi + s.Ve) * abs(rho)

    s.alpha = 1/2 * (1 + s.assoc - (VXe - VXi))
    s.assoc = s.Ve - s.Vi

def hybridStep(s: Individual, g: Group, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    rho = lamda - (sigmaE - sigmaI)

    NVe = s.Ve
    NVi = s.Vi
    if rho >= 0:
        NVe = s.alpha_mack * s.Ve + g.betap * s.alpha_hall * lamda
    else:
        NVi = s.alpha_mack * s.Vi + g.betan * s.alpha_hall * abs(rho)

    VXe = sigmaE - s.Ve
    VXi = sigmaI - s.Vi
    if rho > 0:
        s.alpha_mack += -g.thetaE * (abs(lamda - s.Ve + s.Vi) - abs(lamda - VXe + VXi))
    elif rho < 0:
        s.alpha_mack += -g.thetaI * (abs(abs(rho) - s.Vi + s.Ve) - abs(abs(rho) - VXi + VXe))

    s.alpha_mack = min(max(s.alpha_mack, 0.05), 1)
    s.alpha_hall = g.gamma * abs(rho) + (1 - g.gamma) * s.alpha_hall

    s.Ve = NVe
    s.Vi = NVi

    s.assoc = s.Ve - s.Vi

StepFuncs = dict(
    linear = linearStep,
    exponential = exponentialStep,
    mack = mackStep,
    hall = hallStep,
    macknhall = macknhallStep,
    dualV = dualVStep,
    lepelley = lepelleyStep,
    dualmack = dualmackStep,
    hybrid = hybridStep,
)
