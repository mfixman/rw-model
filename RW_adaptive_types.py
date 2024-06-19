StepFuncs = dict(
    linear = linearStep
    exponential = exponentialStep
    mack = mackStep
    hall = hallStep
    macknhall = macknhallStep
    dualV = dualVStep
    lepelley = lepelleyStep
    dualmack = dualmackStep
    hybrid = hybridStep
)

def linearStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    self.s[cs].alpha *= 1 + sign * 0.05
    self.s[cs].assoc += self.s[cs].alpha * beta * (self.prev_lamda - sigma)

def exponentialStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    if sign == 1:
        self.s[cs].alpha *= (self.s[cs].alpha ** 0.05) ** sign
    self.s[cs].assoc += self.s[cs].alpha * beta * (self.prev_lamda - sigma)

def mackStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    self.s[cs].alpha_mack = self.get_alpha_mack(cs, sigma)
    self.s[cs].alpha = self.s[cs].alpha_mack
    self.s[cs].assoc += self.s[cs].alpha * beta * (self.prev_lamda - sigma)

def hallStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    self.s[cs].alpha_hall = self.get_alpha_hall(cs, sigma, self.prev_lamda)
    self.s[cs].alpha = self.s[cs].alpha_hall
    self.s[cs].assoc += .5 * self.s[cs].alpha * abs(self.prev_lamda)

def macknhallStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    self.s[cs].alpha_mack = self.get_alpha_mack(cs, sigma)
    self.s[cs].alpha_hall = self.get_alpha_hall(cs, sigma, self.prev_lamda)
    self.s[cs].alpha = (1 - abs(self.prev_lamda - sigma)) * self.s[cs].alpha_mack + self.s[cs].alpha_hall
    self.s[cs].assoc += self.s[cs].alpha * beta * (self.prev_lamda - sigma)

def dualVStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    # Ask Esther whether this is lamda^{n + 1) or lamda^n.
    rho = lamda - (sigmaE - sigmaI)

    if rho >= 0:
        self.s[cs].Ve += self.betap * self.s[cs].alpha * lamda
    else:
        self.s[cs].Vi += self.betan * self.s[cs].alpha * abs(rho)

    self.s[cs].alpha = self.gamma * abs(rho) + (1 - self.gamma) * self.s[cs].alpha
    self.s[cs].assoc = self.s[cs].Ve - self.s[cs].Vi

    print(f'{cs}:\t𝛒 = {rho: .3f}; Ve = {self.s[cs].Ve:.3f}; Vi = {self.s[cs].Vi:.3f}')

def lepelleyStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    rho = lamda - (sigmaE - sigmaI)

    VXe = sigmaE - self.s[cs].Ve
    VXi = sigmaI - self.s[cs].Vi

    DVe = 0.
    DVi = 0.
    if rho >= 0:
        DVe = self.s[cs].alpha * self.betap * (1 - self.s[cs].Ve + self.s[cs].Vi) * abs(rho)

        if rho > 0:
            self.s[cs].alpha += -self.thetaE * (abs(lamda - self.s[cs].Ve + self.s[cs].Vi) - abs(lamda - VXe + VXi))
    else:
        DVi = self.s[cs].alpha * self.betan * (1 - self.s[cs].Vi + self.s[cs].Ve) * abs(rho)
        self.s[cs].alpha += -self.thetaI * (abs(abs(rho) - self.s[cs].Vi + self.s[cs].Ve) - abs(abs(rho) - VXi + VXe))

    self.s[cs].alpha = min(max(self.s[cs].alpha, 0.05), 1)
    self.s[cs].Ve += DVe
    self.s[cs].Vi += DVi

    self.s[cs].assoc = self.s[cs].Ve - self.s[cs].Vi

def dualmackStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    rho = lamda - (sigmaE - sigmaI)

    VXe = sigmaE - self.s[cs].Ve
    VXi = sigmaI - self.s[cs].Vi

    if rho >= 0:
        self.s[cs].Ve += self.s[cs].alpha * self.betap * (1 - self.s[cs].Ve + self.s[cs].Vi) * abs(rho)
    else:
        self.s[cs].Vi += self.s[cs].alpha * self.betan * (1 - self.s[cs].Vi + self.s[cs].Ve) * abs(rho)

    self.s[cs].alpha = 1/2 * (1 + self.s[cs].assoc - (VXe - VXi))
    self.s[cs].assoc = self.s[cs].Ve - self.s[cs].Vi

def hybridStep(cs: str, beta: float, lamda: float, sign: int, sigma: float, sigmaE: float, sigmaI: float):
    rho = lamda - (sigmaE - sigmaI)

    NVe = self.s[cs].Ve
    NVi = self.s[cs].Vi
    if rho >= 0:
        NVe = self.s[cs].alpha_mack * self.s[cs].Ve + self.betap * self.s[cs].alpha_hall * lamda
    else:
        NVi = self.s[cs].alpha_mack * self.s[cs].Vi + self.betan * self.s[cs].alpha_hall * abs(rho)

    VXe = sigmaE - self.s[cs].Ve
    VXi = sigmaI - self.s[cs].Vi
    if rho > 0:
        self.s[cs].alpha_mack += -self.thetaE * (abs(lamda - self.s[cs].Ve + self.s[cs].Vi) - abs(lamda - VXe + VXi))
    elif rho < 0:
        self.s[cs].alpha_mack += -self.thetaI * (abs(abs(rho) - self.s[cs].Vi + self.s[cs].Ve) - abs(abs(rho) - VXi + VXe))

    self.s[cs].alpha_mack = min(max(self.s[cs].alpha_mack, 0.05), 1)
    self.s[cs].alpha_hall = self.gamma * abs(rho) + (1 - self.gamma) * self.s[cs].alpha_hall

    self.s[cs].Ve = NVe
    self.s[cs].Vi = NVi

    self.s[cs].assoc = self.s[cs].Ve - self.s[cs].Vi

    print(f'{cs}:\t𝛒 = {rho: .3f}; Ve = {self.s[cs].Ve:.3f}; Vi = {self.s[cs].Vi:.3f}')
