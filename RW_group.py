import math
from itertools import combinations
from RW_strengths import Strengths, History, Individual
from RW_adaptive_types import StepFuncs

class Group:
    name : str

    cs : list[str]
    s : Strengths

    betan : float
    betap : float
    lamda : float

    prev_lamda : float

    use_configurals : bool

    window_size : None | int
    xi_hall : None | float

    def __init__(self, name : str, alphas : dict[str, float], default_alpha : float, betan : float, betap : float, lamda : float, gamma : float, thetaE : float, thetaI : float, cs : None | set[str] = None, use_configurals : bool = False, adaptive_type : None | str = None, window_size : None | int = None, xi_hall : None | float = None):
        if cs is not None:
            alphas = {k: alphas.get(k, default_alpha) for k in cs | alphas.keys()}

        self.name = name

        self.s = Strengths(
            s = {
                k: Individual(assoc = 0, alpha = alphas[k])
                for k in alphas.keys()
            }
        )

        self.xi_hall = xi_hall

        self.betan = betan
        self.betap = betap
        self.lamda = lamda
        self.gamma = gamma
        self.thetaE = thetaE
        self.thetaI = thetaI

        self.step = StepFuncs[adaptive_type]

        self.use_configurals = use_configurals
        self.window_size = window_size

        self.prev_lamda = 0

        # For simplicity, if we use_configurals and some compound stimuli don't
        # have a corresponding \alpha, then we calculate it as the product
        # of their simple stimuli.
        self.cs = [x for x in alphas.keys() if len(x) == 1]

        # use_configurals is not working with strengths - add it later.
        '''
        if use_configurals:
            simples = {k: v for k, v in alphas.items() if len(k) == 1}
            compounds = {
                ''.join(p[0] for p in x): math.prod(p[1] for p in x)
                for x in sum(
                    [
                        list(combinations(simples.items(), t))
                        for t in range(1, len(simples) + 1)
                    ],
                    []
                )
            }
            self.alphas = compounds | alphas
            self.cs = self.alphas.keys()
        '''

    # compounds should probably be moved to Strengths.
    def compounds(self, part : str) -> set[str]:
        compounds = set(part)
        if self.use_configurals:
            compounds.add(part)

        return compounds

    # runPhase runs a single trial of a phase, in order, and returns a list of the Strength values
    # of its CS at every step.
    # It also modifies `self.s` to account for all the strengths modified in this phase.
    def runPhase(self, parts : list[tuple[str, str]], phase_lamda : None | float) -> list[Strengths]:
        hist = dict()

        for part, plus in parts:
            if plus == '+':
                beta, lamda, sign = self.betap, phase_lamda or self.lamda, 1
            else:
                beta, lamda, sign = self.betan, 0., -1

            compounds = self.compounds(part)

            sigma = sum(self.s[x].assoc for x in compounds)
            sigmaE = sum(self.s[x].Ve for x in compounds)
            sigmaI = sum(self.s[x].Vi for x in compounds)

            for cs in compounds:
                if cs not in hist:
                    hist[cs] = History()
                    hist[cs].add(self.s[cs])

                self.step(self.s[cs], self, beta, lamda, sign, sigma, sigmaE, sigmaI)

                if self.window_size is not None:
                    if len(self.s[cs].window) >= self.window_size:
                        self.s[cs].window.popleft()

                    self.s[cs].window.append(self.s[cs].assoc)
                    window_avg = sum(self.s[cs].window) / len(self.s[cs].window)

                    # delta_ma_hall is modified using the previous associated value.
                    self.s[cs].delta_ma_hall = window_avg - hist[cs].assoc[-1]

                hist[cs].add(self.s[cs])
            self.prev_lamda = lamda

        return Strengths.fromHistories(hist)
