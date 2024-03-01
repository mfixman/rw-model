import math
from collections import deque, defaultdict
import numpy as np
from itertools import combinations
from RW_strengths import Strengths

class Group:
    def __init__(self, name, alphas, betan, betap, lamda, cs = None, use_configurals = False, adaptive_type = None, window_size = None, xi_hall = None):
        # If some alphas don't appear, set their alpha to 0.
        if cs is not None:
            initial_alpha = 0.2 if adaptive_type is not 'macknhall' else 0
            alphas = {k: alphas.get(k, initial_alpha) for k in cs | alphas.keys()}

        self.name = name
        # self.alphas_copy = alphas.copy()
        # self.alphas = self.alphas_copy

        self.s = Strengths(cs)

        # self.delta_ma_hall is the difference of the moving averages used for Hall.
        # self.delta_ma_hall^n = V^n_{MA} - V^{n - 1}_{MA}
        # self.delta_ma_hall = {k: None for k in cs}
        self.xi_hall = xi_hall

        self.betan = betan
        self.betap = betap
        self.lamda = lamda
        self.use_configurals = use_configurals
        self.adaptive_type = adaptive_type
        self.window_size = window_size

        self.window = defaultdict(lambda: deque())

        # For simplicity, if we use_configurals and some compound stimuli don't
        # have a corresponding \alpha, then we calculate it as the product
        # of their simple stimuli.
        self.cs = [x for x in alphas.keys() if len(x) == 1]

        # Do something with strengths and use_configurals.
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

    def get_alpha_mack(self, cs):
        # This overflows -- ask Esther.
        # big_mack_error = sum(self.assoc.values()) / len(self.assoc) - self.assoc[cs]
        # return self.alpha_mack[cs] - big_mack_error
        #return 0
        """if len(self.assoc) > 1:
            big_mack_error = ( sum(self.assoc[x] for x in self.assoc if x!=cs)/ (len(self.assoc)-1) ) - self.assoc[cs]
        else:
            big_mack_error = - self.assoc[cs]
        return self.alpha_mack[cs] - 0.2*big_mack_error"""
        
        #return self.alpha_mack[cs] + 0.01*self.assoc[cs]

        return self.s[cs].alpha_mack + 0.01*(2*self.s[cs].assoc - self.s.Sigma())

    def get_alpha_hall(self, cs : str) -> None | float:
        if self.window_size is None:
            return None

        delta_ma_hall = self.s[cs].delta_ma_hall
        if delta_ma_hall is None:
            delta_ma_hall = 0

        try:
            error = self.xi_hall * self.s[cs].alpha_hall * math.exp(- delta_ma_hall**2 / 2)
        except:
            error = 0.0000001
        # print(f"alpha_hall: {error}")
        return error

    def compounds(self, part : str) -> set[str]:
        compounds = set(part)
        if self.use_configurals:
            compounds.add(part)

        return compounds

    def runPhase(self, parts : list[tuple[str, str]], phase_lamda : None | float) -> list[Strengths]:
        s_hist = [self.s.copy()]
        last_hist = {cs: 0 for cs in self.s.cs}

        for part, plus in parts:
            if plus == '+':
                beta, lamda, sign = self.betap, phase_lamda or self.lamda, 1
            else:
                beta, lamda, sign = self.betan, 0, -1

            compounds = self.compounds(part)
            sigma = sum(self.s[x].assoc for x in compounds)

            for cs in compounds:
                self.s[cs].alpha_mack = self.get_alpha_mack(cs)
                self.s[cs].alpha_hall = self.get_alpha_hall(cs)
                
                match self.adaptive_type:
                    case 'linear':
                        self.s[cs].alpha *= 1 + sign * 0.05
                    case 'exponential':
                        if sign == 1:
                            self.s[cs].alpha *= (self.s[cs].alpha ** 0.05) ** sign
                    case 'macknhall':
                        #print(lamda)
                        self.s[cs].alpha = (1-lamda+self.s[cs].assoc) * self.s[cs].alpha_mack + (lamda-self.s[cs].assoc) * self.s[cs].alpha_hall
                        """delta_v_n = self.s[cs].alpha_hall * beta * (lamda - sigma)
                        v_n = self.s[cs].alpha_mack * self.s[cs].assoc
                        self.s[cs].assoc = v_n + delta_v_n"""

                print(f"DV: {self.s[cs].alpha * beta * (lamda - sigma)}, Alpha: {self.s[cs].alpha}, sigma: {sigma}")
                self.s[cs].assoc += self.s[cs].alpha * beta * (lamda - sigma)

                # Redo this after being done with strengths.
                """
                if self.window_size is not None:
                    if len(self.s[cs].window) >= self.window_size:
                        self.s[cs].window.popleft()

                    self.s[cs].window.append(self.s[cs].assoc)

                    window_avg = sum(self.s[cs].window) / len(self.s[cs].window)

                    self.s[cs].delta_ma_hall = window_avg - V[cs][-1]

                    # I need to somehow adapt this to the separate Strengths
                    # V[cs].append(window_avg)
                    """

                last_hist[cs] += 1
                print(last_hist[cs])

                if len(s_hist) <= last_hist[cs]:
                    s_hist.append(Strengths({cs}, {cs: self.s[cs]}))
                elif cs in s_hist[last_hist[cs]].s:
                    s_hist[last_hist[cs]].s[cs] += self.s[cs]
                else:
                    s_hist[last_hist[cs]].s[cs] = self.s[cs]

            # s_hist.append(self.s.copy())

        return s_hist
