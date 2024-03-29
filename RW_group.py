import math
import numpy as np
from collections import deque, defaultdict
from itertools import combinations
from RW_strengths import Strengths, History, Individual
from typing import List

class Group:
    def __init__(self, name, alphas, betan, betap, lamda, cs = None, use_configurals = False, adaptive_type = None, window_size = None, xi_hall = None):
        # If some alphas don't appear, set their alpha to 0.

        # Initially, alpha_mack and alpha_hall are identical.
        if cs is not None:
            initial_alpha = 0.5
            alphas = {k: alphas.get(k, initial_alpha) for k in cs | alphas.keys()}

        self.name = name

        self.s = Strengths(s = {
                k: Individual(assoc = 0, alpha = alphas[k])
                for k in alphas.keys()
            }
        )

        self.xi_hall = xi_hall

        self.betan = betan
        self.betap = betap
        self.lamda = lamda
        self.use_configurals = use_configurals
        self.adaptive_type = adaptive_type
        self.window_size = window_size

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

    def get_alpha_mack(self, cs, sigma):
        """if len(self.assoc) > 1:
            big_mack_error = ( sum(self.assoc[x] for x in self.assoc if x!=cs)/ (len(self.assoc)-1) ) - self.assoc[cs]
        else:
            big_mack_error = - self.assoc[cs]
        return self.alpha_mack[cs] - 0.2*big_mack_error"""
        
        return ((1 + (2*self.s[cs].assoc - sigma))/2)

    def get_alpha_hall(self, cs, sigma, lamda):
        # assert self.window_size is not None
        delta_ma_hall = self.s[cs].delta_ma_hall
        if delta_ma_hall is None:
            delta_ma_hall = 0

        error = 1/2 * (1 - abs(lamda - sigma) * self.s[cs].alpha_hall * (1 - self.xi_hall * math.exp(1/2 * -delta_ma_hall ** 2) + abs(lamda - sigma)))
        # try:
        # error = (((1-abs(lamda-sigma)) * self.s[cs].alpha_hall * (1-self.xi_hall*math.exp(-delta_ma_hall**2 / 2)))+abs(lamda-sigma))/2
        # except:
        # error = ((1-abs(lamda-sigma))*self.s[cs].alpha_hall + abs(lamda-sigma))/2

        return error

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
                beta, lamda, sign = self.betan, 0, -1

            compounds = self.compounds(part)
            sigma = sum(self.s[x].assoc for x in compounds)

            for cs in compounds:
                if cs not in hist:
                    hist[cs] = History()
                    hist[cs].add(self.s[cs])

                self.s[cs].alpha_mack = self.get_alpha_mack(cs, sigma)
                self.s[cs].alpha_hall = self.get_alpha_hall(cs, sigma, lamda)
                
                match self.adaptive_type:
                    case 'linear':
                        self.s[cs].alpha *= 1 + sign * 0.05
                    case 'exponential':
                        if sign == 1:
                            self.s[cs].alpha *= (self.s[cs].alpha ** 0.05) ** sign
                    case 'macknhall':
                        self.s[cs].alpha = (lamda - sigma) * self.s[cs].alpha_hall

                self.s[cs].assoc = self.s[cs].assoc * self.s[cs].alpha_mack + self.s[cs].alpha

                if self.window_size is not None:
                    if len(self.s[cs].window) >= self.window_size:
                        self.s[cs].window.popleft()

                    self.s[cs].window.append(self.s[cs].assoc)
                    window_avg = sum(self.s[cs].window) / len(self.s[cs].window)

                    # delta_ma_hall is modified using the previous associated value.
                    self.s[cs].delta_ma_hall = window_avg - hist[cs].assoc[-1]

                hist[cs].add(self.s[cs])

        return Strengths.fromHistories(hist)
