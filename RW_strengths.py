from __future__ import annotations
from itertools import combinations
from functools import reduce

class Individual:
    assoc : float
    alpha : float
    alpha_mack : float
    alpha_hall : float

    window : list[float]
    delta_ma_hall : float

    def __init__(self, assoc = 0., alpha = .2, alpha_mack = .2, alpha_hall = .2, delta_ma_hall = .2, window = []):
        self.assoc = min(1., assoc)

        self.alpha = alpha
        self.alpha_mack = alpha_mack
        self.alpha_hall = alpha_hall

        self.window = window
        self.delta_ma_hall = delta_ma_hall

    def join(self, other : Individual, op) -> Individual:
        ret = {}
        for prop in self.__dict__.keys():
            this = getattr(self, prop)
            that = getattr(other, prop)

            if type(this) is float or type(this) is int:
                ret[prop] = op(this, that)
            elif type(this) is list:
                assert len(this) == len(that)
                ret[prop] = [op(a, b) for a, b in zip(this, that)]
            else:
                raise ValueError(f'Unknown type {type(this)} for {prop}, which is equal to {this} and {that}')

        return Individual(**ret)

    def __add__(self, other : Individual) -> Individual:
        return self.join(other, lambda a, b: a + b)

    def __mul__(self, other : Individual) -> Individual:
        return self.join(other, lambda a, b: a * b)

    def __repr__(self):
        return f'{self.assoc:g}'

    def copy(self) -> Individual:
        return Individual(**self.__dict__)

class History:
    hist : list[Individual]

    def __init__(self, ind : None | Individual = None):
        self.hist = []

        if ind is not None:
            self.add(ind)

    def add(self, ind : Individual):
        self.hist.append(ind.copy())

    def __getattr__(self, key):
        return [getattr(p, key) for p in self.hist]

class Strengths:
    cs : set[str]
    s : dict[str, Individual]

    def __init__(self, cs : None | set[str] = None, s : None | dict[str, Individual] = None):
        if s is None and cs is not None:
            s = {k: Individual() for k in cs}

        if cs is None and s is not None:
            cs = set(s.keys())

        if cs is None or s is None:
            raise ValueError('Either cs or s have to have a value')

        self.cs = set(cs)
        self.s = dict(s)

    @staticmethod
    def fromHistories(histories : dict[str, History]) -> list[Strengths]:
        longest = max(len(x.hist) for x in histories.values())
        return [
            Strengths(
                s = {
                    cs: h.hist[i]
                    for cs, h in histories.items()
                    if len(h.hist) > i
                }
            )
            for i in range(longest)
        ]

    def combined_cs(self) -> set[str]:
        h = set()

        simples = {k: v for k, v in self.s.items() if len(k) == 1}
        compounds = {k: v for k, v in self.s.items() if len(k) > 1}
        for size in range(1, len(self.s) + 1):
            for comb in combinations(simples.items(), size):
                names = [x[0] for x in comb]
                assocs = [x[1] for x in comb]
                h.add(''.join(names))

        for k, v in compounds.items():
            h.add(k)
            h.add(f'c({k})')

        return h

    # Get the individual values of either a single key (len(key) == 1), or
    # the combined values of a combination of keys (product of values).
    def __getitem__(self, key : str) -> Individual:
        assert len(set(key)) == len(key)
        return reduce(lambda a, b: a + b, [self.s[k] for k in key])

    def __add__(self, other : Strengths) -> Strengths:
        cs = self.cs | other.cs
        return Strengths(cs, {k: self[k] + other[k] for k in cs})

    def __repr__(self):
        return repr(self.s)

    def copy(self) -> Strengths:
        return Strengths(self.cs.copy(), {k: v.copy() for k, v in self.s.items()})

    def Sigma(self):
        return sum(x.assoc for x in self.s.values())
