from __future__ import annotations
from itertools import combinations
from functools import reduce

class Individual:
    def __init__(self, assoc = 0, alpha = .2, alpha_mack = .2, alpha_hall = .2, delta_ma_hall = .2):
        self.assoc = assoc
        self.assoc = min(1., self.assoc)

        self.alpha = alpha
        self.alpha_mack = alpha_mack
        self.alpha_hall = alpha_hall

        self.delta_ma_hall = delta_ma_hall

    def __add__(self, other : Individual) -> Individual:
        ret = Individual()
        for prop in self.__dict__.keys():
            setattr(ret, prop, self.__dict__[prop] + other.__dict__[prop])

        ret.assoc = min(1., self.assoc)

        return ret

    def __mul__(self, other : Individual) -> Individual:
        ret = Individual()
        for prop in self.__dict__.keys():
            setattr(ret, prop, self.__dict__[prop] + other.__dict__[prop])

        ret.assoc = min(1., self.assoc)

        return ret

    def __repr__(self):
        return f'{self.assoc:g}'

    def copy(self) -> Individual:
        ret = Individual()
        for key, val in self.__dict__.items():
            setattr(ret, key, val)

        return ret

class Strengths:
    # Add a way to initialise this with set values.
    def __init__(self, cs : set[str], s : None | dict[str, Individual] = None):
        if s is None:
            s = {k: Individual() for k in cs}

        self.cs = cs
        self.s = s

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
    def __getitem__(self, key):
        assert len(set(key)) == len(key)
        return reduce(lambda a, b: a * b, [self.s.get(k, 0) for k in key])

    def __add__(self, other):
        cs = self.cs | other.cs
        return Strengths(cs, {k: self[k] + other[k] for k in cs})

    def __repr__(self):
        return repr(self.s)

    def copy(self):
        return Strengths(self.cs.copy(), {k: v.copy() for k, v in self.s.items()})

    def Sigma(self):
        return sum(x.assoc for x in self.s.values())
