

from abc import ABC, abstractmethod
from math import sqrt, inf, exp, pi

import numpy as np


GASCONST = 8.31446261815324e-3  # kJ per mol Kelvin


class Species(ABC):
    name = None                 # string containing species name

    def __init__(self, symb, name):
        self.symb, self.name = symb, name

    def __str__(self):
        return "{} ({})".format(self.name, self.symb)

    @abstractmethod
    def effective_diffusion(self, T):
        pass

class KnudsenSpecies(Species):

    # put the species-specific parameters last so we can curry the constructor
    def __init__(self, Dpore, epspore, tau, symb, name, Mi):
        super().__init__(symb, name)
        self.Dpore, self.epspore, self.tau, self.Mi = Dpore, epspore, tau, Mi

    def __repr__(self):
        return (
            "KnudsenSpecies(Dpore = {}, epspore = {}, tau = {}, symb = {}, "
            "name = {}, Mi = {})"
        ).format(
            self.Dpore, self.epspore, self.tau, self.symb, self.name, self.Mi
        )

    def effective_diffusion(self, T):
        # input units are nm for Dpore, kJ per mol K for GASCONST, K for T, and
        # g per mod for Mi; output units are square micrometers per second
        return 1e6 / 3 * self.Dpore * sqrt(
            8 * GASCONST * T / pi / self.Mi
        ) * self.epspore / self.tau

class Reaction(ABC):
    src, dst = None, None       # names of species

    def __init__(self, src, dst):
        self.src, self.dst = src, dst

    def __str__(self):
        return "{} -> {}".format(self.src, self.dst)

    @abstractmethod
    def kij(self, T):
        pass

class ArrheniusReaction(Reaction):

    def __init__(self, src, dst, A, b, Ea, T0):
        super().__init__(src, dst)
        self.A, self.b, self.Ea, self.T0 = A, b, Ea, T0

    def __repr__(self):
        return (
            "ArrheniusReaction(src = {}, dst = {}, A = {}, b = {}, Ea = {}, "
            "T0 = {})"
        ).format(
            self.src, self.dst, self.A, self.b, self.Ea, self.T0
        )

    def kij(self, T):
        A = self.A
        if self.b is not None and self.b != 0.0:
            #TODO we need to be more carful about this; this only works if A is primitive
            A *= T**self.b
        if self.T0 is None or self.T0 == inf or self.T0 == -inf:
            T0inv = 0.0
        else:
            T0inv = self.T0**(-1)
        return A * exp(-self.Ea / GASCONST * (1 / T - T0inv))

class Mechanism(object):

    spcs = None             # list of species objects
    rxns = None             # list of reaction objects
    numactive = None        # ignore products

    # providing an explicit constructor prevents the common error of
    # overwriting default values
    def __init__(self, spcs, rxns):
        self.spcs, self.rxns = spcs, rxns

    def reversible(self):
        edges = [(rxn.src, rxn.dst) for rxn in self.rxns]
        def findroot(es):
            ns = list(set([n for e in es for n in e]))
            for n in ns:
                ss = [s for s, d in es if n == d]
                if len(ss) == 0:
                    return n
            return None
        root = findroot(edges)
        while root is not None:
            edges = [e for e in edges if root not in e]
            root = findroot(edges)
        return len(edges) == 0

    def findnonreacting(self, verb):
        nonreacting = [
            spc.symb for spc in self.spcs
            if not any([
                spc.symb == rxn.src or spc.symb == rxn.dst
                for rxn in self.rxns
            ])
        ]
        if verb:
            if len(nonreacting):
                print("Some species are inert, they are: {}".format(
                    ", ".join(nonreacting)
                ))
            else:
                print("All species are involved in at least one reaction.")
        return nonreacting

# TODO this should be re-enabled after we deal with how sorting is going to
# work
#   def findproducts(self, verb):
#       products = [
#           spc.symb for spc in self.spcs
#           if not any([spc.symb == rxn.src for rxn in self.rxns])
#       ]
#       self.spcs.sort(key=lambda spc: spc.symb in products)
#       self.numactive = len(self.spcs) - len(products)
#       if verb:
#           print((
#               "Found {} product species:\n\t{}"
#           ).format(len(products), ", ".join(products)))
#       return products

    def getnumactive(self):
        if self.numactive is None:
            self.findproducts(False)
        return self.numactive

    def findundefined(self, verb):
        spcnames = [spc.symb for spc in self.spcs]
        undef = []
        for rxn in self.rxns:
            if rxn.src not in spcnames:
                undef.append(rxn.src)
                if verb:
                    print("Source '{}' not in species list.".format(rxn.src))
            if rxn.dst not in spcnames:
                undef.append(rxn.dst)
                if verb:
                    print(
                        "Destination '{}' not in species list.".format(rxn.dst)
                    )
        return undef

    def getkijs(self, T):
        symbols = [spc.symb for spc in self.spcs]
        kijs = np.zeros((len(self.spcs), len(self.spcs)))
        for rxn in self.rxns:
            i, j = symbols.index(rxn.src), symbols.index(rxn.dst)
            kijs[i, j] = rxn.kij(T)
        return kijs

    def getDis(self, T):
        Dis = np.empty((len(self.spcs),))
        for i, spcs in enumerate(self.spcs):
            Dis[i] = spcs.effective_diffusion(T)
        return Dis

    def getmatrix(self, T):
        kijs, Dis = self.getkijs(T), self.getDis(T)
        consumption = np.diag(np.sum(
            kijs / np.tile(Dis.reshape((-1,1)), len(self.spcs)),
            axis=1
        ))
        production = kijs.T / np.tile(Dis.reshape((-1,1)), len(self.spcs))
        B = consumption - production  # yes, this is a bad sign convention
        return B


