import numpy

from amuse.ext.fresco.blackbody import B_lambda

from amuse.units import (
    units,
)


base1 = [
    lambda x: x**0,
    lambda x: x,
    lambda x: x**2,
    lambda x: x**3,
    lambda x: x**4,
]

base2 = [
    lambda x: 1/(x+0.01)**3,
    lambda x: x**0.,
    lambda x: 1/(1.01-x)**3,
    lambda x: x**3,
    lambda x: x**4,
]

base3 = [
    lambda x: B_lambda(x, 5000. | units.K).number,
    lambda x: B_lambda(x, 20000. | units.K).number,
    lambda x: B_lambda(x, 2000. | units.K).number,
    lambda x: B_lambda(x, 3000. | units.K).number,
    lambda x: B_lambda(x, 10000. | units.K).number,
]

default_base = base1


class ColorConverter(object):
    def __init__(
            self, source, target,
            base=None, lmin=None, lmax=None, N=1000,
    ):
        if base is None:
            base = default_base
        if len(source) > len(base):
            raise Exception("provide enough base functions")
        if lmin is None:
            lmin = min(
                [min(x['wavelength']) for x in source]
                + [min(x['wavelength']) for x in target]
            )
        if lmax is None:
            lmax = max(
                [max(x['wavelength']) for x in source]
                + [max(x['wavelength']) for x in target]
            )
        self.lmin = lmin
        self.lmax = lmax
        self.dim = len(source)
        self.base = base[0:self.dim]
        self.source = source
        self.target = target
        self.dim2 = len(target)
        self.N = N

        larray = lmin+(lmax-lmin)*numpy.array(list(range(N+1)))/N

        Amatrix = numpy.zeros((self.dim, self.dim))

        for j, src in enumerate(source):
            xp = src['wavelength']
            fp = src['throughput']
            f = numpy.interp(
                larray.number,
                xp=xp.value_in(larray.unit),
                fp=fp,
                left=0.,
                right=0.,
            )
            for i, b in enumerate(self.base):
                bint = f*self.fbase(b)(larray)
                Amatrix[j, i] = numpy.trapz(bint, x=larray.number)

        Bmatrix = numpy.zeros((self.dim2, self.dim))

        for j, src in enumerate(target):
            xp = src['wavelength']
            fp = src['throughput']
            f = numpy.interp(
                larray.number,
                xp=xp.value_in(larray.unit),
                fp=fp,
                left=0.,
                right=0.,
            )
            for i, b in enumerate(self.base):
                bint = f*self.fbase(b)(larray)
                Bmatrix[j, i] = numpy.trapz(
                    bint,
                    x=larray.number,
                )

        self.larray = larray
        self.Amatrix = Amatrix
        self.Bmatrix = Bmatrix
        self.Ainv = numpy.linalg.inv(Amatrix)
        self.conversion_matrix = Bmatrix.dot(self.Ainv)

    def fbase(self, b):
        return lambda x: b((x-self.lmin)/(self.lmax-self.lmin))

    def convert(self, x):
        return self.conversion_matrix.dot(x)


class XYZ_to_sRGB_linear(object):
    def __init__(self):
        self.conversion_matrix = numpy.array(
            [
                [3.2406, -1.5372, -0.4986],
                [-0.9689, 1.8758, 0.0415],
                [0.0557, -0.2040, 1.0570],
            ],
        )

    def convert(self, x):
        return self.conversion_matrix.dot(x)


class sRGB_linear_to_sRGB(object):
    def __init__(self):
        self.a = 0.055

    def clip(self, x):
        return x.clip(0., 1.)

    def convert(self, x):
        x = self.clip(x)
        a = numpy.where(x <= 0.0031308)
        b = numpy.where(x > 0.0031308)
        x[a] = 12.92*x[a]
        x[b] = (1+self.a) * x[b]**(1/2.4) - self.a
        return x
