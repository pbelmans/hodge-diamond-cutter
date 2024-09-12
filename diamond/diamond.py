r"""
A tool to work with Hodge diamonds, comes with many varieties and constructions
built into it.

* `repository <https://github.com/pbelmans/hodge-diamond-cutter>`_
* `documentation <https://cutter.ncag.info>`_

Hodge diamonds encode the Hodge numbers of a variety, and provide interesting
information about its structure. They provide a numerical incarnation of many
operations one can perform in algebraic geometry, such as blowups, projective
bundles, products. They are also computed for many more specific constructions
such as certain moduli spaces of sheaves, or flag varieties, ...

These Hodge numbers are defined as the dimensions of the sheaf cohomology of
exterior powers of the cotangent bundle, i.e.

.. MATH::

    \mathrm{h}^{p,q}(X)=\dim\mathrm{H}^q(X,\Omega_X^p)

Here $p$ and $q$ range from $0$ to $n=\\dim X$. These numbers satisfy additional
symmetry properties:

    * Hodge symmetry: $\\mathrm{h}^{p,q}(X)=\\mathrm{h}^{q,p}(X)$
    * Serre duality: $\\mathrm{h}^{p,q}(X)=\\mathrm{h}^{n-p,n-q}(X)$

Because of these symmetries they are usually displayed as a diamond (it's
really just a square tilted 45 degrees), so that for a surface it would be::

                        h^{2,2}
                h^{2,1}         h^{1,2}
        h^{2,0}         h^{1,1}         h^{0,2}
                h^{1,0}         h^{0,1}
                        h^{0,0}

One of their famous applications is the mirror symmetry prediction that every
Calabi-Yau 3-fold has a mirror Calabi-Yau threefold, which should imply that
their Hodge diamonds are transpositions. The first instance of this is the
quintic 3-fold and its mirror, whose Hodge diamonds are::

                    1
               0         0
          0         1         0
      1       101       101       1
          0         1         0
               0         0
                    1

and::

                   1
              0         0
          0       101       0
      1       1         1       1
          0       101       0
              0         0
                   1


The following are some very basic examples of operations and constructions one
can use within the Hodge diamond cutter. To get started we do::

    sage: from diamond import *

after starting Sage.

Pretty print the Hodge diamond of a genus 2 curve::

    sage: X = HodgeDiamond.from_matrix([[1, 2], [2, 1]])
    sage: print(X)
          1
      2       2
          1

Compute the Euler characteristic of the product of `X` with itself::

    sage: print((X*X).euler())
    4

Pretty print the Hodge diamond of the Hilbert square of a K3 surface::

    sage: S = HodgeDiamond.from_matrix([[1, 0, 1], [0, 20, 0], [1, 0, 1]])
    sage: print(hilbn(S, 2))
                        1
                   0         0
              1        21        1
          0        0         0        0
      1       21       232       21       1
          0        0         0        0
              1        21        1
                   0         0
                        1

It also possible to generate LaTeX code::

    sage: latex(K3().pprint())
    \begin{tabular}{ccccc}
     &  & $1$ &  &  \\
     & $0$ &  & $0$ &  \\
    $1$ &  & $20$ &  & $1$ \\
     & $0$ &  & $0$ &  \\
     &  & $1$ &  &  \\
    \end{tabular}

There are many varieties built in, e.g. the previously defined K3 surface can
be compared to the built-in one::

    sage: print(S == K3())
    True

Check out the `documentation <https://cutter.ncag.info>`_
for all the available functionality.

If you use the software in your work, please cite it as explained on
`Zenodo <https://doi.org/10.5281/zenodo.3893509>`_, or in
`the README file <https://github.com/pbelmans/hodge-diamond-cutter?tab=readme-ov-file#how-to-cite>`_.

AUTHORS:

- Pieter Belmans (2019-01-27): initial version
- Pieter Belmans (2020-06-16): the version which got assigned a DOI
- Pieter Belmans (2021-08-04): various additions, added unit tests and proper
  documentation
"""

from itertools import groupby

# ****************************************************************************
#       Copyright (C) 2021 Pieter Belmans <pieterbelmans@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.arith.misc import binomial, factorial, gcd
from sage.categories.cartesian_product import cartesian_product
from sage.categories.rings import Rings
from sage.combinat.composition import Compositions
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.partition import Partitions
from sage.combinat.q_analogues import q_binomial
from sage.combinat.root_system.dynkin_diagram import DynkinDiagram
from sage.combinat.subset import Subsets
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.misc.cachefunc import cached_function
from sage.misc.fast_methods import Singleton
from sage.misc.misc_c import prod
from sage.misc.table import table
from sage.modules.free_module_element import vector
from sage.rings.function_field.constructor import FunctionField
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.element import Element
from sage.structure.parent import Parent


def _to_matrix(f):
    r"""
    Convert Hodge--Poincaré polynomial to matrix representation

    EXAMPLES::

        sage: from diamond import *
        sage: x,y = polygens(ZZ,'x,y')
        sage: f = 1+x**2+4*x*y+x**2*y**2+y**2
        sage: diamond._to_matrix(f)
        [1 0 1]
        [0 4 0]
        [1 0 1]
    """
    assert f in HodgeDiamond.R

    if f.is_zero():
        return matrix(ZZ, 1, 1, [0])

    # deal with the size of the diamond in this way because of the following example:
    # X = complete_intersection(5, 3)
    # X*X - hilbtwo(X)
    d = max(max(e) for e in f.exponents()) + 1
    return matrix(ZZ, d, d, lambda i, j: f[i, j])


class HodgeDiamond(Element):
    r"""
    This class implements some methods to work with Hodge diamonds.
    """

    #: polynomial ring used internally for the Hodge-Poincaré polynomial
    #: and can be used externally to create new polynomials (and thus diamonds)
    R = PolynomialRing(ZZ, ("x", "y"))
    #: variables in the polynomial ring for Hodge-Poincaré polynomials
    x, y = R.gens()

    #: static configuration variable for pretty-printing, see :func:`HodgeDiamond.pprint`
    hide_zeroes = None
    #: static configuration variable for pretty-printing, see :func:`HodgeDiamond.pprint`
    quarter = None

    def __init__(self, parent, m):
        r"""
        Constructor for a Hodge diamond if you know what you are doing

        This just uses the given matrix. It is probably advised to use the
        class methods

        * :meth:`HodgeDiamond.from_matrix`
        * :meth:`HodgeDiamond.from_polynomial`

        in most cases though.
        """
        # matrix representation of the Hodge diamond is used internally
        self._m = m
        self._size = m.ncols() - 1
        Element.__init__(self, parent)

    @classmethod
    def from_matrix(cls, m, from_variety=False):
        r"""
        Construct a Hodge diamond from a matrix

        INPUT:

        - ``m`` -- square integer matrix representing a Hodge diamond

        - ``from_variety`` (default: False) -- whether a check should be
          performed that it comes from a variety

        EXAMPLES:

        Hodge diamond of a K3 surface::

            sage: from diamond import *
            sage: S = HodgeDiamond.from_matrix([[1, 0, 1], [0, 20, 0], [1, 0, 1]])
            sage: S == K3()
            True

        The following fails as the lack of symmetry prevents a geometric origin::

            sage: HodgeDiamond.from_matrix([[1, 2], [0, 1]], from_variety=True)
            Traceback (most recent call last):
            ...
            AssertionError: The matrix does not
                            satisfy the conditions satisfied by the Hodge diamond of a
                            smooth projective variety.
        """
        return HodgeDiamondRing()(m, from_variety=from_variety)

    @classmethod
    def from_polynomial(cls, f, from_variety=False):
        r"""
        Construct a Hodge diamond from a Hodge--Poincaré polynomial

        INPUT:

        - ``f`` -- an integer polynomial in the ring ``HodgeDiamond.R``
          representing the Hodge--Poincaré polynomial

        - ``from_variety`` (default: False) -- whether a check should be
          performed that it comes from a variety

        EXAMPLES:

        Hodge diamond of a K3 surface::

            sage: from diamond import *
            sage: x, y = (HodgeDiamond.x, HodgeDiamond.y)
            sage: S = HodgeDiamond.from_polynomial(1 + x**2 + 20*x*y + y**2 + x**2 * y**2)
            sage: S == K3()
            True

        The following fails as the lack of symmetry prevents a geometric origin::

            sage: HodgeDiamond.from_polynomial(1 + x, from_variety=True)
            Traceback (most recent call last):
            ...
            AssertionError: The matrix does not
                            satisfy the conditions satisfied by the Hodge diamond of a
                            smooth projective variety.
        """
        return HodgeDiamondRing()(f, from_variety=from_variety)

    @property
    def polynomial(self):
        r"""The Hodge--Poincaré polynomial describing the Hodge diamond

        :getter: returns the Hodge--Poincaré polynomial
        :setter: sets the Hodge diamond using a Hodge--Poincaré polynomial
        :type: element of :attr:`HodgeDiamond.R`

        EXAMPLES:

        The Hodge--Poincaré polynomial of a K3 surface::

            sage: from diamond import *
            sage: print(K3().polynomial)
            x^2*y^2 + x^2 + 20*x*y + y^2 + 1

        Modifying the Hodge diamond of the projective plane::

            sage: X = Pn(2)
            sage: X.polynomial = X.polynomial + X.x * X.y
            sage: print(X)
                      1
                  0       0
              0       2       0
                  0       0
                      1
        """
        M = self.matrix
        return self.R(
            {(i, j): M[i, j] for i in range(M.nrows()) for j in range(M.ncols())}
        )

    @polynomial.setter
    def polynomial(self, f):
        r"""Setter for the Hodge--Poincaré polynomial"""
        self.matrix = _to_matrix(f)

    @property
    def matrix(self):
        r"""The matrix describing the Hodge diamond

        :getter: returns the matrix
        :setter: sets the Hodge diamond using a matrix
        :type: square matrix of integers
        """
        return self._m

    @matrix.setter
    def matrix(self, m):
        r"""Setter for the Hodge diamond as a matrix"""
        m = matrix(m)

        assert m.base_ring() == ZZ, "Entries need to be integers"
        assert m.is_square()

        self._m = m
        self.__normalise()

    def size(self):
        r"""Internal method to determine the (relevant) size of the Hodge diamond"""
        return self._size

    def __normalise(self):
        r"""Internal method to get rid of trailing zeros"""
        self._m = _to_matrix(self.polynomial)
        self._size = self._m.ncols() - 1

    def __eq__(self, other):
        r"""Check whether two Hodge diamonds are equal

        This compares the Hodge polynomials, not the possibly oversized
        matrices describing the Hodge diamond.

        EXAMPLES:

        A quartic surface is a K3 surface::

            sage: from diamond import *
            sage: K3() == hypersurface(4, 2)
            True
        """
        if not isinstance(other, HodgeDiamond):
            return False
        return self.polynomial == other.polynomial

    def __ne__(self, other):
        r"""Check whether two Hodge diamonds are not equal

        EXAMPLES:

        The projective line is not a genus 2 curve::

            sage: from diamond import *
            sage: Pn(1) != curve(2)
            True

        The point is not the Lefschetz class::

            sage: point() != lefschetz()
            True
        """
        return not self == other

    def _add_(self, other):
        r"""Add two Hodge diamonds together

        This corresponds to taking the disjoint union of varieties, or the
        direct sum of the Hodge structure.

        EXAMPLES:

        Hodge diamond of the projective line is the sum of that of a point
        and the Lefschetz diamond::

            sage: from diamond import *
            sage: Pn(1) == 1 + lefschetz()
            True

        Adding zero doesn't do anything::

            sage: K3() + zero() == K3()
            True
        """
        return self.parent()(self.polynomial + other.polynomial)

    def _sub_(self, other):
        r"""Subtract two Hodge diamonds

        EXAMPLES:

        Hodge diamond of the projective line is the sum of that of a point
        and the Lefschetz diamond, but now we check it the other way around::

            sage: from diamond import *
            sage: Pn(1) - 1 == lefschetz()
            True
        """
        return self.parent()(self.polynomial - other.polynomial)

    def _mul_(self, other):
        r"""Multiply two Hodge diamonds

        This corresponds to taking the product of two varieties.

        EXAMPLES:

        The quadric surface is the product of two projective lines::

            sage: from diamond import *
            sage: Pn(1) * Pn(1) == hypersurface(2, 2)
            True

        The product is commutative::

            sage: K3() * curve(5) == curve(5) * K3()
            True

        The point is the unit::

            sage: K3() * point() == point() * K3() == K3()
            True

        TESTS::

           sage: 2 * K3() == K3() + K3()
           True
        """
        if not isinstance(other, HodgeDiamond):
            # in the rare case someone does X*3 instead of 3*X
            return other * self

        return self.parent()(self.polynomial * other.polynomial)

    def __pow__(self, power):
        r"""Raise a Hodge diamond to a power

        This corresponds to iterated multiplication.

        INPUT:

        - ``power`` -- exponent for the iterated multiplication

        EXAMPLES:

        The product of 2 K3 surfaces in two ways::

            sage: from diamond import *
            sage: K3()**2 == K3()*K3()
            True
        """
        return self.parent()(self.polynomial**power)

    def __call__(self, i, y=None):
        r"""
        The calling operator either does a Lefschetz twist, or an evaluation

        If one parameter is present, then twist by a power of the Lefschetz
        Hodge diamond. If two parameters are present, then evaluate the
        Hodge-Poincaré polynomial

        Negative values are allowed to untwist, up to the appropriate power.

        INPUT:

        - ``i`` -- integer denoting the power of the Lefschetz class, or
          value for the first variable

        - ``y`` -- value of the second variable (default: ``None``), if it is
          non-zero then ``i`` is reinterpreted as the value of the first variable

        EXAMPLES:

        The Lefschetz class is by definition the twist of the point::

            sage: from diamond import *
            sage: lefschetz() == point()(1)
            True

        We can reconstruct projective space as a sum of twists of the point::

            sage: Pn(10) == sum(point()(i) for i in range(11))
            True

        If we supply two parameters we are evaluation the Hodge-Poincaré
        polynomial, e.g. to find the Euler characteristic::

            sage: Pn(10)(1, 1) == 11
            True
        """
        if y is None:
            assert i >= -self.lefschetz_power()

            return self.parent()(self.R(self.polynomial * self.x**i * self.y**i))
        x = i
        return self.polynomial(x, y)

    def __getitem__(self, index):
        r"""Get (p, q)th entry of Hodge diamond or the ith row of the Hodge diamond"""
        # first try it as (p, q)
        try:
            p, q = index
            if p < 0 or q < 0:
                return 0
            else:
                try:
                    return self.matrix[p, q]
                except IndexError:
                    return 0
        # now we assume it's an integer: this is equivalent to HodgeDiamond.row(p)
        except TypeError:
            # we could do something smarter, but this is it for now
            return [self.matrix[p, index - p] for p in range(index + 1)]

    def _repr_(self):
        r"""Output diagnostic information

        This is a one-line string giving some basic information about the Hodge
        diamond. You'll see this when you just evaluate something which returns
        a Hodge diamond. To see something more useful, you'll likely want to use

        * :meth:`HodgeDiamond.__str__` via `print`
        * :meth:`HodgeDiamond.pprint`
        * :meth:`HodgeDiamond.polynomial`

        It is also possible to override this output by using the built-in
        functionality for parents and renaming.

        EXAMPLES:

        The projective line::

            sage: from diamond import *
            sage: Pn(1)
            Hodge diamond of size 2 and dimension 1

        We give it a more descriptive name::

            sage: P1 = Pn(1)
            sage: P1.rename("The projective line")
            sage: P1
            The projective line

        """
        return (
            f"Hodge diamond of size {self._size + 1} and dimension {self.dimension()}"
        )

    @property
    def name(self):
        return self.__repr__()

    def __str__(self):
        r"""Pretty print Hodge diamond

        This gets called when you specifically print the object.

        EXAMPLES:

        The projective line::

            sage: from diamond import *
            sage: print(Pn(1))
                  1
              0       0
                  1
        """
        return str(self.pprint())

    def __table(self, hide_zeroes=None, quarter=None):
        r"""Generate a table object for the Hodge diamond"""
        # take the default values from the static variables
        if hide_zeroes is None:
            if HodgeDiamond.hide_zeroes is None:
                hide_zeroes = False
            else:
                hide_zeroes = HodgeDiamond.hide_zeroes

        if quarter is None:
            if HodgeDiamond.quarter is None:
                quarter = False
            else:
                quarter = HodgeDiamond.quarter

        d = self._size
        T = []

        if self.is_zero():
            T = [[0]]
        else:
            for i in range(2 * d + 1):
                row = [""] * (abs(d - i))

                for j in range(max(0, i - d), min(i, d) + 1):
                    row.extend([self.matrix[j, i - j], ""])

                T.append(row)

        # making sure all rows have same length
        for i in range(len(T)):
            T[i].extend([""] * (2 * d - len(T[i]) + 1))
            T[i] = T[i][: 2 * d + 1]

        # replace zeroes by spaces, if requested
        if hide_zeroes:
            T = [[t if t != 0 else "" for t in row] for row in T]

        # only print the top-left quarter, if requested
        if quarter:
            T = [[T[i][j] for j in range(d + 1)] for i in range(d + 1)]

        # determine the minimum number of leading and trailing empty strings
        # and remove them to align better to the left
        if hide_zeroes:
            empty = []
            for row in T:
                (leading, trailing) = (0, 0)
                groups = [(k, len(list(g))) for k, g in groupby(row)]
                if groups[0][0] == "":
                    leading = groups[0][1]
                if groups[-1][0] == "":
                    trailing = groups[-1][1]
                empty.append(
                    (leading, trailing),
                )

            leading = min(a for (a, _) in empty)
            trailing = min(b for (_, b) in empty)

            for i in range(len(T)):
                T[i] = T[i][leading : len(T) - trailing]

        return table(T, align="center")

    def pprint(self, format="table", hide_zeroes=None, quarter=None):
        r"""Pretty print the Hodge diamond

        INPUT:

        - ``format`` -- output format (default: `"table"`), if table it pretty prints
          a Hodge diamond; all else defaults to the polynomial

        - ``hide_zeroes`` (default: False) -- whether to hide the zeroes if `"table"`
          is used as format

        - ``quarter`` (default: False) -- whether to only print the top-left quarter
          if `"table"` is used for the format

        The parameters ``hide_zeroes`` and ``quarter`` can be set using a static
        variable, in which case providing them will override this value (if it is set).

        EXAMPLES:

        The projective line::

            sage: from diamond import *
            sage: Pn(1).pprint()
                  1
              0       0
                  1
            sage: Pn(1).pprint(format="polynomial")
            x*y + 1

        Don't print the zeroes::

            sage: from diamond import *
            sage: (Pn(2) * curve(3)).pprint(hide_zeroes=True)
                  1
              3       3
                  2
              3       3
                  2
              3       3
                  1

        Only print the top-left quarter::

            sage: from diamond import *
            sage: (Pn(2) * curve(3)).pprint(quarter=True)
                          1
                      3
                  0       2
              0       3

        Only print the top-left quarter whilst hiding zeroes::

            sage: from diamond import *
            sage: (Pn(2) * curve(3)).pprint(hide_zeroes=True, quarter=True)
                  1
              3
                  2
              3

        """
        if format == "table":
            return self.__table(hide_zeroes=hide_zeroes, quarter=quarter)
        else:
            return self.polynomial

    def __is_positive(self):
        r"""Check whether all entries are positive integers"""
        return all(hpq >= 0 for hpq in self.matrix.coefficients())

    def is_hodge_symmetric(self):
        r"""Check whether the Hodge diamond satisfies Hodge symmetry

        This checks the equality

        .. MATH::

            \mathrm{h}^{p,q}(X)=\mathrm{h}^{q,p}(X)

        for $p,q=0,\\ldots,\\dim X$.

        Almost all of the constructions provided with the library satisfy Hodge
        symmetry, because we (somewhat implicitly) work with things which are
        (or behave like) smooth projective varieties over a field of
        characteristic zero.

        Over the complex numbers this can fail for non-Kähler manifolds, such
        as the Hopf surface.

        In positive characteristic this can fail too, with an example given by
        classical and singular Enriques surfaces in characteristic 2, see [MR0491720]
        and Proposition 1.4.2 in [MR0986969]

        * [MR0491720] Bombieri--Mumford, Enriques' classification of surfaces in char. p. III.
        * [MR0986969] Cossec--Dolgachev, Enriques surfaces I, Progress in Mathematics, 1989

        EXAMPLES:

        Constructions satisfy this property::

            sage: from diamond import *
            sage: Pn(5).is_hodge_symmetric()
            True

        The Hopf surface over the complex numbers::

            sage: S = HodgeDiamond.from_matrix([[1, 0, 0], [1, 0, 1], [0, 0, 1]])
            sage: print(S)
                      1
                  0       1
              0       0       0
                  1       0
                      1
            sage: S.is_hodge_symmetric()
            False

        Classical and singular Enriques surfaces in characteristic 2
        (which are smooth, despite their name) also have a Hodge diamond
        violating Hodge symmetry::

            sage: enriques(two="classical").is_hodge_symmetric()
            False
            sage: enriques(two="singular").is_hodge_symmetric()
            False
            sage: enriques(two="supersingular").is_hodge_symmetric()
            True
        """
        return self.matrix.is_symmetric()

    def is_serre_symmetric(self):
        r"""Check whether the Hodge diamond satisfies Serre symmetry

        This checks the equality

        .. MATH::

            \mathrm{h}^{p,q}(X)=\mathrm{h}^{\dim X-p,\dim X-q}(X)

        for $p,q=0,\\ldots,\\dim X$.

        Because Serre duality holds for all smooth projective varieties,
        independent of the characteristic, and also for non-Kähler varieties
        there are no examples where this condition fails. It can of course fail
        for motivic pieces, for silly reasons.

        EXAMPLES:

        The Hilbert scheme of 4 points on a K3 surface satisfies the symmetry::

            sage: from diamond import *
            sage: hilbn(K3(), 4).is_serre_symmetric()
            True

        The Lefschetz diamond fails it for silly reasons::

            sage: lefschetz().is_serre_symmetric()
            False

        """
        d = self._size
        return all(
            self.matrix[p, q] == self.matrix[d - p, d - q]
            for p in range(d + 1)
            for q in range(d + 1)
        )

    def betti(self):
        r"""Betti numbers of the Hodge diamond

        This gives an integer vector.

        EXAMPLES:

        Betti numbers of a K3 surface::

            sage: from diamond import *
            sage: K3().betti()
            [1, 0, 22, 0, 1]

        The second Betti number of the Hilbert scheme of points on a K3 surface
        is 23, not 22::

            sage: [hilbn(K3(), n).betti()[2] for n in range(2, 10)]
            [23, 23, 23, 23, 23, 23, 23, 23]

        """
        d = self._size
        return [
            ZZ.sum(self.matrix[j, i - j] for j in range(max(0, i - d), min(i, d) + 1))
            for i in range(2 * d + 1)
        ]

    def middle(self):
        r"""Middle cohomology of the Hodge diamond

        For smooth projective varieties the middle cohomology sits in degree
        equal to the dimension.

        EXAMPLES:

        There is an interesting link between K3 surfaces and cubic fourfolds
        which can be seen on the level of middle cohomology::

            sage: from diamond import *
            sage: (hypersurface(3, 4) - lefschetz()**2).middle()
            [0, 1, 20, 1, 0]
            sage: K3().middle()
            [1, 20, 1]

        """
        d = self._size
        return [self.matrix[i, d - i] for i in range(d + 1)]

    def row(self, i, truncate=False):
        r"""Get the ith row of the Hodge diamond

        For smooth projective varieties these are the Hodge numbers of the
        Hodge structure on the cohomology in degree `i`.

        Alternatively, you can use ``HodgeDiamond.__getitem__`` with a single
        index `i` (but then you have to truncate yourself).

        INPUT:

        - ``i`` -- the row of the Hodge diamond

        - ``truncate`` (default: False) -- whether you want to omit leading
          and trailing zeroes

        EXAMPLES:

        For a smooth projective variety the middle cohomology is the row
        sitting in the middle dimension::

            sage: from diamond import *
            sage: hypersurface(3, 4).middle() == hypersurface(3, 4).row(4)
            True

        If you don't want to truncate, ``HodgeDiamond.__getitem__`` gives the
        same functionality::

            sage: from diamond import *
            sage: hypersurface(3, 4).row(4) == hypersurface(3, 4)[4]
            True

        For the moduli space of vector bundles on a curve, the cohomology
        in degree 3 is the same as the cohomology of the curve in degree 1::

            sage: from diamond import *
            sage: moduli_vector_bundles(3, 1, 9).row(3, truncate=True)
            [9, 9]
        """
        row = [self.matrix[j, i - j] for j in range(i + 1)]

        if truncate:
            while row[0] == 0 and row[-1] == 0:
                row = row[1:-1]

        return row

    def signature(self):
        r"""The signature of the Hodge diamond

        This is the index of the intersection form on middle cohomology
        taken with real coefficients. By the Hodge index theorem it is given
        by the formula in Theorem 6.33 of Voisin's first book on Hodge theory.

        .. MATH::

            \sigma=\sum_{p,q=0}^{\dim X}(-1)^p\mathrm{h}^{p,q}

        This of course only makes sense if the diamond comes from a compact
        Kähler manifold.

        EXAMPLES:

            sage: from diamond import *
            sage: K3().signature()
            -16
        """
        assert self.arises_from_variety()

        d = self._size
        return sum((-1) ** p * self[p, q] for p in range(d + 1) for q in range(d + 1))

    def euler(self):
        r"""The topological Euler characteristic of the Hodge diamond

        This is the alternating sum of the Betti numbers, so that

        .. MATH::

            \chi_{\mathrm{top}}=\sum_{p,q=0}^{\dim X}(-1)^{p+q}\mathrm{h}^{p,q}

        EXAMPLES:

        The Euler characteristic of projective space grows linearly::

            sage: from diamond import *
            sage: [Pn(n).euler() for n in range(10)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        For Hilbert schemes of points of K3 surfaces these are the
        coefficients of the series expansion of the Dedekind eta-function,
        see A006922 in the OEIS::

            sage: [hilbn(K3(), n).euler() for n in range(10)]
            [1, 24, 324, 3200, 25650, 176256, 1073720, 5930496, 30178575, 143184000]
        """
        return ZZ.sum((-1) ** i * bi for i, bi in enumerate(self.betti()))

    def holomorphic_euler(self):
        r"""Holomorphic Euler characteristic

        This is the Euler characteristic of the structure sheaf, so

        .. MATH::

            \chi(X)=\sum_{i=0}^{\dim X}(-1)^i\mathrm{h}^{0,i}(X)

        EXAMPLES:

        For projective space it is 1::

            sage: from diamond import *
            sage: all(Pn(n).holomorphic_euler() == 1 for n in range(10))
            True

        For a hyperkähler variety of dimension $2n$ this number is $n+1$::

            sage: all(K3n(n).holomorphic_euler() == n+1 for n in range(5))
            True
        """
        return ZZ.sum((-1) ** i * self.matrix[i, 0] for i in range(self.matrix.nrows()))

    def hirzebruch(self):
        r"""Hirzebruch's \chi_y genus

        For a smooth projective variety $X$ Hirzebruch's $\\chi_y$-genus is
        defined as

        .. MATH::

            \chi_y(X)=\sum_{p,q=0}^{\dim X}(-1)^{p+q}\mathrm{h}^{p,q}(X)y^p

        which shows it is the specialisation of the Hodge-Poincaré polynomial
        for $x=-1$. A further specialisation to $y=-1$ gives the Euler characteristic.

        EXAMPLES:

        For a K3 surface we have::

            sage: from diamond import *
            sage: K3().hirzebruch()
            2*y^2 - 20*y + 2
            sage: K3().hirzebruch().subs(y=-1) == K3().euler()
            True

        For the Hilbert square of a K3 surface we get::

            sage: hilbn(K3(), 2).hirzebruch()
            3*y^4 - 42*y^3 + 234*y^2 - 42*y + 3
            sage: hilbn(K3(), 2).hirzebruch().subs(y=-1) == hilbn(K3(), 2).euler()
            True
        """
        return self.polynomial.subs(x=-1)

    def homological_unit(self):
        r"""Dimensions of $\\mathrm{H}^\\bullet(X,\\mathcal{O}_X)$

        A notion introduced by Abuaf.
        """
        return self.matrix.row(0)

    def hochschild(self):
        r"""Dimensions of the Hochschild homology

        Columns of the Hodge diamond are Hochschild homology, by the Hochschild-
        Kostant-Rosenberg theorem.
        """
        d = self._size
        return HochschildHomology.from_list(
            [
                ZZ.sum(
                    self.matrix[d - i + j, j]
                    for j in range(max(0, i - d), min(i, d) + 1)
                )
                for i in range(2 * d + 1)
            ]
        )

    def hh(self):
        r"""Shorthand for :meth:`HodgeDiamond.hochschild`"""
        return self.hochschild()

    def arises_from_variety(self):
        r"""Check whether the Hodge diamond can arise from a smooth projective variety

        The constraints are:

        - satisfy Hodge symmetry
        - satisfy Serre symmetry
        - there is no Lefschetz twist
        """
        return (
            self.is_hodge_symmetric()
            and self.is_serre_symmetric()
            and self.lefschetz_power() == 0
        )

    def is_zero(self):
        r"""Check whether the Hodge diamond is identically zero"""
        return self.matrix.is_zero()

    def lefschetz_power(self):
        r"""
        Return the twist by the Lefschetz motive that is present

        In other words, we see how divisible the Hodge--Poincaré polynomial is
        with respect to the monomial $x^iy^i$
        """
        if self.is_zero():
            return 0
        i = 0
        while (self.x**i * self.y**i).divides(self.polynomial):
            i += 1
        return i - 1

    def dimension(self):
        r"""Dimension of the Hodge diamond

        This takes twists by the Lefschetz class into account: we untwist by
        the maximal power and only then determine how big the diamond is.

        EXAMPLES:

        A point is 0-dimensional::

            sage: from diamond import *
            sage: point().dimension()
            0

        The Lefschetz diamond is also 0-dimensional::

            sage: print(lefschetz())
                  0
              0       0
                  1
            sage: lefschetz().dimension()
            0
        """
        assert self.is_hodge_symmetric()
        if self.is_zero():
            return -1
        return (
            max(
                [
                    i
                    for i in range(self.matrix.ncols())
                    if not self.matrix.column(i).is_zero()
                ]
            )
            - self.lefschetz_power()
        )

    def level(self):
        r"""Compute the level (or complexity) of the Hodge diamond

        This is a measure of the width of the non-zero part
        of the Hodge diamond.

        EXAMPLES:

        The simplest case is projective space, with level zero::

            sage: from diamond import *
            sage: all(Pn(n).level() == 0 for n in range(10))
            True

        For intersections of 2 quadrics it alternates between zero and one::

            sage: all(complete_intersection([2,2], 2*n).level() == 0 for n in range(5))
            True
            sage: all(complete_intersection([2,2], 2*n+1).level() == 1 for n in range(5))
            True

        A Calabi-Yau variety (e.g. a hypersurface of degree $n+1$ in $\\mathbb{P}^n$) has maximal level::

            sage: all(hypersurface(n+2, n).level() == n for n in range(10))
            True
        """
        return max(
            abs(m.degrees()[0] - m.degrees()[1]) for m in self.polynomial.monomials()
        )

    def blowup(self, other, codim=None):
        r"""Compute Hodge diamond of blowup

        No consistency checks are performed, this just naively applies the blowup
        formula from Hodge theory.

        INPUT:

        - ``other`` -- Hodge diamond of the center of the blowup

        - ``codim`` -- codimension of the center (optional), in case it is not
          the Hodge diamond of an honest variety

        EXAMPLES:

        A cubic surface is the blowup of $\\mathbb{P}^2$ in 6 points::

            sage: from diamond import *
            sage: Pn(2).blowup(6*point()) == hypersurface(3, 2)
            True
        """
        # let's guess the codimension
        if codim is None:
            codim = self.dimension() - other.dimension()

        return self + sum(other(i) for i in range(1, codim))

    def bundle(self, rank):
        r"""Compute the Hodge diamond of a projective bundle

        This applies the bundle formula from Hodge theory without any consistency checks.

        INPUT:

        - ``rank``: rank of the vector bundle on ``self``

        EXAMPLES:

        A projective bundle on a point is a projective space::

            sage: from diamond import *
            sage: point().bundle(3) == Pn(2)
            True

        A quadric surface is a $\\mathbb{P}^1$-bundle on $\\mathbb{P}^1$::

            sage: Pn(1).bundle(2) == hypersurface(2, 2)
            True
        """
        return sum(self(i) for i in range(rank))

    def mirror(self):
        r"""Compute the mirror Hodge diamond

        EXAMPLES:

        The mirror to a quintic 3-fold is the following::

            sage: from diamond import *
            sage: print(hypersurface(5, 3).mirror())
                             1
                        0         0
                    0       101       0
                1       1         1       1
                    0       101       0
                        0         0
                             1
        """
        assert self.arises_from_variety()
        n = self.dimension()
        x, y = self.x, self.y
        return self.parent()(
            sum(
                cf * x ** (n - exp[0]) * y ** exp[1]
                for exp, cf in self.polynomial.dict().items()
            )
        )


class HodgeDiamondRing(Singleton, Parent):
    def __init__(self):
        """
        TESTS::

            sage: from diamond import *
            sage: H = HodgeDiamondRing()
            sage: TestSuite(H).run()   # not tested (works only in the console)
        """
        Parent.__init__(self, category=Rings().Commutative())

    def _repr_(self) -> str:
        """ """
        return "Ring of Hodge diamonds"

    def _element_constructor_(self, *args, **keywords):
        m = args[0]
        if m in ZZ:
            m = matrix(ZZ, 1, 1, [m])
        elif isinstance(m, MPolynomial):
            m = _to_matrix(m)
        elif isinstance(m, (list, tuple)):
            m = matrix(m)
        elt = self.element_class(self, m)
        if keywords.get("from_variety", False):
            assert elt.arises_from_variety(), """The matrix does not satisfy the conditions satisfied by the
                Hodge diamond of a smooth projective variety."""
        return elt

    def from_matrix(self, m, from_variety=False):
        diamond = self.element_class(self, matrix(m))

        # get rid of trailing zeroes from the diamond
        diamond.matrix = _to_matrix(diamond.polynomial)

        return diamond

    def one(self):
        return point()

    def zero(self):
        return zero()

    def an_element(self):
        return K3()

    def _coerce_map_from_(self, R):
        return R is ZZ

    Element = HodgeDiamond


class HochschildHomology(Element):
    r"""
    This class implements some methods to work with (the dimensions of)
    Hochschild homology spaces, associated to the :class:`HodgeDiamond` class.

    The documentation is not intended to be complete, as this is mostly
    for my own sake.
    """

    # exposes R and t for external use
    R = LaurentPolynomialRing(ZZ, "t")
    t = R.gen()

    def __init__(self, parent, L):
        r"""
        Constructor for Hochschild homology dimensions of smooth and proper dg categories, so that Serre duality holds.

        INPUT:

        - ``L`` -- a list of integers of length 2n+1 representing $\\mathrm{HH}_{-n}$ to $\\mathrm{HH}_n$, such that ``L[i] == L[2n - i]``

        EXAMPLES::

            sage: from diamond import *
            sage: HochschildHomology.from_list([1,0,22,0,1])
            Hochschild homology vector of dimension 2

        TESTS::

            sage: from diamond import *
            sage: k = K3().hh()
            sage: k + k
            Hochschild homology vector of dimension 2
            sage: k*k
            Hochschild homology vector of dimension 4
            sage: k**2
            Hochschild homology vector of dimension 4
            sage: k.sym(2)
            Hochschild homology vector of dimension 4
            sage: 1 + k
            Hochschild homology vector of dimension 2
            sage: sum(k for i in range(2))
            Hochschild homology vector of dimension 2
            sage: 3*k
            Hochschild homology vector of dimension 2
        """
        assert len(L) % 2, "length needs to be odd, to reflect Serre duality"
        assert all(
            L[i] == L[len(L) - i - 1] for i in range(len(L))
        ), "Serre duality is not satisfied"
        self._L = L
        Element.__init__(self, parent)

    @classmethod
    def from_list(cls, L):
        r"""
        Constructor for Hochschild homology dimensions from a list.

        INPUT:

        - ``L`` -- a list of integers representing $\\mathrm{HH}_{-n}$ to $\\mathrm{HH}_n$

        EXAMPLES::

            sage: from diamond import *
            sage: HochschildHomology.from_list([1,0,22,0,1])
            Hochschild homology vector of dimension 2
        """
        return HochschildHomologies()(L)

    @classmethod
    def from_positive(cls, L):
        r"""
        Constructor for Hochschild homology dimensions from a list when only the positive part is given.

        INPUT:

        - ``L`` -- a list of integers representing $\\mathrm{HH}_0$ to $\\mathrm{HH}_n$

        EXAMPLES::

            sage: from diamond import *
            sage: HochschildHomology.from_positive([22,0,1])
            Hochschild homology vector of dimension 2
        """
        return HochschildHomologies()(L, positive=True)

    @classmethod
    def from_polynomial(cls, f):
        """
        Constructor for Hochschild homology dimensions from Hochschild--Poincaré Laurent polynomial

        INPUT

        - ``f`` -- the Hochschild--Poincaré Laurent polynomial

        EXAMPLES::

            sage: from diamond import *
            sage: x = LaurentPolynomialRing(ZZ, 'x').gen()
            sage: HochschildHomology.from_polynomial(x**-2+20+x**2)
            Hochschild homology vector of dimension 2
        """
        return HochschildHomologies()(f)

    @property
    def polynomial(self):
        """
        EXAMPLES::

            sage: from diamond import *
            sage: h = HochschildHomology.from_list([1,0,22,0,1])
            sage: h.polynomial
            t^-2 + 22 + t^2
        """
        return HochschildHomology.R(
            {i: self[i] for i in range(-self.dimension(), self.dimension() + 1)}
        )

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: from diamond import *
            sage: HochschildHomology.from_list([1,0,22,0,1])
            Hochschild homology vector of dimension 2
        """
        return f"Hochschild homology vector of dimension {self.dimension()}"

    def __str__(self) -> str:
        """
        EXAMPLES::

            sage: from diamond import *
            sage: print(HochschildHomology.from_list([1,0,22,0,1]))
              -2   -1   0    1   2
              1    0    22   0   1
        """
        return str(self.pprint())

    def __table(self):
        if self.is_zero():
            return table([[0], [0]])

        indices = list(range(-self.dimension(), self.dimension() + 1))

        return table([indices, [self[i] for i in indices]])

    def pprint(self, output="table"):
        if output == "table":
            return self.__table()
        return self._L

    def dimension(self):
        r"""Largest index ``i`` such that $\\mathrm{HH}_i\\neq 0$

        EXAMPLES::

            sage: from diamond import *
            sage: h = HochschildHomology.from_list([1,0,22,0,1])
            sage: h.dimension()
            2

        """
        if self.is_zero():
            return -1

        return (len(self._L) // 2) - min([i for i, d in enumerate(self._L) if d != 0])

    def is_zero(self):
        """
        EXAMPLES::

            sage: from diamond import *
            sage: h = HochschildHomology.from_list([1,0,22,0,1])
            sage: h.is_zero()
            False
        """
        return all(cf == 0 for cf in self._L)

    def euler(self):
        """
        Euler characteristic of Hochschild homology

        EXAMPLES::

            sage: from diamond import *
            sage: h = HochschildHomology.from_list([1,0,22,0,1])
            sage: h.euler()
            24
        """
        return self.polynomial(-ZZ.one())

    def _add_(self, other):
        return HochschildHomology.from_polynomial(self.polynomial + other.polynomial)

    def _sub_(self, other):
        return HochschildHomology.from_polynomial(self.polynomial - other.polynomial)

    def _mul_(self, other):
        return HochschildHomology.from_polynomial(self.polynomial * other.polynomial)

    def __pow__(self, i):
        return HochschildHomology.from_polynomial(self.polynomial**i)

    def __eq__(self, other):
        if not isinstance(other, HochschildHomology):
            return False
        return self.polynomial == other.polynomial

    def __ne__(self, other):
        return not self == other

    def __getitem__(self, i):
        if i > len(self._L) // 2:
            return 0
        return self._L[len(self._L) // 2 - i]

    def __iter__(self):
        return self._L.__iter__()

    def symmetric_power(self, k):
        r"""
        Hochschild homology of the Ganter--Kapranov symmetric power of a smooth and proper dg category

        This is possibly only a heuristic (I didn't check for proofs
        in the literature) based on the decomposition of Hochschild
        homology for a quotient stack, as discussed in the paper of
        Polishchuk--Van den Bergh.

        EXAMPLES::

            sage: from diamond import *
            sage: k = K3().hh()
            sage: k.symmetric_power(2)
            Hochschild homology vector of dimension 4
            sage: print(_)
              -4   -3   -2   -1   0     1   2    3   4
              1    0    23   0    276   0   23   0   1
        """

        def summand(f, k):
            assert all(c > 0 for c in f.coefficients())

            t = HochschildHomology.t

            # trivial case
            if f.number_of_terms() == 0:
                return f
            # base case
            elif f.number_of_terms() == 1:
                i = f.exponents()[0]
                a = f.coefficients()[0]

                if i % 2 == 0:
                    return binomial(a + k - 1, k) * t ** (k * i)
                else:
                    return binomial(a, k) * t ** (k * i)
            # general case
            else:
                # splitting f into monomial and the difference
                g = f.coefficients()[0] * t ** (f.exponents()[0])
                h = f - g

                return sum(summand(g, j) * summand(h, k - j) for j in range(k + 1))

        # see the object C^{(\lambda)} in the Polishchuk--Van den Bergh paper
        return HochschildHomology.from_polynomial(
            sum(
                prod(summand(self.polynomial, ri) for ri in P.to_exp())
                for P in Partitions(k)
            )
        )

    def sym(self, k):
        """Shorthand for ```HochschildHomology.symmetric_power```"""
        return self.symmetric_power(k)


class HochschildHomologies(Singleton, Parent):
    def __init__(self):
        """
        TESTS::

            sage: from diamond import *
            sage: H = HochschildHomologies()
            sage: TestSuite(H).run()   # not tested (works only in the console)
        """
        Parent.__init__(self, category=Rings().Commutative())

    def _repr_(self) -> str:
        """
        TESTS::

            sage: from diamond import *
            sage: HochschildHomologies()
            Ring of Hochschild homology vectors
        """
        return "Ring of Hochschild homology vectors"

    def _element_constructor_(self, *args, **keywords):
        if len(args) == 1:
            args = args[0]
        if args in ZZ:
            args = [args]
        if isinstance(args, (tuple, list)):
            if keywords.get("positive", False):
                # right half of the list
                L = list(reversed(args))[:-1] + args
            else:
                # full list
                L = args
        else:
            # a Laurent polynomial
            L = list(args)

        return self.element_class(self, L)

    def one(self):
        return self.element_class(self, [1])

    def zero(self):
        return self.element_class(self, [0])

    def _coerce_map_from_(self, R):
        return R is ZZ

    def an_element(self):
        return K3().hochschild()

    Element = HochschildHomology


# ==== Now a collection of diamonds ====


def zero():
    r"""Hodge diamond for the empty space

    EXAMPLES:

    Zero::

        sage: from diamond import *
        sage: print(zero())
          0
    """
    return HodgeDiamond.from_matrix(matrix([[0]]), from_variety=True)


def point():
    r"""Hodge diamond for the point

    EXAMPLES:

    The point::

        sage: from diamond import *
        sage: print(point())
          1
    """
    return HodgeDiamond.from_matrix(matrix([[1]]), from_variety=True)


def lefschetz():
    r"""Hodge diamond for the Lefschetz motive

    This is the Hodge-Poincaré polynomial of the affine line.

    EXAMPLES:

    The affine line::

        sage: from diamond import *
        sage: print(lefschetz())
              0
          0       0
              1

    We can take powers of it, to get the Hodge-Poincaré polynomial for higher-
    dimensional affine spaces::

        sage: print(lefschetz()**3)
                      0
                  0       0
              0       0       0
          0       0       0       0
              0       0       0
                  0       0
                      1
    """
    return point()(1)


def Pn(n):
    r"""
    Hodge diamond for projective space of dimension $n$

    INPUT:

    - ``n``: dimension, non-negative integer

    EXAMPLES:

    The zero-dimensional case is a point::

        sage: from diamond import *
        sage: Pn(0) == point()
        True

    In general projective space is the sum of powers of the Lefschetz class::

        sage: all(Pn(n) == sum([lefschetz()**i for i in range(n + 1)]) for n in range(1, 10))
        True
    """
    assert n >= 0
    return HodgeDiamond.from_matrix(matrix.identity(n + 1), from_variety=True)


def curve(genus):
    """
    Hodge diamond for a curve of a given genus.

    INPUT:

    - ``genus``: the genus of the curve, non-negative integer

    EXAMPLE:

    A curve of genus 0 is the 1-dimensional projective space::

        sage: from diamond import *
        sage: curve(0) == Pn(1)
        True

    A curve of genus 1 is an abelian variety of dimension 1::

        sage: curve(1) == abelian(1)
        True

    A curve of genus 2::

        sage: print(curve(2))
              1
          2       2
              1
    """
    assert genus >= 0
    return HodgeDiamond.from_matrix(matrix([[1, genus], [genus, 1]]), from_variety=True)


def surface(genus, irregularity, h11):
    r"""
    Hodge diamond for a surface $S$ with given invariants

    These invariants are the geometric genus, the irregularity and the middle
    Hodge numbers ``h11``.

    INPUT:

    - ``genus`` -- geometric genus of the surface, $\\dim\\mathrm{H}^2(S,\\mathcal{O}_S)$,
        a non-negative integer

    - ``irregularity`` -- irregularity of the surface, $\\dim\\mathrm{H}^1(S,\\mathcal{O}_S)$,
        a non-negative integer

    - ``h11`` -- middle Hodge number, $\\dim\\mathrm{H}^1(S,\\Omega_S^1)$,
        a non-negative integer

    EXAMPLES:

    The projective plane::

        sage: from diamond import *
        sage: Pn(2) == surface(0, 0, 1)
        True

    A K3 surface::

        sage: K3() == surface(1, 0, 20)
        True

    """
    assert genus >= 0
    assert irregularity >= 0
    assert h11 >= 0

    pg = genus
    q = irregularity

    return HodgeDiamond.from_matrix(
        matrix([[1, q, pg], [q, h11, q], [pg, q, 1]]), from_variety=True
    )


def symmetric_power(n, genus):
    r"""
    Hodge diamond for the nth symmetric power of a curve of given genus

    For the proof, see Example 1.1(1) of [MR2777820]. An earlier reference, probably in Macdonald, should exist.

    * [MR2777820] Laurentiu--Schuermann, Hirzebruch invariants of symmetric products. Topology of algebraic varieties and singularities, 163–177, Contemp. Math., 538, Amer. Math. Soc., 2011.

    INPUT:

    - ``n`` -- exponent of the symmetric power

    - ``genus`` -- genus of the curve, a non-negative integer

    EXAMPLES:

    The symmetric square of a genus 3 curve::

        sage: from diamond import *
        sage: print(symmetric_power(2, 3))
                  1
              3        3
          3       10       3
              3        3
                  1

    If $n=1$ we get the curve back::

        sage: all(symmetric_power(1, g) == curve(g) for g in range(10))
        True

    If $n=0$ we get the point::

        sage: symmetric_power(0, 4) == point()
        True

    If $n<0$ we have the empty space::

        sage: symmetric_power(-1, 4) == zero()
        True

    """

    def hpq(g, n, p, q):
        assert p <= n
        assert q <= n

        if p <= q and p + q <= n:
            return sum([binomial(g, p - k) * binomial(g, q - k) for k in range(p + 1)])
        if p > q:
            return hpq(g, n, q, p)
        if p + q > n:
            return hpq(g, n, n - p, n - q)

    if n < 0:
        return zero()

    M = matrix(n + 1)
    for i in range(n + 1):
        for j in range(n + 1):
            M[i, j] = hpq(genus, n, i, j)

    return HodgeDiamond.from_matrix(M, from_variety=True)


def jacobian(genus):
    """
    Hodge diamond for the Jacobian of a genus $g$ curve

    This is an abelian variety of dimension `genus`, so we call :func:`abelian`

    INPUT:

    - ``genus`` -- genus of the curve

    EXAMPLES:

    The Jacobian of a genus 3 curve::

        sage: from diamond import *
        sage: print(jacobian(3))
                      1
                  3       3
              3       9       3
          1       9       9       1
              3       9       3
                  3       3
                      1

    For the projective line we get a point::

        sage: jacobian(0) == point()
        True

    The Jacobian of an elliptic curve is isomorphic to it::

        sage: jacobian(1) == curve(1)
        True

    """
    return abelian(genus)


def abelian(dimension):
    r"""
    Hodge diamond for an abelian variety of a given dimension.

    The $g$th power of an elliptic curve is an abelian variety of the given
    dimension, so we just use this computation method.

    INPUT:

    - ``dimension`` -- dimension of the abelian variety

    EXAMPLES:

    A 1-dimensional abelian variety is an elliptic curve::

        sage: from diamond import *
        sage: abelian(1) == curve(1)
        True

    A 2-dimensional abelian variety is a surface with known Hodge numbers::

        sage: abelian(2) == surface(1, 2, 4)
        True

    """
    return curve(1) ** dimension


def kummer_resolution(dimension):
    """
    Hodge diamond for the standard resolution of the Kummer variety of an abelian variety of a given dimension.

    There's an invariant part (Hodge numbers of even degree) and the resolution of the $2^2g$ singularities is added.

    INPUT:

    - ``dimension`` -- dimension of the abelian variety taken as input

    EXAMPLES:

    The Kummer resolution of the involution on an abelian surface is a K3 surface::

        sage: from diamond import *
        sage: kummer_resolution(2) == K3()
        True

    """
    g = dimension

    invariant = sum(
        [
            jacobian(g).polynomial.monomial_coefficient(m) * m
            for m in jacobian(g).polynomial.monomials()
            if m.degree() % 2 == 0
        ]
    )
    return HodgeDiamond.from_polynomial(invariant) + sum(
        [2 ** (2 * g) * point()(i) for i in range(1, g)]
    )


def moduli_vector_bundles(rank, degree, genus):
    r"""
    Hodge diamond for the moduli space of vector bundles of given rank and
    fixed determinant of given degree on a curve of a given genus.

    For the proof, see Corollary 5.1 of [MR1817504].

    * [MR1817504] del Baño, On the Chow motive of some moduli spaces. J. Reine Angew. Math. 532 (2001), 105–132.

    If the Hodge diamond for the moduli space with non-fixed determinant of
    degree `d` is required, this can be obtained by::

        jacobian(g) * moduli_vector_bundles(r, d, g)

    INPUT:

    - `rank` -- rank of the bundles, at least 2

    - `degree` -- degree of the fixed determinant, coprime to rank

    - `genus` -- genus of the curve, at least 2

    EXAMPLES:

    The case of rank 2, degree 1 and genus 2 is famously the intersection of
    2 quadrics in $\\mathbb{P}^5$::

        sage: from diamond import *
        sage: moduli_vector_bundles(2, 1, 2) == complete_intersection([2, 2], 3)
        True

    """
    r = rank
    d = degree
    g = genus

    assert r >= 2, "rank needs to be at least 2"
    assert g >= 2, "genus needs to be at least 2"
    assert gcd(r, d) == 1, "rank and degree need to be coprime"

    R = HodgeDiamond.R
    x = HodgeDiamond.x
    y = HodgeDiamond.y

    def bracket(x):
        """Return the decimal part of a fraction."""
        return x - x.floor()

    def one(C, g):
        """Return the first factor in del Bano's formula."""
        return (
            (-1) ** (len(C) - 1)
            * ((1 + x) ** g * (1 + y) ** g) ** (len(C) - 1)
            / (1 - x * y) ** (len(C) - 1)
        )  # already corrected the factor from the Jacobian

    def two(C, g):
        """Return the second factor in del Bano's formula."""
        return prod(
            [
                prod(
                    [
                        (1 + x**i * y ** (i + 1)) ** g
                        * (1 + x ** (i + 1) * y**i) ** g
                        / ((1 - x**i * y**i) * (1 - x ** (i + 1) * y ** (i + 1)))
                        for i in range(1, C[j])
                    ]
                )
                for j in range(len(C))
            ]
        )

    def three(C, g):
        """Return the third factor in del Bano's formula."""
        return prod([1 / (1 - (x * y) ** (C[j] + C[j + 1])) for j in range(len(C) - 1)])

    def four(C, d, g):
        """Return the fourth factor in del Bano's formula."""
        exponent = sum(
            [sum([(g - 1) * C[i] * C[j] for i in range(j)]) for j in range(len(C))]
        ) + sum(
            [
                (C[i] + C[i + 1]) * bracket(-(sum(C[0 : i + 1]) * d / C.size()))
                for i in range(len(C) - 1)
            ]
        )
        return (x * y) ** exponent

    return HodgeDiamond.from_polynomial(
        R(
            sum(
                [
                    one(C, g) * two(C, g) * three(C, g) * four(C, d, g)
                    for C in Compositions(r)
                ]
            )
        ),
        from_variety=True,
    )


def seshadris_desingularisation(genus):
    r"""
    Hodge diamond for Seshadri's desingularisation of the moduli space of
    rank 2 bundles with trivial determinant on a curve of a given genus $g$.

    For the statement, see Corollary 3.18 of [MR1895918].

    * [MR1895918] del Baño, On the motive of moduli spaces of rank two vector bundles over a curve. Compositio Math. 131 (2002), 1-30.

    INPUT:

    - ``genus`` -- the genus $g$, at least 2

    EXAMPLES:

    For $g=2$ nothing needs to be desingularised, and the answer is $\\mathbb{P}^3$::

        sage: from diamond import *
        sage: seshadris_desingularisation(2) == Pn(3)
        True

    Already for $g=3$ the result is not a familiar variety, so we just check the
    Euler characteristic::

        sage: seshadris_desingularisation(3).euler() == 112
        True

    """
    g = genus

    assert g >= 2, "genus needs to be at least 2"

    R = HodgeDiamond.R
    x = HodgeDiamond.x
    y = HodgeDiamond.y

    L = x * y
    A = (1 + x) * (1 + y)
    B = (1 - x) * (1 - y)

    one = ((1 + x * L) ** g * (1 + y * L) ** g - L**g * A**g) / ((1 - L) * (1 - L**2))
    two = (A**g * (L - L**g) / (1 - L) + (A**g + B**g) / 2) / (
        1 + L
    )  # fixed typo in del Bano: compare to 3.12, motive of P^(g-2)
    three = (A**g - 2 ** (2 * g)) / 2 * ((1 - L ** (g - 1)) / (1 - L)) ** 2
    four = (
        (B**g / 2 - 2 ** (2 * g - 1)) * (1 - L ** (2 * g - 2)) / (1 - L**2)
    )  # fixed typo in del Bano: compare to 3.14, power of L at the end
    five = (
        (1 - L**g)
        * (1 - L ** (g - 1))
        * (1 - L ** (g - 2))
        / ((1 - L) * (1 - L**2) * (1 - L**3))
        + (1 - L**g) * (1 - L ** (g - 1)) / ((1 - L) * (1 - L**2)) * L ** (g - 2)
    ) * 2 ** (2 * g)

    return HodgeDiamond.from_polynomial(R(one - two + three + four + five))


def moduli_parabolic_vector_bundles_rank_two(genus, alpha):
    r"""
    Hodge diamond for the moduli space of parabolic rank 2 bundles with
    fixed determinant of odd degree on a curve of genus $g$.

    See Corollary 5.34 of [2011.14872].

    * [2011.14872] Fu--Hoskins--Pepin Lehalleur, Motives of moduli spaces of
      bundles on curves via variation of stability and flips

    This is not a proof of the formula we implemented per se, but it should be correct.
    Also, it could be that the choice of weights give something singular / stacky. Then it'll
    give bad output without warning. You have been warned.

    INPUT:

    - ``genus`` -- the genus of the curve

    - ``alpha`` -- the weights of the parabolic bundles

    """
    total = sum(alpha)
    N = len(alpha)
    rN = range(N)

    def d(j, alpha):
        return len(
            [
                1
                for I in Subsets(rN)
                if (len(I) - j) % 2 == 0
                and j - 1 < (len(I) + total - 2 * sum(alpha[i] for i in I)) < j + 1
            ]
        )

    def c(j, alpha):
        return binomial(N, j) - d(j, alpha)

    def b(j, alpha):
        return sum(((i + 2) // 2) * c(j - i, alpha) for i in range(j + 1))

    N = len(alpha)

    if genus == 0:
        M = zero()
    elif genus == 1:
        M = curve(1)
    elif genus >= 2:
        M = moduli_vector_bundles(2, 1, genus)

    result = M * (Pn(1) ** N) + sum(
        [b(j, alpha) * jacobian(genus)(genus + j) for j in range(N - 2)]
    )
    assert result.arises_from_variety()

    return result


def fano_variety_intersection_quadrics_odd(g, k):
    r"""
    Hodge diamond for the Fano variety of $k$-planes on the intersection of
    two quadrics in $\\mathbb{P}^{2g+1}$, using [MR3689749].

    We have that for $k=g-2$ we get M_C(2,L) as above, for deg L odd.

    * [MR3689749] Chen--Vilonen--Xue, On the cohomology of Fano varieties and the Springer correspondence, Adv. Math. 318 (2017), 515–533.

    EXAMPLES:

    For $k=0$ we have the intersection of two quadrics::

        sage: from diamond import *
        sage: fano_variety_intersection_quadrics_odd(2, 0) == complete_intersection([2, 2], 3)
        True
        sage: fano_variety_intersection_quadrics_odd(5, 0) == complete_intersection([2, 2], 9)
        True

    For $k=g-2$ we recover the moduli space of rank 2 bundles with odd determinant
    on a curve of genus $g$::

        sage: from diamond import *
        sage: fano_variety_intersection_quadrics_odd(11, 9) == moduli_vector_bundles(2, 1, 11)
        True

    For $k=g-1$ it is the Jacobian of $C$::

        sage: from diamond import *
        sage: fano_variety_intersection_quadrics_odd(12, 11) == jacobian(12)
        True

    For other $k$ it is an interesting variety::

        sage: from diamond import *
        sage: fano_variety_intersection_quadrics_odd(12, 9)
        Hodge diamond of size 51 and dimension 50

    """
    assert g >= 2, "genus needs to be at least 2"
    assert k in range(g), "non-empty only from 0 to g-1"

    if k == g - 1:
        return jacobian(g)

    # go back to the notation of Chen--Vilonen--Xue
    i = g - k
    # dimension of the Fano variety
    d = (g - i + 1) * (2 * i - 1)

    # multiplicity N_i(k, j) as in Theorem 1.1 of [MR3689749]
    def N(i, k, j):
        R = PowerSeriesRing(ZZ, default_prec=2 * d + 1)
        q = R.gen(0)
        return (
            q ** (-(j - i + 1) * (2 * i - 1))
            * (1 - q ** (4 * j))
            * prod([1 - q ** (2 * l) for l in range(j - i + 2, i + j - 1)])
            / prod([1 - q ** (2 * l) for l in range(1, 2 * i - 1)])
        )[k]

    x, y = (HodgeDiamond.x, HodgeDiamond.y)
    polynomial = 0

    for k in range(2 * d + 1):
        for j in range(i - 1, g + 1):
            if N(i, d - k, j) == 0:
                continue

            # the `g-j`th exterior power of the first cohomology of the curve
            # is the `g-j`th cohomology of the Jacobian
            dimensions = jacobian(g).row(g - j)
            # turn the dimensions into a polynomial
            piece = sum(
                [dimensions[m] * x**m * y ** ((g - j) - m) for m in range(g - j + 1)]
            )
            # the appropriate Lefschetz twist to put it in the right cohomological degree
            twist = x ** ((k - (g - j)) / 2) * y ** ((k - (g - j)) / 2)

            # they reindex using `d-k`
            polynomial = polynomial + N(i, d - k, j) * piece * twist

    return HodgeDiamond.from_polynomial(polynomial, from_variety=True)


def fano_variety_intersection_quadrics_even(g, k):
    r"""
    Hodge diamond for the Fano variety of $i-1$-planes on the intersection of
    two quadrics in $\\mathbb{P}^{2g+1}$, using [1510.05986v3].

    * [1510.05986v3] Chen--Vilonen--Xue, Springer correspondence, hyperelliptic curves, and cohomology of Fano varieties

    INPUT:

    - ``g`` -- half of the dimension of the quadrics

    - ``k`` -- dimension of the linear subspaces on the intersection of quadrics, at most $g-1$

    EXAMPLES:

    For $k=0$ we have the intersection of two quadrics::

        sage: from diamond import *
        sage: fano_variety_intersection_quadrics_even(2, 0) == complete_intersection([2, 2], 2)
        True
        sage: fano_variety_intersection_quadrics_even(5, 0) == complete_intersection([2, 2], 8)
        True

    We have that for $k = g-2$ we get the moduli space of parabolic bundles on
    $\\mathbb{P}^1$ with weight $1/2$ in $2g+3$ points::

        sage: from diamond import *
        sage: moduli_parabolic_vector_bundles_rank_two(0, [1/2]*5) == fano_variety_intersection_quadrics_even(2, 0)
        True
        sage: moduli_parabolic_vector_bundles_rank_two(0, [1/2]*9) == fano_variety_intersection_quadrics_even(4, 2)
        True

    For $k=g-1$ we get a finite reduced scheme of length $4^g$::

        sage: from diamond import *
        sage: print(fano_variety_intersection_quadrics_even(4, 3))
          256

    """
    assert g >= 2, "genus needs to be at least 2"
    assert k in range(g), "non-empty only from 0 to g-1"

    def M(k, j):
        index = k - j * (g - i)
        if index < 0:
            return 0
        else:
            return q_binomial(2 * g - i - j, i - j).padded_list(index + 1)[index]

    i = k + 1

    x, _ = (HodgeDiamond.x, HodgeDiamond.y)
    R = x.parent()
    polynomial = R(
        {
            (k, k): sum(M(k, j) * binomial(2 * g + 1, j) for j in range(i + 1))
            for k in range(i * (2 * g - 2 * i) + 1)
        }
    )

    return HodgeDiamond.from_polynomial(polynomial, from_variety=True)


def quot_scheme_curve(genus, length, rank):
    """
    Hodge diamond for the Quot scheme of zero-dimensional quotients of given
    length of a vector bundle of given rank on a curve of given genus.

    For the proof, see Proposition 4.5 of [1907.00826] (or rather, the
    reference [Bif89] in there)

    * [1907.00826] Bagnarol--Fantechi--Perroni, On the motive of zero-dimensional Quot schemes on a curve
    """

    def dn(P):
        # shift in indexing because we start at 0
        return sum(i * ni for i, ni in enumerate(P))

    return sum(
        [
            prod([symmetric_power(ni, genus) for ni in P])(dn(P))
            for P in IntegerVectors(length, rank)
        ]
    )


def hilbtwo(X):
    """
    Hodge diamond for the Hilbert square of any smooth projective variety

    For the proof, see e.g. lemma 2.6 of [MR2506383], or the corollary on page 507 of [MR1382733] (but it can be said to be classical).

    * [MR2506383] Muñoz--Ortega--Vázquez-Gallo, Hodge polynomials of the moduli spaces of triples of rank (2,2). Q. J. Math. 60 (2009), no. 2, 235–272.
    * [MR1382733] Cheah, On the cohomology of Hilbert schemes of points, J. Algebraic Geom. 5 (1996), no. 3, 479-511

    INPUT:

    -  ``X`` - Hodge diamond of the smooth projective variety

    EXAMPLES:

    We recover the Hilbert square of a surface::

        sage: from diamond import *
        sage: hilbtwo(K3()) == hilbn(K3(), 2)
        True

    """
    assert X.arises_from_variety()

    d = X.dimension()

    return HodgeDiamond.from_polynomial(
        X.R(
            ((X.polynomial) ** 2 + X.polynomial(-X.x**2, -X.y**2)) / 2
            + sum([X.x**i * X.y**i * X.polynomial for i in range(1, d)])
        ),
        from_variety=True,
    )


def hilbthree(X):
    """
    Hodge diamond of the Hilbert cube of any smooth projective variety

    The corollary on page 507 of [MR1382733].

    * [MR1382733] Cheah, On the cohomology of Hilbert schemes of points, J. Algebraic Geom. 5 (1996), no. 3, 479-511

    INPUT:

    -  ``X`` - Hodge diamond of the smooth projective variety

    EXAMPLES:

    We recover the Hilbert cube of a surface::

        sage: from diamond import *
        sage: hilbthree(K3()) == hilbn(K3(), 3)
        True

    """
    assert X.arises_from_variety()

    d = X.dimension()
    X2 = X**2
    R = HodgeDiamond.R

    return HodgeDiamond.from_polynomial(
        (X**3).polynomial / 6
        + X.polynomial * X.polynomial(-X.x**2, -X.y**2) / 2
        + X.polynomial(X.x**3, X.y**3) / 3
        + R.sum(X2(i).polynomial for i in range(1, d))
        + R.sum(X(i + j).polynomial for i in range(1, d) for j in range(i, d)),
        from_variety=True,
    )


def K3n(n):
    r"""
    Hodge diamond of the Hilbert scheme of $n$ points on a K3 surface

    This is the first family of hyperkähler varieties, constructed by Beauville.

    INPUT:

    - ``n`` -- number of points

    EXAMPLES:

    For $n=1$ we have a K3 surface::

        sage: from diamond import *
        sage: K3n(1) == K3()
        True

    For $n\\geq 2$ we have second Betti number 23::

        sage: all(K3n(n).betti()[2] == 23 for n in range(2, 5))
        True
    """
    return hilbn(K3(), n)


def generalised_kummer(n):
    r"""
    Hodge diamond of the $n$th generalised Kummer variety

    For the proof, see Corollary 1 of [MR1219901].

    * [MR1219901] Göttsche--Soergel, Perverse sheaves and the cohomology of Hilbert schemes of smooth algebraic surfaces. Math. Ann. 296 (1993), no. 2, 235–245.

    EXAMPLES:

    The first generalised Kummer is just a point::

        sage: from diamond import *
        sage: generalised_kummer(1) == point()
        True

    The second generalised Kummer is the Kummer K3 surface::

        sage: generalised_kummer(2) == K3()
        True

    The higher generalised Kummers are hyperkähler varieties with second Betti number 7::

        sage: all(generalised_kummer(n).betti()[2] == 7 for n in range(3, 10))
        True
    """
    x = HodgeDiamond.x
    y = HodgeDiamond.y

    def product(n):
        hd = sum(
            gcd(a := ap.to_exp_dict()) ** 4
            * (x * y) ** (n - sum(a.values()))
            * prod(
                [
                    sum(
                        [
                            prod(
                                [
                                    ~((j**bj) * factorial(bj))
                                    * ((1 - x**j) * (1 - y**j)) ** (2 * bj)
                                    for j, bj in b.to_exp_dict().items()
                                ]
                            )
                            for b in Partitions(ai)
                        ]
                    )
                    for ai in a.values()
                ]
            )
            for ap in Partitions(n)
        )
        return HodgeDiamond.R(hd(-x, -y))

    # Göttsche--Soergel gives the polynomial for A\times Kum^n A, so we quotient out A
    return HodgeDiamond.from_polynomial(product(n) // product(1), from_variety=True)


def ogrady6():
    r"""
    Hodge diamond for O'Grady's exceptional 6-dimensional hyperkähler variety

    For the proof, see Theorem 1.1 of [MR3798592].

    * [MR3798592] Mongardi--Rapagnetta--Saccà, The Hodge diamond of O'Grady's six-dimensional example. Compos. Math. 154 (2018), no. 5, 984–1013.

    EXAMPLES:

    The second Betti number is 8::

        sage: from diamond import *
        sage: ogrady6().betti()[2] == 8
        True
    """
    H = HodgeDiamond.from_matrix
    return H(
        [
            [1, 0, 1, 0, 1, 0, 1],
            [0, 6, 0, 12, 0, 6, 0],
            [1, 0, 173, 0, 173, 0, 1],
            [0, 12, 0, 1144, 0, 12, 0],
            [1, 0, 173, 0, 173, 0, 1],
            [0, 6, 0, 12, 0, 6, 0],
            [1, 0, 1, 0, 1, 0, 1],
        ],
        from_variety=True,
    )


def ogrady10():
    """
    Hodge diamond for O'Grady's exceptional 10-dimensional hyperkähler variety

    For the proof, see theorem A of [1905.03217]

    * [1905.03217] de Cataldo--Rapagnetta--Saccà, The Hodge numbers of O'Grady 10 via Ngô strings

    EXAMPLES:

    The second Betti number is 24::

        sage: from diamond import *
        sage: ogrady10().betti()[2] == 24
        True
    """
    H = HodgeDiamond.from_matrix
    return H(
        [
            [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
            [0, 22, 0, 22, 0, 23, 0, 22, 0, 22, 0],
            [1, 0, 254, 0, 276, 0, 276, 0, 254, 0, 1],
            [0, 22, 0, 2299, 0, 2531, 0, 2299, 0, 22, 0],
            [1, 0, 276, 0, 16490, 0, 16490, 0, 276, 0, 1],
            [0, 23, 0, 2531, 0, 88024, 0, 2531, 0, 23, 0],
            [1, 0, 276, 0, 16490, 0, 16490, 0, 276, 0, 1],
            [0, 22, 0, 2299, 0, 2531, 0, 2299, 0, 22, 0],
            [1, 0, 254, 0, 276, 0, 276, 0, 254, 0, 1],
            [0, 22, 0, 22, 0, 23, 0, 22, 0, 22, 0],
            [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
        ],
        from_variety=True,
    )


def hilbn(surface, n):
    """
    Hodge diamond for Hilbert scheme of ``n`` points on a smooth projective surface ``S``

    For the proof, see Theorem 2.3.14 of [MR1312161].

    * [MR1312161] Göttsche, Hilbert schemes of zero-dimensional subschemes of smooth varieties. Lecture Notes in Mathematics, 1572. Springer-Verlag, Berlin, 1994. x+196 pp.

    INPUT:

    - ``S`` -- Hodge diamond for smooth projective surface

    - ``n`` -- number of points
    """
    assert surface.arises_from_variety()
    assert surface.dimension() == 2

    ring_ab = PolynomialRing(ZZ, "a,b")
    a, b = ring_ab.gens()
    R = PowerSeriesRing(ring_ab, "t", default_prec=n + 1)
    t = R.gen().O(n + 1)

    series = R.one().O(n + 1)

    # theorem 2.3.14 of Göttsche's book
    for k in range(1, n + 1):
        for p in range(3):
            for q in range(3):
                eps_pq = (-1) ** (p + q + 1)
                term = (1 + eps_pq * a ** (p + k - 1) * b ** (q + k - 1) * t**k).O(
                    n + 1
                )
                series *= term ** (eps_pq * surface[p, q])

    coeff_n = series[n]

    # read off Hodge diamond from the (truncated) series
    M = matrix(2 * n + 1)
    for p in range(2 * n + 1):
        for q in range(2 * n + 1):
            M[p, q] = coeff_n.coefficient(
                [min(p, 2 * n - q), min(q, 2 * n - p)]
            )  # use Serre duality

    return HodgeDiamond.from_matrix(M, from_variety=True)


def nestedhilbn(surface, n):
    """
    Hodge diamond for the nested Hilbert scheme ``S^{[n-1,n]}``

    This is the unique nested Hilbert scheme of a smooth projective
    surface ``S`` which is itself smooth (of dimension $2n$)
    """
    assert surface.arises_from_variety()
    assert surface.dimension() == 2

    ring_xy = PolynomialRing(ZZ, "x,y")
    x, y = ring_xy.gens()
    R = PowerSeriesRing(ring_xy, "t", default_prec=n + 1)
    t = R.gen().O(n + 1)

    series = R.one().O(n + 1)
    for k in range(1, n + 1):
        for p in range(3):
            for q in range(3):
                s_pq = surface[p, q]
                if p + q % 2:
                    term = (1 + x ** (p + k - 1) * y ** (q + k - 1) * t**k).O(n + 1)
                    series *= term**s_pq
                else:
                    term = (1 - x ** (p + k - 1) * y ** (q + k - 1) * t**k).O(n + 1)
                    series *= term ** (-s_pq)

    series = series * R(surface.polynomial) * t / (1 - x * y * t)
    top_poly = series[n]

    # read off Hodge diamond from the (truncated) series
    M = matrix(2 * n + 1)
    for p in range(2 * n + 1):
        for q in range(2 * n + 1):
            M[p, q] = top_poly.coefficient([min(p, 2 * n - q), min(q, 2 * n - p)])
            # use Serre duality

    return HodgeDiamond.from_matrix(M, from_variety=True)


def complete_intersection(degrees, dimension):
    r"""
    Hodge diamond for a complete intersection of multidegree $(d_1,\\ldots,d_k)$ in $\\mathbb{P}^{n+k}$

    For a proof, see théorème 2.3 of exposé XI in SGA7.

    INPUT:

    - ``degrees`` -- the multidegree, if it is an integer we interpret it as a hypersurface

    - ``dimension`` -- the dimension of the complete intersection (not of the ambient space)

    EXAMPLES:

    For multidegrees $(1,\\ldots,1) we get a lower-dimension projective space::

        sage: from diamond import *
        sage: complete_intersection(1, 2) == Pn(2)
        True
        sage: complete_intersection([1, 1], 2) == Pn(2)
        True
        sage: complete_intersection([1, 1], 5) == Pn(5)
        True

    The Euler characteristics of cubic hypersurfaces::

        sage: [complete_intersection(3, n).euler() for n in range(10)]
        [3, 0, 9, -6, 27, -36, 93, -162, 351, -672]

    The Euler characteristics of intersections of 2 quadrics::

        sage: [complete_intersection([2, 2], n).euler() for n in range(10)]
        [4, 0, 8, 0, 12, 0, 16, 0, 20, 0]

    """
    # hypersurface as complete intersection
    try:
        degrees = list(degrees)
    except TypeError:
        degrees = [degrees]

    R = PowerSeriesRing(ZZ, ("a", "b"), default_prec=dimension + 2)
    a, b = R.gens()
    H = ~((1 + a) * (1 + b)) * (
        prod(
            [
                ((1 + a) ** di - (1 + b) ** di)
                / (a * (1 + b) ** di - b * (1 + a) ** di)
                for di in degrees
            ]
        )
        - 1
    ) + ~(1 - a * b)

    H = H.polynomial()

    M = matrix.identity(dimension + 1)
    for i in range(dimension + 1):
        M[i, dimension - i] = H.coefficient([i, dimension - i])

    return HodgeDiamond.from_matrix(M, from_variety=True)


def hypersurface(degree, dimension):
    r"""Shorthand for a complete intersection of the given dimension where $k=1$"""
    return complete_intersection(degree, dimension)


def K3():
    """
    Hodge diamond for a K3 surface

    EXAMPLES:

    The K3 surface::

        sage: from diamond import *
        sage: print(K3())
                  1
              0        0
          1       20       1
              0        0
                  1
    """
    dia = complete_intersection(4, 2)
    dia.rename("Hodge diamond of K3 surface")
    return dia


def enriques(two=None):
    r"""Hodge diamond for an Enriques surface

    It is possible to ask for the Hodge diamond of a classical or (super)singular
    Enriques surface in characteristic 2.

    In characteristic 2 the invariants are given in Proposition 1.4.2 of [MR0986969].

    * [MR0986969] Cossec--Dolgachev, Enriques surfaces I, Progress in Mathematics, 1989

    INPUT:

    - ``two`` -- optional parameter to indicate the type of surface in characteristic 2
        possible values are `"classical"`, `"singular"`, `"supersingular"`

    EXAMPLES:

    An ordinary Enriques surface::

        sage: from diamond import *
        sage: print(enriques())
                  1
              0        0
          0       10       0
              0        0
                  1

    Enriques surfaces in characteristic 2::

        sage: print(enriques(two="classical"))
                  1
              0        1
          0       12       0
              1        0
                  1
        sage: print(enriques(two="singular"))
                  1
              1        0
          1       10       1
              0        1
                  1
        sage: print(enriques(two="supersingular"))
                  1
              1        1
          1       12       1
              1        1
                  1
    """
    if two:
        if two == "classical":
            return HodgeDiamond.from_matrix([[1, 0, 0], [1, 12, 1], [0, 0, 1]])
        if two == "singular":
            return HodgeDiamond.from_matrix([[1, 1, 1], [0, 10, 0], [1, 1, 1]])
        if two == "supersingular":
            return HodgeDiamond.from_matrix([[1, 1, 1], [1, 12, 1], [1, 1, 1]])
        raise ValueError("invalid choice for characteristic 2")
    return surface(0, 0, 10)


def ruled(genus):
    r"""Hodge diamond for a ruled surface

    These are $\\mathbb{P}^1$-bundles over a curve of given genus.

    INPUT:

    - ``genus`` -- genus of the base curve

    EXAMPLES:

    For genus 0 we get Hirzebruch surfaces, whose Hodge diamond is that of
    the quadric surface::

        sage: from diamond import *
        sage: ruled(0) == hypersurface(2, 2)
        True

    For higher genus the Hodge diamond looks as follows::

        sage: print(ruled(5))
                  1
              5       5
          0       2       0
              5       5
                  1
    """
    return surface(0, genus, 2)


def weighted_hypersurface(degree, weights):
    """
    Hodge diamond for a weighted hypersurface of degree ``d`` in ``P(w_0,...,w_n)``

    This implements Theorem 7.2 of Fletcher's notes, Working with weighted complete
    intersections.

    INPUT:

    - ``degree`` -- degree of the hypersurface

    - ``weights`` -- the weights of the weighted projective space, if it is an
      integer we interpret it as the dimension of P^n

    EXAMPLES:

    Elliptic curves can be realised as hypersurfaces in 3 ways::

        sage: from diamond import *
        sage: weighted_hypersurface(3, 2) == weighted_hypersurface(3, [1, 1, 1])
        True
        sage: weighted_hypersurface(3, 2) == curve(1)
        True
        sage: weighted_hypersurface(4, [1, 1, 2]) == curve(1)
        True
        sage: weighted_hypersurface(6, [1, 2, 3]) == curve(1)
        True

    The Fano 3-fold 1.1 is a weighted hypersurface::

        sage: fano_threefold(1, 1) == weighted_hypersurface(6, [1,1,1,1,3])
        True

    If the variety is only quasismooth, not smooth, then we have to interpret
    the Hodge numbers accordingly. For instance, number 2 on Reid's list of 95 K3s
    has middle Hodge number 19, because the surface is has one node.

        sage: weighted_hypersurface(5, [1, 1, 1, 2]).middle()
        [1, 19, 1]

    """
    # weights should be interpreted as dimension of unweighted P^n
    if isinstance(weights, Integer):
        weights = [1] * (weights + 1)

    n = len(weights) - 1

    def hij(d, W, i, j):
        if i + j != n - 1 and i != j:
            return 0
        if i + j != n - 1 and i == j:
            return 1

        w = sum(W)

        R = PowerSeriesRing(ZZ, default_prec=j * d + d + 1)
        t = R.gen(0)

        H = prod((1 - t ** (d - wi)) / (1 - t**wi) for wi in W)

        if i + j == n - 1 and i != j:
            return H[j * d + d - w]
        if i + j == n - 1 and i == j:
            return H[j * d + d - w] + 1

    M = matrix.identity(n)
    for i in range(n):
        M[i, n - i - 1] = hij(degree, weights, i, n - i - 1)

    return HodgeDiamond.from_matrix(M, from_variety=True)


def cyclic_cover(ramification_degree, cover_degree, weights):
    r"""
    Hodge diamond of a cyclic cover of weighted projective space

    Implementation taken from
    https://github.com/jxxcarlson/math_research/blob/master/hodge.sage.

    INPUT:

    - ``ramification_degree`` -- degree of the ramification divisor

    - ``cover_degree`` -- size of the cover

    - ``weights`` -- the weights of the weighted projective space, if it is an
      integer we interpret it as the dimension of ``P^n``

    EXAMPLES:

    Some K3 surfaces are double covers of $\\mathbb{P}^2$ in a sextic curve::

        sage: from diamond import *
        sage: cyclic_cover(6, 2, 2) == K3()
        True

    The Fano 3-fold 1.1 is a double cover of $\\mathbb{P}^3$ in a sextic::

        sage: cyclic_cover(6, 2, 3) == fano_threefold(1, 1)
        True
    """
    # weights should be interpreted as dimension of unweighted P^n
    if isinstance(weights, Integer):
        weights = [1] * (weights + 1)

    weights.append(ramification_degree // cover_degree)

    return weighted_hypersurface(ramification_degree, weights)


def partial_flag_variety(D, I):
    r"""
    Hodge diamond of a partial flag variety G/P

    This is computed by counting the number of Schubert cells in the
    appropriate dimension.

    INPUT:

    - ``D`` -- Dynkin type

    - ``I`` -- indices of vertices to be omitted in defining the parabolic
      subgroup

    EXAMPLES:

    An absolute baby case is projective space::

        sage: from diamond import *
        sage: partial_flag_variety("A5", [2,3,4,5]) == Pn(5)
        True

    The next easiest case are quadrics::

        sage: partial_flag_variety("B5", [2,3,4,5]) == hypersurface(2, 9)
        True
        sage: partial_flag_variety("D5", [2,3,4,5]) == hypersurface(2, 8)
        True

    """
    R = PolynomialRing(ZZ, "x")
    x = R.gen()

    D = DynkinDiagram(D)

    WG = D.root_system().root_lattice().weyl_group()
    P = prod(sum(x**i for i in range(n)) for n in WG.degrees())

    if not I:
        Q = 1
    else:
        WL = D.subtype(I).root_system().root_lattice().weyl_group()
        Q = prod(sum(x**i for i in range(n)) for n in WL.degrees())

    return HodgeDiamond.from_matrix(
        diagonal_matrix((P // Q).coefficients()), from_variety=True
    )


def generalised_grassmannian(D, k):
    r"""
    Hodge diamond of the generalised Grassmannian of type D modulo the maximal parabolic subgroup $P_k$.

    This is just shorthand for :func:`partial_flag_variety(D, I)` where ``I`` is the complement of a singleton.

    INPUT:

    - ``D`` -- Dynkin type

    - ``k`` -- the vertex in the Dynkin diagram defining the maximal parabolic

    """
    return partial_flag_variety(D, [i for i in range(1, int(D[1:]) + 1) if i != k])


def grassmannian(k, n):
    r"""
    Hodge diamond of the Grassmannian $\\operatorname{Gr}(k,n)$ of $k$-dimensional subspaces in
    an $n$-dimensional vector space

    INPUT:

    - ``k`` -- dimension of the subspaces

    - ``n`` -- dimension of the ambient vector space

    EXAMPLES:

    Grassmannians are projective spaces if `k` is one or `n-1`::

        sage: from diamond import *
        sage: grassmannian(1, 5) == Pn(4)
        True
        sage: grassmannian(7, 8) == Pn(7)
        True

    The Grassmannian of 2-planes in a 4-dimensional vector space is the Kleiin quadric::

        sage: grassmannian(2, 4) == hypersurface(2, 4)
        True

    """
    assert 0 <= k <= n
    if n in [0, k]:
        return point()

    x, y = HodgeDiamond.x, HodgeDiamond.y
    return HodgeDiamond.from_polynomial(q_binomial(n, k)(q=x * y))


def orthogonal_grassmannian(k, n):
    r"""
    Hodge diamond of the orthogonal Grassmannian $\\operatorname{OGr}(k, n)$ of $k$-dimensional
    subspaces in an $n$-dimensional vector space isotropic with respect to a
    non-degenerate symmetric bilinear form

    INPUT:

    - ``k`` -- dimension of the subspaces

    - ``n`` -- dimension of the ambient vector space

    """
    if n % 2 == 0:
        assert k < n // 2
    else:
        assert k <= n // 2

    if n % 2 == 0:
        D = "D" + str(n // 2)
        if k - 1 == n // 2:
            # exceptional case: need submaximal parabolic
            I = list(range(1, n // 2 - 1))
        else:
            I = [i for i in range(1, n // 2 + 1) if i != k]
        return partial_flag_variety(D, I)
    else:
        D = "B" + str(n // 2)
        I = [i for i in range(1, n // 2 + 1) if i != k]
        return partial_flag_variety(D, I)


def symplectic_grassmannian(k, n):
    r"""
    Hodge diamond of the symplectic Grassmannian $\\operatorname{SGr}(k, n)$ of $k$-dimensional
    subspaces in an $n$-dimensional vector space isotropic with respect to a
    non-degenerate skew-symmetric bilinear form

    INPUT:

    - ``k`` -- dimension of the subspaces

    - ``n`` -- dimension of the ambient vector space

    """
    assert n % 2 == 0

    D = "C" + str(n // 2)
    I = [i for i in range(1, n // 2 + 1) if i != k]
    return partial_flag_variety(D, I)


def lagrangian_grassmannian(n):
    """Shorthand for the symplectic Grassmannian of Lagrangian subspaces"""
    return symplectic_grassmannian(n, 2 * n)


def horospherical(D, y=0, z=0):
    r"""
    Horospherical varieties as discussed in [1803.05063], with labelling
    and notation as in op. cit.

    INPUT:

    - ``D``: either a Dynkin type from the (small) list of allowed types
             in the classification
             _or_ a plaintext label from X1(n), X2, X3(n,m), X4, X5

    - ``y``: index for the parabolic subgroup for Y, see classification

    - ``z``: index for the parabolic subgroup for the closed orbit Z

    ``y`` and ``z`` must be omitted if a plaintext description is given.

    * [1803.05063] Gonzales--Pech--Perrin--Samokhin, Geometry of horospherical
      varieties of Picard rank one
    """
    # treat D as the plaintext description of a horospherical variety
    # not supposed to be 100% robust
    if y == 0 and z == 0:
        X = D  # rename for less confusion

        i = int(X[1])

        if i == 1:
            n = X[3:-1]
            return horospherical("B" + n, int(n) - 1, int(n))
        if i == 2:
            return horospherical("B3", 1, 3)
        if i == 3:
            n = X[3:].split(",")[0]
            m = int(X[:-1].split(",")[1])
            return horospherical("C" + n, m, m - 1)
        if i == 4:
            return horospherical("F4", 2, 3)
        if i == 5:
            return horospherical("G2", 1, 2)

        # didn't recognise it so far, so must be wrong
        raise Exception

    # determine the Hodge diamond from the blowup description
    Y = generalised_grassmannian(D, y)
    Z = generalised_grassmannian(D, z)

    n = int(D[1:])

    if D[0] == "B":
        assert (n == 3 and y == 1 and z == 3) or (n >= 3 and y == n - 1 and z == n)
        if n == 3 and y == 1:
            dimension = 9
        else:
            dimension = n * (n + 3) / 2
    elif D[0] == "C":
        assert n >= 2 and y in range(2, n + 1) and z == y - 1
        dimension = y * (2 * n + 1 - y) - y * (y - 1) / 2
    elif D == "F4":
        assert y == 2 and z == 3
        dimension = 23
    elif D == "G2":
        assert y == 1 and z == 2
        dimension = 7

    codimXY = dimension - Y.dimension()
    codimXZ = dimension - Z.dimension()

    return Y.bundle(codimXY + 1) + Z - Z.bundle(codimXZ)


def odd_symplectic_grassmannian(k, n):
    r"""
    Hodge diamond of the odd symplectic Grassmannian $\\operatorname{SGr}(k,n)$

    Here $n$ is odd. This is just shorthand for a call to :func:`horospherical_variety`
    for type C, with parameters $\lfloor n/2\rlfloor$, and $Y$ and $Z$ determined
    by $k$ and $k - 1$.
    """
    assert n % 2 == 1

    return horospherical("C" + str(n // 2), k, k - 1)


def gushel_mukai(n):
    r"""
    Hodge diamond for a smooth $n$-dimensional Gushel--Mukai variety

    See Proposition 3.1 of [1605.05648v3].

    * [1605.05648v3] Debarre--Kuznetsov, Gushel-Mukai varieties: linear spaces and periods

    INPUT:

    - ``n`` - the dimension, where $n=1,\\ldots,6$
    """

    assert n in range(
        1, 7
    ), """There is no Gushel--Mukai variety of this
        dimension"""

    if n == 1:
        return curve(6)
    if n == 2:
        return K3()
    if n == 3:
        return curve(10)(1) + lefschetz() ** 0 + lefschetz() ** 3
    if n == 4:
        return K3()(1) + lefschetz() ** 0 + 2 * lefschetz() ** 2 + lefschetz() ** 4
    if n == 5:
        return curve(10)(2) + Pn(5)
    return K3()(2) + lefschetz() ** 3 + Pn(6)


def fano_variety_lines_cubic(n):
    r"""
    Hodge diamond for the Fano variety of lines on a smooth $n$-dimensional
    cubic hypersurface.

    This follows from the "beautiful formula" or X-F(X)-relation due to
    Galkin--Shinder, Theorem 5.1 of [1405.5154v2].

    * [1405.5154v2] Galkin--Shinder, The Fano variety of lines and rationality
      problem for a cubic hypersurface

    INPUT:

    - ``n`` -- the dimension, where ``n`` is at least 2

    EXAMPLES:

    There are 27 lines on a cubic surface::

        sage: from diamond import *
        sage: fano_variety_lines_cubic(2) == 27*point()
        True

    The Fano surface of lines on a cubic threefold is a surface of general type::

        sage: fano_variety_lines_cubic(3) == surface(10, 5, 25)
        True

    The Fano fourfold of lines on a cubic fourfold is deformation equivalent
    to the Hilbert square on a K3 surface::

        sage: fano_variety_lines_cubic(4) == hilbn(K3(), 2)
        True

    """
    assert n >= 2

    X = hypersurface(3, n)
    return (hilbtwo(X) - Pn(n) * X)(-2)


def Mzeronbar(n):
    r"""
    Hodge diamond for the moduli space of $n$-pointed stable curves of genus 0

    Taken from (0.12) in [MR1363064]. Keel's original paper has a recursion
    on page 550, but that seems to not work.

    * [MR1363064] Manin, Generating functions in algebraic geometry and sums
      over trees

    EXAMPLES:

    The first few cases are a point, the projective line, and the blowup
    of $\\mathbb{P}^2$ in 4 points::

        sage: from diamond import *
        sage: Mzeronbar(3) == point()
        True
        sage: Mzeronbar(4) == Pn(1)
        True
        sage: Mzeronbar(5) == Pn(2).blowup(4*point())
        True

    """
    assert n >= 2

    x = HodgeDiamond.x
    y = HodgeDiamond.y

    def Manin(n):
        if n in [2, 3]:
            return HodgeDiamond.R.one()
        else:
            return Manin(n - 1) + x * y * sum(
                binomial(n - 2, i) * Manin(i + 1) * Manin(n - i)
                for i in range(2, n - 1)
            )

    return HodgeDiamond.from_polynomial(Manin(n))


# make Theta into slope stability
def slope(Theta):
    r"""Helper function to turn stability for quiver representations in slope stability"""
    return lambda d: sum(a * b for (a, b) in zip(Theta, d)) / sum(d)


def quiver_moduli(Q, d, **kwargs):
    r"""
    Hodge diamond for the moduli space of semistable quiver representations
    for a quiver Q, dimension vector d, and slope-stability condition mu.

    Taken from Corollary 6.9 of [MR1974891]

    * [MR1974891] Reineke, The Harder-Narasimhan system in quantum groups and cohomology of quiver moduli.

    INPUT:

    - ``Q`` -- adjacency matrix of an acyclic quiver

    - ``d`` -- dimension vector

    - ``mu`` -- stability condition, these can be produced using :func:`slope`,
      if left unspecified then the canonical stability condition is used

    The canonical stability condition for dimension vector `d` is given by
    the antisymmetrised Euler form pairing with `d`.

    EXAMPLES:

    Let's consider moduli spaces for the Kronecker quiver::

        sage: from diamond import *
        sage: def kronecker(d): return matrix([[0, d], [0, 0]])

    For the 2-Kronecker quiver and dimension vector `(1,1)` a representation
    is given by 2 scalars, and the stability condition `(1,-1)` encodes that
    they are not both zero. This way we obtain the projective line::

        sage: quiver_moduli(kronecker(2), (1, 1), mu=slope((1, -1))) == Pn(1)
        True

    There is only one relevant stability chamber here, which contains the
    canonical stability condition. Thus we get the same result not specifying
    the stability function::

        sage: quiver_moduli(kronecker(2), (1, 1)) == Pn(1)
        True

    Similar to the first example, the $d$-Kronecker quiver gives rise to
    projective spaces::

        sage: all(quiver_moduli(kronecker(d), (1, 1)) == Pn(d - 1) for d in range(3, 10))
        True

    We can also realise Grassmannians using the $d$-Kronecker quiver, for
    dimension vector $(1,k)$ and stability condition $(k,1)$ (or the canonical
    stability condition) we get the Grassmannian $\\operatorname{Gr}(k,d)$::

        sage: quiver_moduli(kronecker(4), (1, 2), mu=slope((2, 1))) == grassmannian(2, 4)
        True
        sage: quiver_moduli(kronecker(7), (1, 3)) == grassmannian(3, 7)
        True

    Any stability function in the same chamber gives the same variety::

        sage: quiver_moduli(kronecker(7), (1, 3)) == quiver_moduli(kronecker(7), (1, 3), mu=slope((1, -1)))
        True

    The following is an example of wall-crossing, and thus different varieties::

        sage: M = matrix([[0, 1, 1], [0, 0, 2], [0, 0, 0]])
        sage: quiver_moduli(M, (1, 1, 1)) == Pn(2).blowup(point())
        True
        sage: quiver_moduli(M, (1, 1, 1), mu=slope((1, 0, 0))) == Pn(2)
        True


    The flag variety $\\operatorname{Fl}(n,r_1,\\ldots,r_s)$ is also a quiver
    moduli space, for $\\mathrm{A}_s$ quiver prefixed with an $n$-Kronecker
    quiver, dimension vector $(1,r_1,\\ldots,r_s)$ with $r_1>\\ldots>r_s$
    and stability condition the indicator function at the first vertex::

        sage: def flags(n, s): return matrix(ZZ, s+1, s+1, lambda i,j: 0 if i != j-1 else (n if i == 0 else 1))
        sage: quiver_moduli(flags(4, 1), (1, 1)) == Pn(3)
        True
        sage: quiver_moduli(flags(3, 2), (1, 2, 1)) == fano_threefold(2, 32)
        True
        sage: quiver_moduli(flags(5, 3), (1, 4, 3, 1)) == partial_flag_variety("A4", [2])
        True
    """
    K = FunctionField(QQ, "v")
    v = K.gen(0)
    # we will use v, not v^2, throughout

    # solve Ax=b for A upper triangular via back substitution
    def solve(A, b):
        assert A.is_square() and A.nrows() == len(b)

        n = len(b) - 1
        x = [0] * (n + 1)

        # start
        x[n] = b[n] / A[n, n]

        # induct
        for i in range(n - 1, -1, -1):
            x[i] = (b[i] - sum([A[i, j] * x[j] for j in range(i + 1, n + 1)])) / A[i, i]

        return x

    @cached_function
    def GL(n):
        r"""Cardinality of general linear group $\\mathrm{GL}_n(\\mathbb{F}_v)$"""
        return prod([v**n - v**k for k in range(n)])

    # not caching this one seems faster
    # @cached_function
    def Rd(Q, d):
        """Cardinality of ``R_d`` from Definition 3.1"""
        return v ** sum([d[i] * d[j] * Q[i, j] for i, j in Q.dict()])

    @cached_function
    def Gd(d):
        """Cardinality of ``G_d`` from Definition 3.1"""
        return prod([GL(di) for di in d])

    def Id(Q, d, mu):
        """Returns the indexing set from Corollary 5.5

        These are the dimension vectors smaller than ``d`` whose slope is bigger
        than ``d``, together with the zero dimension vector and ``d`` itself.
        """
        # all possible dimension vectors e <= d
        E = cartesian_product([range(di + 1) for di in d])

        # predicate from Corollary 5.5, E[0] is the zero dimension vector
        return [E[0]] + list(filter(lambda e: mu(e) > mu(d), E[1:])) + [d]

    def Td(Q, d, mu):
        """Returns the upper triangular transfer matrix from Corollary 5.5"""
        # Euler form
        chi = matrix.identity(len(d)) - Q
        # indexing set for the transfer matrix
        I = Id(Q, d, mu)
        # make them vectors now so that we only do it once
        I = list(map(vector, I))

        def entry(Q, e, f):
            """Entry of the transfer matrix, as per Corollary 6.9"""
            fe = f - e

            if all(fei >= 0 for fei in fe):
                return v ** (-fe * chi * e) * Rd(Q, fe) / Gd(fe)
            return 0

        T = matrix(K, len(I), len(I))

        for i, Ii in enumerate(I):
            for j in range(i, len(I)):  # upper triangular
                T[i, j] = entry(Q, Ii, I[j])

        return T

    Q = DiGraph(matrix(Q))
    # see Section 2
    assert Q.is_directed_acyclic(), "Q needs to be acyclic"

    # see Definition 6.3 and following lemmas
    assert gcd(d) == 1, "dimension vector is not coprime"

    # if mu is not provided we resort to the canonical stability condition
    mu = kwargs.get("mu", None)
    if mu is None:
        A = matrix.identity(len(d)) - Q.adjacency_matrix()
        mu = slope(vector(d) * A - A * vector(d))

    # (0,d)-entry of the inverse of the transfer matrix, as per Corollary 6.9
    T = Td(Q.adjacency_matrix(), d, mu)
    # doing a back-substitution is faster
    result = solve(T, [0] * (T.nrows() - 1) + [1])[0] * (1 - v)
    # result = T.inverse()[0,-1] * (1 - v)

    assert result.denominator() == 1, "result needs to be a polynomial"

    x = HodgeDiamond.x
    y = HodgeDiamond.y

    result = HodgeDiamond.R(result.numerator().subs(v=x * y))

    return HodgeDiamond.from_polynomial(result)


def fano_threefold(rho, ID):
    r"""
    Hodge diamond of a Fano threefold

    INPUT:

    - ``rho`` - Picard rank

    - ``ID`` - numbering from the Mori-Mukai classification

    EXAMPLES:

    The 17th Fano 3-fold of rank 1 is projective threespace::

        sage: from diamond import *
        sage: print(fano_threefold(1, 17))
                      1
                  0       0
              0       1       0
          0       0       0       0
              0       1       0
                  0       0
                      1

    The 4th Fano 3-fold of rank 1 is an intersection of 3 quadrics::

        sage: fano_threefold(1, 4) == complete_intersection((2, 2, 2), 3)
        True

    The 27th Fano 3-fold of rank 3 is the triple product of projective lines::

        sage: fano_threefold(3, 27) == Pn(1)**3
        True

    """
    h12 = {
        (1, 1): 52,
        (1, 2): 30,
        (1, 3): 20,
        (1, 4): 14,
        (1, 5): 10,
        (1, 6): 7,
        (1, 7): 5,
        (1, 8): 3,
        (1, 9): 2,
        (1, 10): 0,
        (1, 11): 21,
        (1, 12): 10,
        (1, 13): 5,
        (1, 14): 2,
        (1, 15): 0,
        (1, 16): 0,
        (1, 17): 0,
        (2, 1): 22,
        (2, 2): 20,
        (2, 3): 11,
        (2, 4): 10,
        (2, 5): 6,
        (2, 6): 9,
        (2, 7): 5,
        (2, 8): 9,
        (2, 9): 5,
        (2, 10): 3,
        (2, 11): 5,
        (2, 12): 3,
        (2, 13): 2,
        (2, 14): 1,
        (2, 15): 4,
        (2, 16): 2,
        (2, 17): 1,
        (2, 18): 2,
        (2, 19): 2,
        (2, 20): 0,
        (2, 21): 0,
        (2, 22): 0,
        (2, 23): 1,
        (2, 24): 0,
        (2, 25): 1,
        (2, 26): 0,
        (2, 27): 0,
        (2, 28): 1,
        (2, 29): 0,
        (2, 30): 0,
        (2, 31): 0,
        (2, 32): 0,
        (2, 33): 0,
        (2, 34): 0,
        (2, 35): 0,
        (2, 36): 0,
        (3, 1): 8,
        (3, 2): 3,
        (3, 3): 3,
        (3, 4): 2,
        (3, 5): 0,
        (3, 6): 1,
        (3, 7): 1,
        (3, 8): 0,
        (3, 9): 3,
        (3, 10): 0,
        (3, 11): 1,
        (3, 12): 0,
        (3, 13): 0,
        (3, 14): 1,
        (3, 15): 0,
        (3, 16): 0,
        (3, 17): 0,
        (3, 18): 0,
        (3, 19): 0,
        (3, 20): 0,
        (3, 21): 0,
        (3, 22): 0,
        (3, 23): 0,
        (3, 24): 0,
        (3, 25): 0,
        (3, 26): 0,
        (3, 27): 0,
        (3, 28): 0,
        (3, 29): 0,
        (3, 30): 0,
        (3, 31): 0,
        (4, 1): 1,
        (4, 2): 1,
        (4, 3): 0,
        (4, 4): 0,
        (4, 5): 0,
        (4, 6): 0,
        (4, 7): 0,
        (4, 8): 0,
        (4, 9): 0,
        (4, 10): 0,
        (4, 11): 0,
        (4, 12): 0,
        (4, 13): 0,
        (5, 1): 0,
        (5, 2): 0,
        (5, 3): 0,
        (6, 1): 0,
        (7, 1): 0,
        (8, 1): 0,
        (9, 1): 0,
        (10, 1): 0,
    }

    M = matrix(
        [
            [1, 0, 0, 0],
            [0, rho, h12[(rho, ID)], 0],
            [0, h12[(rho, ID)], rho, 0],
            [0, 0, 0, 1],
        ]
    )

    return HodgeDiamond.from_matrix(M, from_variety=True)
