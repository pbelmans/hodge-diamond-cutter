class HodgeDiamond:
  """
  This class implements some methods to work with Hodge diamonds in an implementation-agnostic way.
  """

  # exposes R, x and y for external use
  R.<x, y> = PolynomialRing(ZZ, 2)

  def __init__(self, m):
    r"""
    Constructor for a Hodge diamond if you know what you are doing: just use the given matrix.

    It is probably advised to use the class methods ``HodgeDiamond.from_matrix`` and ``HodgeDiamond.from_polynomial``.

    INPUT:

    - ``m`` -- an integer matrix representing a Hodge diamond
    """
    # matrix representation of the Hodge diamond is used internally
    self._m = m


  @classmethod
  def from_matrix(cls, m, from_variety=False):
    """
    Constructor when given a polynomial

    INPUT:

    - ``m`` -- a square symmetric integer matrix representing a Hodge diamond

    - ``from_variety`` (default: False) -- whether a check should be performed that it comes from a variety
    """
    diamond = cls(matrix(m))
    diamond.matrix = matrix(m)
    diamond.__normalise() # get rid trailing zeroes from the diamond

    if from_variety:
      assert diamond.arises_from_variety(), "The matrix does not satisfy the conditions satisfied by the Hodge diamond of a smooth projective variety."

    return diamond


  @classmethod
  def from_polynomial(cls, f, from_variety=False):
    """
    Constructor when given a Hodge--Poincaré polynomial

    INPUT:

    - ``f`` -- an integer polynomial in the ring ``HodgeDiamond.R`` representing the Hodge--Poincaré polynomial

    - ``from_variety`` (default: False) -- whether a check should be performed that it comes from a variety
    """
    diamond = cls(cls.__to_matrix(f))

    if from_variety:
      assert diamond.arises_from_variety(), "The matrix does not satisfy the conditions satisfied by the Hodge diamond of a smooth projective variety."

    return diamond


  @property
  def polynomial(self):
    """Getter for the Hodge--Poincaré polynomial"""
    return sum([self.matrix[i,j] * self.x^i * self.y^j for (i, j) in cartesian_product([range(self.matrix.nrows()), range(self.matrix.ncols())])])


  @polynomial.setter
  def polynomial(self, f):
    """Setter for the Hodge--Poincaré polynomial"""
    self.matrix = HodgeDiamond.__to_matrix(f)


  @property
  def matrix(self):
    """Getter for the Hodge diamond as a matrix"""
    return self._m


  @matrix.setter
  def matrix(self, m):
    """Setter for the Hodge diamond as a matrix

    Performs the following consistency checks: it needs to be an square, symmetric integer matrix
    """
    m = matrix(m)

    assert m.base_ring() == ZZ, "Entries need to be integers"
    assert m.is_symmetric()
    assert m.is_square()

    self._m = matrix(m)
    self.__normalise()


  @staticmethod
  def __to_matrix(f):
    """Convert Hodge-Poincaré polynomial to matrix representation"""
    assert f in HodgeDiamond.R

    if f.is_zero():
      m = matrix([[0]])
    else:
      # deal with the size of the diamond in this way because of the following example:
      # X = complete_intersection(5, 3)
      # X*X - hilbtwo(X)
      d = max([max(e) for e in f.exponents()]) + 1
      m = matrix(d)

      for (i, j) in cartesian_product([range(d), range(d)]):
        m[i, j] = f.monomial_coefficient(HodgeDiamond.x^i * HodgeDiamond.y^j)

    return m


  def __size(self):
    """Internal method to determine the (relevant) size of the Hodge diamond"""
    return self.matrix.ncols() - 1


  def __normalise(self):
    """Internal method to get rid of trailing zeros"""
    self.matrix = HodgeDiamond.__to_matrix(self.polynomial)


  """
  Output methods
  """
  def __repr__(self):
    """Output diagnostic information"""
    return "Hodge diamond of size {} and dimension {}".format(self.__size() + 1, self.dimension())


  def __str__(self):
    """Pretty print Hodge diamond"""
    return str(self.pprint())


  def __table(self):
    """Generate a table object for the Hodge diamond"""
    d = self.__size()
    T = []

    if self.is_zero():
      T = [[0]]

    for i in range(2*d + 1):
      row = [""]*(abs(d - i))

      for j in range(max(0, i - d), min(i, d) + 1):
        row.extend([self.matrix[j, i - j], ""])

      T.append(row)

    # padding all rows to full length
    for i in range(len(T)):
      T[i].extend([""]*(2*d - len(T[i]) + 1))

    return table(T, align="center")


  def pprint(self, output="table"):
    """Pretty print the Hodge diamond"""
    if output == "table":
      return self.__table()
    else:
      return self.polynomial


  """
  Overloaded operators
  """
  def __add__(self, other):
    """Add Hodge diamonds"""
    return HodgeDiamond.from_polynomial(self.polynomial + other.polynomial)


  def __radd__(self, other):
    """Add Hodge diamonds"""
    # to make sum() work as intended, which by default starts with 0
    if other == 0: return self

    return HodgeDiamond.from_polynomial(self.polynomial + other.polynomial)


  def __sub__(self, other):
    """Add Hodge diamonds"""
    return HodgeDiamond.from_polynomial(self.polynomial - other.polynomial)


  def __mul__(self, other):
    """Multiply two Hodge diamonds"""
    if not isinstance(other, HodgeDiamond): # in the rare case someone does X*3 instead of 3*X
      return other * self

    return HodgeDiamond.from_polynomial(self.polynomial * other.polynomial)


  def __rmul__(self, factor):
    """Multiply a Hodge diamond with a factor"""
    return HodgeDiamond.from_polynomial(factor * self.polynomial)


  def __pow__(self, i):
    """Raise a Hodge diamond to a power"""
    return HodgeDiamond.from_polynomial(self.polynomial ** i)


  def __eq__(self, other):
    """Check whether two Hodge diamonds are equal"""
    return self.polynomial == other.polynomial


  def __ne__(self, other):
    """Check whether two Hodge diamonds are not equal"""
    return not self == other


  def __call__(self, i, y=None):
    """
    Either
    - twist by a power of the Lefschetz Hodge diamond, when a single parameter is present
    - evaluate the Hodge-Poincaré polynomial

    Negative values are allowed to untwist, up to the appropriate power.
    """
    if y == None:
      assert i >= -self.lefschetz_power()

      return HodgeDiamond.from_polynomial(self.R(self.polynomial * self.x^i * self.y^i))
    else:
      x = i

      return self.polynomial(x, y)


  def __getitem__(self, index):
    """Get (p, q)th entry of Hodge diamond or the ith row of the Hodge diamond"""
    # first try it as (p, q)
    try:
      (p, q) = index
      return self.matrix[p, q]
    # now we assume it's an integer
    except TypeError:
      # we could do something smarter, but this is it for now
      return [self.matrix[p, index - p] for p in range(index + 1)]


  """
  Consistency checks
  """
  def __is_positive(self):
    """Check whether all entries are positive integers"""
    return all([hpq >= 0 for hpq in self.matrix.coefficients()])


  def is_symmetric(self):
    """Check whether the Hodge diamond is symmetric"""
    return self.matrix.is_symmetric()


  def is_hodge_symmetric(self):
    """Check whether the Hodge diamond satisfies Hodge symmetry: h^{p,q}(X)=h^{q,p}(X)"""
    d = self.__size()
    return all([self.matrix[p, q] == self.matrix[q, p] for (p, q) in cartesian_product([range(d + 1), range(d + 1)])])


  def is_serre_symmetric(self):
    """Check whether the Hodge diamond satisfies Serre symmetry: h^{p,q}(X)=h^{n-p,n-q}(X)"""
    d = self.__size()
    return all([self.matrix[p, q] == self.matrix[d - p, d - q] for (p, q) in cartesian_product([range(d + 1), range(d + 1)])])


  """
  Associated invariants
  """
  def betti(self):
    """Betti numbers of the Hodge diamond"""
    d = self.__size()
    return [sum([self.matrix[j, i - j] for j in range(max(0, i - d), min(i, d) + 1)]) for i in range(2*d + 1)]


  def euler(self):
    """Euler characteristic of the Hodge diamond"""
    return sum([(-1)^i * bi for i, bi in enumerate(self.betti())])


  def hirzebruch(self):
    r"""Hirzebruch's \chi_y genus"""
    return self.polynomial.subs(x=-1)


  def homological_unit(self):
    """Dimensions of H^*(X,O_X)"""
    return self.matrix.row(0)


  def hochschild(self):
    """Dimensions of the Hochschild homology"""
    d = self.__size()
    return HochschildHomology([sum([self.matrix[d - i + j, j] for j in range(max(0, i - d), min(i, d) + 1)]) for i in range(2*d + 1)])


  def hh(self):
    """Shorthand for ```HodgeDiamond.hochschild()```"""
    return self.hochschild()




  """
  Properties
  """
  def arises_from_variety(self):
    """Check whether the Hodge diamond can arise from a smooth projective variety

    The constraints are:
    - satisfy Hodge symmetry
    - satisfy Serre symmetry
    - there is no Lefschetz twist
    """
    return self.is_hodge_symmetric() and self.is_serre_symmetric() and self.lefschetz_power() == 0


  def is_zero(self):
    """Check whether the Hodge diamond is identically zero"""
    return self.matrix.is_zero()


  def lefschetz_power(self):
    """Return the twist by the Lefschetz motive that is present

    In other words, we see how divisible the Hodge--Poincaré polynomial is with respect to the monomial x^iy^i"""
    if self.is_zero(): return 0

    i = 0
    while (self.x^i * self.y^i).divides(self.polynomial):
      i = i + 1

    return i - 1


  def dimension(self):
    """Dimension of the Hodge diamond

    This is not just naively the size of the Hodge diamond: zero padding is taken into account"""
    assert self.is_symmetric()
    if self.is_zero():
      return -1
    else:
      return max([i for i in range(self.matrix.ncols()) if not self.matrix.column(i).is_zero()]) - self.lefschetz_power()


  def level(self):
    """Compute the level (or complexity) of the Hodge diamond"""
    return max([abs(p - q) for (p, q) in [m.degrees() for m in self.polynomial.monomials()]])


  """
  Constructions
  """
  def blowup(self, other, codim=None):
    """Compute Hodge diamond of blowup"""
    # let's guess the codimension
    if codim == None:
      codim = self.dimension() - other.dimension()

    return self + sum([other(i) for i in range(1, codim)])

  def bundle(self, rank):
    """Compute the Hodge diamond of a projective bundle"""
    return sum([self(i) for i in range(rank)])


class HochschildHomology:
  """
  This class implements some methods to work with (the dimensions of) Hochschild homology spaces, associated to the `HodgeDiamond` class.
  """

  # exposes R and t for external use
  R.<t> = LaurentPolynomialRing(ZZ)

  def __init__(self, L):
    r"""
    Constructor for Hochschild homology dimensions of smooth and proper dg categories, so that Serre duality holds.

    INPUT:

    - ``L`` -- a list of integers of length 2n+1 representing $\mathrm{HH}_{-n}$ to $\mathrm{HH}_n$, such that ``L[i] == L[2n - i]``
    """
    assert len(L) % 2 == 1, "length needs to be odd, to reflect Serre duality"
    assert all([L[i] == L[len(L) - i - 1] for i in range(len(L))]), "Serre duality is not satisfied"

    self._L = L


  @classmethod
  def from_list(cls, L):
    r"""
    Constructor for Hochschild homology dimensions from a list.

    INPUT:

    - ``L`` -- a list of integers representing $\mathrm{HH}_{-n}$ to $\mathrm{HH}_n$
    """
    return cls(L)


  @classmethod
  def from_positive(cls, L):
    """
    Constructor for Hochschild homology dimensions from a list when only the positive part is given.

    INPUT:

    - ``L`` -- a list of integers representing ``HH_0`` to ``HH_n``
    """
    double = list(reversed(self._L))[:-1] + L

    return cls(double)


  @classmethod
  def from_polynomial(cls, f):
    """
    Constructor for Hochschild homology dimensions from Hochschild--Poincaré Laurent polynomial

    INPUT

    - ``f`` -- the Hochschild--Poincaré Laurent polynomial
    """
    if f.is_zero(): return cls([0])

    L = [f.dict()[i] if i in f.exponents() else 0 for i in range(-f.degree(), f.degree() + 1)]

    return cls(L)


  @property
  def polynomial(self):
    return HochschildHomology.R(sum([self[i] * (self.t)^i for i in range(-self.dimension(), self.dimension() + 1)]))


  """
  Output methods
  """
  def __repr__(self):
    return "Hochschild homology vector of dimension {}".format(self.dimension())


  def __str__(self):
    return str(self.pprint())


  def __table(self):
    if self.is_zero():
      return table([[0], [0]], header_row=True)

    indices = range(-self.dimension(), self.dimension() + 1)

    return table([indices, [self[i] for i in indices]], header_row=True)


  def pprint(self, output="table"):
    if output == "table":
      return self.__table()
    else:
      return self._L


  """
  Properties
  """
  def dimension(self):
    r"""Largest index ``i`` such that $\mathrm{HH}_i\neq 0$"""
    if self.is_zero(): return -1

    return (len(self._L) // 2) - min([i for (i, d) in enumerate(self._L) if d != 0])


  def is_zero(self):
    return set(self._L) == set([0])



  """
  Overloaded operators
  """
  def __add__(self, other):
    return HochschildHomology.from_polynomial(self.polynomial + other.polynomial)


  def __radd__(self, other):
    # to make sum() work as intended, which by default starts with 0
    if other == 0: return self

    return HochschildHomology.from_polynomial(self.polynomial + other.polynomial)


  def __sub__(self, other):
    return HochschildHomology.from_polynomial(self.polynomial - other.polynomial)


  def __mul__(self, other):
    if not isinstance(other, HochschildHomology): # in the rare case someone does X*3 instead of 3*X
      return other * self

    return HochschildHomology.from_polynomial(self.polynomial * other.polynomial)


  def __rmul__(self, factor):
    return HochschildHomology.from_polynomial(factor * self.polynomial)


  def __pow__(self, i):
    return HochschildHomology.from_polynomial(self.polynomial ** i)


  def __eq__(self, other):
    return self.polynomial == other.polynomial


  def __ne__(self, other):
    return not self == other


  def __getitem__(self, i):
    if i > len(self._L) // 2:
      return 0
    else:
      return self._L[len(self._L) // 2 - i]


  def __iter__(self):
    return self._L.__iter__()


  def symmetric_power(self, k):
    """
    Hochschild homology of the Ganter--Kapranov symmetric power of a smooth and proper dg category

    This is possibly only a heuristic (I didn't check for proofs in the literature) based on the decomposition of Hochschild homology for a quotient stack, as discussed in the paper of Polishchuk--Van den Bergh
    """
    def summand(f, k):
      assert [c > 0 for c in f.coefficients()]

      t = HochschildHomology.t

      # trivial case
      if f.number_of_terms() == 0:
        return f
      # base case
      elif f.number_of_terms() == 1:
        i = f.exponents()[0]
        a = f.coefficients()[0]

        if i % 2 == 0:
          return binomial(a + k - 1, k) * t^(k * i)
        else:
          return binomial(a, k) * t^(k * i)
      # general case
      else:
        # splitting f into monomial and the difference
        g = f.coefficients()[0] * t^(f.exponents()[0])
        h = f - g

        return sum([summand(g, j) * summand(h, k - j) for j in range(k + 1)])

    # see the object C^{(\lambda)} in the Polishchuk--Van den Bergh paper
    return HochschildHomology.from_polynomial(sum([product([summand(self.polynomial, ri) for ri in P.to_exp()]) for P in Partitions(k)]))


  def sym(self, k):
    """Shorthand for ```HochschildHomology.symmetric_power```"""
    return self.symmetric_power(k)


"""Hodge diamond for the point"""
point = HodgeDiamond.from_matrix(matrix([[1]]), from_variety=True)


"""Hodge diamond for the Lefschetz motive"""
lefschetz = point(1)


def Pn(n):
  """
  Hodge diamond for projective space of dimension ``n``
  """
  return HodgeDiamond.from_matrix(matrix.identity(n + 1), from_variety=True)


def curve(g):
  """
  Hodge diamond for a curve of genus ``g``
  """
  return HodgeDiamond.from_matrix(matrix([[1, g], [g, 1]]), from_variety=True)



def surface(pg, q, h11):
  """
  Hodge diamond for a surface with geometric genus ``pg``, irregularity ``q`` and middle Hodge numbers ``h11``
  """
  return HodgeDiamond.from_matrix(matrix([[1, q, pg], [q, h11, q], [pg, q, 1]]), from_variety=True)



def symmetric_power(n, g):
  """
  Hodge diamond for the ``n``th symmetric power of a genus ``g`` curve

  For the proof, see example 1.1(1) of [MR2777820]. An earlier reference, probably in Macdonald, should exist.

  * [MR2777820] Laurentiu--Schuermann, Hirzebruch invariants of symmetric products. Topology of algebraic varieties and singularities, 163–177, Contemp. Math., 538, Amer. Math. Soc., 2011.
  """
  def hpq(g, n, p, q):
    assert(p <= n)
    assert(q <= n)

    if p <= q and p + q <= n:
      return sum([binomial(g, p - k) * binomial(g, q - k) for k in range(p + 1)])
    if p > q:
      return hpq(g, n, q, p)
    if p + q > n:
      return hpq(g, n, n - p, n - q)

  M = matrix(n + 1)
  for (i,j) in cartesian_product([range(n+1), range(n+1)]):
    M[i,j] = hpq(g, n, i, j)

  return HodgeDiamond.from_matrix(M, from_variety=True)


def jacobian(g):
  """
  Hodge diamond for the Jacobian of a genus $g$ curve

  This description is standard, and follows from the fact that the cohomology is the exterior power of the H^1, as graded bialgebras. See e.g. proposition 7.27 in the Edixhoven--van der Geer--Moonen book in progress on abelian varieties.
  """
  M = matrix(g + 1)
  for (i,j) in cartesian_product([range(g+1), range(g+1)]):
    M[i,j] = binomial(g, i) * binomial(g, j)

  return HodgeDiamond.from_matrix(M, from_variety=True)


def abelian(g):
  """
  Hodge diamond for an abelian variety of dimension $g$

  This is just an alias for ``jacobian(g)``.
  """
  return jacobian(g, from_variety=True)


def moduli_vector_bundles(r, d, g):
  """
  Hodge diamond for the moduli space of vector bundles of rank ``r`` and fixed determinant of degree ``d`` on a curve of genus ``g``

  For the proof, see corollary 5.1 of [MR1817504].

  * [MR1817504] del Baño, On the Chow motive of some moduli spaces. J. Reine Angew. Math. 532 (2001), 105–132.

  If the Hodge diamond for the moduli space with non-fixed determinant of degree d is required, this can be obtained by

      sage: jacobian(g) * moduli_vector_bundles(r, d, g)

  """
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
    return (-1)^(len(C)-1) * ((1+x)^g * (1+y)^g)^(len(C)-1) / (1-x*y)^(len(C)-1) # already corrected the factor from the Jacobian

  def two(C, g):
    """Return the second factor in del Bano's formula."""
    return prod([prod([(1 + x^i*y^(i+1))^g * (1 + x^(i+1)*y^i)^g / ((1 - x^i*y^i) * (1 - x^(i+1)*y^(i+1))) for i in range(1, C[j])]) for j in range(len(C))])

  def three(C, g):
    """Return the third factor in del Bano's formula."""
    return prod([1 / (1-(x*y)^(C[j] + C[j+1])) for j in range(len(C) - 1)])

  def four(C, d, g):
    """Return the fourth factor in del Bano's formula."""
    exponent = sum([sum([(g-1) * C[i] * C[j] for i in range(j)]) for j in range(len(C))]) + sum([(C[i] + C[i+1]) * bracket(-(sum(C[0:i+1]) * d / C.size())) for i in range(len(C) - 1)])
    return (x*y)^exponent

  return HodgeDiamond.from_polynomial(R(sum([one(C, g) * two(C, g) * three(C, g) * four(C, d, g) for C in Compositions(r)])), from_variety=True)


def quot_scheme_curve(g, n, r):
  """
  Hodge diamond for the Quot scheme of zero-dimensional quotients of length ``r`` of a vector bundle of rank ``r`` on a curve of genus ``g``

  For the proof, see proposition 4.5 of [1907.00826] (or rather, the reference [Bif89] in there)

  * [1907.00826] Bagnarol--Fantechi--Perroni, On the motive of zero-dimensional Quot schemes on a curve
  """
  def dn(P):
    # shift in indexing because we start at 0
    return sum([i * ni for (i, ni) in enumerate(P)])

  return sum([product([symmetric_power(ni, g) for ni in P])(dn(P)) for P in IntegerVectors(n, r)])


def hilbtwo(X):
  """
  Hodge diamond for the Hilbert square of any smooth projective variety

  For the proof, see e.g. lemma 2.6 of [MR2506383], or the corollary on page 507 of [MR1382733] (but it can be said to be classical).

  * [MR2506383] Muñoz--Ortega--Vázquez-Gallo, Hodge polynomials of the moduli spaces of triples of rank (2,2). Q. J. Math. 60 (2009), no. 2, 235–272.
  * [MR1382733] Cheah, On the cohomology of Hilbert schemes of points, J. Algebraic Geom. 5 (1996), no. 3, 479-511

  INPUT:

  -  ``X`` - Hodge diamond of the smooth projective variety
  """
  assert X.arises_from_variety()

  d = X.dimension()

  return HodgeDiamond.from_polynomial(X.R(1/2 * ((X.polynomial)^2 + X.polynomial(-X.x^2, -X.y^2)) + sum([X.x^(i+1) * X.y^(i+1) * X.polynomial for i in range(d - 1)])), from_variety=True)


def hilbthree(X):
  """
  Hodge diamond for the Hilbert cube of any smooth projective variety

  The corollary on page 507 of [MR1382733].

  * [MR1382733] Cheah, On the cohomology of Hilbert schemes of points, J. Algebraic Geom. 5 (1996), no. 3, 479-511

  INPUT:

  -  ``X`` - Hodge diamond of the smooth projective variety
  """
  assert X.arises_from_variety()

  d = X.dimension()

  return HodgeDiamond.from_polynomial(X.R(
    1/6 * (X^3).polynomial
    + 1/2 * X.polynomial * X.polynomial(-X.x^2, -X.y^2)
    + 1/3 * X.polynomial(X.x^3, X.y^3)
    + sum([(X^2)(i).polynomial for i in range(1, d)])
    + sum([X(i+j).polynomial for i in range(1, d) for j in range(i, d)])
    ), from_variety=True)


def generalisedkummer(n):
  """
  Hodge diamond of the ``n``th generalised Kummer variety

  The first generalised Kummer (i.e. when ``n`` is 1) is the point, the second the Kummer K3, etc.

  For the proof, see corollary 1 of [MR1219901].

  * [MR1219901] Göttsche--Soergel, Perverse sheaves and the cohomology of Hilbert schemes of smooth algebraic surfaces. Math. Ann. 296 (1993), no. 2, 235–245.
  """
  def product(n):
    x = HodgeDiamond.x
    y = HodgeDiamond.y

    return sum([gcd(a.keys())^4 * (x*y)^(n - sum(a.values())) * prod([sum([prod([1 / (j^b[j] * factorial(b[j])) * ((1-x^j) * (1-y^j))^(2*b[j]) for j in b.keys()]) for b in [b.to_exp_dict() for b in Partitions(a[i])]]) for i in a.keys()]) for a in [a.to_exp_dict() for a in Partitions(n)]])(-x, -y)

  return HodgeDiamond.from_polynomial(HodgeDiamond.R(product(n) / product(1)), from_variety=True)


"""
Hodge diamond for O'Grady's exceptional 6-dimensional hyperkähler variety

For the proof, see theorem 1.1 of [MR3798592].

* [MR3798592] Mongardi--Rapagnetta--Saccà, The Hodge diamond of O'Grady's six-dimensional example. Compos. Math. 154 (2018), no. 5, 984–1013.
"""
ogrady6 = HodgeDiamond.from_matrix(
    [[1, 0, 1, 0, 1, 0, 1],
    [0, 6, 0, 12, 0, 6, 0],
    [1, 0, 173, 0, 173, 0, 1],
    [0, 12, 0, 1144, 0, 12, 0],
    [1, 0, 173, 0, 173, 0, 1],
    [0, 6, 0, 12, 0, 6, 0],
    [1, 0, 1, 0, 1, 0, 1]], from_variety=True)


"""
Hodge diamond for O'Grady's exceptional 10-dimensional hyperkähler variety

For the proof, see theorem A of [1905.03217]

* [1905.03217] de Cataldo--Rapagnetta--Saccà, The Hodge numbers of O'Grady 10 via Ngô strings
"""
ogrady10 = HodgeDiamond.from_matrix(
    [[1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
    [0, 22, 0, 22, 0, 23, 0, 22, 0, 22, 0],
    [1, 0, 254, 0, 276, 0, 276, 0, 254, 0, 1],
    [0, 22, 0, 2299, 0, 2531, 0, 2299, 0, 22, 0],
    [1, 0, 276, 0, 16490, 0, 16490, 0, 276, 0, 1],
    [0, 23, 0, 2531, 0, 88024, 0, 2531, 0, 23, 0],
    [1, 0, 276, 0, 16490, 0, 16490, 0, 276, 0, 1],
    [0, 22, 0, 2299, 0, 2531, 0, 2299, 0, 22, 0],
    [1, 0, 254, 0, 276, 0, 276, 0, 254, 0, 1],
    [0, 22, 0, 22, 0, 23, 0, 22, 0, 22, 0],
    [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]], from_variety=True)


def hilbn(S, n):
  """
  Hodge diamond for Hilbert scheme of ``n`` points on a smooth projective surface ``S``

  For the proof, see theorem 2.3.14 of [MR1312161].

  * [MR1312161] Göttsche, Hilbert schemes of zero-dimensional subschemes of smooth varieties. Lecture Notes in Mathematics, 1572. Springer-Verlag, Berlin, 1994. x+196 pp.

  INPUT:

  - ``S`` -- Hodge diamond for smooth projective surface

  - ``n`` -- number of points
  """
  assert S.arises_from_variety()
  assert S.dimension() == 2

  R.<a, b, t> = PowerSeriesRing(ZZ, default_prec = (3*n + 1))
  series = R(1)

  # theorem 2.3.14 of Göttsche's book
  for k in range(1, n+1):
    for (p,q) in cartesian_product([range(3), range(3)]):
      series = series * (1 + (-1)**(p+q+1) * a**(p+k-1) * b**(q+k-1) * t**k)**((-1)**(p+q+1) * S[p,q])

  # read off Hodge diamond from the (truncated) series
  M = matrix([[0 for _ in range(2*n + 1)] for _ in range(2*n + 1)])
  for (p, q) in cartesian_product([range(2*n + 1), range(2*n + 1)]):
    monomial = a**min(p, 2*n-q)*b**min(q, 2*n-p)*t**n # use Serre duality
    if monomial in series.coefficients():
      M[p,q] = series.coefficients()[monomial]

  return HodgeDiamond.from_matrix(M, from_variety=True)


def nestedhilbn(S, n):
  """
  Hodge diamond for the nested Hilbert scheme ``S^[n-1,n]`` of a smooth projective surface ``S``

  This is a smooth projective variety of dimension 2n.
  """
  assert S.arises_from_variety()
  assert S.dimension() == 2

  R.<x, y, t> = PowerSeriesRing(ZZ, default_prec = (5*n + 1)) # TODO lower this?

  series = R(1)
  for k in range(1, 3*n): # TODO lower this?
    for (p, q) in cartesian_product([range(3), range(3)]):
      if p + q % 2 == 1:
        series = series * (1 + x**(p+k-1) * y^(q+k-1) * t^k)^(S[p, q])
      else:
        series = series * (1 - x**(p+k-1) * y^(q+k-1) * t^k)^(-S[p, q])

  series = series * R(S.polynomial) * t / (1 - x*y*t)

  # read off Hodge diamond from the (truncated) series
  M = matrix(2*n + 1)
  for (p, q) in cartesian_product([range(2*n + 1), range(2*n + 1)]):
    monomial = x^min(p, 2*n-q) * y^min(q, 2*n-p) * t^n # use Serre duality
    if monomial in series.coefficients():
      M[p,q] = series.coefficients()[monomial]

  return HodgeDiamond.from_matrix(M, from_variety=True)



def complete_intersection(degrees, dimension):
  """
  Hodge diamond for a complete intersection of multidegree ``(d_1,...,d_k)`` in ``P^{n+k}``

  For a proof, see théorème 2.3 of exposé XI in SGA7.

  INPUT:

  - ``degrees`` -- the multidegree, if it is an integer we interpret it as a hypersurface

  - ``dimension`` -- the dimension of the complete intersection (not of the ambient space)
  """
  # hypersurface as complete intersection
  if type(degrees) == sage.rings.integer.Integer:
    degrees = [degrees]

  R.<a, b> = PowerSeriesRing(ZZ, default_prec=dimension+2)
  H = 1/((1+a)*(1+b)) * (prod([((1+a)^di - (1+b)^di) / (a*(1+b)^di - b*(1+a)^di) for di in degrees]) - 1) + 1/(1-a*b)

  middle = [H.coefficients()[a^i * b^(dimension-i)] if a^i * b^(dimension-i) in H.coefficients() else 0 for i in range(dimension + 1)]

  M = matrix.identity(dimension + 1)
  for i in range(dimension + 1):
    M[i, dimension - i] = middle[i]

  return HodgeDiamond.from_matrix(M, from_variety=True)


def hypersurface(degree, dimension):
  """Shorthand for a complete intersection of the given dimension where ``k=1``"""
  return complete_intersection(degree, dimension)


"""Hodge diamond for a K3 surface"""
K3 = complete_intersection(4, 2)


def weighted_hypersurface(degree, weights):
  """
  Hodge diamond for a weighted hypersurface of degree ``d`` in ``P(w_0,...,w_n)``

  Implementation from https://github.com/jxxcarlson/math_research/blob/master/hodge.sage. If I'm not mistaken, the generalisation to weighted complete intersections depends on which polynomials are used, not just the degrees, even if the result is smooth?

  INPUT:

  - ``degree`` -- degree of the hypersurface

  - ``weights`` -- the weights of the weighted projective space, if it is an integer we interpret it as the dimension of P^n
  """
  var('t')

  def poincare_function(L):
    return expand(product([1] + [simplify((1 - t^a) / (1 - t^b)) for (a, b) in L]))

  def jacobian_weights(degree, weights):
    return [(degree - weight, weight) for weight in weights]

  def poincare(degree, weights):
    return poincare_function(jacobian_weights(degree, weights))


  # weights should be interpreted as dimension of unweighted P^n
  if type(weights) == sage.rings.integer.Integer:
    weights = [1]*(weights + 2)

  P = poincare(degree, weights)
  n = len(weights)
  vdeg = sum(weights)
  tdeg = (n + 1) * degree - vdeg + 1
  T = P.taylor(t, 0, tdeg)
  middle = []

  for q in range(0, n - 1):
    middle.append(T.coefficient(t, (q+1)*degree - vdeg))

  # adding in non-primitive cohomology if necessary
  if len(middle) % 2 == 1: middle[len(middle) // 2] += 1

  M = matrix.identity(n - 1)
  for i in range(n - 1):
    M[i, n - i - 2] = middle[i]

  return HodgeDiamond.from_matrix(M, from_variety=True)


def cyclic_cover(ramification_degree, cover_degree, weights):
  """
  Hodge diamond of a cyclic cover of weighted projective space

  Implementation from https://github.com/jxxcarlson/math_research/blob/master/hodge.sage.

  INPUT:

  - ``ramification_degree`` -- degree of the ramification divisor

  - ``cover_degree`` -- size of the cover

  - ``weights`` -- the weights of the weighted projective space, if it is an integer we interpret it as the dimension of ``P^n``
  """
  # weights should be interpreted as dimension of unweighted P^n
  if type(weights) == sage.rings.integer.Integer:
    weights = [1]*(weights + 1)

  weights.append(ramification_degree // cover_degree)

  return weighted_hypersurface(ramification_degree, weights)


def partial_flag_variety(D, I):
  """
  Hodge diamond of a partial flag variety G/P

  This is computed by counting the number of Schubert cells in the appropriate dimension.

  INPUT:

  - ``D`` -- Dynkin type

  - ``I`` -- indices of vertices to be omitted in defining the parabolic subgroup
  """
  W = WeylGroup(D)

  cosets = set([w.coset_representative(index_set=I) for w in W])
  lengths = [len(w.reduced_word()) for w in cosets]
  d = max(lengths)

  diagonal = [0]*(d + 1)

  for i in range(d + 1):
    diagonal[i] = len([l for l in lengths if l == i])

  return HodgeDiamond.from_matrix(diagonal_matrix(diagonal), from_variety=True)


def grassmannian(k, n):
  """
  Hodge diamond of the Grassmannian Gr(k, n) of k-dimensional subspaces in an n-dimensional vector space


  INPUT:

  - ``k`` -- dimension of the subspaces

  - ``n`` -- dimension of the ambient vector space
  """
  return partial_flag_variety("A" + str(n - 1), [i for i in range(1, n) if i != k])


def orthogonal_grassmannian(k, n):
  """
  Hodge diamond of the orthogonal Grassmannian OGr(k, n) of k-dimensional subspaces in an n-dimensional vector space isotropic with respect to a non-degenerate symmetric bilinear form


  INPUT:

  - ``k`` -- dimension of the subspaces

  - ``n`` -- dimension of the ambient vector space
  """
  if n % 2 == 0:
    assert k < n // 2
  else:
    assert k <= n // 2

  if n % 2 == 0:
    # exceptional case: need submaximal parabolic associated to the 2 tails
    if k - 1 == n // 2:
      return partial_flag_variety("D" + str(n // 2), [i for i in range(1, n // 2 - 1)])
    else:
      return partial_flag_variety("D" + str(n // 2), [i for i in range(1, n // 2 + 1) if i != k])
  else:
    return partial_flag_variety("B" + str(n // 2), [i for i in range(1, n // 2 + 1) if i != k])


def symplectic_grassmannian(k, n):
  """
  Hodge diamond of the symplectic Grassmannian SGr(k, n) of k-dimensional subspaces in an n-dimensional vector space isotropic with respect to a non-degenerate skew-symmetric bilinear form


  INPUT:

  - ``k`` -- dimension of the subspaces

  - ``n`` -- dimension of the ambient vector space
  """
  assert n % 2 == 0

  return partial_flag_variety("C" + str(n // 2), [i for i in range(1, n // 2 + 1) if i != k])


def gushel_mukai(n):
  r"""
  Hodge diamond for a smooth $n$-dimensional Gushel--Mukai variety.

  See proposition 3.1 of [1605.05648v3].

  * [1605.05648v3] Debarre--Kuznetsov, Gushel-Mukai varieties: linear spaces and periods

  INPUT:

  - ``n`` -- the dimension, where $n=1,\ldots,6
  """

  assert n in range(1, 7), "There is no Gushel--Mukai variety of this dimension"

  if n == 1:   return curve(6)
  elif n == 2: return K3
  elif n == 3: return curve(10)(1) + point + point(3)
  elif n == 4: return K3(1) + point + 2*point(2) + point(4)
  elif n == 5: return curve(10)(2) + Pn(5)
  else:        return K3(2) + point(3) + Pn(6)


def fano_variety_lines_cubic(n):
  r"""
  Hodge diamond for the Fano variety of lines on an ``n``-dimensional cubic hypersurface.

  This follows from the "beautiful formula" or X-F(X)-relation due to Galkin--Shinder, theorem 5.1 of [1405.5154v2].

  * [1405.5154v2] Galkin--Shinder, The Fano variety of lines and rationality problem for a cubic hypersurface

  INPUT:

  - ``n`` -- the dimension, where ``n`` is at least 2
  """
  assert n >= 2

  X = hypersurface(3, n)
  return (hilbtwo(X) - Pn(n) * X)(-2)

