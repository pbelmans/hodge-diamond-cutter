[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3893509.svg)](https://doi.org/10.5281/zenodo.3893509)

# Hodge diamond cutter

A collection of Python classes and functions in Sage to deal with Hodge diamonds (and Hochschild homology) of smooth projective varieties, together with many constructions.

If you have used this code in any way (including the interactive versions on my blog), please consider citing it as explained on [Zenodo](https://doi.org/10.5281/zenodo.3893509). You can choose to cite a specific version, or the library in general.


## Getting started

It suffices to load ``diamond.sage`` in Sage to get started.


# Functionality

## Varieties and constructions

Below is the list of currently implemented constructions. See their respective documentation if their usage is not clear. They are roughly grouped into themes.
```
zero
point
lefschetz
Pn(n)

curve(genus)
symmetric_power(n, genus)
jacobian(genus)
abelian(dimension)
moduli_vector_bundles(rank, degree, genus)
quot_scheme_curve(genus, length, rank)

fano_variety_intersection_quadrics_odd(g, i)
fano_variety_intersection_quadrics_even(g, i)
fano_variety_lines_cubic(n)

hilbtwo(X)
hilbthree(X)

generalisedkummer(n)
ogrady6
ogrady10

surface(genus, irregularity, h11)
hilbn(surface, n)
nestedhilbn(surface, n)

complete_intersection(degrees, dimension)
hypersurface(degree, dimension)

weighted_hypersurface(degree, weights)
cyclic_cover(ramification, cover, weights)

partial_flag_variety(D, I)
grassmannian(k, n)
orthogonal_grassmannian(k, n)
symplectic_grassmannian(k, n)

gushel_mukai(n)

Mzeronbar(n)
```


## The library
Of course, Hodge diamonds are little more than collections of numbers. But to make manipulating them easy, there is a convenient `HodgeDiamond` class which has many convenient methods, so that manipulating them becomes easy. Describing the whole interface is a bit tedious, let me just give the ones which are not obvious operations.

The `HodgeDiamond` class:
```
HodgeDiamond.betti()
HodgeDiamond.middle()
HodgeDiamond.euler
HodgeDiamond.hirzebruch
HodgeDiamond.homological_unit()
HodgeDiamond.hochschild()

HodgeDiamond.level()

HodgeDiamond.blowup(other, codim)
HodgeDiamond.bundle(rank)
```

The `HochschildHomology` class:
```
HochschildHomology.symmetric_power(k)
```


## Examples

As an example of why it is interesting to consider Hodge diamonds and , let us consider the cubic fourfold. There is an intricate connection to K3 surfaces. Let us try and see this well-known connection.

To see the Hodge diamond of a K3 surface, it suffices to do

```
print(K3)
```

because K3 surfaces are hardcoded (well, all of them have the same Hodge diamond, which is just that of a quartic surface), which is

```
          1
      0        0
  1       20       1
      0        0
          1
```

Now for the cubic fourfold we do

```
print(hypersurface(3, 4))
```

which gives

```
                  1
              0        0
          0       1        0
      0       0        0       0
  0       1       21       1       0
      0       0        0       0
          0       1        0
              0        0
                  1
```

Removing the primitive part of the middle cohomology gives you back the Hodge diamond of a K3 surface, and this is the first glimpse at a very interesting story relating the two.

```
print(hypersurface(3, 4) - K3(1))
```

gives

```
                  1
              0       0
          0       0       0
      0       0       0       0
  0       0       1       0       0
      0       0       0       0
          0       0       0
              0       0
                  1
```
suggesting that there is a "boring" and an "interesting" part of the cohomology of a cubic fourfold.


## Contributing

Please feel free to make suggestions for more examples of Hodge diamonds. Preferably with a link to a closed formula, generating series or method of computation.

Feature requests are also very welcome. And suggestions on improving the documentation are also welcome.
