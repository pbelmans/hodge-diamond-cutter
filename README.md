# Hodge diamond cutter

A collection of Python classes and functions to deal with Hodge diamonds (and Hochschild homology) of smooth projective varieties, together with many constructions.

## Getting started

It suffices to load ``diamond.sage`` in Sage to get started. It was written using Sage 8.3, but presumably it works without modification in more recent and not too old versions, as nothing fancy is used.

## Examples

1. As a first example, let us consider the cubic fourfold. There is an intricate connection to K3 surfaces. To see the Hodge diamond of a K3 surface, it suffices to do

   ```
   print K3
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
   print complete_intersection(3, 4)
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

   Removing the primitive part of the cohomology gives you back the Hodge diamond of a K3 surface, and this is the first glimpse at a very interesting story relating the two.

2. As a second example, let us discuss a conjectural semiorthogonal decomposition for the moduli space of rank 2 bundles with fixed determinant of degree 1 on a curve $C$ of genus $g$. As observed by [Kyoung-Seog Lee](https://arxiv.org/abs/1806.11101) and [myself](https://www.mfo.de/document/1822/preliminary_OWR_2018_24.pdf) the Hodge diamond of this variety can be decomposed in terms of symmetric powers of the curve $C$, which can be checked as follows

   ```
   for g in range(2, 10):
     assert moduli_vector_bundles(2, 1, g) == sum([symmetric_power(i, g)(i) for i in range(g)]) + sum([symmetric_power(i, g)(3*g - 3 - 2*i) for i in range(g - 1)])
   ```

## Contributing

Please feel free to make suggestions for more examples of Hodge diamonds. Preferably with a link to a closed formula, generating series or method of computation.

Feature requests are also very welcome. And suggestions on improving the documentation are also welcome.
