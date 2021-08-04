[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3893509.svg)](https://doi.org/10.5281/zenodo.3893509)

# Hodge diamond cutter

A collection of Python classes and functions in Sage to deal with Hodge diamonds (and Hochschild homology) of smooth projective varieties, together with many constructions.

If you have used this code in any way (including the interactive versions on my blog), please consider citing it as explained on [Zenodo](https://doi.org/10.5281/zenodo.3893509). You can choose to cite a specific version, or the library in general.


## Getting started

It suffices to load ``diamond.py`` in Sage to get started. The documentation with lots of examples can be [read online](https://pbelmans.ncag.info/hodge-diamond-cutter/) or as [a pdf](https://pbelmans.ncag.info/hodge-diamond-cutter/hodgediamondcutter.pdf).


## Contributing

Please feel free to make suggestions for more examples of Hodge diamonds. Preferably with a link to a closed formula, generating series or method of computation.

Feature requests are also very welcome.

## Instructions to myself

To build the documentation:

```
sage -sh -c "make html"
cp -r _build/html/ docs
sage -sh -c "make latexpdf"
cp _build/latex/hodgediamondcutter.pdf docs
```

To perform the unit tests:

```
sage -t diamond.py
```

And suggestions on improving the documentation are also welcome.

