[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3893509.svg)](https://doi.org/10.5281/zenodo.3893509)
[![releases](https://badgen.net/github/release/pbelmans/hodge-diamond-cutter?color=green)](https://github.com/pbelmans/hodge-diamond-cutter/releases)
[![license](https://badgen.net/github/license/pbelmans/hodge-diamond-cutter)](https://github.com/pbelmans/hodge-diamond-cutter/blob/master/LICENSE)

[![tests](https://github.com/pbelmans/hodge-diamond-cutter/actions)](https://github.com/pbelmans/hodge-diamond-cutter/actions/workflows/main.yml/badge.svg)

# Hodge diamond cutter

A collection of Python classes and functions in Sage to deal with Hodge diamonds (and Hochschild homology) of smooth projective varieties, together with many constructions.




## Getting started

It suffices to put ``diamond/diamond.py`` in your directory and load it using ``load("diamond.py")`` in Sage to get started.

Alternatively you can install it as follows:

``sage --pip install git+https://github.com/pbelmans/hodge-diamond-cutter.git``

and then you can use

``from diamond import *``

to use it.

The documentation with lots of examples can be [read online](https://pbelmans.ncag.info/hodge-diamond-cutter/) or as [a pdf](https://pbelmans.ncag.info/hodge-diamond-cutter/hodgediamondcutter.pdf).


## How to cite

If you have used this code in any way (including the interactive versions on my blog), please consider citing it as explained on [Zenodo](https://doi.org/10.5281/zenodo.3893509). You can choose to cite a specific version, or always the latest version. For the latter you can use `doi:10.5281/zenodo.3893509`.

The following BibTeX entry is a good starting point:

```bibtex
@software{hodge-diamond-cutter,
  author = {Belmans, Pieter},
  title = {Hodge diamond cutter},
  url = {https://github.com/pbelmans/hodge-diamond-cutter},
  doi = {10.5281/zenodo.3893509},
}
```

which leads to something like

> Pieter Belmans. _Hodge diamond cutter_. doi:10.5281/zenodo.3893509. url: ht<span>tps://github.com/pbelmans/hodge-diamond-cutter.


## Contributing

Please feel free to make suggestions for more examples of Hodge diamonds. Preferably with a link to a closed formula, generating series or method of computation.

Feature requests are also very welcome.

## A warning

In the [words of Simon Pepin Lehalleur on Twitter](https://twitter.com/plain_simon/status/1355599647893549056),
in order to use the Hodge diamond cutter, you must always answer the question

> You are solving a PDE; do you want to proceed?

positively.

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

