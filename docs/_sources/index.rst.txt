Documentation for the Hodge diamond cutter
==========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: diamond
   :no-members:

Hodge diamonds
==============

.. autoclass:: diamond.HodgeDiamond
    :members:
    :special-members:
    :member-order: bysource

Hochschild homology
===================

.. autoclass:: diamond.HochschildHomology
    :members:
    :member-order: bysource

Constructions
=============

Elementary constructions
------------------------

.. autofunction:: diamond.zero
.. autofunction:: diamond.point
.. autofunction:: diamond.lefschetz
.. autofunction:: diamond.Pn
.. autofunction:: diamond.hypersurface
.. autofunction:: diamond.weighted_hypersurface
.. autofunction:: diamond.cyclic_cover
.. autofunction:: diamond.complete_intersection

Curves and moduli spaces of sheaves on them
-------------------------------------------

.. autofunction:: diamond.curve
.. autofunction:: diamond.symmetric_power
.. autofunction:: diamond.jacobian
.. autofunction:: diamond.moduli_vector_bundles
.. autofunction:: diamond.seshadris_desingularisation
.. autofunction:: diamond.moduli_parabolic_vector_bundles_rank_two
.. autofunction:: diamond.quot_scheme_curve


Surfaces and moduli spaces of sheaves on them
---------------------------------------------

.. autofunction:: diamond.surface
.. autofunction:: diamond.ruled
.. autofunction:: diamond.K3
.. autofunction:: diamond.enriques
.. autofunction:: diamond.hilbn
.. autofunction:: diamond.nestedhilbn


Abelian varieties and related objects
-------------------------------------

.. autofunction:: diamond.abelian
.. autofunction:: diamond.kummer_resolution


Fano varieties
--------------
These are Hodge diamonds of Fano varieties in the sense that their anticanonical bundle is ample.

The term "Fano variety" can also mean a variety parametrising linear subspaces on another variety. Some of these are Fano in the first sense, others are not (always). See e.g. :func:`diamond.fano_variety_lines_cubic`.

.. autofunction:: diamond.fano_threefold
.. autofunction:: diamond.gushel_mukai
.. autofunction:: diamond.fano_variety_intersection_quadrics_even
.. autofunction:: diamond.fano_variety_intersection_quadrics_odd


Homogeneous varieties and closely related constructions
-------------------------------------------------------
These are also all Fano varieties, but they are grouped together because of their similar origin.

.. autofunction:: diamond.partial_flag_variety
.. autofunction:: diamond.generalised_grassmannian
.. autofunction:: diamond.grassmannian
.. autofunction:: diamond.orthogonal_grassmannian
.. autofunction:: diamond.symplectic_grassmannian
.. autofunction:: diamond.lagrangian_grassmannian
.. autofunction:: diamond.horospherical
.. autofunction:: diamond.odd_symplectic_grassmannian


Moduli spaces attached to quivers
---------------------------------

.. autofunction:: diamond.quiver_moduli


Hyperkähler varieties
---------------------

.. autofunction:: diamond.K3n
.. autofunction:: diamond.generalised_kummer
.. autofunction:: diamond.ogrady6
.. autofunction:: diamond.ogrady10


Other
-----

.. autofunction:: diamond.Mzeronbar
.. autofunction:: diamond.fano_variety_lines_cubic
