PEN-ASSEMBLE v0.5.1
===================

**Programmable Enzyme Networks — Automated Strategy and Scoring Engine for Molecular Bridge-recombinase Library Engineering**

Part of `PEN-STACK <https://github.com/ahmedanees-m>`_ (Paper 4 of 5):
`genome-atlas <https://github.com/ahmedanees-m/genome-atlas>`_ →
`mech-class <https://github.com/ahmedanees-m/mech-class>`_ →
`pen-score <https://github.com/ahmedanees-m/pen-score>`_ →
**pen-assemble** → PEN-COMPARE *(in prep)*

----

PEN-ASSEMBLE is the computational nomination pipeline for IS110-family bridge recombinase
design. It generates and evaluates **1,029 candidate designs** across four orthogonal engineering
strategies — domain-swap chimeras, IS110 ortholog discovery, Monte Carlo deimmunization, and
ProteinMPNN backbone redesign — using the **8-axis PenScore** composite metric (pen-score v0.1.2).
The pipeline is pre-registered, fully reproducible, and produces wet-lab synthesis sheets
ready for experimental validation.

.. note::
   All five pre-registered predictions **PASS** (5/5). Publication policy: *"PUBLISH with
   strong claim"*. 16 designs beat IS621 verbatim lockpoint (PenScore > 0.929). See :doc:`predictions`.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   quickstart
   api
   predictions
   deviations
   changelog

Quick links
-----------

* :ref:`genindex`
* :ref:`modindex`
* :doc:`api`
* `GitHub <https://github.com/ahmedanees-m/pen-assemble>`_
* `genome-atlas <https://github.com/ahmedanees-m/genome-atlas>`_
* `mech-class <https://github.com/ahmedanees-m/mech-class>`_
* `pen-score <https://github.com/ahmedanees-m/pen-score>`_
