Quick Start
===========

Installation
------------

Install from the repository root::

    pip install -e ".[dev,docs]"

Or without extras::

    pip install -e .

Requirements: Python ≥ 3.10, pandas ≥ 2.0, pyarrow ≥ 14.0, numpy ≥ 1.24.

Generating the catalog
----------------------

Run the Part F scripts in order::

    cd scripts/
    py 50_assemble_catalog.py   # Step 23 — catalog CSV/Parquet/FASTA/JSON
    py 51_build_browser.py      # Step 24 — self-contained HTML browser
    py 52_generate_wetlab_reference.py  # Step 25 — wet-lab Markdown files

Outputs appear in ``catalog/release_v0.5.0/``.

Loading the catalog in Python
------------------------------

.. code-block:: python

    from pen_assemble.catalog import load_catalog, load_p1_beaters, load_top5

    # Full 1,029-design scorecard
    df = load_catalog()
    print(df.shape)  # (1029, 28)

    # 16 P1-beating designs (pen_score > 0.929)
    p1 = load_p1_beaters()
    print(p1[["design_id", "strategy", "pen_score"]].head())

    # P5-compliant top-5
    t5 = load_top5()
    print(t5[["design_id", "strategy", "pen_score"]])

Computing PenScore
------------------

.. code-block:: python

    from pen_assemble.pen_score import pen_score, PenScoreAxes, beats_is621

    ax = PenScoreAxes(
        S_DSB=1.0, S_Spec=1.0, S_Cargo=1.0, S_Deliv=1.0,
        S_Immuno=0.8777, S_Prog=1.0, S_Mature=1.0,
    )
    s = pen_score(ax)
    print(f"PenScore = {s:.4f}")       # 0.9678
    print(f"Beats IS621: {beats_is621(s)}")  # True

Codon optimisation
------------------

.. code-block:: python

    from pen_assemble.codon import build_expression_orf, check_restriction_sites

    aa_seq = "MDRFFPVIRICKVGFTME..."  # your design sequence
    orf = build_expression_orf(aa_seq, kozak=True, stop=True)

    # Check for restriction sites before ordering synthesis
    hits = check_restriction_sites(orf)
    if hits:
        print(f"Warning: {hits} sites present — silent-mutate before cloning")

Running the pre-registered prediction tests
-------------------------------------------

::

    cd scripts/
    py 40_test_pred_P1_beat_is621.py
    py 41_test_pred_P2_cargo_deliv.py
    py 42_test_pred_P3_deimmunization.py
    py 43_test_pred_P4_orthologs.py
    py 44_test_pred_P5_diversity.py
    py 45_summarise_predictions.py

Results are written to ``pipeline_results_local_test/validation/``.

Running the test suite
----------------------

::

    py -m pytest tests/ -v

Or with coverage::

    py -m pytest tests/ -v --cov=pen_assemble --cov-report=term-missing
