.. algorithm::

.. warning::

    This algorithm is a prototype for comparing performance with :ref:`algm-AlignAndFocusPowderSlim`
    on VULCAN data. It implements only part of that algorithm, and its properties may change.


.. summary::

.. relatedalgorithms::

.. properties::

Description
-----------

This is a Python version of :ref:`algm-AlignAndFocusPowderSlim`, written to measure how fast the reduction can be
when the events never become Mantid objects. For the options it supports, it produces the same histograms as
:ref:`algm-AlignAndFocusPowderSlim`, bin for bin.

Each bank's detector ids and times-of-flight are read in chunks into reused arrays. A compiled `numba
<https://numba.pydata.org/>`_ kernel applies the calibration and increments the output bin of every event, using all
available cores. Only the final focused histograms are made into a workspace. The instrument definition is never
loaded; the output gets its geometry from :ref:`algm-EditInstrumentGeometry`, using ``L1``, ``L2``, ``Polar`` and
``Azimuthal``.

When the event data are stored uncompressed, the algorithm reads the file directly rather than through HDF5, with
several reads in flight at once. Compressed data are read through HDF5.

The output is in time-of-flight. ``XMin``, ``XDelta`` and ``XMax`` are converted to it using each output spectrum's
focused geometry, as :ref:`algm-AlignAndFocusPowderSlim` does.

Supported:

- a calibration file from :ref:`algm-SaveDiffCal`, giving the calibration, the grouping and the mask
- one set of ``XMin``, ``XDelta`` and ``XMax`` for all spectra, in d-spacing, time-of-flight or momentum transfer
- logarithmic or linear binning
- loading the logs, with an allow or block list

Not supported, compared to :ref:`algm-AlignAndFocusPowderSlim`:

- running without a calibration file
- grouping and calibration workspaces, and separate grouping files
- different binning for each spectrum
- filtering by time, filtering bad pulses, and splitting
- setting the sample name from the file

Usage
-----

.. code-block:: python

    from mantid.simpleapi import AlignAndFocusPowderSlimNumba

    ws = AlignAndFocusPowderSlimNumba(
        "VULCAN_218062.nxs.h5",
        CalFileName="VULCAN_calibration.h5",
        XMin=0.3,
        XDelta=0.0016,
        XMax=3.0,
        BinningUnits="dSpacing",
        BinningMode="Logarithmic",
        L1=43.755,
        L2=[2.296, 2.296, 2.070, 2.070, 2.070, 2.530],
        Polar=[90.0, 90.0, 120.0, 150.0, 157.0, 65.5],
    )

.. categories::

.. sourcelink::
