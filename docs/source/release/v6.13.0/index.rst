.. _v6.13.0:

===========================
Mantid 6.13.0 Release Notes
===========================

.. figure:: ../../images/6_13_release/mantid-6_13-summary.png
   :class: screenshot
   :width: 600px
   :align: right

.. contents:: Table of Contents
   :local:

We are proud to announce version 6.13.0 of Mantid.

This release showcases a wide range of new algorithms and behind-the-scenes improvements, including:

- Availability of a native Apple Silicon package for macOS users. With Rosetta no longer being required, users of this
  new package should notice a significant performance improvement on Apple Silicon based machines. `Installing via
  conda <https://www.mantidproject.org/installation/index#conda>`_ will automatically install the best version of the
  package for your machine.
- A :ref:`complete overhaul of the NeXus API <6_13_dependencies>` for loading HDF5-based files.
- A :ref:`set of brand new algorithms <6_13_sans>` for producing and saving multiple polarized SANS files into a single
  output file.
- Significant improvements to the :ref:`PEARLTransFit <6_13_powder_diffraction>` algorithm.
- An :ref:`upgrade to Python 3.11 <6_13_dependencies>`.

These are just some of the many improvements in this release, so please take a
look at the release notes, which are filled with details of the
important changes and improvements in many areas. The development team
has put a great effort into making all of these improvements within
Mantid, and we would like to thank all of our beta testers for their
time and effort helping us to make this another reliable version of Mantid.

Throughout the Mantid project we put a lot of effort into ensuring
Mantid is a robust and reliable product. Thank you to everyone that has
reported any issues to us. Please keep on reporting any problems you
have, or crashes that occur on our `forum`_.

Installation packages can be found on our `download page`_
which now links to the assets on our `GitHub release page`_, where you can also
access the source code for the release.

Citation
--------

Please cite any usage of Mantid as follows:

- *Mantid 6.13.0: Manipulation and Analysis Toolkit for Instrument Data.; Mantid Project*. `doi: 10.5286/SOFTWARE/MANTID6.13 <https://dx.doi.org/10.5286/SOFTWARE/MANTID6.13>`_

- Arnold, O. et al. *Mantid-Data Analysis and Visualization Package for Neutron Scattering and mu-SR Experiments.* Nuclear Instruments
  and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment 764 (2014): 156-166
  `doi: 10.1016/j.nima.2014.07.029 <https://doi.org/10.1016/j.nima.2014.07.029>`_
  (`download bibtex <https://raw.githubusercontent.com/mantidproject/mantid/master/docs/source/mantid.bib>`_)


Changes
-------

.. toctree::
   :hidden:
   :glob:

   *

- :doc:`Framework <framework>`
- :doc:`Mantid Workbench <mantidworkbench>`
- :doc:`Diffraction <diffraction>`
- :doc:`Muon Analysis <muon>`
- Low Q

  - :doc:`Reflectometry <reflectometry>`

  - :doc:`SANS <sans>`
- Spectroscopy

  - :doc:`Direct Geometry <direct_geometry>`

  - :doc:`Indirect Geometry <indirect_geometry>`

  - :doc:`Inelastic <inelastic>`

Full Change Listings
--------------------

For a full list of all issues addressed during this release please see the `GitHub milestone`_.

.. _download page: https://download.mantidproject.org

.. _forum: https://forum.mantidproject.org

.. _GitHub milestone: https://github.com/mantidproject/mantid/pulls?utf8=%E2%9C%93&q=is%3Apr+milestone%3A%22Release+6.13%22+is%3Amerged

.. _GitHub release page: https://github.com/mantidproject/mantid/releases/tag/v6.13.0
