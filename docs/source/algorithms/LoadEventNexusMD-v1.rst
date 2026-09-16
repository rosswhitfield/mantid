.. algorithm::

.. summary::

.. relatedalgorithms::

.. properties::

Description
-----------

Loads single crystal diffraction event data from a NeXus file straight into an
:py:obj:`MDEventWorkspace <mantid.api.IMDWorkspace>` in momentum space. It is the equivalent of
:ref:`algm-LoadEventNexus` followed by :ref:`algm-ConvertToMD`, but the intermediate
:ref:`EventWorkspace <EventWorkspace>` is never built.

That intermediate workspace is what makes the usual route expensive: it holds one event list per pixel and a
time-of-flight plus a pulse time for every neutron, all of which is discarded one step later once each event has
become a *Q* vector. This algorithm reads ``event_id`` and ``event_time_offset`` out of each bank in chunks, converts
them, and adds them to the box tree directly, so the event lists are never built. It follows the same approach as
:ref:`algm-AlignAndFocusPowderSlim` does for powder data.

On a 128 million event TOPAZ run this takes about 13 seconds and 7.7 GiB against about 18 seconds and 10.5 GiB for
:ref:`algm-LoadEventNexus` followed by :ref:`algm-ConvertToMD`, so roughly a third faster in about three quarters of
the memory. Asking for ``MDLeanEvent`` output saves a further third or so of the output workspace, at the cost of the
per-event detector id.

Most of the remaining time goes on filling the MD box tree rather than on reading or converting, so ``SplitInto`` is
the property worth tuning if the algorithm is slower than expected: every event walks the tree from the root, and
the shallower tree that a larger ``SplitInto`` gives is markedly quicker to fill. Splitting into 2 rather than the
default 5 makes that TOPAZ run about 70% slower.

For the same ``MinValues`` and ``MaxValues`` the output is identical, event for event, to
:ref:`algm-LoadEventNexus` followed by :ref:`algm-ConvertToMD` with ``QDimensions=Q3D`` and
``dEAnalysisMode=Elastic``.

Momentum transfer
#################

Scattering is assumed to be elastic. For a detector seen from the sample along the unit vector
:math:`\hat{e}`, with the beam along :math:`\hat{z}`, an event of time-of-flight :math:`t` is placed at

.. math:: \vec{Q} = |k| \left( \hat{z} - \hat{e} \right), \qquad |k| = \frac{2\pi}{\lambda} = \frac{2 \pi m_n (L_1 + L_2)}{h\,t}

The sign follows the ``Q.convention`` setting: the default ``Inelastic`` gives :math:`k_i - k_f` as above, while
``Crystallography`` gives :math:`k_f - k_i`.

The instrument's ``T0`` parameter is added to every time-of-flight before the conversion, exactly as
:ref:`algm-LoadEventNexus` does, and is recorded as the ``T0`` sample log on the output. This matters on the single
crystal instruments, which all define it.

Choosing the frame
##################

``Q (lab frame)`` gives the momentum transfer in the laboratory frame. ``Q (sample frame)`` takes the goniometer
rotation out by applying the inverse of the goniometer matrix, which is read from the ``omega``, ``chi`` and ``phi``
sample logs while the file is loaded. If the run has no goniometer the two frames are the same and a warning is
logged.

Extents
#######

``MinValues`` and ``MaxValues`` define the box the output workspace covers; events outside it are not added. If they
are left empty the enclosing box is calculated from the detector positions and the range of time-of-flight in the
file, which costs one extra pass over ``event_time_offset``. Supplying them explicitly skips that pass, and is worth
doing for data whose shortest times of flight would otherwise stretch the box over mostly empty space.

Limitations
###########

This first version deliberately does no masking, no filtering, no calibration and no corrections. In particular
there is no time-of-flight or pulse-time filtering, no bank selection, no detector masking and no Lorentz
correction. Events whose ``event_id`` matches no detector in the instrument definition are discarded, and their
number is reported, which is also what :ref:`algm-LoadEventNexus` does with them.

Usage
-----

**Example - loading a CORELLI run into Q_lab**

.. testcode:: ExLoadEventNexusMD

   md = LoadEventNexusMD(Filename='CORELLI_2100.nxs.h5',
                         QFrame='Q (lab frame)',
                         MinValues=[-10, -10, -10],
                         MaxValues=[10, 10, 10])

   print('Workspace type is: {}'.format(md.id()))
   print('It has {} events in {} dimensions'.format(md.getNEvents(), md.getNumDims()))
   print('Dimensions: {}'.format(', '.join(md.getDimension(d).name for d in range(3))))

**Output:**

.. testoutput:: ExLoadEventNexusMD

   Workspace type is: MDEventWorkspace<MDEvent,3>
   It has 2502 events in 3 dimensions
   Dimensions: Q_lab_x, Q_lab_y, Q_lab_z

.. categories::

.. sourcelink::
