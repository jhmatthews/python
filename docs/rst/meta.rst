==================
Meta-documentation
==================
How to use the auto-documentation tools
---------------------------------------

*python* has two documentation tools, stored in the *py_progs* directory, and a common file
format for documenting the python functions.

YAML documentation
------------------

The documentation for *python* variables takes the form of *YAML* input files stored in
*docs/parameters*. These are functionally more or less just *Python* dictionaries, and have a very
flexible structure. A full example outline is:

.. literalinclude:: reference_yaml.txt
   :language: yaml
   :lines: 1-12


The **required** keys that must be present in any file are ``name``, ``description``, ``type``
and ``file``.

For enumerators, a more complex setup can be used for the ``values`` key outlining all the
possible choices:

.. literalinclude:: reference_yaml.txt
   :language: yaml
   :lines: 14-25

The ``parent`` key is used to link between parameters, so you can see which options depend on others
(for example, :ref:`Reverb.matom_lines` depends on :ref:`Line_transfer` being a macro-atom mode).
Parent should be a list of parameters, followed by one (or a list of) values of the parent that result
in this parameter being used, for example:

.. literalinclude:: reference_yaml.txt
   :language: yaml
   :lines: 27-28

Parent will automatically link to the page for each parent. Parent is also used when figuring out the structure
within each file, so put any local parents at the top of the list. 
If you'd like to reference one parameter elsewhere within another, you can use the following format

.. literalinclude:: reference_yaml.txt
   :language: yaml
   :lines: 30-30


TODO: Make structure of parent explicit.

autogenerate_parameter_docs.py
------------------------------

This tool takes the *python* .c input files, and scans them for calls to input parameters from file.
It is run from the command line as::

  autogenerate_parameter_docs.py

This mode will tell you which parameters are **new**, and undocumented, and which files are
**deprecated** and refer to parameters that no longer exist. In many cases, deprecated parameters
will simply be renamed. You can use this as a reminder to rename the parameter and tweak the *YAML*
documentation file.

Once this has been done, it can be run using::

  autogenerate_parameter_docs.py -w

This will move any deprecated documentation to *docs/parameters/old*, and generate new *YAML*
documentation files for the new parameters in *docs/parameters*. You can now proceed to fill
in the documentation skeletons produced by the utility.


----------------------------------------

autogenerate_rtd_pages.py
-------------------------

This tool takes the `YAML documentation`_ input files, and converts them into *RST* files and from there a *HTML*
documentation page using *Sphinx*. You don't need to know anything about the formats in order to use
it- you simply need *Sphinx* installed, and it is available as a MacPorts package or through apt
on 'nix systems.

It will read the files in the *docs/parameters* directory and build a *HTML* page,
then automatically open it in your browser. You may see some errors during the creation::

   while scanning a block scalar
    in "/Users/amsys/python/docs/parameters/sv.diskmin.yaml", line 5, column 7
   expected chomping or indentation indicators, but found 'c'
    in "/Users/amsys/python/docs/parameters/sv.diskmin.yaml", line 5, column 8

This occurs during the *YAML* step when a file has an invalid entry, and that parameter file will be
skipped as a consequence. Typically it is cause when a text line is malformed. For inline text
(e.g. ``values: condition``) it occurs when the string starts with a non-alphanumeric character like
``>`` (so ```values: >0`` is invalid).

You may also see it when using block text (e.g. ``description:|``) if the following lines are not
correctly indented e.g.::

   description:|
   The text for values has not been indented correctly
   and so will throw an error.

   description:|
     This text has been indented
     and so will not throw an error

A common warning is::

   WARNING: Unknown target name:

This occurs during the *Sphinx* step when it cannot link to the parent for a variable. If a
parameter has a parent provided as a ``key:value`` pair e.g. ``reverb.type: Anything above 3``.

*Sphinx* will automatically try to link to ``reverb.type`` as its parent. If the name is misspelled
(e.g. ``rverb.type``) or the parameter simply doesn't exist any more, it will throw a warning and
generate a dead link.
