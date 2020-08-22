Developer's Guide and Roadmap
=============================

The `chemicals` project has the following core ideas as its goals:

* Efficient
    * Functions do only the work required.
    * Caching various values, precomputing others.
    * Using various macros and automated expressions to run code with Numba at its optimal speed.
    * Not using Numpy/SciPy most of the time, allowing PyPy or Numba to speed code up when speed is important.
* Capable of vectorized computation
    * Wrapped with numpy's np.vectorize.
    * Wrapped with numba's ufunc machinery.
* Comprehensive
    * Most correlations taught at the undergrad level included.
* Capable of handling units
    * Pint interface
    * All docstrings/code in base SI units

This is a hobby project primarily by Caleb Bell. Contributors are welcome! Chemicals properties is huge big field.

Scope and Future Features
-------------------------
The library has a focus on pure-component properties. Lot sof data files are included and there is ample room for more files. However, it is not feasible to add data on an individual chemical alone - a compilation of data which includes that chemicals must be found instead.




The following properties have been looked at a little but may have too much data fit on PyPi. If you are interested in implementing one of them it may take multiple data packages uploaded to PyPi alongside `chemicals`:

* IR, NMR, MS, and UV-Vis spectra. Files are typically 4-40 KB. A collection of ~2000 UV files from NIST takes 6 MB space, so the space issue would not be ran into right away. Some databases like NIST and `NMRShiftDB <https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdb/media-type/html/user/anon/page/default.psml/js_pane/P-Home>`_ are open.


Contributing
------------
`Chemicals` has a lot of infrastructure that makes it attractive to add code to the project. Adding functionality to `chemicals` may also make your work more visible to more people.

Adding new functionality is possible without compromising load speed, RAM usage or maintainability. If you have a chemical property you are interested in adding, a utility function, or a new data source, please feel free to open a PR and we can make any changes needed. There is no template - just do your best.

In an ideal world, new contributions would come with unit tests, docstrings, an addition to the tutorial if relevant.

Running Tests
-------------
From the root directory of the project you downloaded with `git clone https://github.com/CalebBell/chemicals.git`, run the following command:

python3 -m pytest .

This will run all of the tests. Additionally pytest can be used to run the doctests:

python3 -m pytest --doctest-modules .

The test suite can take some time to run; tests are marked with various markers to allow a fast subset of tests to run.

python3 -m pytest -m "not slow" .

This should only take a few seconds, and show red output if a test is broken. To keep the main test suite fast, pytest allows a flag which shows how long each test takes.

python3 -m pytest -m "not slow" --durations=100

If a test you added appears in this list, consider splitting it into a fast portion and a slow portion.

Docstrings
----------
The docstrings follow Pep8, most of the numpydoc standard,
More information about numpydoc can be found `here <https://numpydoc.readthedocs.io/en/latest/format.html>`_

In addition to being documentation, the docstrings in `chemicals` serve the following purposes:

* Contain LaTeX math formulas for implemented formulas. This makes it easy for the reader and authors to follow code. This is especially important when the code can be optimized by hand significantly, and end up not looking like the math formulas.
* Contain doctests for every public method. These examples often make debugging really easy since they can just be copy-pasted into Jupyter or an IDE/debugger.
* Contain type information for each variable, which is automatically parsed by the unit handling framework around `pint`.
* Contain the units of each argument, which is used by the unit handling framework around `pint`.
* Contain docstrings for every argument - these are checked by the unit tests programatically to avoid forgetting to add a description, which the author did often before the checker was added.

No automated style tool is ran on the docstrings at present, but the following command
was used once to format the docstrings with the tool `docformatter <https://github.com/myint/docformatter>`_

python3 -m docformatter --wrap-summaries=80 --wrap-descriptions=80 --in-place --recursive .

This does not quite match numpydoc's recommended 75 character limit.

Doctest
-------
As anyone who has used doctest before knows, floating-point calculations have trivially different results across platforms. An example cause of this is that most compilers have different sin/cos implementations which are not identical. However, docstrings are checked bit-for-bit, so consistent output is important. Python is better than most languages at maintaining the same results between versions but it is still an issue.

The docstrings are not considered sufficiently consistent to be part of the automated CI infrastructure. All functionality tested by docstrings should also be tested as part of the unit tests.

CPython 3.7 64 bit on Linux compiled with gcc 9.2 is currently the platform used to generate the final/official results of docstring examples. Docstrings are should be added by contributors for new functionality, but **don't worry about this floating point issue**. The principal author will make any necessary changes before each release.

Type Hints
----------
The Python ecosystem is gradually adding support for type information, which may allow static analyzers to help find bugs in code even before it is ran. The author has not found these helpful in Python yet - the tools are too slow, missing features, and most libraries do not contain type information. However, type hints might still be useful for your program that uses `chemicals`!

For that reason `chemicals` includes a set of type hints as stub files (.pyi extension). These are not generated by hand - they use the cool `MonkeyType <https://github.com/Instagram/MonkeyType/>`_ library.
An included script `make_test_stubs` interfaces with this library, which runs the test suite and at the end generates the type hints including the types of every argument to every function seen in the test suite. This is another reason comprehensive test suite coverage is required.

Monkeytype on the `chemicals` test suite takes ~5 minutes to run, and generates a ~1 GB database file which is deleted at the end of the run. Some manipulation of the result by hand may be required in the future, or MonkeyType may be replaced by making the type hints by hand. It is planned to incorporate the type stubs into the main file at some point in the future when the tooling is better.

**If you are contributing, the main developer can do this step for your contribution.**

Supported Python Versions
-------------------------
It is strongly recommended to use Chemicals with Python 3.6 or later.

Numpy's latest release supports Python 3.6 or later as of August 2020.
Pint, the unit interface, supports Python 3.6 or later as of August 2020.
If using the Numba interface of Chemicals, the latest version of Numba is required. New features added to Numba may be added to Chemicals quite quickly.

Chemicals may import in an earlier Python but that is not an indication you should use it in that case.

Packaging
---------
The most up to date chemicals can be obtained on GitHub, and new releases are pushed to PyPi whenever a new release is made.


Code Formatting
---------------
Pep8 is loosely followed. Do your best to follow it if possible, otherwise don't worry about it. Please don't submit a PR for just style changes.

Documentation
-------------
Sphinx is used with readthedocs. Readthedocs is configured to build whatever is on the release branch. From the root of the `chemicals` project, the documentation can be built with the following command, which will output html files into a "_build" folder:

sphinx-build -b html docs _build

Sample Notebooks
----------------
The `nbval <https://pypi.org/project/nbval/>`_ pytest plugin can be used to check the results of running the notebooks against the stored notebooks.

On UNIX/Mac OS/WSL, the notebook results can be regenerated with the following shell command, from the directory with the notebooks:

for i in *.ipynb ; do python3 -m jupyter nbconvert --to notebook --inplace --execute "$i" ; done

Continuous Integration
----------------------
Travis and Appveyor are presently used. They test only code in the `release` branch. Some tests, like those that download data from the internet, are not ran by design on their platforms. The same goes for testing `numba` online - getting an up to date version of numba is challenging.

Load Speed
----------
On CPython, `chemicals` will load Numpy on load if it is available and `SciPy` and `Pandas` when they are needed. Numpy loads in ~150 ms. chemicals alone loads in ~15 ms. It is intended for this to increase only slowly.

RAM Usage
---------
Adding new data and methods is well worth the price of increasing RAM, but it is intended to keep RAM consumption small via lazy-loading all data sets.

Additional Material
-------------------
More information about the interfaces with PyPy and Numba can be found `here <https://fluids.readthedocs.io/developers.html>`_.