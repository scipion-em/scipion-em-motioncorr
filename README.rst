=================
Motioncorr plugin
=================

This plugin allows to use motioncor3 program within the Scipion framework.

Motioncor3 is a GPU-accelerated program for correction of electron beam-induced sample motion. It is developed by `Shawn Zheng <https://github.com/czimaginginstitute/MotionCor3>`_.

.. image:: https://img.shields.io/pypi/v/scipion-em-motioncorr.svg
        :target: https://pypi.python.org/pypi/scipion-em-motioncorr
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-motioncorr.svg
        :target: https://pypi.python.org/pypi/scipion-em-motioncorr
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-motioncorr.svg
        :target: https://pypi.python.org/pypi/scipion-em-motioncorr
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-motioncorr?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-motioncorr
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-motioncorr
        :target: https://pypi.python.org/pypi/scipion-em-motioncorr
        :alt: Downloads

**IMPORTANT!**

    1. If you have imported movies with a gain file in **DM4** format, you need to **flip the gain reference upside-down** in the motioncor protocol! (`see details <https://github.com/I2PC/xmippCore/issues/39>`_)
    2. When importing EER movies, you should specify dose per single EER frame during import step.
    3. When importing EER tilt-series movies, you should specify dose per tilt during import step. It is not relevant for motion-correction (no dose-weighting is done here), but will be used later on.
    4. If you are processing EER movies and providing \*.gain reference file camera defects will be automatically extracted from the gain file header and converted to Motioncor format. This step is omitted if you provide a defects file yourself.

Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-motioncorr

b) Developer's version

   * download repository 
   
   .. code-block::
   
      git clone -b devel https://github.com/scipion-em/scipion-em-motioncorr.git

   * install
   
   .. code-block::

      scipion installp -p /path/to/scipion-em-motioncorr --devel

**Important:** Starting from plugin v3.15, the config variables have been renamed. See: `scipion3 config -p motioncorr`

- Motioncor binaries will be installed automatically with the plugin, but you can also link an existing installation.
- Default installation path assumed is ``software/em/motioncor3-1.1.2``, if you want to change it, set *MOTIONCOR_HOME* in ``scipion.conf`` file to the folder where the Motioncor3 is installed.
- Depending on your CUDA version this plugin will guess the right default binary from ``MotionCor3_1.1.2_CudaXY_06-11-2024`` (X is for cuda major version, Y for the minor). You can always set a different one by explicitly setting *MOTIONCOR_BIN* variable.
- If you need to use CUDA different from the one used during Scipion installation (defined by CUDA_LIB), you can add *MOTIONCOR_CUDA_LIB* variable to the config file. Various binaries can be downloaded from the official UCSF website.

For an automatically updated installation of motioncor binaries, do not define neither *MOTIONCOR_HOME* nor *MOTIONCOR_BIN*.

To check the installation, simply run the following Scipion tests:

.. code-block::

    scipion test motioncorr.tests.test_protocols_motioncor.TestMotioncorAlignMovies
    scipion test motioncorr.tests.test_protocols_tomo.TestMotioncorTiltSeriesAlignMovies

Licensing
---------

Motioncor3 is provided under BSD-3-Clause license

Supported versions
------------------

* motioncor3 1.1.1, 1.1.2
* motioncor2 1.6.4

Protocols
---------

* movie alignment
* align tilt-series movies

References
----------

1.  Shawn Q Zheng, Eugene Palovcak, Jean-Paul Armache, Kliment A Verba, Yifan Cheng & David A Agard. MotionCor2: anisotropic correction of beam-induced motion for improved cryo-electron microscopy. Nature Methods volume 14, pages 331–332 (2017).
