=================
Motioncorr plugin
=================

This plugin allows to use motioncor2 program within the Scipion framework.

Motioncor2 is a GPU-accelerated program for correction of electron beam-induced sample motion. It is developed by `Shawn Zheng <https://emcore.ucsf.edu/ucsf-motioncor2>`_.

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

    1. If you have imported movies with a gain file in **DM4** format, you need to **flip the gain reference upside-down** in the motioncor2 protocol! (`see details <https://github.com/I2PC/xmippCore/issues/39>`_)
    2. When importing EER movies, you should specify dose per single EER frame during import step.
    3. If you are processing EER movies and providing \*.gain reference file camera defects will be automatically extracted from the gain file header and converted to Motioncor2 format. This step is omitted if you provide a defects file yourself.

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

- Motioncor2 binaries will be installed automatically with the plugin, but you can also link an existing installation.
- Default installation path assumed is ``software/em/motioncor2-1.5.0``, if you want to change it, set *MOTIONCOR2_HOME* in ``scipion.conf`` file to the folder where the Motioncor2 is installed.
- Depending on your CUDA version this plugin will guess the right default binary from ``MotionCor2_1.5.0_CudaXY_05-31-2022`` (X is for cuda major version, Y for the minor). You can always set a different one by explicitly setting *MOTIONCOR2_BIN* variable.
- If you need to use CUDA different from the one used during Scipion installation (defined by CUDA_LIB), you can add *MOTIONCOR2_CUDA_LIB* variable to the config file. Various binaries can be downloaded from the official UCSF website.

For an automatically updated installation of motioncor2 binaries, do not define neither *MOTIONCOR2_HOME* nor *MOTIONCOR2_BIN*.

To check the installation, simply run the following Scipion test:

``scipion test motioncorr.tests.test_protocols_motioncor2.TestMotioncor2AlignMovies``

Licensing
---------

Motioncor2 is free for academic use only. For commercial use, please contact David Agard or Yifan Cheng for licensing:

    * agard@msg.ucsf.edu
    * Yifan.Cheng@ucsf.edu

Supported versions
------------------

1.4.0, 1.4.2, 1.4.4, 1.4.5, 1.4.7, 1.5.0

Protocols
---------

* movie alignment
* align tilt-series movies

References
----------

1.  Shawn Q Zheng, Eugene Palovcak, Jean-Paul Armache, Kliment A Verba, Yifan Cheng & David A Agard. MotionCor2: anisotropic correction of beam-induced motion for improved cryo-electron microscopy. Nature Methods volume 14, pages 331–332 (2017).
