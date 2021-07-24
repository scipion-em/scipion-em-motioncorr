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


+--------------+----------------+--------------------+
| prod: |prod| | devel: |devel| | support: |support| |
+--------------+----------------+--------------------+

.. |prod| image:: http://scipion-test.cnb.csic.es:9980/badges/motioncorr_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/motioncorr_devel.svg
.. |support| image:: http://scipion-test.cnb.csic.es:9980/badges/motioncorr_support.svg

**IMPORTANT!**

    If you have imported movies with a gain file in **DM4** format, you need to **flip the gain reference upside-down** in the motioncor2 protocol! (`bug details <https://github.com/I2PC/xmippCore/issues/39>`_)

Installation
------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-motioncorr

b) Developer's version

   * download repository 
   
   .. code-block::
   
      git clone https://github.com/scipion-em/scipion-em-motioncorr.git

   * install
   
   .. code-block::

      scipion installp -p path_to_scipion-em-motioncorr --devel

Motioncor2 binaries will be installed automatically with the plugin, but you can also link an existing installation. 
Default installation path assumed is ``software/em/motioncor2-1.4.2``, if you want to change it, set *MOTIONCOR2_HOME* in ``scipion.conf`` file to
the folder where the Motioncor2 is installed. Depending on your CUDA version you might want to change the default binary from ``MotionCor2_1.4.2_Cuda101-02-15-2020``
to a different one by explicitly setting *MOTIONCOR2_BIN* variable. If you need to use CUDA different from the one used during Scipion installation
(defined by CUDA_LIB), you can add *MOTIONCOR2_CUDA_LIB* variable to the config file. Various binaries can be downloaded from the official UCSF website.

To check the installation, simply run the following Scipion test: 

``scipion test motioncorr.tests.test_protocols_motioncor2.TestMotioncor2AlignMovies``

Supported versions
------------------

1.2.6, 1.3.0, 1.3.1, 1.3.2, 1.4.0, 1.4.2

Protocols
---------

* movie alignment
* align tilt-series movies

References
----------

1.  Shawn Q Zheng, Eugene Palovcak, Jean-Paul Armache, Kliment A Verba, Yifan Cheng & David A Agard. MotionCor2: anisotropic correction of beam-induced motion for improved cryo-electron microscopy. Nature Methods volume 14, pages 331–332 (2017).
