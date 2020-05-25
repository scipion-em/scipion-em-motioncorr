=================
Motioncorr plugin
=================

This plugin allows to use motioncor2 program within the Scipion framework.

Motioncor2 is a GPU-accelerated program for correction of electron beam-induced sample motion. It is developed by `Shawn Zheng <https://emcore.ucsf.edu/ucsf-motioncor2>`_.

+--------------+----------------+--------------------+
| prod: |prod| | devel: |devel| | support: |support| |
+--------------+----------------+--------------------+

.. |prod| image:: http://scipion-test.cnb.csic.es:9980/badges/motioncorr_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/motioncorr_devel.svg
.. |support| image:: http://scipion-test.cnb.csic.es:9980/badges/motioncorr_support.svg


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
Default installation path assumed is ``software/em/motioncor2-1.3.1``, if you want to change it, set *MOTIONCOR2_HOME* in ``scipion.conf`` file to the folder where the Motioncor2 is installed. Depending on your CUDA version you might want to change the default binary from ``MotionCor2_v1.3.1-Cuda92`` to a different one by explicitly setting *MOTIONCOR2* variable. If you need to use CUDA different from the one used during Scipion installation (defined by CUDA_LIB), you can add *MOTIONCOR2_CUDA_LIB* variable to the config file. Various binaries can be downloaded from the official UCSF website.

To check the installation, simply run the following Scipion test: 

``scipion test motioncorr.tests.test_protocols_motioncor2.TestMotioncor2AlignMovies``

Supported versions
------------------

1.2.3, 1.2.6, 1.3.0, 1.3.1

Protocols
---------

* movie alignment

References
----------

1.  Shawn Q Zheng, Eugene Palovcak, Jean-Paul Armache, Kliment A Verba, Yifan Cheng & David A Agard. MotionCor2: anisotropic correction of beam-induced motion for improved cryo-electron microscopy. Nature Methods volume 14, pages 331–332 (2017).
