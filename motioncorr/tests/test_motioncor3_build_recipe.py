"""MotionCor3 binary-build recipe regression tests for ``scipion3 test``.

These are pure unit tests: the install command is inspected but never executed.
Run with:
    scipion3 test motioncorr.tests.test_motioncor3_build_recipe.TestMotioncor3BuildRecipe
"""

import unittest
from unittest.mock import patch

import pwem
from motioncorr import Plugin
from motioncorr.constants import MOTIONCOR_BIN


class TestMotioncor3BuildRecipe(unittest.TestCase):
    """Ensure the compiler toolchain and Conda prefix are consistent."""

    def _build_command(self):
        packages = []

        class FakeEnvironment:
            def addPackage(self, name, **kwargs):
                packages.append((name, kwargs))

        # Do not depend on or alter the user's configured CUDA or plugin vars.
        with patch.object(pwem.Config, "CUDA_LIB", "/usr/local/cuda-11.8/lib64"), \
                patch.object(Plugin, "_supportedVersions", ["1.2.4"]), \
                patch.object(
                    Plugin,
                    "getVar",
                    new=classmethod(
                        lambda cls, name, default=None: (
                            "MotionCor3" if name == MOTIONCOR_BIN else default
                        )
                    ),
                ):
            Plugin.defineBinaries(FakeEnvironment())

        self.assertEqual(len(packages), 1)
        package_name, options = packages[0]
        self.assertEqual(package_name, "motioncor3")
        commands = options["commands"]
        self.assertEqual(len(commands), 1)
        parts, marker = commands[0]
        self.assertEqual(marker, "MOTIONCOR_BIN_installed")
        self.assertIsInstance(parts, list)
        return " ".join(parts)

    def test_MotionCor3BinaryBuildUsesCompatibleHostCompilerForCudaAndCpp(self):
        """Use a CUDA-compatible host compiler for NVCC and generated C++."""
        command = self._build_command()

        self.assertIn("g++-11", command)
        self.assertIn("MOTIONCOR_HOST_CXX", command)
        self.assertIn("CC=", command)
        self.assertIn("NVCC=", command)
        self.assertTrue("-ccbin" in command or "--compiler-bindir" in command)

    def test_MotionCor3BinaryBuildUsesActualCondaEnvironment(self):
        """Never use a hard-coded, unrelated Miniforge prefix."""
        command = self._build_command()

        self.assertNotIn("$(HOME)/miniforge3", command)
        self.assertIn("CONDA_PREFIX", command)
