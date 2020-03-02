import sys
import platform
import subprocess

import distutils
import setuptools

from distutils.core import Extension
from setuptools.command.build_ext import build_ext
"""
setup.py support multiple commands. The default for building C/C++ code is
build_ext, and we need to overload it in order to call the build chain

There are some so-called Extensions to handle Makefiles etc., but setuptools
is such a voodoo piece of code difficult to debug, understand and use, that for
our case it is easier and safer to call the automated chain from processes.

In order to not confuse the installation procedures between usual repository 
clone and "pip install", the goal is to keep the toolchain as much as possible
unchanged. For that reason the AutoconfiCommand class is just a way to step into
the extensions building process of setup.py and follow afterwards step-by-step
the usual installation instructions.
"""

with open("README.md", "r") as fh:
    long_description = fh.read()

# Read in requirements
requirements = open('scripts/requirements.txt').readlines()
requirements = [r.strip() for r in requirements]


class AutoconfigCommand(build_ext):
    """A custom command to compile qFlex"""

    description = 'Compile qFlex python binary'

    def run(self):
        # Compile before install
        self.autoconf()

        # This fake call will generate a so lib,
        # but the next method will overwrite it
        build_ext.run(self)

        # Make the pybind11 interface
        self.make()

    def autoconf(self):
        try:
            out = subprocess.check_output(['autoconf', '--version'])
        except OSError:
            raise RuntimeError(
                'Autoconf must be installed to build the following extensions: '
                + ', '.join(e.name for e in self.extensions))
        if platform.system() == 'Windows':
            raise RuntimeError('For Windows use Docker image of QFlex.')
        """
          Run autoreconf.
        """
        self.runcommand(['autoreconf', '-i'])
        """
          Run autoconf.
        """
        self.runcommand(['autoconf'])
        """
          Run configure.
        """
        self.runcommand(['./configure'])
        """
        Git
        """
        try:
            # After pip install the folder is not a git repo
            self.runcommand(['git', 'init', '.'])
            # Add the docopt submodule
            self.runcommand([
                'git', 'submodule', 'add',
                'https://github.com/docopt/docopt.cpp.git', 'src/docopt'
            ])
        except:
            e = sys.exc_info()[0]
            print("SUBMODULE ERROR...", e)

    def make(self):
        # The destination should not be the one specified in the Makefile
        # but the one necessary to get the .so in the wheel
        # CIRQ_IF_DIR is defined in src/Makefile.in
        destination_lib = "CIRQ_IF_DIR=\'../" + self.build_lib + "/qflexcirq/\'"

        # Override the destination path
        # http://www.gnu.org/software/make/manual/make.html#Overriding
        command = ['make', '-j8', destination_lib, 'pybind']

        self.runcommand(command)

    def runcommand(self, command):
        self.announce('Running command: %s' % str(command),
                      level=distutils.log.INFO)
        subprocess.check_call(command)


setuptools.setup(
    name="qflexcirq",
    version="0.0.5",
    author="The QFlex Contributors",
    description="A Google Cirq interface for QFlex",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ngnrsaa/qflex",
    python_requires=('>=3.6'),
    install_requires=requirements,
    license='Apache 2',

    # Consider these directories as parts of the package
    # Will be copied to the wheel
    packages=[
        'qflexcirq', 'qflexcirq/interface', 'qflexcirq/circuits',
        'qflexcirq/ordering'
    ],

    # This is just to fool setuptools to start build_ext
    # The compiled library will be called 'qflex'
    # and will be copied in 'qflexcirq'
    ext_modules=[Extension('qflexcirq.qflex', sources=[])],

    # A customised command to start the C++ compilation
    cmdclass={
        'build_ext': AutoconfigCommand,
    },
)
