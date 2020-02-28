import os
import sys
import platform
import subprocess

import distutils
import setuptools

from distutils.core import Extension

from setuptools.command.build_ext import build_ext

with open("README.md", "r") as fh:
    long_description = fh.read()

# Read in requirements
requirements = open('scripts/requirements.txt').readlines()
requirements = [r.strip() for r in requirements]


# HOWTO and why
# https://jichu4n.com/posts/how-to-add-custom-build-steps-and-commands-to-setuppy/

class AutoconfigCommand(build_ext):
    """A custom command to compile qFlex"""

    description = 'Compile qFlex python binary'

    def run(self):

        # Compile before install
        self.autoconf()

        # This fake call will generate a lib, but the next method will overwrite it
        build_ext.run(self)

        # Make the pybind11 interface
        self.make()

    def autoconf(self):
        try:
            out = subprocess.check_output(['autoconf', '--version'])
        except OSError:
            raise RuntimeError(
                'Autoconf must be installed to build the following extensions: ' +
                ', '.join(e.name for e in self.extensions))
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

    def make(self):
        # The destination should not be the one specified in the Makefile
        # but the one necessary to get the SO in the wheel
        destination_lib = "CIRQ_IF_DIR=\'../"+self.build_lib+"/qflexcirq/\'"

        # Override the destination path
        # http://www.gnu.org/software/make/manual/make.html#Overriding
        command = ['make',
                   destination_lib,
                   'pybind']

        self.runcommand(command)

    def runcommand(self, command):
        self.announce(
            'Running command: %s' % str(command),
            level=distutils.log.INFO)
        subprocess.check_call(command)


setuptools.setup(
    name="qflexcirq",
    version="0.0.1",
    author="The QFlex Contributors",
    description="A Google Cirq interface for QFlex",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ngnrsaa/qflex",
    packages=['qflexcirq',
              'qflexcirq/interface',
              'qflexcirq/circuits',
              'qflexcirq/ordering'],
    python_requires=('>=3.6'),
    install_requires=requirements,
    license='Apache 2',

    # This is just to fool setuptools to start build_ext
    ext_modules=[Extension('qflexcirq.qflex', sources=[])],

    # A customised command to start the compilation
    cmdclass={
        'build_ext': AutoconfigCommand,
    },
)
