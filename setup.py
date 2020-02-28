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

# HOWTO
# https://jichu4n.com/posts/how-to-add-custom-build-steps-and-commands-to-setuppy/

# class AutoconfigCommand(distutils.cmd.Command):
class AutoconfigCommand(build_ext):
  """A custom command to compile qflex"""

  description = 'Compile qflex python binary'

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
      command = ['autoreconf', '-i']
      self.announce(
          'Running command: %s' % str(command),
          level=distutils.log.INFO)
      subprocess.check_call(command)
      """
        Run autoconf.
      """
      command = ['autoconf']
      self.announce(
          'Running command: %s' % str(command),
          level=distutils.log.INFO)
      subprocess.check_call(command)
      """
        Run configure.
      """
      command = ['./configure']
      self.announce(
          'Running command: %s' % str(command),
          level=distutils.log.INFO)
      subprocess.check_call(command)


  def make(self):
      # The destination should not be the one specified in the Makefile
      # but the one necessary to get the SO in the wheel
      destination_lib = 'QFLEXLIB=\'../{}\''.format(
          "build/lib.linux-x86_64-3.6/qflexcirq/qflex.cpython-36m-x86_64-linux-gnu.so"
      )

      command = ['make', destination_lib , 'pybind']

      # command.append(os.getcwd())
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

    # This is just to mock setuptools to start build_ext
    ext_modules=[Extension('qflexcirq.qflex', sources = [])],

    # A customised command to start the compilation
    cmdclass={
        'build_ext': AutoconfigCommand,
    },
)