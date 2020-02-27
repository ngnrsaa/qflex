import setuptools

import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

with open("README.md", "r") as fh:
    long_description = fh.read()

# Read in requirements
requirements = open('scripts/requirements.txt').readlines()
requirements = [r.strip() for r in requirements]

class CMakeExtension(Extension):

  def __init__(self, name, sourcedir=''):
    Extension.__init__(self, name, sources=[])
    self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):

  def run(self):
    try:
      out = subprocess.check_output(['autoconf', '--version'])
    except OSError:
      raise RuntimeError(
          'Autoconf must be installed to build the following extensions: ' +
          ', '.join(e.name for e in self.extensions))

    if platform.system() == 'Windows':
        raise RuntimeError('For Windows use Docker image of QFlex.')

    for ext in self.extensions:
      self.build_extension(ext)

  def build_extension(self, ext):
    extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
    # cmake_args = [
    #     '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
    #     '-DPYTHON_EXECUTABLE=' + sys.executable
    # ]
    #
    # cfg = 'Debug' if self.debug else 'Release'
    # build_args = ['--config', cfg]
    #
    # if platform.system() == 'Windows':
    #   cmake_args += [
    #       '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)
    #   ]
    #   if sys.maxsize > 2**32:
    #     cmake_args += ['-A', 'x64']
    #   build_args += ['--', '/m']
    # else:
    #   cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
    #   build_args += ['--', '-j2']
    #
    env = os.environ.copy()
    # env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
    #     env.get('CXXFLAGS', ''), self.distribution.get_version())


    if not os.path.exists(self.build_temp):
      os.makedirs(self.build_temp)

    # subprocess.check_call(
    #     ['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
    # subprocess.check_call(
    #     ['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

    print("SELF SELF", self.build_temp)


    subprocess.check_call(
        ['autoreconf', '-i'],
        cwd=self.build_temp, env=env
    )
    subprocess.check_call(
        ['autoconf'],
        cwd=self.build_temp, env=env
    )
    subprocess.check_call(
        ['./configure'],
        cwd=self.build_temp, env=env
    )
    subprocess.check_call(
        ['make'], cwd=self.build_temp
    )

setuptools.setup(
    name="qflexcirq",
    version="0.0.1",
    author="The QFlex Contributors",
    description="A Google Cirq interface for QFlex",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ngnrsaa/qflex",
    packages=['qflexcirq'],
    python_requires=('>=3.6'),
    install_requires=requirements,
    license='Apache 2',

    ext_modules=[CMakeExtension('qflexcirq')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)