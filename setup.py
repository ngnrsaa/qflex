import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

# Read in requirements
requirements = open('scripts/requirements.txt').readlines()
requirements = [r.strip() for r in requirements]

setuptools.setup(
    name="qflex-cirq", # Replace with your own username
    version="0.0.1",
    author="The QFlex Contributors",
    description="A Google Cirq interface for QFlex",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ngnrsaa/qflex",
    packages=setuptools.find_packages(),
    python_requires=('>=3.6.0'),
    install_requires=requirements,
    license='Apache 2',
)
