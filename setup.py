try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='AqPy',
    version='0.1.0',
    author='Michael Skarbek',
    author_email='maskarb@gmail.com',
    packages=['aqpy',],
    license='MIT license',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=[
        "matplotlib",
        "numpy",
        "gurobipy",
    ],
)
