from distutils.core import setup

setup(
    name='AqPy',
    version='0.1.0',
    author='Michael Skarbek',
    author_email='maskarb@gmail.com',
    packages=['aqpy',],
    license='MIT license',
    long_description=open('README.md').read(),
    install_requires=[
        "matplotlib",
        "numpy",
        "gurobipy",
    ],
)