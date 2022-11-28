from setuptools import setup, find_packages

from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

setup(
    name='openrxn',
    version='0.0.1',

    author="Alex Dickson",
    author_email="alexrd@msu.edu",
    description="OpenRXN",
    license="MIT",
    url="https://gitlab.com/ADicksonLab/OpenRXN.git",
    classifiers=[
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
        'Programming Language :: Python :: 3'
    ],

  # package
    packages=find_packages(where='src'),

    package_dir={'' : 'src'},

    # if this is true then the package_data won't be included in the
    # dist. Use MANIFEST.in for this
    include_package_data=True,


    # SNIPPET: this is a general way to do this
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],

    entry_points = {},

    install_requires=[
        'numpy',
        'scipy',
        'networkx',
        'pint',
        'pandas',
        'matplotlib'
    ],
)
