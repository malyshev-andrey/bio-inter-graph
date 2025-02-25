from os import path
from setuptools import setup, find_packages


__version__ = None
version_path = path.join('biointergraph', '__version__.py')
exec(open(version_path).read())
assert __version__ is not None, f'Something is wrong with {version_path} file!'

setup(
    name="bio-inter-graph",
    version=__version__,
    packages=find_packages(),
    install_requires=[
        'numpy>=1.21.0',
        'pandas>=1.3.0',
        'tqdm>=4.42.0',
        'requests>=2.20.0',
        'scipy>=1.0.0,<2.0.0',
        'PyRanges>=0.0.114',
        'joblib>=1.2.0,<2.0.0',
        'networkx>=2.6'
    ],
    python_requires=">=3.10",
    description="Graph Algorithms for Biological Interactomes",
    author="Andrey Malyshev",
    author_email="malyshew.a.d@yandex.ru",
    url="https://github.com/malyshev-andrey/bio-inter-graph",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_data={
        'bio-inter-graph': ['../static/*']  # Adjust to your static files path
    },
    include_package_data=True
)
