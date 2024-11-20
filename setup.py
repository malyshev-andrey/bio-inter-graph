from setuptools import setup, find_packages

setup(
    name="bio-inter-graph",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'numpy>=1.21.0',
        'pandas>=1.3.0'
    ],
    python_requires=">=3.7",
    description="Graph Algorithms for Biological Interactomes",
    author="Andrey Malyshev",
    author_email="malyshew.a.d@yandex.ru",
    url="https://github.com/malyshev-andrey/bio-inter-graph",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
