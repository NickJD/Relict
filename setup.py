from setuptools import setup, find_packages

setup(
    name="phylo16s",
    version="0.1",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "phylo16s=phylo16s.cli:main"
        ]
    },
)
