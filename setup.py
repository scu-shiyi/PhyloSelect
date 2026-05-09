from setuptools import setup, find_packages


setup(
    name="phyloselect",
    version="1.1.0",
    description="PhyloSelect: A comprehensive toolkit for phylogenetic analysis and evolutionary selection.",
    author="scu-shiyi",
    author_email="shi@stu.scu.edu.cn",
    url="https://github.com/scu-shiyi/PhyloSelect",
    packages=find_packages(include=["phyloselect", "phyloselect.*"]),
    package_data={"phyloselect": ["utils/*.ttf"]},
    python_requires=">=3.10",
    install_requires=[
        "biopython",
        "matplotlib",
        "numpy",
        "pandas",
        "scipy",
        "seaborn",
        "psutil",
        "ete3",
        "requests",
    ],
    entry_points={
        "console_scripts": [
            "phyloselect=phyloselect.cli:main",
        ]
    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
)