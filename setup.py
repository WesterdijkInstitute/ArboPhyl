import setuptools

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name="ArboPhyl",
    version="v1.1.0",
    author="Tim Verschuren",
    author_email="t.verschuren@wi.knaw.nl",
    description= "ArboPhyl is a BUSCO based pipeline for the construction of\
                Phylogenetic trees",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WesterdijkInstitute/ArboPhyl",
    project_urls={
        "Bug Tracker":"https://github.com/WesterdijkInstitute/ArboPhyl/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=["ArboPhyl_src/ArboPhyl.sh"],
    install_requires=["Bio"],
    packages=setuptools.find_packages(),
    python_requires=">=3.10",
    entry_points={
        "console_scripts": [
            "arbophyl = ArboPhyl_src.cli:main",
        ]
    }
)
