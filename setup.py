from setuptools import find_packages, setup

VERSION = "0.0.1"
DESCRIPTION = "abaqus rve simulator "
LONG_DESCRIPTION = "A package for running rve simulation via abaqus "

# Setting up
setup(
    name="rvesimulator",
    version=VERSION,
    author="Jiaxiang Yi (Delft University of Technology)",
    author_email="<J.Yi@tudelft.nl>",
    url="https://github.com/bessagroup/rvesimulator.git",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages("src"),
    package_dir={'': "src"},
    readme="README.md",
    install_requires=["numpy", "matplotlib", "scipy"],
    keywords=["python", "rve simulation", "abaqus", "multiscale simulation"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ],
)
