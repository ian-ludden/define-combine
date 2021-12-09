from setuptools import find_packages, setup

requirements = [
    "gerrychain", 
    "networkx", 
    "scipy"
]

setup(
    name="definecombine", 
    description="Analyze the define-combine procedure for political redistricting", 
    author="Ian Ludden", 
    author_email="iludden2@illinois.edu", 
    maintainer="Ian Ludden", 
    maintainer_email="iludden2@illinois.edu", 
    url="https://github.com/ian-ludden/define-combine", 
    install_requires=requirements, 
    keywords="Define-combine"
)
