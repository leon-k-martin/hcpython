from setuptools import setup, find_packages

setup(
    name="hcpy",
    version="1.0.0",
    description="A Python package for HCP data processing",
    author="Leon Martin",
    author_email="leon.martin@bih-charite.de",
    url="https://github.com/yourusername/hcpype",
    packages=find_packages(include=["hcpy", "hcpy.*"]),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "nipype",
        "nibabel",
        "nilearn",
        "pybids",
        # Add any other dependencies here
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
    ],
)
