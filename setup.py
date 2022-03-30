from setuptools import setup, find_packages


setup(
    name = "zebra-finch-memory-lesions",
    version = "1.0",
    packages = ["zf_data"],
    include_package_data = True,
    zip_safe = False,
    description = "Code and notebooks for zebra finch vocalizer memory paper",
    author = "Theunissen Lab",
    author_email = "theunissen@berkeley.edu",
    url = "https://github.com/theunissenlab/zebra-finch-memory-lesions",
    keywords = "",
    classifiers = [
        "Programming Language :: Python :: 3.0",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        "numpy",
        "matplotlib",
        "pandas",
        "scipy",
        "SoundFile",
        "statsmodels",
        "jupyterlab",
    ],
)
