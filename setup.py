from setuptools import setup

version = "0.5.5"
name = 'amuse-fresco'
author = 'Steven Rieder and Inti Pelupessy'
author_email = 'steven@rieder.nl'
license_ = "MIT"
url = 'http://amusecode.org'

classifiers=[
    # Python versions supported by amuse-fresco
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",

    # License
    "License :: OSI Approved :: MIT License",

    # OS support
    "Operating System :: OS Independent",

    # Maturity of amuse-fresco
    "Development Status :: 4 - Beta",

    # Intended audience
    "Intended Audience :: Science/Research",
]

install_requires = [
    'wheel>=0.32',
    'amuse-framework>=12.0.0',
    'scipy',
    'matplotlib',
    'astropy',
]
description = 'Make a realistic visualisation of a star cluster'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

extensions = []

all_data_files = [
]

packages = ['amuse.ext.fresco']

package_dir = {
    'amuse.ext.fresco': 'src/amuse/ext/fresco'
}

package_data = {
}

setup(
    name=name,
    version=version,
    classifiers=classifiers,
    url=url,
    author_email=author_email,
    author=author,
    license=license_,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    install_requires=install_requires,
    ext_modules=extensions,
    package_dir=package_dir,
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
    include_package_data=True,
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    scripts=["fresco.py"],
)
