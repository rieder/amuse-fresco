import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='amuse_fresco',
    version='0.2.0',
    description='Make a visualisation of a star cluster',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/rieder/fresco',
    author='Steven Rieder',
    author_email='steven@rieder.nl',
    license='MIT',
    packages=setuptools.find_packages(),
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    classifiers=[
        # Python versions supported by amuse_masc
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",

        # License
        "License :: OSI Approved :: MIT License",

        # OS support
        "Operating System :: OS Independent",

        # Maturity of amuse_masc
        "Development Status :: 4 - Beta",

        # Intended audience
        "Intended Audience :: Science/Research",
    ],
)

