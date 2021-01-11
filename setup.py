'''
Setup for the ssv package
'''
import setuptools
import versioneer

setuptools.setup(
    name='ssv1',
    version='0.0.1',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
        "altair",
        "astropy",
        "pandas",
        "specutils",
    ],
    python_requires='>=3.7',
    author='Thomas Reichardt, Ray Seikel, James Tocknell',
    author_email='thomas.reichardt@students.mq.edu.au, rseikel@bigpond.com, james.tocknell@mq.edu.au',
    description='Simple Spectra Viewer for Data Central',
    license='MIT',
    keywords=[
        'spectrum',
        'viewer',
        'marz'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3 :: Only',
    ],
)
