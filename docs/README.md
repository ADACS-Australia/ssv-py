# Python Simple Spectrum Viewer
The Simple Spectrum Viewer package has been written to allow a simple, quick method for visualising spectra.

## Developer hints
### Requirements
Python 3.7+
virtualenv

### Install
First obtain a clone of this repository with
`git clone https://github.com/ADACS-Australia/ssv-py.git`

We recommend starting with a fresh virtual environment
`virtualenv venv`
`source venv/bin/activate`

Then, change to the ssv-py directory, and install the required python packages in your virtualenv, along with the ssv package
`cd ssv-py`
`pip install -r requirements.txt`
`pip install .`


### Running tests
Within the above virtual environment, run the command `tox` to run the test suite.
