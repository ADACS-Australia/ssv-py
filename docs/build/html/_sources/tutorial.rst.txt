Simple Spectrum Viewer Tutorial
===============================

Here we'll just cover the very basics of creating Charts with the ``ssv`` python package.

Loading spectrum data
---------------------

We make use of the Astropy :doc:`Registry <astropy:io/registry>` to help load spectrum data into a format usable by SSV.
We have provided the function :func:`~ssv.utils.read_spectra_file` to help with loading spectrum data.
This function returns a SpectrumList object, which is able to be read by our :class:`~ssv.viewer.SimpleSpectrum` class.

.. code-block:: python

   >>> from ssv.utils import read_spectra_file
   >>> spectrum_data = read_spectra_file('path/to/file.fits')

It's worth noting that you can also directly pass a :class:`~astropy.io.fits.HDUList` object into :func:`~ssv.utils.read_spectra_file` as well.

SimpleSpectrum
--------------

The :class:`~ssv.viewer.SimpleSpectrum` class is the basic container class to hold and manipulate the data for the spectrum charts.
We can create a new SimpleSpectrum containing the spectrum data as follows:

.. code-block:: python

   >>> from ssv.viewer import SimpleSpectrum
   >>> spectrum = SimpleSpectrum('Tutorial', spectrum_data)

We can immediately build an Altair :class:`~altair:altair.Chart` by calling the :meth:`~ssv.viewer.SimpleSpectrum.SimpleSpectrum.build_chart` method.
This method returns a Chart object, which should be viewed in a web browser.

SimpleSpectrumViewer
--------------------

The SimpleSpectrum object itself implements no features to be able to edit the Chart. To fill in this role, we use the :class:`~ssv.viewer.SimpleSpectrumViewer` object.
We can create the SimpleSpectrumViewer and add a SimpleSpectrum object to it, giving us much more control over the Chart.

.. code-block:: python

   >>> from ssv.viewer import SimpleSpectrumViewer
   >>> viewer = SimpleSpectrumViewer()
   >>> viewer.add_spectrum(spectrum)

We are able to add multiple SimpleSpectrum objects to the viewer using the :meth:`~ssv.viewer.SimpleSpectrumViewer.SimpleSpectrumViewer.add_spectrum` method. Calling the :meth:`~ssv.viewer.SimpleSpectrumViewer.SimpleSpectrumViewer.build_chart` method of the SimpleSpectrumViewer then calls the ``build_chart()`` method from the component spectra and plots them together.
Further, we can begin to control some aspects of the Chart visuals, such as toggling a grid and setting the axis titles.

SimpleSpectralLines
-------------------

One common feature of spectrum plots is showing spectral lines. We have implemented another class, :class:`~ssv.viewer.SimpleSpectralLines`, for this purpose.
By default, the SimpleSpectralLines class contains a reasonably comprehensive set of common spectral lines. This set can be specified when a new instance of the class is defined.
SimpleSpectralLines objects can be included in the Chart by using the :meth:`~ssv.viewer.SimpleSpectrumViewer.SimpleSpectrumViewer.add_lines` method of the SimpleSpectrumViewer class.

