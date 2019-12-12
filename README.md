# Fluorescence Viewer
*FluorescenceViewer.m* is a short script for viewing fluorescence spectra. It requires the folder containing fluorescent spectra to be in the same directory as the .m file.



To add spectra, save them as a .csv file in the Spectra Raw Data folder with the file name being the name of the fluorophore and the columns being: Wavelength, excitation, wavelength, emission, wavelength, two-photon excitation.



You can save the currently plotted list of fluorophores/lasers/detectors as a .csv file, or load a new set of settings from an existing .csv file.



To get started run *FluorescenceViewer* from the MATLAB command line, this will populate the plots with some defaults. Next, either edit these settings directly in the user interface, or save the defaults as a .csv by clicking the save button then edit the .csv file and load the new settings by clicking the load button.

![Main User Interface](https://github.com/DrStoltzfus/fluorescenceviewer/blob/master/ExampleOutputs/MainUserInterface.png?raw=true)

![Plot Output](https://github.com/DrStoltzfus/fluorescenceviewer/blob/master/ExampleOutputs/FluorescencePlots.png?raw=true)





