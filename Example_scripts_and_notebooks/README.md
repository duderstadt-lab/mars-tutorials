This examples database contains examples that are briefly discussed in the [mars docs]() as well as some example scripts for functions also accessible through the user interfaces.

### Discussed Example Scripts
[1. Plot the Kinetic Change Point Rates with Python](https://github.com/duderstadt-lab/mars-tutorials/blob/master/Example_scripts_and_notebooks/KCP_widget_and_jupyter_plot.ipynb)  
Jupyter notebook showing how to retrieve the segment tables in a Python environment and how to plot the calculated rates with seaborn as well as in the **Rover** scriptable widgets.

[2. Generate a Table of Parameter Values from the Archive](https://github.com/duderstadt-lab/mars-tutorials/blob/master/Example_scripts_and_notebooks/Generate_a_table_of_parameter_values.groovy)  
Groovy script to make a MarsTable showing the value for the specified parameter for each molecule entry. Run the script and a pop-up dialog will be created in which the parameter of interest can be selected.

[3. Plot the Tracking Results as Color Coded Overlay in the Image](https://github.com/duderstadt-lab/mars-tutorials/blob/master/Example_scripts_and_notebooks/Color_coded_tracks_overlay)  
Script to show the tracking results in the MoleculeArchive as color coded ROIs on the original video. Also possible to show the overlay for tagged molecules only (see (repository example)[https://github.com/duderstadt-lab/mars-tutorials/blob/master/Example_scripts_and_notebooks/Color_coded_tracks_overlay_tagged.groovy]).


### Example Scripts
1. Add flowrates from a table to the metadata
2. Add flowrates from the metadata to molecule records
3. Calculate the slope between the start and the end of a trace
4. Create color coded track overlays for all molecules
5. Create color coded track overlays for tagged molecules
6. Convert from pixels to biologically relevant length
7. Create a MarsTable
8. Filter molecules by parameter values and tag
9. Generate a table of parameter values
10. Get a molecule DataTable for a specified UUID
11. Get the intensity for the first time point of all molecules
12. Get the number of molecules with a specified tag
13. Kinetic Change Point results analysis in a Jupyter notebook
