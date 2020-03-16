## Metformin Resistance Project

Description of the scripts for the Metformin Resistance Project


**GC_plotter.R**:

This script allows to plot growth curves from processed Biolog plate reader files. To process raw files, you will need to use the Growth.py script from [Pov's repository](https://github.com/PNorvaisas/Growth_analysis). 

As an example, some files have been included in the folder _biolog_example_files_, where you can find 8 .txt files, the design file, and all pattern media files with the conditions corresponding to each plate. 

To use this R script, you will need to have installed some R packages in your computer. Just open an interactive R session, and copy and paste this:

```
install.packages(c('optparse', 'tidyverse', 'readxl', 'here', 'viridis'))
```

After installing all these R libraries, go to the terminal and type, in the same folder as you have your script, this command for help:

```
Rscript GC_plotter.R -h
```

The use of the script is quite simple, just write:

```
Rscript GC_plotter.R -i path/Timeseries.csv -o path/output_folder
```

If the output folder does not exist, it will be created. At the end of the processing, you will have as many pdf files as plates you are analysing. You might have a file named _Rplots.pdf_, which is a bug from the ggplot2 library that has not been solved at the moment. You can just delete it.
