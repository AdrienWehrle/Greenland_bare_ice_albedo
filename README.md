# Greenland_bare_ice_albedo

+ These tools have been developed using a conda virtual environment that can be identically recreated. To this end, create a new      
  environment using [[./PROMICE.yml]] as below:
  ```bash
  conda env create -f PROMICE.yml
  ```
  ```ResolvePackageNotFound``` error can be raised. In that case, run ```bash conda env export --no-builds > PROMICE.yml``` instead. \\
  Then, run =conda activate PROMICE= to activate this new virtual environment.

+ The repository can also be run interactively on Binder, without any download:

  [[https://mybinder.org/v2/gh/AdrienWehrle/Greenland_bare_ice_albedo/master][https://mybinder.org/badge_logo.svg]]
  
 
+ The three functions presented in [[#PROMICE_processing_toolspy][PROMICE_processing_tools.py]] can be used as follows:

  ```python
  import sys
  sys.path.append(path/to/module/)
  import PROMICE_processing_tools as ppt

  #load PROMICE dataset for a given station, all available years
  ds=ppt.load_data(file='path/to/dataset.txt', year='all')

  #process ice ablation and albedo time series around the onset of bare ice conditions 
  ds_proc=ppt.BIC_processing(ds, visualisation=True)

  #compute multi-year and multi-site composites for air temperature, snow height, ice ablation and albedo time series
  #centered on bare ice appearance and spanning ± dt days
  composite=ppt.BIC_composite('path/to/dataset/folder/', dt=45)
  ```
