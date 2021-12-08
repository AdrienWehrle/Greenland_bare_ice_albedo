[![DOI](https://zenodo.org/badge/273865759.svg)](https://zenodo.org/badge/latestdoi/273865759)
# Greenland bare ice albedo

Post-processing tools for a dataset of PROMICE daily ice-ablation, albedo, snow height and temperature measurements (225 station years from 26 stations, 2007-2019) including corrections for signal shifts caused by station movement, sensor reinstallation and measurement failure. Using the resulting processed ice-ablation measurements in combination with seasonal snow layer thickness and air temperature, the seasonal timing of bare-ice onset is then determined for each station year. The albedo value at the timing of bare-ice onset is finally extracted to compute an average Greenland bare-ice albedo at ice-ablation onset (called bare-ice-onset albedo; [Wehrlé et al, 2021](https://geusbulletin.org/index.php/geusb/article/view/5284)). Accurate definition of this variable has applications in classifying the bare-ice area over large areas of the ice sheet, in constraining polar regional climate models used to estimate the surface mass balance of the Greenland ice sheet and in climate monitoring.

![](https://geusbulletin.org/index.php/geusb/article/download/5284/12394/41081)

*Fig. 1 Multi-site and multi-year composite surface conditions synchronised to bare-ice onset from ice ablation (black vertical dashed lines). a: Air temperature. b: Snow height. c: Ice ablation. d: Albedo where the red horizontal dashed line indicates the bare-ice-onset albedo. Grey shading corresponds to ± one standard deviation around daily averages. Figure 2 of [Wehrlé et al, 2021](https://geusbulletin.org/index.php/geusb/article/view/5284).*

+ These tools have been developed using a conda virtual environment that can be identically recreated. To this end, 
create a new environment using [PROMICE.yml](https://github.com/AdrienWehrle/Greenland_bare_ice_albedo/blob/master/PROMICE.yml) as below:  
  ```bash
  conda env create -f PROMICE.yml
  ```

+ The repository can also be run interactively on Binder using [PROMICE_processing_tools.ipynb](https://github.com/AdrienWehrle/Greenland_bare_ice_albedo/blob/master/PROMICE_processing_tools.ipynb), without any download:

  [<img src="https://mybinder.org/badge_logo.svg">](https://mybinder.org/v2/gh/AdrienWehrle/Greenland_bare_ice_albedo/master)

 
+ The three functions presented in [PROMICE_processing_tools.py](https://github.com/AdrienWehrle/Greenland_bare_ice_albedo/blob/master/PROMICE_processing_tools.py) can be used as follows:

  ```python
  import sys
  import PROMICE_processing_tools as ppt
  sys.path.append(path/to/module/)

  # load PROMICE dataset for a given station, all available years
  ds = ppt.load_data(file='path/to/dataset.txt', year='all')

  # process ice ablation and albedo time series around the onset of bare ice conditions 
  ds_proc = ppt.BIC_processing(ds, visualisation=True)

  # compute multi-year and multi-site composites for air temperature, snow height, ice ablation  
  # and albedo time series centered on bare ice appearance and spanning ± dt days
  composite = ppt.BIC_composite('path/to/dataset/folder/', dt=45)
  ```
