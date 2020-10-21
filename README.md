# Zebra Finch Auditory Memory: Behavior, Lesions, and Ephys

## Data

Summaries of behavioral data and ephys data are located in [BehavioralDataSummary.md](BehavioralDataSummary.md) and [EphysDataSummary.md]( EphysDataSummary.md).

Data directory should be placed in the top level of this project (i.e. `zebra-finch-memory-lesions/data`) and can be downloaded from [Google Drive](https://drive.google.com/drive/folders/1M76aCU6dXOHVGbm1duyboV_wILjh0Ovq?usp=sharing). Last count shows it as 356MB zipped. I excluded 20GB of unit waveforms data which would be located in `data/ephys/unit_waveforms`, so the function `load_data.load_unit_waveforms(unit_id)` will not be usable.

## Code

This project implements functions for loading behavioral, lesion, and ephys data for our auditory memory project, as well as code for data analysis.

### Dependencies

Code is run and tested on Python3.6, using requirements listed in requirements.txt (install with `pip install -r requirements.txt`). A local installation of R is required for statistical functions and installation of `rpy2` library.
