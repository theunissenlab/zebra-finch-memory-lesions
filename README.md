# Zebra Finch Auditory Memory: Behavior, Lesions, and Ephys

## Data

Summaries of behavioral data and ephys data are located in [BehavioralDataSummary.md](BehavioralDataSummary.md) and [EphysDataSummary.md]( EphysDataSummary.md).

Data directory should be placed in the top level of this project (i.e. `zebra-finch-memory-lesions/data`) and can be downloaded from [Google Drive](https://drive.google.com/drive/folders/1M76aCU6dXOHVGbm1duyboV_wILjh0Ovq?usp=sharing). Last count shows it as 356MB zipped. I excluded 20GB of unit waveforms data which would be located in `data/ephys/unit_waveforms`, so the function `load_data.load_unit_waveforms(unit_id)` will not be usable.

### Subject Data

| Subject Name | Sex |
|--------------|---|
|GreBla5671F   | F |
|GreBla7410M   | M |
|WhiRed9510F   | F |
|RedBla0907M   | M |
|XXXOra0037F   | F |
|HpiBlu6194F   | F |
|YelPur7906M   | M |
|WhiWhi2526M   | M |
|BluYel2571F   | F |
|YelRed3010F   | F |
|GraWhi4040F   | F |
|BlaGre1349M   | M |
|XXXHpi0038M   | M |
|GreBlu5039F   | F |
|GreBla3404M   | M |
|XXXRed0088M   | M |
|XXXOra0039F   | F |
|XXXBla0054F   | F |
|XXXBla0055M   | M |
|XXXBla0081M   | M |

## Code

This project implements functions for loading behavioral, lesion, and ephys data for our auditory memory project, as well as code for data analysis.

### Dependencies

Code is run and tested on Python3.6, using requirements listed in requirements.txt (install with `pip install -r requirements.txt`). A local installation of R is required for statistical functions and installation of `rpy2` library.
