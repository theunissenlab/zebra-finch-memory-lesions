# Description of Ephys Data Files

## Data Files

Data is organized into pandas dataframes or summarized below. Helper functions for loading the data are located in the ffile "code/load_data.py". Data files are in (data/ephys/)

* Spike Data: Neural responses indexed by unit and stim, with spike times for all trials (data/ephys/UnitData.pkl)

* Stimulus Info: Information about each stimulus file played during anesthetized recordings as well as time and amplitude arrays (data/ephys/StimData.pkl)

* Stimulus Files: WAV files played as stimuli during the recordings (data/ephys/stimuli)

* Unit Waveforms: The spike waveform snippets for every unit, saved in separate npz files in (data/ephys/unit_waveforms)


