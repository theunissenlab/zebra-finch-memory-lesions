# Description of Behavioral Data Files

## Data Files

Data is organized into pandas dataframes or summarized below

* Subject Data: Information about each subject who took the test

* Trial Data: Subject response data from all trials in the test across all days (data/behavior/TrialData.csv)

* Test Contexts: The stages of the song learning ladder to cross reference with the Trial Data

* Stimulus Files: WAV files played as stimuli during the test, referenced by Trial Data (data/behavior/stimuli)

* Stimulus File Metadata: Sex, Age, Lab recorded from?

### Subject Data

### Trial Data

Trial data represents the stimulus playback and behavioral response data for all trials and subjects. Data has been cleaned to remove trials triggered erroneously due to hardware issues (e.g. peck double registering due to button sensitivity, button getting stuck in the down position).

| Column Name | Data Type | Description |
|-------------|-----------|-------------|
|Subject      |String     |             |
|Subject Group|String     |Experimental group (see Subject Group section below) |
|Trial        |Integer    |Trial number (within day) |
|Time         |Datetime   |System time of trial start|
|Date         |Date       |             |
|Interrupt    |Boolean    |True if subject pecked to interrupt playback|
|RT           |Float      |Response time (in seconds)|
|Stimulus File |String           |Name of stimulus wav file|
|Stimulus Vocalizer | String    | Name of vocalizing subject |
|Stimulus Call Type | String    | Call type (SO or DC) |
|Stimulus Class        |String     |"Rewarded" or "Unrewarded"|
|Rewarded     |Boolean    |True if subject received food reward (derived from Class and Interrupt)|
|Informative Trials Seen |Integer |Number of times that a stimulus from this vocalizer had previously been uninterrupted| 
|Test Context |String     | Test context (references Test Context table) |
|Condition    |String     | Normally NaN, but for some tests that occured after a month without reinforcement, indicated by "MonthLater" |
|Ladder Group |String     | Which stimulus set / ladder is being tested. One of PrelesionSet1, PostlesionSet1 (retest of previously learned data) and PostlesionSet2 (introduction of new vocalizers after lesion) |

### Test Context

|Ladder| Test Name   |  # Rewarded Vocalizers |  # Non-rewarded Vocalizers  | Description | Ladder Group  |
|------|-------------|-----------|-------------|---|---|
|Week 1 / Week 2      |SovsSo_1v1       | 1 song | 1 song |  | PrelesionSet1 |
|      |SovsSo_4v4       | 4 songs (3 new) | 4 songs (3 new) | | PrelesionSet1 |
|      |SovsSo_8v8_d1    | 8 songs (4 new) | 8 songs (4 new) | New vocalizers played twice as frequently | PrelesionSet1 |
|      |SovsSo_8v8_d2    | 8 songs | 8 songs | All vocalizers played at equal frequency | PrelesionSet1 |
|      |DCvsDC_1v1       | 1 dc | 1 dc | | PrelesionSet1 |
|      |DCvsDC_4v4       | 4 dcs (3 new) | 4 dcs (3 new) | | PrelesionSet1 |
|      |DCvsDC_6v6_d1    | 6 dcs (2 new) | 6 dc (2 new) | New vocalizers played twice as frequently | PrelesionSet1 |
|      |DCvsDC_6v6_d2    | 6 dcs | 6 dcs | All vocalizers played at equal frequency | PrelesionSet1 |

#### a. Lesion and sham lesion subjects

After Week 1 / 2, most subjects undergo procedure outlined by their Subject Group (lesion or sham). They are then retested on the first set (Weeks 3/4) and then tested on a brand new set of stimuli (Weeks 5/6). Subject groups that do not follow this pattern are `Subject Group == EPHYS` (did not undergo lesions), and those with `Subject Group == DEAD` (did not survive surgery and so do not have postlesion data).

|Ladder| Test Name   |  # Rewarded Vocalizers |  # Non-rewarded Vocalizers  | Description | Ladder Group  |
|------|-------------|-----------|-------------|---|---|
|Week 3      |SovsSo_1v1       | 1 song | 1 song | Refresher day | PostlesionSet1 |
|      |SovsSo_8v8_d2    | 8 songs | 8 songs | All vocalizers played at equal frequency | PostlesionSet1 |
|      |DCvsDC_1v1       | 1 dc | 1 dc | Refresher day | PostlesionSet1 |
|      |DCvsDC_6v6_d2    | 6 dcs | 6 dcs | All vocalizers played at equal frequency | PostlesionSet1 |
|Week 4 / Week 5|SovsSo_1v1_S2       | 1 song | 1 song |  | PostlesionSet2 |
|      |SovsSo_4v4_S2       | 4 songs (3 new) | 4 songs (3 new) | | PostlesionSet2 |
|      |SovsSo_8v8_d1_S2    | 8 songs (4 new) | 8 songs (4 new) | New vocalizers played twice as frequently | PostlesionSet2 |
|      |SovsSo_8v8_d2_S2    | 8 songs | 8 songs | All vocalizers played at equal frequency | PostlesionSet2 |
|      |DCvsDC_1v1_S2       | 1 dc | 1 dc | | PostlesionSet2 |
|      |DCvsDC_4v4_S2       | 4 dcs (3 new) | 4 dcs (3 new) | | PostlesionSet2 |
|      |DCvsDC_6v6_d1_S2    | 6 dcs (2 new) | 6 dc (2 new) | New vocalizers played twice as frequently | PostlesionSet2 |
|      |DCvsDC_6v6_d2_S2    | 6 dcs | 6 dcs | All vocalizers played at equal frequency | PostlesionSet2 |

#### b. No lesion subjects

In four subjects who did not undergo a lesion or a sham lesion, instead of breaking up the stimulus sets, we extended the ladder to include all stimuli. Note that though they see the same vocalizers as the lesioned subjects, their Ladder Groups are all labeled "PrelesionSet" since they do not undergo a lesion and all the stimuli are accumulated over days. These ladders proceed in the following way:

|Ladder| Test Name   |  # Rewarded Vocalizers |  # Non-rewarded Vocalizers  | Description | Ladder Group  |
|------|-------------|-----------|-------------|---|---|
|Week 3 / Week 4|SovsSo_1v1_S2       | 1 song | 1 song |  | PrelesionSet1 |
|      |SovsSo_4v4_S2       | 4 songs (3 new) | 4 songs (3 new) | | PrelesionSet1 |
|      |SovsSo_8v8_d1_S2    | 8 songs (4 new) | 8 songs (4 new) | New vocalizers played twice as frequently | PrelesionSet1 |
|      |SovsSo_8v8_d2_S2    | 8 songs | 8 songs | All vocalizers played at equal frequency | PrelesionSet1 |
|      |DCvsDC_1v1_S2       | 1 dc | 1 dc | | PrelesionSet1 |
|      |DCvsDC_4v4_S2       | 4 dcs (3 new) | 4 dcs (3 new) | | PrelesionSet1 |
|      |DCvsDC_6v6_d1_S2    | 6 dcs (2 new) | 6 dc (2 new) | New vocalizers played twice as frequently | PrelesionSet1 |
|      |DCvsDC_6v6_d2_S2    | 6 dcs | 6 dcs | All vocalizers played at equal frequency | PrelesionSet1 |
|Week 5|DCvsDC_12v12     | 12 dcs | 12 dcs | Combined stimuli from DCvsDC_6v6_d2 and DCvsDC_6v6_d2_S2 | PrelesionSet1 |
|      |SovsSo_16v16     | 16 songs | 16 songs | Combined stimuli from SovsSo_8v8_d2 and SovsSo_8v8_d2_S2 | PrelesionSet1 |
|Week 6|AllvsAll_4v4     | 2 songs + 2 dcs | 2 songs + 2 dcs | All vocalizers previously learned in earlier sets (refresher set) | PrelesionSet1 |
|      |AllvsAll_28v28     | 16 songs + 12 dcs | 16 songs + 12 dcs | Combined stimuli from DCvsDC_12v12 and SovsSo_16v16 | PrelesionSet1 |

### Subject Group

Subject group column of TrialData.csv indicates the experimental group the subject belongs to.

|Subject Group| Description  |
|---|---|
|EPHYS | Subjects used in ephys experiments and did not undergo lesion/sham lesion operation. Follows tests outlined in (b) above. |
|NCM | Subjects who underwent focal lesions of NCM via injections of NMA. Follows tests in (a) above |
|CM | Subjects who underwent focal lesions of CM via injections of NMA. Follows tests in (a) above |
|DEAD | Subjects who died during lesion operation, thus not reaching tests in (a) |
|CTRL | Subjects who underwent sham lesions of NCM or CM via injections of saline. Follows tests in (a) |


### Stimulus Files

Stimulus .wav files are saved in a zip folder. Paths in the Trial Data reference directories relative to the top level of the `stimuli` folder.

Filenames are formatted with the following convention:
`CALLTYPE_Stim_RENDITION_VOCALIZER_RENDITIONID_norm[_COPY].wav`
