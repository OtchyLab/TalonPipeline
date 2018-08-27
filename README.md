# Talon Pipeline

This is a collection of the software used in the Otchy Lab to process and analyze bird songs

## Important/Recurring File Types

This pipeline is designed to work with three main file types:
1. **.WAV**: .WAV is the standard format for uncompressed audio files.
2. **.DAT**: these files are more complex, containing audio and coinciding electrophysiology readings (allowing observation of bird song and brain activity in tandem).
3. **.STM**: files containing audio and information about the type of stimulation applied to the animal.

**.DAT Format**: .DAT files have either five or seven independent channels of data. When a .DAT file is opened, the first two pieces of information that are read in are the sample rate of the data (this will be consistent across all channels, typically 44150Hz) and the number of channels. After reading this, it will move on to the actual data channels themselves. The first channel is always the audio of the birdsong that was recorded. All of the subsequent channels are electrophysiology readings from different electrodes in the bird's brain. The file is formatted such that one full second of each of these channels is read in first (so if there were 7 channels, there would be seven seconds in which each individual second was the data from a different channel). This results in 44150 readings for each individual channel. After this, the remaining data in the file is "multiplexed," which simply means they are interleaved so there is one reading from each channel and then it switches to the next one. 

**.STM Format**: .STM files follow a very similar format to .DAT files. The only difference is that instead of the channels after the audio channel being electrophysiology readings, they contain data for current injection, voltage command, target flag, and STM flag (in that order). 

## Typical Workflow
The analysis of birdsong data typically follows a standard sequence of steps that convert large number of audio and electrophysiology data files into a format from which interesting conclusions can be drawn. This sequence is as follows:

(Note: each process in this workflow can be launched by running the MATLAB file whose name matches that of the folder it is stored in. For example, to use the SongBlaster package, you would run SongBlaster.m in the MATLAB)

For more detailed instructions on how to use the following tools, see the documentation within the folders.

### **1. Scrubbing:** SongBlaster
The first step in the analysis of the raw data files (described above) is throwing out "garbage" files. The song-triggered recording software is designed to record whenever it hears bird songs, but sometimes it picks up background noise that is not useful in song analysis. So, the audio files must be manually scrubbed to remove all of these unwanted files. SongBlaster is a GUI that allows you to browse and delete unwanted .WAV, .DAT, and .STM files.

**Input**: Folder of all raw audio files from song-triggered recording software

**Output**: Folder of audio files containing only relevant data (bird song recordings)

### **2. Spike Sorting:** Wave_clus
The Wave_clus package is used to sort spikes in electrophysiology data and try to identify individual neurons based on the shape of the neural activity. This step is only necessary for .DAT files since .WAV files do not have any neural data. The application works on the principle that the firing of different neurons have consistent and unique shapes.  

**Input**: Folder of .DAT files.

**Output**: Spike time files.

*Note*: This step is usually done in tandem with the annotation step so the spike time files can be evaluated for accuracy in TweetVision.

### **3. Annotating:** TweetVision and TweetVisionLite
These two packages are used to identify different recurring syllables in the audio files and take note of _where_ they occur. They have tools (both manual and automated) for identifying the beginnings and ends of syllables and classifying them as types. For the built in tools to be effective, the recordings must be _high quality_ and the syllables must be _sufficiently unique_.

(Note: this is the most critical step for producing high quality analysis; be conservative and careful by marking a syllable as 'unknown' if you are unsure of how to classify it or it looks like it has noise superimposed on it)

#### TweetVisionLite: Use for .WAV files
TweetVisionLite does not have any functionality for working with .DAT files, so it should be used for .WAV files. This package has many convenient tools for autiomated syllable identification and labeling that makes annotating audio files much easier. 

**Input**: Folder of .WAV files.

**Output**: 
1. Annotation File: a list of syllable types and start/stop times. There will be one annotation file per bird per day. Naming conventions are as follows: Birdname_YYMMDD_annotation.mat
2. Syllable Pattern File: A graphical list of all the syllable types for a bird (typically excludes calls and low probability syllables). There is typically one of these files per bird. Naming conventions: Birdname_YYMMDD_SylPattern.mat
3. Neural Network File: if you use the network classifier for automating your annotations, you can save the neural network for later use (it will work well as long as the recording conditions don't change too much!). Naming conventions: Birdname_YYMMDD_network.mat 

#### TweetVision: Use for .DAT files
TweetVision has similar functionality to TweetVisionLite, but it is designed to work with .DAT files and has features to analyze electrophysiology data as well as audio data. It is missing several of the automated features that TweetVisionLite has. Use it for .DAT files. The application is able to do annotations on audio files as well as create and edit cell files that contain information about when specific neurons fire during the bird songs. 

**Input**: 
1. Folder of .DAT files.
2. Any spike times files corresponding to the .DAT files

**Output**: 
1. Annotation File: A list of syllable types and start/stop times. There will be one annotation file per bird per day. Naming conventions are as follows: Birdname_YYMMDD_annotation.mat
2. Syllable Pattern File: A graphical list of all the syllable types for a bird (typically excludes calls and low probability syllables). There is typically one of these files per bird. Naming conventions: Birdname_YYMMDD_SylPattern.mat
3. Cell File: Cell files contain information about how different spikes in electrophysiology data, assumed to be action potentials of different neurons, are clustered. The names of these files adhere to the following conventions: Birdname_YYYY_MM_DD_cells.mat (if there are multiple cell files for one bird, add a number to the file name after cells).

### **4. Aligning:** StretchEm and StretchEmLite
These two packages aid in the temporal aligning of syllables and syllable sequences identified using the TweetVision functions.  They provides a better way to find the beginnings and endings of syllables, and their goal is to find corresponding points of different songs that match up. 

#### StretchEmLite: Use for .WAV files
Like TweetVisionLite, StretchEmLite lacks the functionality to deal with electrophysiology data. As such, it should be used for pure audio data only. 

**Input**: 
1. Annotation file: generated during the annotation stage of analysis
2. Syllable pattern file: also from the annotation stage
3. Folder containing song recordings (.WAV files)

**Output**: 
1. Talon file: contains alignment paths defining the temporal structural of each rendition and raw recording data for each song rendition. There is one talon file per bird per day, named according to the following convention: Birdname__YYMMDD_dataset.mat.
2. Song Template file: includes template spectrogram and identified sullable start and stop times. There should be one of these files per bird, named according to the convention: Birdname_YYMMDD_template[#s].mat, where [#s] identifies the syllables in the motif.

#### StretchEm: Use for .DAT files
StretchEm is almost exactly the same as StretchEmLite with just a few added features that allow you to deal with electrophysiology data. Any analysis of .DAT files should be done using StretchEm. 

**Input**:
1. Annotation file: generated during the annotation stage of analysis
2. Syllable pattern file: also from annotation stage
3. Folder containing song/E-phys recordings (.DAT files) 

**Output**: 
1. Talon file: contains alignment paths defining the temporal structural of each rendition and raw recording data for each song rendition. There is one talon file per bird per day, named according to the following convention: Birdname__YYMMDD_dataset.mat.
2. Song Template file: includes template spectrogram and identified sullable start and stop times. There should be one of these files per bird, named according to the convention: Birdname_YYMMDD_template[#s].mat, where [#s] identifies the syllables in the motif.

### 5. **Analyzing:** Metermeter2
Metermeter2 is concerned with exploring the way that bird songs change over time, both temporally and acoustically. This package is a major work in progress as it is rather counterintuitive and difficult to use, so it is likely to be changed significantly in the coming months. Nonetheless, we can document the main functions and goals of the application which will not change.

The use of Metermeter2 is much more open to user preference/objective than any of the other steps in analysis. It provides a variety of tools for assessing the variability of bird songs over time and you can choose the one that is most interesting to you. The main tools that are provided are variance decomposition, in which the variability in the bird songs is fit to a mathematical model of how sequences are structured in the brain and separates all the differences in the recordings into four parameters, spectral analysis, which is used  to analyze the acoustic variance in bird song, and temporal analysis, which deals with differences in the temporal structure of the motifs. You can use this application to examine any of these types of variance and export data and graphs. 

**Input**: Several talon files from different days for the same bird and the same syllable motif.

**Output**: Whatever you want! Export whichever data/graph you deem useful to your project.

### 6. **Electrophysiology Analysis**: Talon
This package is used to align electrophysiology data for song motifs and obtain statistics about the neural spikes and obtain spike train statistics.

**Input**: 
1. One or more cell files. 
2. One or more talon files

**Output**: Spike train statistics for each cell file stored in one data structure.






