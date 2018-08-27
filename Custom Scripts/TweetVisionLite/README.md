# TweetVisionLite: Basic Use Instructions
TweetVisionLite is used to annotate audio recordings, i.e., break them down into syllables and document which syllables occur when.
## Interface Description:
A view of the basic GUI:
![tweetvisionlitegui](https://user-images.githubusercontent.com/18174572/40993075-8f060074-68c6-11e8-841a-772aef6210c5.png)
An example syllable pattern:
![syllablepattern](https://user-images.githubusercontent.com/18174572/40992922-1aa074f8-68c6-11e8-8b71-9fad54787a59.png)
* **File Control**: The file control panel, located in the bottom left corner of the GUI, is used to load folders containing audio files into the program so they can be annotated.
    * Load Folder: this button opens a finder window that allows you to choose the folder you want to annotate
    * Previous/Next: let you scroll through the audio files in the loaded folder
    * List View: shows the names of all the files in the loaded folder; files that are already annotated turn green
    * Quarentine/Empty Quarentine: a now outdated tool (due to the existence of SongBlaster), these buttons were used to flag/unflag certain files for deletion
    * Delete: delete quarentined files
    * Export Spec: open spectrogram in a new figure window
    * Gain/High Pass/Low Pass: display settings for audio envelope and spectrogram
    * Load WAVs/Load DATs/Load STMs: allow you to choose which file types you want to load into the program
* **Graphical Components**: Located in the upper left, the graphical components include a folder display, audio envelope, and spectrogram. This is where the majority of the manual work that goes into file annotation occurs. 
    * Top Box: contains the address of the folder that the TweetVisionLite program loaded
    * Audio Envelope: a graphical representation of the amplitude (loudness) of the bird's calls; the box with the arrows at the upper left controls the aspect ratio of the y-axis; click and drag in the audio envelope to zoom in on the spectrogram below
    * Spectrogram: a visual representation of the frequency and amplitude of the bird's calls; the x-axis represents time, the y-axis represents frequency, and the color denotes amplitude (with red being the loudest); you will add new syllables by clicking and dragging on this view and then numbering your selections
    * Help: opens a new window displaying a list of keyboard shortcuts and the legends for plots
* **Annotation**:
    * New Annotation: creates a new annotation file; you must creat a new annotation for each folder you load before doing anything else; name according to conventions described in the use instructions
    * Add Syllable: press this button prior to clicking and dragging to select sections on the spectrogram to create syllables
    * Delete Syllable: click a selection on the spectrogram and then press this button in order to delete an unwanted syllable
    * Edge Finder: OBSOLETE
    * Delete Pause: joins two separate syllables
    * Load Annotation: load a previously generated annotation file to add to it
    * View Template: view a previously generated syllable pattern
    * Infused Drugs/Directed-Undirected: provide more information about the annotations you are creating
* **Annotation Progress**: this panel (located in the center on the bottom of the GUI) allows you to see how much you've annotated and how many of each pattern of syllables you have identified 
* **Segmenting**: The tools in this panel (which is located in the top right corner of the application) eliminate the need to manually segment each file into syllables; with the help of the panel this process can be completely automated.
    * Multi-Segment: auto-segments entire file
    * Frame Segment: auto-segments all syllables in the currently selected box in the audio envelope
    * Clear All: removes segmentation on spectrogram
    * Reset Params: resets all parameters to their default values
    * Power Env/Dynamic/Seg Guide: segmentation settings
    * Syllables/Bouts: different parameters that can be adjusted to improve the accuracy of the auto-segmenting functions
* **Automated Syllable ID**: This panel (in the bottom right) allows you to automate the process of identifying the syllables in each audio file. In order to implement this functionality, some syllables must be annotated manually and a neural network must be trained.
    * Train Network: train a network to identify syllables for you
    * Load Network: load a previously trained network for use in this folder; calls should be fairly consistent for the same bird, so once you train one neural network for a particular bird you should be able to reuse it for any recordings from that bird
    * Save network: allows you to save a trained network for later use (adhere to naming conventions described in use instructions)
    * ID Syllables: use the loaded neural network to classify all segmented syllables in the current file
    * Batch Process: automates segmentation AND syllable ID's across multiple flies
    * Overhang/SylDurWeighting: neural network settings
    * Decision Pnt: confidence threshold for neural network assigning syllables or "unknown"
    * Time: assigns time frame for batch process
## Instructions:
### Manual Use
1. Load folder using file control panel, select a folder containing .wav files, and click open
2. An "Input Experiment Parameters" window will appear; all you have to do is enter the bird's name
3. A warning will appear saying that you have not yet created an annotation file; this is fine
4. Click "new annotation", name according to the convention Birdname_YYMMDD_annotation.mat
5. Double click a file name in the file control panel. This will bring up the audio envelope and spectrogram for that file.
6. Click and drag on the audio envelope to zoom in
7. Click **add syllable**, click and drag around syllable in spectrogram
8. Click on newly created syllable and then number it
9. Add syllables to the syllable pattern by doubleclicking them
10. To save a syllable pattern, use "save template" button in the annotation panel
11. Get annotating! :)

* Note that annotated files turn green after you've worked with them to indicate that they have already been annotated! 

### Automated Use
First, follow steps 1-5 in the manual use instructions. 

6. Click and drag on the audio envelope to select a subset of the recording. Then click "frame segment" to auto-segment this subset into syllables. Alternatively, use "multi-segment" to auto-segment the entire file at once
7. If you have done enough manual IDing of different syllables, you can train a neural network to do the rest for you! Click the "Train Network" button in the Automated Syllable ID panel
8. Once the "Train Network" button turns green, click the "ID syllables" button to automatically classify each syllable in the selected file
----------------------------
The **Batch Process** can be used to automatically segment and annotate multiple files at once. Note that a network must be trained for this function to be used. Select the start and end times for the files that you want annotated and then press the button. 

**** TEMPLATE = SYLLABLE PATTERN ****

