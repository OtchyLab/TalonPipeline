# Talon: Basic Use Instructions
This application is used to align neural activity with song motifs and extract spike train statistics. 

## Interface Description
![talongui.png](https://user-images.githubusercontent.com/18174572/41473191-8ff31b8a-7086-11e8-9a33-64382b8e920d.png)
* **File Control**: This panel is used to load in all of the data needed in order to use this package.
    * Add Cell File: Add the cell files you wish to analyze
    * Clear Cells: Clear all currently added cell files
    * Add Dataset: Add the talon file(s) you want to analyze
    * Clear Datasets: Clear loaded talon file(s)
    * Sync Data: Match up talon file and cell file
    * Ex Flag'd: Exclude all flagged renditions in the talon file
    * Lag: Adjust the time difference between the neural activity and the audio data
* **Spotlight**: This panel is not currently functional but what follows is a description of what it *should* do.
    * Arrows: Scroll through the spike rasters in a new window
    * Trash: Delete currently selected spike raster
    * Single: If this is not selected, every time you move to a new spike raster a new window will pop up 
* **Record Selection/Sorting**: This panel is for filtering *which* files you want to use in your analysis. 
    * From Record/To Record: allow you to filter for specific times within the folder of recordings based on the universal time stamp; when you change this it will change the entries in From Time/To Time to the corresponding values
    * From Time/To Time: see above
    * Directed Status Dropdown/Drug Status Dropdown: provide more information about which recordings you would like to include
    * Sort Drugs/Sort Direct/Sort Chrono: Sorts recordings based on selected criteria
* **Alignment**: 
    * Dropdown Menu: Choose desired method of alignment
        * Local Linear: align by start and end of each syllable
        * Global Linear: align just based on start and end of the recordings
        * Align Start: only align to onset of each motif
        * Full DTW: uses dynamic time warping algorithm 
    * ALIGN!: Align files
* **Analysis**: 
    * Spike Train Stats: Displays a variety of statistics for each cell file:
        * Burst: Fraction of spikes that appear in bursts (patterns/chunks > 150Hz)
        * FR/100: Average # of spikes per second
        * Sparse: Measures how disperesed the spike activity is
        * Corr: Measures how well aligned the rows are
    * Corr Matrix: Compares each row with every other row to 
    * ISI Distribution: Interspike interval distribution (just like in Wave_Clus!)
    * Test 4: Creates a mini window with containing the graphical components on the left of the GUI
    * New Figure: Not functional
* **Graphical Components**: A visual representation of analysis performed in this application. 
    * Top graph: Alignment template
    * Large Middle Graph: Displays spike rasters for all cell files colorcoded by file
    * Peristimulus Time Histogram: Graph on the bottom left of the GUI; displays average "instantaneous" firing rate of each spike raster
    * Smooth: Not functional
    * PSTH Bin Size: Set the bin size (in milleseconds) of the PSTH graph
* **Export Stats**: This button exports the spike train statistics for each cell file and puts it into one data structure.

## Instructions
1. Add cell file(s) 
2. Add talon file from the corresponding day/bird
3. Sync data
4. Select time interval and sorting method in the Record Selection/Sorting panel 
5. Select alignment type in Alignment panel
6. Align!
7. Use the Analysis panel as necessary to explore your data
8. Export spike train stats with the "Export Stats" button  
