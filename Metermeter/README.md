# Metermeter2: Basic Use Instructions
Metermeter2 is a tool that can be used to analyze the way that bird songs change over time. The GUI is messy and needs work, so expect changes in the near future. The basic functionality should remain the same. 

## Interface Description
![metermeter2gui.png](https://user-images.githubusercontent.com/18174572/41473080-4724e97e-7086-11e8-89f4-a646c99b16ad.png)
* **Status Box**: This black box in the upper left of the GUI provides updates on the state of adding, clearing, and loading data sets. 
* **Help Button**: This button opens a second figure with keyboard shortcuts and a legend for the plots in the GUI. 
* **Input/Output**: The Input/Output panel is where all of the file control for this application happens. 
    * Add Dataset: Add the talon files that you are interested in analyzing; they must be the same motif from the same bird
    * Clear Datasets: Clear the currently loaded datasets 
    * Load Datasets: Load the data sets you added into memory; this may take a while depending on how many you chose to load
    * Temporal Only: If you are only interested in temporal analysis, select this box and the raw audio data will not be loaded; this will speed up the loading/analysis processes significantly
    * Exclude Flagged: do not load files that were flagged in TweetVision for improper alignment
* **Interval Selection**: In this panel, you select which times intervals you are interested in studying
    * Parse Time: Here you put the start and end times of the day that you want to look at; for example, if you are interested in 8am to 10am, you would enter "8,10" in this box; if you want to do multiple time intervals, separate them with a semicolon
    * Max # Motifs: Number of motifs you care about in each interval; it will analyze the first x motifs in each interval; enter -1 to analyze all of the motifs in the selected intervals
    * Zero Day (MM/DD/YY): Set the first day on the axis of the large graph in the GUI
    * Normalize Output: Normalize data based on day(s) selected
* **Variance Decomposition**: In this panel, the temporal data is fit to a model of how sequences are structured in the brain. It takes the temporal structure amd fits it to a model of how the intervals should be related, extracting four parameters (described below). 
    * Variability Decomp: Initiate variability decomposition for data
    * Dropdown Menu: Choose which parameter you would like to view on the graph
        * 1st Global: Uniform stretch of entire motif
        * 2nd Global: Stretch on whole song as it unfolds over the course of the day (track circadian rhythm of the bird)
        * Independent: Models how much each interval (syllable or gap) can vary on its own without affecting the ones next to it
        * Jitter: Boundaries that shift because of measurement error
    * Show Ints: Plots individual intervals (means/standard deviations for all selected comparisons)
    * Mean Syls: show mean similarity across all syllables +/- standard deviation
    * Mean Gaps: show mean similarity across all gaps
* **Spectral Analysis**: The Spectral Analysis panel is used to explore the acoustic changes in bird songs over time. 
    * Syllable Similarity: Looks at all pairs of syllables and makes comparisons between them at each millesecond; only compares within the same day and interval; tells you about the amount of variability within a single interval
    * Syllable Recovery: Compares all instances of a single syllable selected in the box to the right
    * Show Ints: Plots individual intervals (mean and standard deviation) for all selected comparisons
    * Mean Syls: Show mean similarity across all syllables +/- standard deviation
    * Mean RMS: Not functional
* **Temporal Analysis**: Analyze changes in the temporal structure of bird songs over time
    * Plot Simple: Plots mean +/- standard deviation of each syllable and gap selected
    * Plot All Rend: Creates a scatterplot of each individual data point (motif) in the intervals selected plotted over time
    * Plot Simple CV: Plots CV (standard deviation divided by mean, a measure of variance) for each interval
    * % Motif Length: Not functional
    * Syls and Gaps: Analyze syllables and gaps as separate entities
    * Syl+Gap: Group each syllable with the gap that follows it, analyze as a unit
    * Sum Intervals: Not functional (sum together and plot)
* Figure View: The graph is used to visually represent the analyses you are performing. 
    * Graph: Displays the most recent analysis performed
    * Export Figure: Exports figure into a separate window so it can be saved
* **Data Export**: The data export functions of this application need a major overhaul! Expect changes in the near future. 
    * Save Final Data: Saves time, interval means and standard deviations, data type, parse time, and number of motifs in a data structure which changes depending on the most recently run analysis
    * Save Parsed Data: Saves almost exact same data as "save raw data" in a different data struct... for some reason
    * Save Raw Data: DELETE
    * Write Syls 2 Disk: DELETE
* **Controls on the Right**: These buttons were used to test Metermeter during its development and are not important.
    * Test 1: No functionality
    * Test 2: No functionality
    * Motif Sim: Compare entire motifs
    * YIN Pitch: Measure pitch

## Instructions: 
1. Add talon files of the same bird and motif for multiple days.
2. Load datasets using the green button.
3. Select the desired interval(s) in the Interval Selection Panel.
4. The rest is up to you! Perform whichever analysis you deem useful to your research.
5. Export data using the "Save Final Data" button; this will save your most recent analysis. 
