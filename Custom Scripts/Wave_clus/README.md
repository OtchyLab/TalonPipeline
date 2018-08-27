# Wave_clus: Basic Use Instructions
The Wave_clus package is used to analyze electrophysiology data and try to identify the firing of specific neurons based on the shape of neural activity. Each neuron has a very specific and consistent shape when it fires, so the spikes can be sorted according to their shapes.

## Interface Description
![wave_clus.png](https://user-images.githubusercontent.com/18174572/41366213-64ed9042-6f09-11e8-80d9-dab335bdf569.png)
* **File Type Selection Dropdown**: Use this menu to select which type of file you are analyzing.  
* **Status Box**: The black box in the upper right of the GUI provides updates about the status of the file selection.
* **Load Filelist**: This button opens a finder window in which you select the folder of data you would like to analyze. 
* **Set Parameters**: This is used to adjust the parameters for the threshold used to find spikes in the electrophysiology data. The three most important parameters for setting this threshold are:
    * par.sr: this sets the sampling rate (constant at 44150 for all of our purposes)
    * par.stdmin: sets the minimum amplitude a spike much reach in order to be counted as a spike
    * par.detection: this can be set as positive, negative, or both; sets on which side the threshold will be placed 
* **Channel Selection Dropdown**: Select the E-phys data channel you would like to analyze using this menu.
* **Preprocessing**: Pre-processing opens a new window in which you can alter the filtering on the electrophysiology data, update applies these new filters to the data. To apply the new settings, check the "PrePro" box.
* **Load Data**: Load the selected files into the application
* **Plot Polytrode**: Not currently functional.
* **Plot All/Plot Average**: Change what is plotted in the large graph directly below the buttons
* **CMS**: Apply common mode subtraction to E-phys data
* **Save Clusters**: Save current clusters to a spike times file
* **Save Screen**: Save an image of the current GUI screen
* **Undo**: Undo previous action
* **Top Graph**: This graph shows either all of the electrophysiology voltage graphs for the selected files and channel or an average of them. It also displays the threshold line that is used to identify spikes. 
* **First Row of Graphs**: 
    * The graph on the far left displays all of the spikes from the loaded files together (colorcoded by cluster) 
        * Plot all projections: this button opens up a graph that breaks down the clusters and tries to find the optimum classifications for the recordings
    * Cluster 1 - Cluster 3: These graphs display all of the spikes in each of the three largest clusters
    * Clustor 0: Displays all of the spikes that were not grouped into a specific cluster
* **Second Row of Graphs**: 
    * Temperature vs. Cluster Size: The x-axis of this graph displays the "temperature" setting on the clustering algorithm, which can be thought of as the amount of variation required from spikes for them to be considered "different". If the temperature is high, it is harder for a shape to be part of a cluster and if it is low it is easy for it to be part of a cluster. The y-axis is the cluster size, which is the minimum number of spikes that a cluster can have. Both of these settings can be adjusted using the "change temperature" button and clicking a spot on the graph.
    * ISI Histograms: These graphs display histograms of the inter-spike intervals for the clusters directly above them. The inter-spike intervals give a good way to judge whether clusters are accurate. Because of the refractory period typical of neurons after they fire, it us unlikely that the same neuron would fire twice within 3 or fewer milleseconds. So, if a cluster's ISI histogram is heavily skewed right, it is unlikely that it is a "good" or accurate cluster. The number of intervals less than 3 milleseconds is displayed at the top of each of the graphs.
        * Accept/Reject: Select "reject" to remove a cluster and put all of its contents in cluster 0
        * Max: Adjust the x-axis scale on the histogram
        * Step: Adjust the bin size of the histogram
* **Force**: Forces as many of the unclassified spikes into clusters as possible
* **Merge**: "Fix" two or more clusters and click this button to merge them together

## Instructions
1. Select file type (you will typically use .DAT files) 
2. Load file list by pressing the "Load Filelist" button and selecting the desired file using in the finder window that appears
3. Select the E-phys data channel you wish to analyze
4. Load data (click one or multiple files)
5. The program will do its best to organize the spikes into clusters for you, but it is your job to refine the work it does until you are satisfied with it; edit the clusters by modifying the temperature or cluster size, merging clusters, rejecting clusters you believe are incorrect, and forcing unclustered spikes into clusters
6. When you are finished spike sorting, Hit the "save clusters" button in the top right of the GUI screen; this will mark the channel whose cluster you saved in the file view and create a new spike times file (format: times_birdname_universaltimestamp_YY_MM_DD_HH_MM_SS.NUMBERchan[#].mat)

