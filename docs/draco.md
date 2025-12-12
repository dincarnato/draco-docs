DRACO uses a sliding window approach to analyze short reads tiling long target transcripts. The default size of the window is 100 nt, slid with an offset equal to 1% the window's size.
<br/><br/>
![DRACO Workflow](http://www.incarnatolab.com/images/datasets/DRACO_Morandi_2021.png)
<br/><br/>
For each window, the following analysis steps are performed:<br/>

1. Reads covering the entire window are used to build a graph. In the graph, each base is a vertex; when two bases are observed to co-mutate within the same read, they are connected by an edge in the graph. Each edge, is weighted according to the co-mutation frequency of the two bases (or essentialy, the number of reads in which both bases are simultaneously mutated).

2. From the adjacency matrix of the graph, the normalized Laplacian matrix is calculated.

3. Spectral clustering is applied to this matrix, to identify the number of co-existing structural conformations (clusters) for the window being analyzed.

4. Graph-Cut is then used to weight each vertex in the graph according to its affinity to each conformation.

5. The window is slid by the chosen offset, and steps 1-4 are repeated, until the entire transcript has been analyzed.

6. Consecutive windows forming the same number of conformations are grouped into *window sets* and merged.

7. Reads are assigned to each *window set*

8. For each *window set*, the relative conformation stoichiometries, and the individual reactivity profiles are reconstructed. Window sets found to form different numbers of conformations, are reported separately.

When analyzing a transcript, DRACO keeps in memory all the reads mapping to that transcript. At extreme sequencing depths (&gt;500,000X), the memory consumption can become prohibitive on some system. If your system does not have sufficient RAM, reads can be randomly downsampled along the transcript, to achieve a lower mean coverage. Usually, a coverage of __10,000-20,000X__ is sufficient for DRACO to efficiently deconvolve the underlying structural heterogeneity.<br/>
<br/>
For additional information concerning downsampling, please consult the documentation for the [__rf-count__](https://rnaframework-docs.readthedocs.io/en/latest/rf-count/#downsampling-reads-for-analysis-with-draco) module of the __RNA Framework__.

<br/>

# Usage
To list the required parameters, simply type:

```bash
$ draco --help
```

Parameter         | Type | Description
----------------: | :--: |:------------
__--mm__ | string | Input [mutation map (MM)](https://rnaframework.readthedocs.io/en/latest/rf-count/#mm-mutation-map-format) file<br/>__Note:__ to pass multiple replicates, simply call the parameter multiple times (e.g., `--mm rep1.mm --mm rep2.mm) 
__--output__ | string | Output JSON file (Default: __draco_deconvolved.json__)
__--processors__ | int | Number of processors to use (Default: __0__)<br/>__Note:__ when set to 0, all available processors will be used
__--whitelist__ | string | A whitelist file, containing the IDs of the transcripts to be analyzed, one per raw (Default: __all transcripts will be analyzed__)
__--shape__ | | Enables spectral analysis on all four bases (Default: __only A/C bases__)<br/>__Note:__ this feature is highly experimental
__--outputRawNClusters__ | string | Outputs a BED-like file containing, for each transcript, the raw number of clusters detected for each window
__--log-level__ | string | Specifies the log level (allowed values: __trace__, __debug__, __info__, __warn__, __error__)
 | | __Mutation filtering options__
__--minBaseCoverage__ | int | Minimum coverage per base (Default: __100__)
__--minBaseMutations__ | int | Minimum mutations per base (Default: __2__)
__--minMutationFreq__ | float | Bases below this mutation frequency will be discarded as noise (Default: __0.005__)
__--minReadMutations__ | float | Reads with fewer than this number of mutations will be discarded as non-informative (Default: __2__)
 | | __Spectral deconvolution options__
__--maxClusters__ | int | Maximum allowed number of clusters/conformations (Default: __unlimited__)
__--minFilteredReads__ | int | Minimum number of reads (post-filtering) to perform spectral analysis (Default: __100__)
__--minPermutations__ | int | Minimum number of permutations performed to build the null model (Default: __10__)
__--maxPermutations__ | int | Maximum number of permutations performed to build the null model (Default: __100__)
__--ignoreFirstEigengap__ | | The first eigengap is ignored
__--eigengapCumRelThresh__ | float | Minimum relative difference between the eigengap and the null model, as a fraction of the cumulative difference between the previous eigengaps and their respective null models (Default: __0.1__)
__--eigengapAbsThresh__ | float | Minimum absolute difference between the eigengap and its null model (Default: __0.03__)
__--alpha__ | float | Below this p-value, the null hypothesis is rejected and the eigengap is marked as informative (Default: __0.01__)
__--beta__ | float | Above this p-value, the alternative hypothesis is rejected and the eigengap is marked as non-informative (Default: __0.2__)
__--minNullStdev__ | float | Minimum standard deviation for the null model (Default: __0.025__)<br/>__Note:__ when this threshold is not met, ``--extraPermutations`` additional permutations will be performed
__--extraPermutations__ | int | Additional permutations to perform when the standard deviation of the null model is &lt; ``--minNullStdev`` (Default: __50__)
__--minWinBases__ | int | Minimum number of bases in window (post-filtering) to perform spectral analysis (Default: __10__)
__--lookaheadEigengaps__ | int | Number of eigengaps to look ahead after a non-informative eigengap is encountered (Default: __0__)
__--saveEigengapData__ | | Saves eigengap data for plotting (Default: __off__)
__--eigengapDataOut__ | string | Eigengap data output folder (Default: __./eigengap_data__)
 | | __Graph-Cut options__
__--minClusterFraction__ | float | Minimum fraction of reads assigned to each cluster/conformation (Default: __0.005__)<br/>__Note:__ if this threshold is not met, the number of clusters is automatically decreased
__--softClusteringInits__ | float | Number of iterations for the initialization process of the graph-cut.<br/>__Note:__ The initialization with the lowest score is picked (Default: __500__)
__--softClusteringIters__ | float | Number of iterations performed on graph-cut<br/>__Note:__ The cut with the lowest score is picked (Default: __50__)
__--softClusteringWeightModule__ | float | The module of the weight that is used to change the cluster weights in order to find the lowest score (Default: __0.005__)
 | | __Windowed analysis options__
__--winLen__ | float | Length of the window. If this value is comprised between 0 and 1, it is interpreted as a fraction of the median read length. If this value is &gt; 1, it is interpreted as the absolute length of the window (Default: __100__)<br/>__Note:__ this parameter and ``--winLenFracRnaLen`` are mutually exclusive
__--winLenFracRnaLen__ | float | Length of the window as fraction of the length of the RNA<br/>__Note:__ this parameter and ``--winLen`` are mutually exclusive
__--winOffset__ | float | Window sliding offset. If this value is comprised between 0 and 1, it is interpreted as a fraction of the window's length. If this value is &ge; 1 , it is interpreted as the number of bases to slide the window by (Default: __0.01__)
__--maxIgnoreWins__ | int | Maximum number of internal windows with a different number of clusters to ignore when merging two external sets of windows (Default: __0__)
__--minExtWins__ | int | Minimum number of external windows, having the same number of clusters, needed to trigger merging (Default: __6__)
__--nonInformativeToSurround__ | int | Non-informative windows (windows with 0 detected clusters) are set to the same number of clusters of surrounding windows (Default: __false__)
__--allNonInformativeToOne__ | int | If *all* windows in the transcript are non-informative (0 clusters), the number of clusters is forced to 1 (Default: __off__)
__--reportNonInformative__ | int | Reports also non-informative windows in the output JSON file
__--assignmentsDumpDir__ | string | When specified, a folder is generated containing an MM file for each set of merged windows, containing a dump of the reads assigned to each conformation
__--skipAmbiguousAssignments__ | | When specified, if the best assignment score for a read is equal across multiple clusters, the read is discarded
__--minWindowsOverlap__ | | Minimum overlap between non-contiguous windows with the same number of conformations in order to merge them in the same set. The value must be comprised between 0 and 1, and it is interpreted as a fraction of the window length (Default: __1__)
<br/>
## Understanding the algorithm

!!! alert "Note"
    Since v1.3 it is possible to follow step-by-step the analysis perfromed by DRACO by setting the parameter `--log-level trace`

While it is advisable for most users to run DRACO with its default parameters, as these are the results of a careful and thorough calibration, it might be useful to adjust the analysis on a case-by-case basis.<br/>
<br/>
DRACO analysis is performed in sliding windows. The size of the window (``--winLen``) is by default 100, and the sliding offset (``--winOffset``) is by default 1% of the length of the window. We observed that a window size of 100 is optimal in most cases, as it allows the identification of both small and large structurally-heterogeneous regions. This length is recommended even if longer reads have been sequenced, as, if the structurally-heterogeneous region in the transcript is much smaller than the read length (for example, 300 bp-long reads and a structurally-heterogeneous region of 90 nt), information might get lost in the spectral deconvolution.<br/>
<br/>
The __spectral deconvolution__ represents the most critical step of the algorithm, as this is responsible for determining the number of conformations populated by each window. During this step, the *eigengaps*, representing the distance between consecutive *eigenvalues*, are compared to a *null model*, built by performing random permutations over the original data matrix. The distribution of the null model is well approximated by a Weibull distribution. The number of permutations performed varies for each eigengap, and it is comprised between ``--minPermutations`` and ``--maxPermutations``. Permutations are performed as long as the standard deviation (*&sigma;*) of the null model is &lt; ``--minNullStdev``, increasing by ``--extraPermutations``.<br/>
<br/>
Starting with the first eigengap, DRACO compares each eigengap to its respective null model. Aim of this comparison is to determine which eigengaps can be considered to be *informative*, as this directly translates into the number of conformations/clusters populated by the analyzed window. For instance, if the first two eigengaps are informative, two conformations (clusters) are present; if the first three eigengaps are informative, three conformations are present, and so on. 
<br/><br/>
![DRACO eigengaps](http://www.incarnatolab.com/images/docs/draco/DRACO_eigengaps.png)
<br/><br/>
In order to be considered informative, each eigengap must fulfill a number of criteria:

1. The distance (*d*) between the mean (*&mu;*) of the null model and the eigengap must be &ge; ``--eigengapAbsThresh``.

2. The distance between the *&mu;* and the eigengap must be greater than ``--egengapCumRelThreshold`` &times; the cumulative sum of the previous eigengaps. For instance, *d<sub>4</sub>*, the distance between the 4<sup>th</sup> eigengap and its null model, must be &ge; ``--egengapCumRelThreshold`` &times; (*d<sub>1</sub>* + *d<sub>2</sub>* + *d<sub>3</sub>*). This, of course, does not apply to the first eigengap

3. A *t-test* is used to assess whether the distance between the eigengap and the null model is significant. If the p-value is &le; ``--alpha``, the null hypothesis is rejected, and the eigengap is marked as __informative__. If the p-value is &ge; ``--beta``, the alternative hypothesis is rejected, and the eigengap is marked as __non-informative__. If the p-value is comprised between these two thresholds, neither of the two hypotheses can be rejected, so additional permutations are performed. If, after having performed ``--maxPermutations`` permutations the p-value is still comprised between the two thresholds, then the alternative hypothesis is rejected and the eigengap is marked as __non-informative__.

The analysis terminates when a non-informative eigengap is encountered as, normally, the subsequent eigengaps are non-informative as well. In some particular cases, however, we have eigengaps to behave unexpectedly. For instance, at times the first eigengap might overlap with its null model, which would cause DRACO to report 0 clusters and terminate the analysis for that window. This can be avoided by enabling the `--ignoreFirstEigengap` parameter. In other cases, although an eigengap is marked as non-informative, one or more eigengaps immediately downstream of it are again informative. To account for these situations, whenever a non-informative eigengap is encountered, DRACO can perform a *lookahead* evaluation of the ``--lookaheadEigengaps`` subsequent eigengaps. If one of these is marked as informative, then the analysis of the eigengaps is allowed to continue. This analysis is disabled by default.<br/>
<br/>
Unless the parameter `--outputRawNClusters` had been enabled (which would cause DRACO to only report the raw number of identified clusters/conformations per window), once the number of conformations has been determined, the algorithm uses a __normalized Graph-Cut__ approach to weight each base in the window, according to its affinity to each of the different conformations. This approach tries to find the best way to partition the previously built graph. It is composed of two steps:

1. __Initialization__. During this step, each base of each cluster is assigned a random weight. The step is repeated ``--softClusteringInits`` times and the solution with the lowest score is picked.
2. __Refinement__. During this step, the weights of the bases are dynamically altered by ``--softClusteringWeightModule`` in an attempt to further minimize the score. This step is repeated ``--softClusteringIters`` and the solution with the lowest score is picked.

!!! alert "Note"
    In the original implementation of DRACO, this step was carried out only once. While being significantly faster, the downside was a higher risk of finding suboptimal graph partitioning solutions, corresponding to local score minima. Lowering the ``--softClusteringWeightModule`` and increasing the ``--softClusteringIters`` can slow down the analysis, but it significantly increases the likelihood that the *true* optimal solution is identified, further ensuring reproducibility of the analysis results.


Final step of the analysis involves merging consecutive windows found to form the same number of conformations, into *window sets*. Let's however consider the following case:
<br/><br/>
![DRACO eigengaps](http://www.incarnatolab.com/images/docs/draco/DRACO_mergewindows_1.png)
<br/><br/>
In this situation, two sets of windows found to form two conformations, are interrupted by one window forming one conformation. By default (see figure above, left panel), DRACO would report three window sets. There are two ways to handle these cases:

1. The first possibility is to make DRACO ignore at most ``--maxIgnoreWins`` internal windows forming a discordant number of conformations. For the merging to occur, the total number of external windows on either sides of the discordant internal windows must be &ge; ``--minExtWins`` (see figure above, right panel).

2. The second possibility (preferred, no information loss) is to make DRACO report every *window set* it did identify, but allowing non-contiguous *window sets* forming the same number of conformations to be merged, provided that the overlap between the last window of one set and the first window of the following set is at least `--minWindowsOverlap` &times; the window length (e.g., ``--minWindowsOverlap 0.5`` requires a minimum overlap of 50% of the window length; see figure below).
<br/><br/>
Following window merging, reads are assigned to each conformation. If the fraction of reads assigned to a conformation is &lt; ``--minClusterFraction`` &times; the number of conformations for that *window set* is decreased by one, and the step is repeated.
<br/><br/>
### Handling of replicates
Since v1.3, DRACO can handle replicate experiments. Multiple MM files can be passed by simply calling the `--mm` parameter multiple times. 

During the spectral analysis DRACO will separately determine the number of conformations populated by corresponding windows across replicates, and the number of conformations for a given window will be set as the median of the number of conformations identified across all replicates. For instance, in an experiment with 3 replicates, if two identified 1 conformation and one identified 2 conformations, the number of conformations for the window will be set to 1. Similarly, if the three replicates respectively identified 1, 2 and 3 conformations, the number of conformations for the window will be set to 2.

During the Graph-Cut part of the algorithm, graph partitioning and score minimization will be simultaneously performed on the graphs of each replicate, hence avoiding experiment-specific local minima and ensuring optimal parititioning.

<br/><br/>
## Output JSON files

DRACO produces a JSON file, with the following structure:

```json
{
   "filenames":["sample.mm"],
   "transcripts":[
      {
         "id":"RNA_1",
         "sequence":"TCTATTCTACATTGATAGAACAACTCT...TCTGTCTATCACCGCTAGAGCACTCGGTGATTGCA",
         "nReads":[20000],
         "windows":[
            [{
               "start":0,
               "end":294,
               "stoichiometries":[ 0.388,0.304,0.308 ],
               "counts":[
                  [ 0,1,0,23,0,0,33,0,54,0,5,0,0,...,301,0,334,0,3,1,0,0,0,0,0,0 ],
                  [ 0,6,0,0,0,0,0,0,0,0,36,0,0,2,...,0,0,0,0,0,0,134,0,86,0,0,62,63 ],
                  [ 0,1,0,17,0,0,57,0,0,0,75,0,11,...,181,0,0,0,0,0,0,0,0,118,0,0,0,0 ]
               ],
               "coverage":[
                  [ 40,69,121,157,194,245,297,337,378,416,465,...,521,482,462,430,404,359,311,287,250,211 ],
                  [ 37,73,108,135,163,198,238,279,317,340,369,...,534,512,481,463,440,414,380,240,204,173 ],
                  [ 35,66,107,141,168,206,250,292,337,380,415,...,428,394,364,336,309,271,240,199,181,111 ]
               ],
               "weights":[
                  [ 0.000,1.000,0.000,0.086,0.000,0.055,0.000,0.000,...,0.110,0.110,0.000,0.000,0.000,0.000,0.000,0.000 ],
                  [ 0.039,0.035,0.000,0.000,0.000,0.000,0.000,0.000,...,0.000,0.922,0.000,0.000,0.035,0.925,0.000,0.000 ],
                  [ 0.165,0.000,0.000,0.855,0.000,0.020,0.000,0.855,...,0.000,0.000,0.000,0.055,0.000,0.055,0.000,0.835 ]
               ],
               "preCoverage":[ 23979,24261,24506,24751,24996,25257,25539,...,25813,26078,26340,26616,26889,27182,27441,27719,28009 ]
            }]
         ]
      }
   ]
}
```
The __transcripts__ array contains all the transcripts processed by DRACO. Each transcript includes the following keys:<br/>

Key    | Description
-----: | :----------
__id__ | Transcript's ID
__sequence__ | Transcript's sequence
__nReads__ | An array containing the total number of reads mapping to the transcript in each replicate
__windows__ | An array of arrays containing the deconvolution results for each window set across each replicate being analyzed

<br/>
The __windows__ array is further structured as follows:<br/>

Key    | Description
-----: | :----------
__start__ | Start position of the window (0-based)
__end__ | End position of the window (1-based)
__stoichiometries__ | An array containing the relative abundances for the reconstructed conformations
__counts__ | An array of arrays, each containing the number of mutations per base for each of the reconstructed conformations
__coverage__ | An array of arrays, each containing the coverage per base for each of the reconstructed conformations, calculated using all the assigned reads post-filtering
__weights__ | An array of arrays, each containing the weight per base from the Graph-Cut analysis
__preCoverage__ | An array containing the total coverage per base along the window, calculated using only the reads pre-filtering, covering the entire window, that have been used to perform the spectral analysis (therefore, it can be &lt; than the sum of the __coverage__ arrays)

<br/>
Since v1.3, in case multiple replicates are being analyzed, an additional __clusterConfidence__ array is present right after the __windows__ array, which is structured as it follows:<br/>

Key    | Description
-----: | :----------
__nClusters__ | Number of detected clusters/conformations
__confidence__ | A value between 0 and 1 representing the fraction of replicates for which the spectral analysis found the above number of conformations (e.g., if __nClusters__ is 2 and out of two replicates one was found to form 2 conformations and the other 3 conformations, __confidence__ will be 0.5)
__start__ | Start position of the window (0-based)
__end__ | End position of the window (1-based)

<br/><br/>
## Post-processing of DRACO output
DRACO JSON output files can be somehow hard to interpret and process. To facilitate this operation, we developed the ``rf-json2rc`` module, that has been included in the [__RNA Framework__](https://rnaframework-docs.readthedocs.io/en/latest/). The module allows converting the reactivity profiles deconvolved by DRACO into RNA Count (RC) files, that can be subsequently normalized with the ``rf-norm`` module, and then fed into the ``rf-fold`` module for secondary structure modeling. For additional details, please refer to the ``rf-json2rc`` [docs](https://rnaframework-docs.readthedocs.io/en/latest/rf-json2rc).