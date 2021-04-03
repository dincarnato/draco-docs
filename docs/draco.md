DRACO uses a sliding window approach to analyze short reads tiling long target transcripts. The default size of the window is equal to 90% the median length of the reads, slid with an offset equal to 5% the window's size.
<br/><br/>
![DRACO Workflow](http://www.incarnatolab.com/images/datasets/DRACO_Morandi_2021.png)
<br/><br/>
For each window, the following analysis steps are performed:<br/>

1. Reads covering the entire window are used to build a graph. In the graph, each base is a vertex; when two bases are observed to co-mutate within the same read, they are connected by an edge in the graph. Each edge, is weighted according to the co-mutation frequency of the two bases (or essentialy, the number of reads in which both bases are simultaneously mutated).

2. From the adjacency matrix of the graph, the normalized Laplacian matrix is calculated.

3. Spectral clustering is applied to this matrix, to identify the number of co-existing structural conformations (clusters) for the window being analyzed.

4. Graph-Cut is then used to weight each vertex in the graph according to its affinity to each conformation.

5. Reads are assigned to each window

5. The window is slid by the chosen offset, and steps 1-4 are repeated, until the entire transcript has been analyzed.

6. Consecutive windows forming the same number of conformations are grouped into *window sets* and merged.

7. For each *window set*, the relative conformation stoichiometries, and the individual reactivity profiles are reconstructed. Window sets found to form different numbers of conformations, are reported separately.

When analyzing a transcript, DRACO keeps in memory all the reads mapping to that transcript. At extreme sequencing depths (&gt;500,000X), the memory consumption can become prohibitive, so it might be beneficial to randomly downsample the reads along the transcript, to achieve a lower mean coverage. Furthermore, as the coverage along the transcript might be unevenly distributed because of library preparation and sequencing biases, downsampling will result in a more uniform coverage. Usually, a coverage of __20,000X__ is sufficient for DRACO to deconvolute the underlying structural heterogeneity.<br/>
<br/>
For additional information concerning downsampling, please consult the documentation for the [__rf-count__](https://rnaframework.readthedocs.io/en/latest/rf-count/#downsampling-reads-for-analysis-with-draco) module of the __RNA Framework__.

<br/>

# Usage
To list the required parameters, simply type:

```bash
$ draco --help
```

Parameter         | Type | Description
----------------: | :--: |:------------
__--mm__ | string | Input [mutation map (MM)](https://rnaframework.readthedocs.io/en/latest/rf-count/#mm-mutation-map-format) file
__--output__ | string | Output JSON file (Default: __draco_deconvoluted.json__)
__--processors__ | int | Number of processors to use (Default: __0__)<br/>__Note #1:__ when set to 0, all available processors will be used<br/>__Note #2:__ only the analysis of multiple transcripts can be parallelized. If analyzing a single transcript, setting this parameter to a value > 1 will not speed up the execution
__--whitelist__ | string | A whitelist file, containing the IDs of the transcripts to be analyzed, one per raw (Default: __all transcripts will be analyzed__)
__--shape__ | | Enables spectral analysis on all four bases (default is only A/C bases)<br/>__Note:__ this feature is highly experimental
 | | __Mutation filtering options__
__--minBaseCoverage__ | int | Minimum coverage per base (Default: __100__)
__--minBaseMutations__ | int | Minimum mutations per base (Default: __2__)
__--minMutationFreq__ | float | Bases below this mutation frequency will be discarded as noise (Default: __0.005__)
__--minReadMutations__ | float | Reads with less than these mutations will be discarded as non-informative (Default: __2__)
 | | __Spectral deconvolution options__
__--maxClusters__ | int | Maximum allowed number of clusters/conformations (Default: __unlimited__)
__--minFilteredReads__ | int | Minimum number of reads (post-filtering) to perform spectral analysis (Default: __5__)
__--minPermutations__ | int | Minimum number of permutations performed to build the null model (Default: __8__)
__--maxPermutations__ | int | Maximum number of permutations performed to build the null model (Default: __400__)
__--firstEigengapThresh__ | float | Threshold to consider the first eigengap (Default: __0.9__)<br/>__Note:__ when this threshold is not met, 0 clusters are reported
__--eigengapCumRelThresh__ | float | Minimum relative difference between the eigengap and the null model, as a fraction of the cumulative difference between the previous eigengaps and their respective null models (Default: __0.1__)<br/>__Note:__ this does not apply to the first eigengap
__--eigengapAbsThresh__ | float | Minimum absolute difference between the eigengap and the null model (Default: __0.03__)
__--alpha__ | float | Below this p-value, the null hypothesis is rejected and the eigengap is marked as informative (Default: __0.01__)
__--beta__ | float | Above this p-value, the alternative hypothesis is rejected and the eigengap is marked as non-informative (Default: __0.2__)<br/>__Note:__ this threshold does not apply to the first eigengap
__--firstEigengapBeta__ | float | Beta p-value threshold for the first eigengap (Default: __0.4__)
__--minNullStdev__ | float | Minimum standard deviation for the null model (Default: __0.025__)<br/>__Note:__ when this threshold is not met, ``--extraPermutations`` additional permutations will be performed
__--extraPermutations__ | int | Additional permutations to perform when the standard deviation of the null model is &lt; ``--minNullStdev`` (Default: __50__)
__--minWinBases__ | int | Minimum number of bases in window (post-filtering) to perform spectral analysis (Default: __10__)
__--lookaheadEigengaps__ | int | Number of eigengaps to look ahead after a non-informative eigengap is encountered (Default: __3__)
__--saveEigengapData__ | | Saves eigengap data for plotting (Default: __off__)
__--eigengapDataOut__ | string | Eigengap data output folder (Default: __./eigengap_data__)
 | | __Graph-Cut options__
__--minClusterFraction__ | float | Minimum fraction of reads assigned to each cluster/conformation (Default: __0.05__)<br/>__Note:__ if this threshold is not met, the number of clusters is automatically decreased
 | | __Windowed analysis options__
__--winLenFraction__ | float | Length of the window as a fraction of the median read length (Default: __0.9__)<br/>__Note:__ this parameter and ``--absWinSize`` are mutually exclusive
__--absWinLen__ | int | Absolute length of the window (Default: __0__)<br/>__Note:__ this parameter and ``--winSizeFraction`` are mutually exclusive
__--winOffsetFraction__ | float | Slide offset as fraction of the size of the window (Default: __0.05__)<br/>__Note:__ this parameter and ``--absWinOffset`` are mutually exclusive
__--absWinOffset__ | int | Absolute slide offset (Default: __0__)<br/>__Note:__ this parameter and ``--winOffsetFraction`` are mutually exclusive
__--maxIgnoreWins__ | int | Maximum number of internal windows with a different number of clusters to ignore when merging two external sets of windows (Default: __0__)
__--minExtWins__ | int | Minimum number of external windows, having the same number of clusters, needed to trigger merging (Default: __6__)
__--nonInformativeToSurround__ | int | Non-informative windows (windows with 0 detected clusters) are set to the same number of clusters of surrounding windows (Default: __false__)
__--allNonInformativeToOne__ | int | If *all* windows in the transcript are non-informative (0 clusters), the number of clusters is forced to 1 (Default: __off__)
__--reportNonInformative__ | int | Reports also non-informative windows in the output JSON file

<br/>
## Understanding the algorithm
While it is advisable for most users to run DRACO with its default parameters, as these are the results of a careful and thorough calibration, it might be useful to adjust the analysis on a case-by-case basis.<br/>
<br/>
DRACO analysis is performed in windows. The size of the window is by default ``winLenFraction`` &times; the median read length; an arbitrary length (in bp) for the window can also be used, by setting ``absWinLen``. Similarly, the slide offset is by default ``winOffsetFraction`` &times; the length of the window, but an arbitrary offset (in bp) can be specified via ``absWinOffset``. When long reads are used, reducing the size of the window might allow the identification of small highly dynamic regions of the target transcript. For instance, if 300 bp-long reads are being used, and the structurally-dynamic region within the target transcript is only 90 nt-long, the spectral analysis might not be able to detect it (especially if probing with DMS, as only ~50% of the bases will be informative); in these situations, reducing the size of the window can increase the sensitivity of the analysis.<br/>
<br/>
The __spectral deconvolution__ represents the most critical step of the algorithm, as this is responsible of determining the number of conformations formed by each window. During this step, the *eigengaps*, representing the distance between consecutive *eigenvalues*, are compared to a *null model*, built by performing random permutations over the original data matrix. The distribution of the null model is well approximated by a Weibull distribution. The number of permutations performed varies for each eigengap, and it is comprised between ``minPermutations`` and ``maxPermutations``. Permutations are performed as long as the standard deviation (*&sigma;*) of the null model is &lt; ``minNullStdev``, increasing by ``extraPermutations``.<br/>
<br/>
Starting with the first eigengap, DRACO compares each eigengap to its respective null model. Aim of this comparison is to determine which eigengaps can be considered to be *informative*, as this directly translates into the number of conformations/clusters formed by the analyzed window. For instance, if the first two eigengaps are informative, two conformations (clusters) are present; if the first three eigengaps are informative, three conformations are present, and so on. 
<br/><br/>
![DRACO eigengaps](http://www.rnaframework.com/images/DRACO_eigengaps.png)
<br/><br/>
In order to be considered informative, each eigengap must fulfill a number of criteria:

1. The distance (*d*) between the mean (*&mu;*) of the null model and the eigengap must be &ge; ``eigengapAbsThresh``.

2. The distance between the *&mu;* and the eigengap must be greater than ``egengapCumRelThreshold`` &times; the cumulative sum of the previous eigengaps. For instance, *d<sub>4</sub>*, the distance between the 4<sup>th</sup> eigengap and its null model, must be &ge; ``egengapCumRelThreshold`` &times; (*d<sub>1</sub>* + *d<sub>2</sub>* + *d<sub>3</sub>*). This, of course, does not apply to the first eigengap

3. A *t-test* is used to assess whether the distance between the eigengap and the null model is significant. If the p-value is &le; ``alpha``, the null hypothesis is rejected, and the eigengap is marked as __informative__. If the p-value is &ge; ``beta`` (or in the case of the first eigengap, &ge; ``firstEigengapBeta``), the alternative hypothesis is rejected, and the eigengap is marked as __non-informative__. If the p-value is comprised between these two thresholds, neither of the two hypotheses can be rejected, so additional permutations are performed. If, after having performed ``maxPermutations`` permutations the p-value is still comprised between the two thresholds, then the alternative hypothesis is rejected and the eigengap is marked as __non-informative__.

The analysis terminates when a non-informative eigengap is encountered as, normally, the subsequent eigengaps are non-informative as well. In some particular cases, however, eigengaps have been observed to behave unexpectedly. In such situations, although an eigengap is marked as non-informative, one or more eigengaps immediately downstream of it are again informative. To account for these situations, whenever a non-informative eigengap is encountered, DRACO can perform a *lookahead* evaluation of the ``lookaheadEigengaps`` subsequent eigengaps. If one of these is marked as informative, then the analysis of the eigengaps is allowed to continue.<br/>
<br/>
Once the number of conformations has been determined, the algorithm uses a __Graph-Cut__ approach to weight each base in the window, according to its affinity to each of the different conformations. Reads are then assigned to each conformation. If the fraction of reads assigned to a conformation is &lt; ``minClusterFraction``, the number of conformations for the window is decreased by one, and the step is repeated.<br/>
<br/>
Final step of the analysis involves merging consecutive windows, found to form the same number of conformations, into *window sets*. Let's however consider the following case:
<br/><br/>
![DRACO eigengaps](http://www.rnaframework.com/images/DRACO_mergewindows1.png)
<br/><br/>
In this situation, two sets of windows found to form two conformations, are interrupted by one window forming one conformation. This often happens with short reads (&lt;100 bp), as the analyzed window might not contain enough information to identify the coexisting conformations. By default (left), DRACO would report three window sets. It is possible to account for these cases, hence allowing DRACO to ignore at most ``maxIgnoreWins`` internal windows forming a discordant number of conformations. For the merging to occur, the total number of external windows on either sides of the discordant internal windows, must be &ge; ``minExtWins``.<br/>
It is however possible for the windows on the left and on the right sides of the discordant internal windows, to form a different number of conformations:
<br/><br/>
![DRACO eigengaps](http://www.rnaframework.com/images/DRACO_mergewindows2.png)
<br/><br/>
In such a case, only if the number of windows on one side is &gt; than the number of windows on the other side, as well as &ge; ``minExtWins``, the internal discordant windows will be merged into a single set with the windows on that side.
<br/><br/>
## Output JSON files

DRACO produces a JSON file, with the following structure:

```json
{
   "filename":"sample.mm",
   "transcripts":[
      {
         "id":"RNA_1",
         "sequence":"TCTATTCTACATTGATAGAACAACTCT...TCTGTCTATCACCGCTAGAGCACTCGGTGATTGCA",
         "nReads":20000,
         "windows":[
            {
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
            }
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
__nReads__ | Total number of reads mapping to the transcript
__windows__ | An array containing the deconvolution results for each window set

<br/>
The __windows__ array is further structured as follows:<br/>

Key    | Description
-----: | :----------
__start__ | Start position of the window (0-based)
__end__ | End position of the window (1-based)
__stoichiometries__ | An array containing the relative abundances for the reconstructed conformations
__counts__ | An array of arrays, each one containing the number of mutations per base for each of the reconstructed conformations
__coverage__ | An array of arrays, each one containing the coverage per base for each of the reconstructed conformations, calculated using all the assigned reads post-filtering
__weights__ | An array of arrays, each one containing the weight per base from the Graph-Cut analysis
__preCoverage__ | An array containing the total coverage per base along the window, calculated using only the reads pre-filtering, covering the entire window, that have been used to perform the spectral analysis (therefore, it can be &lt; than the sum of the __coverage__ arrays)

<br/>