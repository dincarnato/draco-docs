The RF Count module is the core component of the framework. It can process any number of SAM/BAM files to calculate per-base RT-stops/mutations and read coverage on each transcript.<br /><br />

# Usage
To list the required parameters, simply type:

```bash
$ rf-count -h
```

Parameter         | Type | Description
----------------: | :--: |:------------
__-p__ *or* __--processors__ | int | Number of processors (threads) to use (Default: __1__)
__-wt__ *or* __--working-threads__ | int | Number of working threads to use for each instance of SAMTools/Bowtie (Default: __1__).<br/>__Note:__ RT Counter executes 1 instance of SAMTools for each processor specified by ``-p``.  At least ``-p <processors>`` * ``-wt <threads>`` processors are required.
__-o__ *or* __--output-dir__ | string | Output directory for writing counts in RC (RNA Count) format (Default: __rf_count/__)
__-ow__ *or* __--overwrite__ | | Overwrites the output directory if already exists
__-t__ *or* __--tmp-dir__ | string | Path to a directory for temporary files creation (Default: __<output-dir>/tmp__)<br/>__Note:__ If the provided directory does not exist, it will be created
__-s__ *or* __--samtools__ | string | Path to ``samtools`` executable (Default: assumes ``samtools`` is in PATH)
__-r__ *or* __--sorted__ | | In case SAM/BAM files are passed, assumes that they are already sorted lexicographically by transcript ID, and numerically by position
__-t5__ *or* __--trim-5prime__ | int[,int] | Comma separated list (no spaces) of values indicating the number of bases trimmed from the 5'-end of reads in the respective sample SAM/BAM files (Default: __0__)<br/>__Note #1:__ Values must be provided in the same order as the input files (e.g. rf-count -t5 0,5 file1.bam file2.bam, will consider 0 bases trimmed from file1 reads, and 5 bases trimmed from file2 reads)<br/>__Note #2:__ If a single value is specified along with multiple SAM/BAM files, it will be used for all files
__-fh__ *or* __--from-header__ | | Instead of providing the number of bases trimmed from 5'-end of reads through the ``-t5`` (or ``--trim-5prime``) parameter, RF Count will try to guess it automatically from the header of the provided SAM/BAM files
__-f__ *or* __--fasta__ | string | Path to a FASTA file containing the reference transcripts<br/>__Note:__ Transcripts in this file must match transcripts in SAM/BAM file headers
__-mf__ *or* __--mask-file__ | string | Path to a mask file
__-pn__ *or* __--primary-only__ | | Considers only primary alignments (SAM bitwise flag != 256)
__-po__ *or* __--paired-only__ | | When processing SAM/BAM files from paired-end experiments, only those reads for which both mates are mapped will be considered
__-pp__ *or* __--properly-paired__ | | When processing SAM/BAM files from paired-end experiments, only those reads mapped in a proper pair will be considered
__-i__ *or* __--include-clipped__ | | Include reads that have been soft/hard-clipped at their 5'-end when calculating RT-stops<br/>__Note:__ The default behavior is to exclude soft/hard-clipped reads. When this option is active, the RT-stop position is considered to be the position preceding the clipped bases. This option has no effect when ``-m`` (or ``--count-mutations``) is enabled.
__-mq__ *or* __--map-quality__ | int | Minimum mapping quality to consider a read (Default: __10__)
__-co__ *or* __--coverage-only__ | | Only calculates per-base coverage (disables RT-stops/mutations count)
__-m__ *or* __--count-mutations__ | | Enables mutations count instead of RT-stops count (for SHAPE-MaP/DMS-MaPseq)
 | | __Mutation count mode options__
__-ds__ *or* __--discard-shorter__ | int | Discards reads shorter than this length (excluding clipped bases, Default: __1__)<br/>__Note:__ when set to *MEDIAN* (case-insensitive), the median read length will be used
__-q__ *or* __--min-quality__ | int | Minimum quality score value to consider a mutation (Phred+33, requires ``-m``, Default: __20__)
__-es__ *or* __--eval-surrounding__ | | When considering a mutation/indel, also evaluates the quality of surrounding bases (&#177;1 nt)<br/>__Note:__ the quality score threshold set by ``-q`` (or ``--min-quality``) also applies to these bases
__-nd__ *or* __--no-deletions__ | | Ignores deletions
__-ni__ *or* __--no-insertions__ | | Ignores insertions
__-na__ *or* __--no-ambiguous__ | | Ignores ambiguously mapped deletions<br/>__Note:__ the default behavior is to re-align them to their right-most valid position (or to their left-most valid position if ``-la`` has been specified)
__-la__ *or* __--left-align__ | | Re-aligns ambiguously mapped deletions to their left-most valid position
__-rd__ *or* __--right-deletion__ | | Only the right-most base in a deletion is marked as mutated
__-ld__ *or* __--left-deletion__ | | Only the left-most base in a deletion is marked as mutated
__-md__ *or* __--max-deletion-len__ | int | Ignores deletions longer than this number of nucleotides (Default: __10__)
__-me__ *or* __--max-edit-distance__ | float | Discards reads with editing distance frequency higher than this threshold (0<m&le;1, Default: __0.15__ [15%])
__-eq__ *or* __--median-quality__ | int | Median quality score threshold for discarding low-quality reads (Phred+33, Default: __20__)
__-dc__ *or* __--discard-consecutive__ | int | Discards consecutive mutations within this distance from eachothers
__-cc__ *or* __--collapse-consecutive__ | | Collapses consecutive mutations/indels toward the 3'-most one (recommended for SHAPE-MaP experiments)
__-mc__ *or* __--max-collapse-distance__ | int | Maximum distance between consecutive mutations/indels to allow collapsing (requires ``-cc``, &ge;0, Default_ __2__)
__-mv__ *or* __--max-coverage__ | int | Downsamples reads to achieve this maximum mean per-base coverage (&ge;1000, Default: __off__)
__-mm__ *or* __--mutation-map__ | | Generates a mutation map (MM) file for alternative structure deconvolution with [DRACO](https://github.com/dincarnato/draco)

<br/>
## Coverage calculation
It is important to note that, by default (so when counting RT-stop events), RF Count will consider the contribution of the RT drop-off event to the coverage of the base on which the drop-off has occurred.<br/>
Take into account the following example:
<br/><br/>
![Coverage calculation](http://www.rnaframework.com/images/coverage.png)
<br/><br/>
In this case, 3 sets of reads have been aligned to the reference, resulting in 3 RT-stop sites, respectively corresponding to 2, 4 and 3 RT drop-off events. When executing RF count without specifying either parameters ``-co`` or ``-m``, these additional counts will be added to the coverage of the base corresponding to the site on which the RT dropped off (-1 with respect to the read start mapping position).
<br/><br/>
## Deletions re-alignment in mutational profiling-based methods
Mutational profiling (MaP) methods for RNA structure analysis are based on the ability of certain reverse transcriptase enzymes to read-through the sites of SHAPE/DMS modification under specific reaction conditions. Some of them (e.g. SuperScript II) can introduce deletions when encountering a SHAPE/DMS-modified residue. When performing reads mapping, the aligner often reports a single possible alignment of the deletion, although many equally-scoring alignments are possible.<br/>
To avoid counting of ambiguously aligned deletions, that can introduce noise in the measured structural signal, RF Count performs a *deletion re-alignment step* to detect and re-align/discard these ambiguously aligned deletions (see also "__Handling of mutations/indels__"):
<br/><br/>
![Ambiguous deletions](http://www.rnaframework.com/images/ambiguous_deletions.png)
<br/><br/>
For more information, please refer to Smola *et al*., 2015 (PMID: [26426499](https://www.ncbi.nlm.nih.gov/pubmed/26426499)).
<br/><br/>
## Handling of mutations/indels

By giving a rapid look to the numerous parameters provided by RF Count, it appears immediately clear that different parameter combinations produce very different outcomes.
Here follows a brief scheme aimed at illustrating the different behaviors of RF Count with different parameter combinations (dots correspond to sites of assigned mutations):
<br/><br/>
![RF Count MaP handling](http://www.rnaframework.com/images/rf-count_MaP.png)
<br/><br/>
## Downsampling reads for analysis with DRACO

One of the possibilities offered by RF Count is the generation of mutation map (MM) files via the ``-mm`` (or ``--mutation-map``) parameter. These files store, for each mapped read, the index of the mutations with respect to the reference transcripts, and they can be used as input for [__DRACO__](https://draco-docs.readthedocs.io/en/latest/draco/).<br/><br/>
DRACO (<u>D</u>econvolution of <u>R</u>NA <u>A</u>lternative <u>CO</u>nformations) is a method that allows deconvoluting the individual structural profiles of coexisting RNA structural conformations, from MaP experiments.<br/>
When analyzing a transcript, DRACO keeps in memory all the reads mapping to that transcript. At extreme sequencing depths (&gt;500,000X), the memory consumption can become prohibitive, so it might be beneficial to randomly downsample the reads along the transcript, to achieve a lower mean coverage. Furthermore, as the coverage along the transcript might be unevenly distributed because of library preparation and sequencing biases, downsampling will result in a more uniform coverage. Usually, a coverage of __20,000X__ is sufficient for DRACO to deconvolute the underlying structural heterogeneity.
<br/><br/>
![DRACO downsampling](http://www.rnaframework.com/images/rf-count_DRACO_downsample.png)
<br/><br/>
Downsampling is controlled via the ``-mv`` (or ``--max-coverage``). In the above example, reads have been downsampled to reach a final maximum coverage per base of 1,000X. To do this, RF Count estimates the median read length, and calculates the theoretical number of reads of that length, starting at each position of the transcript, needed to achieve that mean coverage per base. However, as a consequence of adapter trimming and soft clipping, read lengths can significantly differ within the same experiment; this results in a less uniform coverage than it would be theoretically expected, as it is possible to observe in the second track (in green).<br/>
To partially compensate for this, it is advisable to set the ``-ds`` (or ``--discard-shorter``) parameter to the value __MEDIAN__ (case insensitive). In this way, reads shorter than the median read length will be discarded, resulting in a much more uniform coverage, as it is possible to observe in the third track (in blue).
<br/><br/>
## RC (RNA Count) format

RF Count produces a RC (RNA Count) file for each analyzed sample. RC files are binary files storing the transcript’s sequence, per-base RT-stop/mutation counts, per-base read coverage, and total number of mapped reads. These files can be indexed for fast random access.<br/>
Each entry in a RC file is structured as follows:

Field             | Description    |  Type
----------------: | :------------: | :----------
__len\_transcript\_id__ | Length of the transcript ID (plus 1, including NULL) | uint32\_t
__transcript\_id__ | Transcript ID (NULL terminated) | char[len\_transcript\_id]
__len\_seq__ | Length of sequence | uint32\_t
__seq__ | 4-bit encoded sequence: 'ACGTN' -> \[0,4] (High nybble first) | uint8\_t\[(len_seq+1)/2]
__counts__ | Transcript's per base RT-stops (or mutations) | uint32\_t[len\_seq]
__coverage__ | Transcript's per base coverage | uint32\_t[len\_seq]
__n<sub>t</sub>__ | Total reads mapping on transcript | unint64\_t

RC files EOF stores the number of total mapped reads, and is structured as follows:

Field             | Description    |  Type
----------------: | :------------: | :----------
__n__ | Total experiment reads (only primary alignments are counted,<br/>SAM bitwise flag != 256) | uint64\_t
__version__ | RC file version | uint16\_t
__marker__ | EOF marker (\\x5b\\x65\\x6f\\x66\\x72\\x63\\x5d) | char[7]

The current RC standard's version is __1__.<br/>
RCI (RC Index) files enable random access to transcript data within RC files.<br/>
The RCI index is structured as follows:

Field             | Description    |  Type
----------------: | :------------: | :----------
__len\_transcript\_id__ | Length of the transcript ID (plus 1, including NULL) | uint32\_t
__transcript\_id__ | Transcript ID (NULL terminated) | char[len\_transcript\_id]
__offset__ | Offset of transcript in the RC file | uint64\_t

!!! note "Information"
    All values are forced to be in little-endian byte-order.

<br/>
## MM (Mutation Map) format

When invoked with the ``-mm`` parameter, RF Count produces a MM (Mutation Map) file for each analyzed sample. These files serve as input for the deconvolution of coexisting RNA conformations with [DRACO](https://github.com/dincarnato/draco). MM files are binary files storing, for each aligned read, the start and end mapping position, the position of the mutations (or indels) with respect to the reference, as well as the reference transcript’s ID and sequence, . These files can be indexed for fast random access.<br/>
Each entry in a MM file is structured as follows:

Field             | Description    |  Type
----------------: | :------------: | :----------
__len\_transcript\_id__ | Length of the transcript ID (plus 1, including NULL) | uint32\_t
__transcript\_id__ | Transcript ID (NULL terminated) | char[len\_transcript\_id]
__len\_seq__ | Length of sequence | uint32\_t
__seq__ | 4-bit encoded sequence: 'ACGTN' -> \[0,4] (High nybble first) | uint8\_t\[(len_seq+1)/2]
__read\_count__ | Number of reads mapping on the transcript | uint32\_t

This is followed by *n* blocks (where *n* is equal to __read\_count_), with the following structure:

Field             | Description    |  Type
----------------: | :------------: | :----------
__start__ | Read's start mapping position | uint32\_t
__end__ | Read's end mapping position | uint32\_t
__mut\_count__ | Number of mutations (or indels) | uint32\_t
__indices__ | Indices of the mutations (or indels) | uint32\_t[mut\_count]

MM files always end with the EOF marker \\x5b\\x6d\\x6d\\x65\\x6f\\x66\\x5d.
MMI (MM Index) files enable random access to transcript data within MM files.<br/>
The MMI has the same structure of the RCI:

Field             | Description    |  Type
----------------: | :------------: | :----------
__len\_transcript\_id__ | Length of the transcript ID (plus 1, including NULL) | uint32\_t
__transcript\_id__ | Transcript ID (NULL terminated) | char[len\_transcript\_id]
__offset__ | Offset of transcript in the RC file | uint64\_t

!!! note "Information"
    All values are forced to be in little-endian byte-order.

<br/>
## Mask file
The mask file allows excluding specific transcript regions from being counted. This is particularly useful when performing targeted MaP analyses, in order to mask the primer pairing regions.<br/>
The mask file is composed of one or more lines, each one reporting the transcript ID and a comma (or semicolon) separated list of either base ranges (0-based, inclusive), or of nucleotide sequences of the regions that need to be masked:<br/>

```
Transcript_1;AGCGTATTAGCGATGCGATGCGA;25-38;504-551
Transcript_2,331-402,AUAUGGAUCGGACG,984-1008
Transcript_3;GUUACAUUCGA,98-123;47-68
```
Transcript regions specified in the mask file will have both 0 counts and coverage in the resulting RC file.
