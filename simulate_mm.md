The ``simulate_mm`` utility provides a simple way to generate simulated mutation map (MM) files for transcripts forming any number of conformations, mixed at arbitrary stoichiometries.<br/><br/>

# Usage
To list the required parameters, simply type:

```bash
$ simulate_mm -h
```

Parameter         | Type | Description
----------------: | :--: |:------------
__-o__ *or* __--outProfiles__ | string | Output file with structure profiles updated according to the simulation
__-p__ *or* __--stoichiometry__ | string | Comma-separated list of % conformation stoichiometries<br/>__Note #1:__ the stoichiometries must sum to approx. 100 (tollerance: 97-103)<br/>__Note #2:__ When no stoichiometry is specified, the conformations are assumed to be equimolar
__-c__ *or* __--meanCoverage__ | int | Mean sequencing depth (coverage) per base<br/>__Note:__ this parameter and ``--nReads`` are mutually exclusive
__-n__ *or* __--nReads__ | int | Number of reads mapping to each transcript<br/>__Note:__ this parameter and ``--meanCoverage`` are mutually exclusive
__--probability__ | float | Sets the *p* value for generation of the binomial distribution of mutations (Default: __0.01927__)<br/>__Note:__ the default value has been learnt empirically from [Homan et al., 2014](https://pubmed.ncbi.nlm.nih.gov/25205807/)
__-s__ *or* __--readLen__ | int | Length (in bp) of the simulated reads
__-t__ *or* __--text__ | string | Output MM file's "human-readable" version

<br/>
## Input structure profile file
RNAs to be generated can be provided in the form of structure profile files.<br/>

```text
TCTATTCTACATTGATAGAAC...ACCGCTAGAGCACTCGGTGATTGCA
x.xxxxxxxx.xxxxx.xxxx...xxxxxx.xxx.xxxx.xxx.xxx.x
xxx.xx.x.xxxxx.xxx.x....xxxxxx.x.xx.xxxxxxxxxxxx.
0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,...,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0
0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0,1,...,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1


TATCTTATCACTTGCTCGCCA...CAAGATCGCGACATAGGTGCTTGAC
(((.)).).(((((.(((.(....)))))).).)).xxxxxxxxxxxx.
(.((((((((.(((((.((((...)))))).))).)))).))).))).)
0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0,1,...,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1
0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,...,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0

```
<br/>
These files contain entries composed of four parts:

1. The transcript's sequence
2. A textual representation of the possible structures formed by the RNA
3. A numeric profile representing the structures indicated in __2__
4. An empty line, marking the end of the entry

The textual representation uses dots ("__.__") to represent unpaired bases, and "__x__" or parantheses ("__(__" and "__)__") to represent paired bases. No check is made on the proper balancing of parentheses in dot-bracket structures.<br/>The numeric profile must match the textual representation of the structure. In these profiles, __0__s indicate paired bases, while any value __&ge;1__ represents an upaired base (the used numeric value is not relevant to the simulation).<br/><br/>
Besides generating an MM file, a new structure profile file will also be generated (controlled via the ``-o`` or ``--outProfiles`` parameter), identical to the one provided as input, but with the numeric profiles updated to the actual mutation counts from the simulation.