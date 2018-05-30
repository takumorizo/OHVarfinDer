OHVarfinDer
======================
OHVarfinDer is a **somatic mutation caller**, which can utilize the information of overlapping paired-end reads, heterozygous germline variants near candidate mutations, strand bias. It takes a tumor bam and a normal bam as inputs, and generates a list of candidate somatic mutations.

How to build
----------
### Preliminary ###
Please make sure developer tools for the OS are available (usually installable in batch, but should include make, gcc/g++ etc.) The version we tested is GCC 4.9.3.

### Prepare build dependencies ###
For convenience, the prerequisites of OHVarfinDer are all included in this repository (You can find them at `libs/`). The prerequisites are:
* Boost 1.45.0 (modified to include only the subset of the library used by HapMuC)
* SAMtools 0.1.19
    * SAMtools has two dependencies (vid. [Samtools - Primer](http://biobits.org/samtools_primer.html#InstallingSAMtools)):
        * [GNU curses library](http://www.gnu.org/software/ncurses/)
        * [ZLib compression library](http://zlib.net/)

To prepare them, please execute the following script in your shell:
```sh
% make dependencies
```

### Build OHVarfinDer ###
Just execute the following script in your shell:
```sh
% make
```

How to run
----------

```sh
sh ./utils/runOHVarfinDer.sh ${referenceSequence} ${tumorBam} ${normalBam} ${outputDir} ${region}
```
referenceSequence  : reference sequence used for generating ${tumorBam} and ${normalBam}
outputDir : output files of ${outputDir}/output.variant, ${outputDir}/output.filt.variant is generated
region : ex) chr1:1-1000, same as samtools mpileup region specification.

Convert to VCF
----------

```sh
python ./utils/toVCF.py ${referenceSequence} ${output} ${outputVCF}
```
referenceSequence  : reference sequence used for generating ${tumorBam} and ${normalBam}
output : An output file of output.variant or output.filt.variant
outputVCF : ex) chr1:1-1000, same as samtools mpileup region specification.

The above script uses pysam. Please make sure that pysam(https://pysam.readthedocs.io/en/latest/) is already installed.
```sh
pip install pysam
```


Publication
----------
Under submission.

License
----------
Copyright &copy; 2017 Takuya Moriyama
Released under the [GNU General Public License, Version 3.0][GPL].

This implementation forked from the program [HapMuC][hapmuc] and [Dindel][dindel], which is licensed under the GPLv3.

[GPL]: http://www.gnu.org/licenses/gpl.html
[dindel]: http://www.sanger.ac.uk/resources/software/dindel/
[hapmuc]: https://github.com/usuyama/hapmuc

Note: installing @ HGC SHIROKANE super computer
----------
If you have problem around preparing build dependencies, i.e., make dependencies, please try to compile dependencies at login node.

To prepare them, please execute the following script in your shell @ login node.:
```sh
% make dependencies
```

### Build OHVarfinDer ###
Just execute the following script in your shell after you qlogin:
```sh
% qlogin
% make
```

