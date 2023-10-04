# SpikChIP

### A novel computational methodology to compare multiple ChIP-seq using spike-in chromatin

### Enrique Blanco, Luciano Di Croce, Sergi Aranda

### What is spikChIP?

spikChIP is a standalone pipeline designed in Perl to perform the normalization
of multiple ChIP-seq experiments with spike-in. spikChIP adopts a wide range
of available normalization strategies (raw, traditional, ChIP-Rx, tag removal).
In addition, it implements a new approach based in local regression that is able
to focus the impact of the normalization over the ChIP signal enriched regions
(previously identified by external peak callers), diminishing secondary effects
of the normalization over the background.

### How to run spikChIP

spikChIP is a command line that runs in Linux and Mac OS-X environments:

    ./spikChIP.pl -rtxgs -vclw -b <bin_size_kbp> -k <chrom_key_spikein> -p <0|1|2|3> <configuration_file> <chrominfo_file>

Users can configure the behavior of the program with the following options:

    --clean|-c:  remove intermediate files to reduce the size of the output folder
    --lessMillion|-l: allow the process of BAM files of < 1 Million reads
    --binsize|-b: bin size (default: 10000 bps)
    --chromkey | -k: Key word used to indicate the spike-in species in the chromInfo file (default: FLY)
    --palette|-p: palette (1 for reds, 2 for greens, 3 for blues and 0 for B&W)
    --help|-h: short help
    --verbose|-v: verbose option
    --overwrite|-w: overwrite existing result files
    --raw|-r: Perform the raw normalization
    --traditional|-t: Perform the traditional normalization
    --chiprx|-x: Perform the ChIPRx normalization
    --tagremoval|-g: Perform the tag removal normalization
    --spikchip|s: Perform the spikchip normalization with loess
    --outputfolder|-o: Path to the result folder (default: results/)

The configuration file of spikChIP is a plain-text file in which each line
contains the information about the files of a particular experimental condition:

    #sample	 bam_sample	bam_spike	peaks_sample	peaks_spike
    condition1	map_files/1_sample.bam	map_files/1_spike.bam	peak_files/1_sample_peaks.bed	peak_files/1_spike_peaks.bed
    condition2	map_files/2_sample.bam	map_files/2_spike.bam	peak_files/2_sample_peaks.bed	peak_files/2_spike_peaks.bed
    ...

BAM files of each sample/spike genome (e.g. human reads and fruit fly reads) must be
generated with a mapping tool prior to running spikChIP by the user. Peak files in
BED format must be also generated with a peak calling tool before too. We employed
Bowtie for mapping (multi-locus reads discarded, http://bowtie-bio.sourceforge.net)
and MACS2 for peak calling (broad peaks, https://pypi.org/project/MACS2).

The ChromInfo file is a plain-text file with the list of chromosomes from the
sample and spike-in genomes. This is an example for such a genome type made
of human (hg19) and fruit fly (dm3):

    chr1	249250621
    chr2	243199373
    chr3	198022430
    chr4	191154276
    chr5	180915260
    chr6	171115067
    chr7	159138663
    chr8	146364022
    chr9	141213431
    chr10	135534747
    chr11	135006516
    chr12	133851895
    chr13	115169878
    chr14	107349540
    chr15	102531392
    chr16	90354753
    chr17	81195210
    chr18	78077248
    chr19	59128983
    chr20	63025520
    chr21	48129895
    chr22	51304566
    chrX	155270560
    chrY	59373566
    chr2L_FLY	23011544
    chr2R_FLY	21146708
    chr3L_FLY	24543557
    chr3R_FLY	27905053
    chr4_FLY	1351857
    chrX_FLY	22422827

Please, notice that we use the "FLY" suffix to denote the chromosomes of the spike-in species.
Moreover, it is mandatory that such chromosome names are identically written to describe the location
of the reads in the BAM files of the sample and spike-in experiments for all the conditions as
formally indicated in the configuration file.

### Software requirements

The following programs must be installed in the computer (PATH variable):

* SeqCode (to count reads of BAM files, overlap bins and peaks of BED files):

http://ldicrocelab.crg.eu/index.php
https://github.com/eblancoga/seqcode

* Samtools (to inform about number of reads of BAM files and for downsampling):

http://www.htslib.org/

* gawk (to manage plain-text files divided into columns):

http://www.gnu.org

* R (to perform local normalization and generate the final boxplots):

https://www.r-project.org/
The affy and MASS libraries must be installed

### Input folders and files

The choices in the names of the folders are optional. Users can decide
which folders to use when running spikChIP in the command line.

* config/:
It contains the configuration files employed for the examples shown in the article.

* map_files/:
Please, notice that we do not provide the BAM files of the examples shown in the
article to save space. However, raw data of each example from NCBI GEO is appropriately
referenced elsewhere to perform the mapping apart. Users should deposit the BAM files
of sample/spike genomes within this folder.

* peak files/:
It contains the BED files employed for the examples shown in the article.
Users should deposit the BED files of peaks of their particular examples.

### Output folders and files

spikChIP will create/recycle the following folders in the working directory.

* results/:
It is the folder to save intermediate results and the final set of normalized values
for each strategy (using the average number of reads or the maximum within each bin) in
the current set of ChIP-seq experiments. With this distribution, we provide the results
for the example of H3K79me2 in 2 conditions (25-75 and 75-25) in the case of the
spikChIP strategy. This is the list of files (average values):

(1) Genome-wide list of bins in the human genome (sample) and in the fruit fly (spike-in)
with the values normalized using the SPIKCHIP (RAW, TRADITIONAL, CHIPRX and TAG
REMOVAL are other normalization approaches also calculated by the program):

    FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_avg_normalized_sample.txt.gz
    FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_avg_normalized_spike.txt.gz

(2) List of non-null bins belonging to peaks or background regions in samples/spike-ins:

    FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_sample_avg_bg.txt.gz
    FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_sample_avg_peaks.txt.gz
    FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_spike_avg_bg.txt.gz
    FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_spike_avg_peaks.txt.gz

* plots/

Boxplots automatically generated with and without labels in PDF to show the distribution
of corrected values by each normalization strategy calculated by the spikChIP script
for the current set of ChIP-seq experiments.

* Rscripts/

Scripts written in R to perform the local regression and the final boxplots.

### Main stages of the analysis

The computational pipeline of spikChIP consists of the following stages:

(Step 0)
Read the options, load the configuration file and check the existence of BAM and BED files.

(Step 1)
Segmentation into bins of the same size of the "sample" genome (e.g. human) and the "spike-in"
genome (e.g. fruit fly) according to the list of chromosomes as declared in the ChromInfo.txt file.

(Step 2)
Perform the counting of reads inside the bins (average and maximum values) and normalize these
initial values with multiple approaches (e.g. raw, traditional, ChIP-Rx, tag removal and spikChIP).

(Step 3)
Classification of the bins of both genomes into peaks or background regions according to the peaks
provided by the user for each experiment as declared in the configuration file.

(Step 4)
Generation of the final boxplots to show the distribution of normalized values for each normalization
strategy and class of region (peaks and background) using average or maximum values.

(Step 5)
Cleaning temporary files and compressing resulting final files of all normalization strategies.

### How to generate Custom tracks from spikChIP final results to be visualized in the UCSC browser:

Users can generate custom tracks in BedGraph format from the final files of genome-wide bins. Below,
we show how to generate the human genome track from the spikChIP results:

    FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_avg_normalized_sample.txt.gz

    chr1*1*1001	0.0997764035926034	0.100224097481314
    chr1*100000001*100001001	0.11228442208177	0.107762054394226
    chr1*10000001*10001001	0.191068011376152	0.15125504155222
    ...

    zcat results/FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_avg_normalized_sample.txt.gz | sed 's/\*/ /g' | gawk 'BEGIN{print "track type=bedGraph name=ALL_AVG_SPIKCHIP_25 description=ALL_AVG_SPIKCHIP_25 visibility=full maxHeightPixels=60 color=200,0,0";OFS="\t"}{print $1,$2,$3,$4}' > UCSCtracks/SPIKCHIP_25.bg; gzip UCSCtracks/SPIKCHIP_25.bg;

    track type=bedGraph name=ALL_AVG_SPIKCHIP_25 description=ALL_AVG_SPIKCHIP_25 visibility=full maxHeightPixels=60 color=200,0,0
    chr1	1	1001	0.0997764035926034
    chr1	100000001	100001001	0.11228442208177
    chr1	10000001	10001001	0.191068011376152
    ...

    zcat results/FINAL_H3K79me2_25-75-75-25_SPIKCHIP_1000_avg_normalized_sample.txt.gz | sed 's/\*/ /g' | gawk 'BEGIN{print "track type=bedGraph name=ALL_AVG_SPIKCHIP_75 description=ALL_AVG_SPIKCHIP_75 visibility=full maxHeightPixels=60 color=255,100,100";OFS="\t"}{print $1,$2,$3,$5}' > UCSCtracks/SPIKCHIP_75.bg; gzip UCSCtracks/SPIKCHIP_75.bg;

    track type=bedGraph name=ALL_AVG_SPIKCHIP_75 description=ALL_AVG_SPIKCHIP_75 visibility=full maxHeightPixels=60 color=255,100,100
    chr1	1	1001	0.100224097481314
    chr1	100000001	100001001	0.107762054394226
    chr1	10000001	10001001	0.15125504155222
    ...

We have created a UCSC session with all tracks generated by spikChIP using all normalization strategies:

https://genome.ucsc.edu/s/DiCroceLab/spikChIP_NAR%2DGB_2021

### Benchmarking datasets

We have used multiple ChIP-seq datasets with spike-in from three different publications:

* H3K79me2: two samples

Orlando DA, Chen MW, Brown VE, Solanki S et al. Quantitative ChIP-Seq normalization reveals global
modulation of the epigenome. Cell Rep 2014 Nov 6;9(3):1163-70. PMID: 25437568

https://www.ncbi.nlm.nih.gov/pubmed/25437568

[Jurkat_K79_25%_R1] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1465005

[Jurkat_K79_75%_R1] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1465007

    ./spikChIP.pl -vcp 1 config/H3K79me2_25-75_config.txt config/ChromInfo.txt 2> logs/H3K79me2_25-75.log;

* H3K79me2: five samples

Orlando DA, Chen MW, Brown VE, Solanki S et al. Quantitative ChIP-Seq normalization reveals global
modulation of the epigenome. Cell Rep 2014 Nov 6;9(3):1163-70. PMID: 25437568

https://www.ncbi.nlm.nih.gov/pubmed/25437568

[Jurkat_K79_0%_R1]   

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1465004

[Jurkat_K79_25%_R1]  

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1465005

[Jurkat_K79_50%_R1]  

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1465006

[Jurkat_K79_75%_R1]  

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1465007

[Jurkat_K79_100%_R1] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1465008

    ./spikChIP.pl -vcp 1 config/H3K79me2_0-25-50-75-100_config.txt config/ChromInfo.txt 2> logs/H3K79me2_0-25-50-75-100.log;

* H3K27me3: two samples

Egan B, Yuan C, Craske ML, Labhart P et al. An Alternative Approach to ChIP-Seq Normalization Enables 
Detection of Genome-Wide Changes in Histone H3 Lysine 27 Trimethylation upon EZH2 Inhibition.
PLoS One. 2016 Nov 22;11(11):e0166438.

https://pubmed.ncbi.nlm.nih.gov/27875550

[PC9_control_H3K27me3_Dmspike] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1890165

[PC9_EZH2inh_H3K27me3_Dmspike] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1890166

    ./spikChIP.pl -vcp 3 config/H3K27me3_control-EZH2inh_config.txt config/ChromInfo.txt 2> logs/H3K27me3_control-EZH2inh.log;

* ER: two samples

Guertin MJ, Cullen AE, Markowetz F, Holding AN. Parallel factor ChIP provides essential internal
control for quantitative differential ChIP-seq. Nucleic Acids Res 2018 Jul 6;46(12):e75. 

https://www.ncbi.nlm.nih.gov/pubmed/29672735

[SLX-8047_1b_ER_none] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2747692

[SLX-8047_1a_ER_Fulvestrant] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2747691

[SLX-8047_2b_ER_none] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2747694

[SLX-8047_2a_ER_Fulvestrant] 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2747693

    ./spikChIP.pl -vcp 0 config/ER_R1R2_config.txt config/ChromInfo.txt 2> logs/ER_R1R2.log;

### Examples of running messages of spikChIP for multiple datasets

Users can find in the logs/ folder the whole list of messages provided by spikChIP for
processing each dataset.

### References to cite spikChIP:

SpikChIP: a novel computational methodology to compare multiple ChIP-seq using spike-in chromatin
Enrique Blanco, Luciano Di Croce, Sergi Aranda.
NAR Genom Bioinform. (2021) Jul 27;3(3):lqab064.

https://pubmed.ncbi.nlm.nih.gov/34327329/

Comparative ChIP-seq (Comp-ChIP-seq): a novel computational methodology for genome-wide analysis.
Enrique Blanco, Luciano Di Croce, Sergi Aranda.
bioRxiv 532622 (2019). 

https://www.biorxiv.org/content/10.1101/532622v2
