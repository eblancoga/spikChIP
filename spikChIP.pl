#!/usr/bin/perl -w                                #
#                                                 #
# Script spikChIP.pl                              #
#                                                 #
#                                                 #
# Input : configuration_file, chrominfo_file      # 
# Output: Normalization files and boxplots        #
#                                                 #
#                                                 #
# by Enrique Blanco @ CRG (2021)                  #
###################################################

use strict;
use Getopt::Std;
use Term::ANSIColor;



#DEFINEs
my $PROGRAM = "spikChIP";
my $TRUE = 1;
my $FALSE = 0;
my $SUCCESS = 0;
my $DEFAULT_BIN_SIZE = 10000;
my $MIN_BIN_SIZE = 1000;
my $MAX_BIN_SIZE = 100000;
my $DEFAULT_PALETTE = 1;
my $MIN_PALETTE = 0;
my $MAX_PALETTE = 3;
my $MEGA = 1000000;
my $PSEUDOCOUNT = 0.1;
my $RAW_TOKEN = "RAW";
my $TRADITIONAL_TOKEN = "TRADITIONAL";
my $CHIPRX_TOKEN = "CHIPRX";
my $TAGREMOVAL_TOKEN = "TAGREMOVAL";
my $SPIKCHIP_TOKEN = "SPIKCHIP";
my $DOWNSAMPLING_TOKEN = "adjusted";
my $PEAKS_TOKEN = "binpeaks";
my $FINAL_TOKEN = "FINAL";
my $SAMPLE_TOKEN = "sample";
my $SPIKE_TOKEN = "spike";
my $AVG_TOKEN = "avg";
my $MAX_TOKEN = "max";
my $RESULTS = "results/";
my $PLOTS = "plots/";
my $RSCRIPTS = "Rscripts/";
my $NOLABELS_TOKEN = "nolabels";
my $N_STRATEGIES = 5;


## Step 0. Reading arguments
my %opt;
my ($start,$date);
my ($configuration_filename,$chrominfo_file);
my $n_files;
my @NAMES;
my (@BAM_SAMPLES,@BAM_SPIKES);
my (@READS_SAMPLES,@READS_SPIKES);
my (@BAM_SAMPLES_DOWN,@BAM_SPIKES_DOWN);
my (@READS_SAMPLES_DOWN,@READS_SPIKES_DOWN);
my (@PEAKS_SAMPLES,@PEAKS_SPIKES);
my (@BINS_PEAKS_SAMPLE,@BINS_PEAKS_SPIKE);
my $bin_size;
my $palette;
my $newdir;
my @CLEAN_PROCEDURE;
my @SAVE_PROCEDURE;

# -c: clean intermediate files of results to save space
# -d: allow the process of BAM files of < 1 Million reads
# -k: bin size (default: 10000 bps)
# -p: palette (1 for reds, 2 for greens, 3 for blues and 0 for B&W)
# -h: short help
# -v: verbose option
# -w: overwrite existing result files

(getopts('whvk:p:cd',\%opt)) or print_error("COMMAND LINE: Problems reading options\n");

print_mess("$PROGRAM.pl by Enrique Blanco @ CRG (2021)\n");
print_help();
print_mess("\n");

# 0.1 Save the starting time
$start = time();
$date = localtime();

# 0.2 Acquire options and arguments
print_mess("[$date] Stage 0.  Configuration of the pipeline\n");
print_mess("Reading options of the command line\n");
$n_files = $#ARGV+1;
($n_files == 2) or print_error("ERROR: Two arguments are required but $n_files are provided!\n$PROGRAM -vhcdk <bin_size_kbp> -p <0|1|2|3> <configuration_file> <chrominfo_file>\nPlease, type spikChIP -h for further assistance or spikChIP -v for verbose mode\n");
($configuration_filename,$chrominfo_file) = @ARGV;
print_mess("Syntax is correct and options are acquired");
print_ok();
print_mess("\n");

# 0.3 Checking binaries of SeqCode and samtools are available
print_mess("Checking for external software availability\n");
#CheckExternalSoftware();
print_mess("Closing external software check...");
print_ok();
print_mess("\n");

# 0.4 Prepare the output folders
print_mess("Generating output folders (if necessary)\n");
# create the results/ folder
$newdir = $RESULTS;
print_mess("Trying the $RESULTS directory\n");
mkdir($newdir) or print_mess("It is already existing\n");
# create the Rscripts/ folder
$newdir = $RSCRIPTS;
print_mess("Trying the $RSCRIPTS directory\n");
mkdir($newdir) or print_mess("It is already existing\n");
# create the plots/ folder
$newdir = $PLOTS;
print_mess("Trying the $PLOTS directory\n");
mkdir($newdir) or print_mess("It is already existing");
print_ok();
print_mess("\n");

# 0.5 Check the chrominfo file 
print_mess("Checking information on chromosome sizes\n");
CheckChromInfoFile($chrominfo_file);
print_mess("Closing ChromInfo file check...");
print_ok();
print_mess("\n");

# 0.6 Checking bin size option (if selected)
print_mess("Establishing the bin size\n");
$bin_size = SettingBinsize();
print_mess("Effective bin size: $bin_size");
print_ok();
print_mess("\n");

# 0.7 Checking palette (if selected)
print_mess("Establishing the palette\n");
$palette = SettingPalette();
print_mess("Effective palette: $palette");
print_ok();
print_mess("\n");

# 0.8 Loading filenames from configuration file
print_mess("Reading configuration file, calculating reads and down-sampling\n");
print_mess("Processing the configuration file\n");
ReadingConfigurationFile($configuration_filename);
print_mess("\n");
print_mess("Starting downsampling of the experiments\n");
DownsamplingOperations();
print_mess("Finishing Stage 0. Configuration...");
print_ok();
print_mess("\n");

# Step 1. Generating the segmentation of the sample and spike genomes
my $command;
my ($spike_bins,$sample_bins);
my ($n_spike_bins,$n_sample_bins);
my $line;


$date = localtime();
print_mess("[$date] Stage 1.  Producing the segmentation of both genomes in bins ($bin_size bps)\n");

# spike segmentation
$spike_bins = $RESULTS.join("-",@NAMES)."_"."spike_".$bin_size.".bed";
if(!(-e $spike_bins) or exists($opt{w}))
{
    $command = "grep FLY $chrominfo_file | gawk 'BEGIN{OFS=\"\\t\";offset=$bin_size;}{for(i=1;i<\$2-offset;i=i+offset) print \$1,i,i+offset;}' > $spike_bins";
    print_mess("$command\n");
    system($command);
}else{
    print_mess("\t The file ", $spike_bins, " already exist. Skipping spike segmentation\n");
}

# count the number of spike bins
$n_spike_bins = 0;
(open(FILEBINS,$spike_bins)) or print_error("SPIKE BINS: FILE $spike_bins file can not be opened");
while($line=<FILEBINS>)
{
    $n_spike_bins++;
}
close(FILEBINS);
print_mess("$n_spike_bins bins generated in the segmentation of the spike genome\n");

# sample genome segmentation
$sample_bins = $RESULTS.join("-",@NAMES)."_"."sample_".$bin_size.".bed";
if(!(-e $sample_bins) or exists($opt{w}))
{
    $command = "grep -v FLY $chrominfo_file | gawk 'BEGIN{OFS=\"\\t\";offset=$bin_size;}{for(i=1;i<\$2-offset;i=i+offset) print \$1,i,i+offset;}' > $sample_bins";
    print_mess("$command\n");
    system($command);
}else{
    print_mess("\t The file ", $sample_bins, " already exist. Skipping sample segmentation\n");
}

# count the number of sample bins
$n_sample_bins = 0;
(open(FILEBINS,$sample_bins)) or print_error("SAMPLE BINS: FILE $sample_bins file can not be opened");
while($line=<FILEBINS>)
{
    $n_sample_bins++;
}
close(FILEBINS);
print_mess("$n_sample_bins bins generated in the segmentation of the sample genome");
CleanFile($spike_bins);
CleanFile($sample_bins);
#
print_ok();
print_mess("Finishing Stage 1. Segmentation...");
print_ok();
print_mess("\n");
exit 42;

# Step 2. Processing each experiment using distinct normalization strategies
my $i;
my $n_experiments;
my ($name,$bam_sample,$bam_spike);


$date = localtime();
print_mess("[$date] Stage 2.  Processing each experiment with different normalization methods\n");

$n_experiments = scalar(@NAMES);
print_mess("Total pairs of sample/spike experiments to be analyzed: $n_experiments\n");

for($i=0; $i<$n_experiments; $i++)
{
    $name = $NAMES[$i];
    print_mess("Processing the $name experiment\n");
    #
    print_mess("Starting RAW normalization...\n");
    NormalizationRaw($i);
    print_ok();
    print_mess("\n");
    #
    print_mess("Starting TRADITIONAL normalization...\n");
    NormalizationTraditional($i);
    print_ok();
    print_mess("\n");
    #
    print_mess("Starting CHIP-RX normalization...\n");
    NormalizationChIPRX($i);
    print_ok();
    print_mess("\n");
    #
    print_mess("Starting TAG_REMOVAL normalization...\n");
    NormalizationTagRemoval($i);
    print_ok();
    print_mess("\n");
}

# join the values of all experiments in a single file per normalization class
print_mess("Join raw data values from all experiments\n");
JoinNormValues($RAW_TOKEN);
print_mess("\n");
print_mess("Join traditional data values from all experiments\n");
JoinNormValues($TRADITIONAL_TOKEN);
print_mess("\n");
print_mess("Join ChIP-RX data values from all experiments\n");
JoinNormValues($CHIPRX_TOKEN);
print_mess("\n");
print_mess("Join tag removal data values from all experiments\n");
JoinNormValues($TAGREMOVAL_TOKEN);
print_mess("\n");

print_mess("Starting $PROGRAM local normalization...\n");
print_mess("Prepare $PROGRAM data values for local normalization\n");
PreparespikChIPValues();
print_mess("Run local regression adjusting sample values with spike values\n");
RunspikChIPValues();
print_ok();
print_mess("Finishing Stage 2. Normalization...");
print_ok();
print_mess("\n");


# Step 3. Identification of bins associated to peaks/bg
$date = localtime();
print_mess("[$date] Stage 3.  Distinguishing bins of spike and sample genomes into peaks or bg according to peak files\n");

print_mess("Generic classification of genome segmentations into peaks and bg\n");
ClassifyBins();
print_ok();
print_mess("\n");

print_mess("Processing values of each normalization strategy into peaks and bg\n");
print_mess("Working with $RAW_TOKEN values\n");
ClassifyNormalizationValues($RAW_TOKEN);
print_ok();
print_mess("\n");
print_mess("Working with $TRADITIONAL_TOKEN values\n");
ClassifyNormalizationValues($TRADITIONAL_TOKEN);
print_ok();
print_mess("\n");
print_mess("Working with $CHIPRX_TOKEN values\n");
ClassifyNormalizationValues($CHIPRX_TOKEN);
print_ok();
print_mess("\n");
print_mess("Working with $TAGREMOVAL_TOKEN values\n");
ClassifyNormalizationValues($TAGREMOVAL_TOKEN);
print_ok();
print_mess("\n");
print_mess("Working with $SPIKCHIP_TOKEN values\n");
ClassifyNormalizationValues($SPIKCHIP_TOKEN);
print_ok();
print_mess("Finishing Stage 3. Segregation...");
print_ok();
print_mess("\n");


# Step 4. Generation of the resulting boxplots of values in bins of peaks and bg
$date = localtime();
print_mess("[$date] Stage 4.  Generating the final boxplots of values\n");

print_mess("Boxplots using the average values\n");
GenerateBoxplot($SPIKE_TOKEN,$AVG_TOKEN);
GenerateBoxplot($SAMPLE_TOKEN,$AVG_TOKEN);
print_ok();
print_mess("\n");

print_mess("Boxplots using the maximum values\n");
GenerateBoxplot($SPIKE_TOKEN,$MAX_TOKEN);
GenerateBoxplot($SAMPLE_TOKEN,$MAX_TOKEN);
print_ok();
print_mess("Finishing Stage 4. Drawing...");
print_ok();
print_mess("\n");


## Step 5. End of the analysis
my $stop;
my ($hours,$mins,$secs);

$date = localtime();
print_mess("[$date] Stage 5.  Successfully finishing the work\n");

# Cleaning intermediate files (option -c)
print_mess("Cleaning temporary files");
CleaningFiles();
print_ok();
print_mess("\n");

# Compressing final files
print_mess("Saving final files (gzip)\n");
SavingFiles();
print_mess("Finishing Stage 5. Cleaning and Saving...");
print_ok();
print_mess("\n");

# Save the ending time
$stop = time();
$hours = int((($stop-$start)/3600)*1000)/1000;
$mins = int((($stop-$start)/60)*1000)/1000;
$secs = int(($stop-$start)*1000)/1000;

print_mess("Total running time (hours):",$hours," hours\n");
print_mess("Total running time (minutes):",$mins," mins\n");
print_mess("Total running time (seconds):",$secs," secs");
print_ok();
exit(0);
##

############################### Subroutines ###############################

sub print_help
{
    if (exists($opt{h}))
    {
	print STDERR color("bold blue"),"spikChIP_v1.0\t\t\t\tUser commands\t\t\t\tspikChIP\n\n";
	print STDERR color("bold blue"),"NAME\n\tspikChIP, a tool to normalize ChIP-seq experiments by spike-in correction\n\n";
	print STDERR color("bold blue"),"SYNOPSIS:\n\t$PROGRAM -vcdk <bin_size_kbp> -p <0|1|2|3> <configuration_file> <chrominfo_file>\n\n";
	print STDERR color("bold blue"),"OPTIONS:\n";
	print STDERR color("bold blue"),"\t-c: clean intermediate files of results to save space\n";
	print STDERR color("bold blue"),"\t-d: allow the process of BAM files of < 1 Million reads\n";
	print STDERR color("bold blue"),"\t-k: bin size (default: 10000 bps)\n";
	print STDERR color("bold blue"),"\t-p: palette (1 for reds, 2 for greens, 3 for blues and 0 for B&W)\n";
	print STDERR color("bold blue"),"\t-h: short help\n";
	print STDERR color("bold blue"),"\t-w: overwrite existing files option\n";
	print STDERR color("bold blue"),"\t-v: verbose option\n\n";
	print STDERR color("bold blue"),"SEE ALSO\n";
	print STDERR color("bold blue"),"\tGitHub source code: https://github.com/eblancoga/spikChIP\n\n";
	print STDERR color("bold blue"),"AUTHORS\n";
	print STDERR color("bold blue"),"\tWritten by Enrique Blanco, Luciano Di Croce and Sergi Aranda (2021)\n\n";
	print STDERR color("bold blue"),"spikChIP_v1.0\t\t\t\tUser commands\t\t\t\tspikChIP\n";
	print STDERR color("reset");

	exit(0);
    }
}

sub print_mess
{
        my @mess = @_;

        print STDERR color("bold green"),"%%%% @mess" if (exists($opt{v}));
	print STDERR color("reset");
}

sub print_error
{
        my @mess = @_;
	my $stop;
	my ($hours,$mins,$secs);
	
        print STDERR color("bold red"),"%%%% @mess\n\n";
	print STDERR color("reset");

	$date = localtime();
	print_mess("[$date] Stage X.  Abruptly finishing the work\n");
	
	# Save the ending time
	$stop = time();
	$hours = int((($stop-$start)/3600)*1000)/1000;
	$mins = int((($stop-$start)/60)*1000)/1000;
	$secs = int(($stop-$start)*1000)/1000;
	
	print_mess("Total running time (hours):",$hours," hours\n");
	print_mess("Total running time (minutes):",$mins," mins\n");
	print_mess("Total running time (seconds):",$secs," secs");
	print_ok();
	exit();
}

sub print_ok
{
    if (exists($opt{v}))
    {
	print STDERR color("bold green"), " [OK]\n";
	print STDERR color("reset");
    }
}

sub CheckChromInfoFile
{
    my $file = $_[0];
    my $command;

    
    # confirm the existence of the file
    if (!(-e $file))
    {
	print_error("CHROM INFO FILE: The file $file does not exist!");
    }
    else
    {
	$command = "cat $file 1>&2";
	system($command);
    }
}

sub CalculateReads
{
    my $bam_file = $_[0];
    my $name = $_[1];
    my $command;
    my $out_file;
    my $n_reads;
    my ($line,@record);

    
    # running the samtools to calculate the number of reads of the BAM file
    $out_file = $RESULTS.$name."_flagstat.txt";
    
    if(!(-e $out_file) or exists($opt{w}))
    {
        $command = "samtools flagstat $bam_file > $out_file";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file ", $name."_flagstat.txt", " already exist. Skipping samtools flagstats\n");
    }

    (open(FILEFLAG,$out_file)) or print_error("SAMTOOLS FLAGSTAT: FILE $out_file file can not be opened");
    while($line=<FILEFLAG>)
    {
	if ($line=~/total/)
	{
	    @record = split(/\s+/,$line);
	    $n_reads = $record[0];

	    if ($n_reads < $MEGA && !exists($opt{d}))
	    {
		print_error("The number of reads is lower than $MEGA reads. Please, run again spikChIP with the option -d");
	    }
	}
    }
    close(FILEFLAG);
    #
    CleanFile($out_file);
    
    return($n_reads);
}

sub ReadingConfigurationFile
{
    my $filename = $_[0];
    my ($line,@record);
    my ($name,$bam_sample,$bam_spike,$peaks_sample,$peaks_spike);
    my $n_reads;
    
    
    (open(FILE,$filename)) or print_error("CONFIG_FILE: FILE $filename configuration file can not be opened");
    while($line=<FILE>)
    {
	# skip comment lines
	next if ($line=~/^#/);

	# extract the information of the current experiment
	@record = split(/\s+/,$line);

	$name = $record[0];
	$bam_sample = $record[1];
	$bam_spike = $record[2];
	$peaks_sample = $record[3];
	$peaks_spike = $record[4];

	# confirm the existence of BAM and BED files
	if (!(-e $bam_sample))
	{
	    print_error("CONFIG_FILE: ($filename, $name) The file $bam_sample does not exist!");
	}
	else
	{
	    push(@BAM_SAMPLES,$bam_sample);
	    $n_reads = CalculateReads($bam_sample,$name."_sample");
	    push(@READS_SAMPLES,$n_reads);
	    print_mess("($name, $n_reads reads) Confirming filename $bam_sample\n");
	}
	#
	if (!(-e $bam_spike))
	{
	    print_error("CONFIG_FILE: ($filename, $name) The file $bam_spike does not exist!");
	}
	else
	{
	    push(@BAM_SPIKES,$bam_spike);
	    $n_reads = CalculateReads($bam_spike,$name."_spike");
	    push(@READS_SPIKES,$n_reads);
	    print_mess("($name, $n_reads reads) Confirming filename $bam_spike\n");
	}
	#
	if (!(-e $peaks_sample))
	{
	    print_error("CONFIG_FILE: ($filename, $name) The file $peaks_sample does not exist!");
	}
	else
	{
	    push(@PEAKS_SAMPLES,$peaks_sample);
	    print_mess("($name) Confirming filename $peaks_sample\n");
	}
	#
	if (!(-e $peaks_spike))
	{
	    print_error("CONFIG_FILE: ($filename, $name) The file $peaks_spike does not exist!");
	}
	else
	{
	    push(@PEAKS_SPIKES,$peaks_spike);
	    print_mess("($name) Confirming filename $peaks_spike\n");
	}

	# save this experiment in the collection of samples
	push(@NAMES,$name);
    }
    close(FILE);
}

sub DownsamplingOperations
{
    my $i;
    my $min_reads;
    my $factor;
    my $command;
    my ($input_file,$output_file);
    my $n_reads;

    
    # first, identify the lowest value of spike reads
    $min_reads = $READS_SPIKES[0];
    for ($i=1; $i<scalar(@READS_SPIKES); $i++)
    {
	if ($READS_SPIKES[$i] < $min_reads)
	{
	    $min_reads = $READS_SPIKES[$i];
	}
    }
    print_mess("Minimum value of spike reads: $min_reads\n");

    # second, apply this factor for downsampling each spike and sample BAM files
    for ($i=0; $i<scalar(@NAMES); $i++)
    {
	$factor = $min_reads / $READS_SPIKES[$i];
	$factor = "1.0" if ($factor == 1);
	print_mess($NAMES[$i].": factor of $factor = $min_reads / ".$READS_SPIKES[$i]."\n");

	# spike bam
	$input_file = $BAM_SPIKES[$i];
	($output_file = $input_file) =~ s/\.bam/\_$DOWNSAMPLING_TOKEN\.bam/g;
	$BAM_SPIKES_DOWN[$i] = $output_file;
	if(!(-e $output_file) or exists($opt{w}))
    {
        $command = "samtools view -h -b -s $factor -o $output_file $input_file";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file ", $output_file, " already exist. Skipping down-sampling\n");
    }
	
	
	$n_reads = CalculateReads($output_file,$NAMES[$i]."_spike_adjusted");
	push(@READS_SPIKES_DOWN,$n_reads); 
	print_mess("($output_file, $n_reads) Final number of reads after downsampling\n");
	
	# sample bam
	$input_file = $BAM_SAMPLES[$i];
	($output_file = $input_file) =~ s/\.bam/\_$DOWNSAMPLING_TOKEN\.bam/g;
	$BAM_SAMPLES_DOWN[$i] = $output_file;
	if(!(-e $output_file) or exists($opt{w}))
    {
        $command = "samtools view -h -b -s $factor -o $output_file $input_file";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file ", $output_file, " already exist. Skipping down-sampling\n");
    }
	
	$n_reads = CalculateReads($output_file,$NAMES[$i]."_sample_adjusted");
	push(@READS_SAMPLES_DOWN,$n_reads);
	print_mess("($output_file, $n_reads) Final number of reads after downsampling\n");
    }
}

sub SettingBinsize
{
    my $bin_size;

    if (exists($opt{k}))
    {
	$bin_size = $opt{"k"};
    }
    else
    {
	$bin_size =$DEFAULT_BIN_SIZE;
    }

    # range of valid bin sizes
    if (($bin_size < $MIN_BIN_SIZE) || ($bin_size > $MAX_BIN_SIZE))
    {
	print_error("Please, set the bin size in the valid range (current $bin_size): $MIN_BIN_SIZE - $MAX_BIN_SIZE");
    }

    return($bin_size);
}

sub SettingPalette
{
    my $palette;

    if (exists($opt{p}))
    {
	$palette = $opt{"p"};
    }
    else
    {
	$palette =$DEFAULT_PALETTE;
    }

    # range of valid paletter
    if (($palette < $MIN_PALETTE) || ($palette > $MAX_PALETTE))
    {
	print_error("Please, set the palette in the valid range (current $palette): $MIN_PALETTE - $MAX_PALETTE");
    }

    if ($palette == 1)
    {
	$palette = "c(\"pink\",\"red\",\"brown\")";
    }
    else
    {
	if ($palette == 2)
	{
	    $palette = "c(\"darkolivegreen1\",\"green4\",\"darkgreen\")";
	}
	else
	{
	    if ($palette == 3)
	    {
		$palette = "c(\"lightskyblue\",\"blue\",\"darkblue\")";
	    }
	    else
	    {
		$palette = "c(\"white\")";
	    }
	}
    }
    
    return($palette);
}

sub CheckExternalSoftware
{
    no warnings 'exec';

    my $command;
    my $check;
    

    print_mess("Searching samtools...");
    $command = "samtools --version > /dev/null";
    $check = system($command);
    if ($check!=$SUCCESS)
    {
	print_error("Please, install samtools in your computer (http://www.htslib.org/download/)");
    }
    else
    {
	print_ok();
    }

    print_mess("Searching gawk...");
    $command = "gawk --version > /dev/null";
    $check = system($command);
    if ($check!=$SUCCESS)
    {
	print_error("Please, install gawk in your computer (http://www.gnu.org)");
    }
    else
    {
	print_ok();
    }

    print_mess("Searching R...");
    $command = "R --version > /dev/null";
    $check = system($command);
    if ($check!=$SUCCESS)
    {
	print_error("Please, install R in your computer (https://www.r-project.org)");
    }
    else
    {
	print_ok();
    }

    print_mess("Searching SeqCode |recoverChIPlevels...");
    $command = "recoverChIPlevels -h 2> /dev/null";
    $check = system($command);
    if ($check!=$SUCCESS)
    {
	print_error("Please, install SeqCode|recoverChIPlevels in your computer (https://github.com/eblancoga/seqcode)");
    }
    else
    {
	print_ok();
    }

    print_mess("Searching SeqCode |matchpeaks...");
    $command = "matchpeaks -h 2> /dev/null";
    $check = system($command);
    if ($check!=$SUCCESS)
    {
	print_error("Please, install SeqCode|matchpeaks in your computer (https://github.com/eblancoga/seqcode)");
    }
    else
    {
	print_ok();
    }
}

sub NormalizationRaw
{
    my $i = $_[0];
    my ($bam_sample,$bam_spike);
    my $out_name;
    my ($file_all,$file_avg,$file_max);

    # input bam files
    $bam_sample = $BAM_SAMPLES[$i];
    $bam_spike = $BAM_SPIKES[$i];

    # Spike bins
    $out_name = $NAMES[$i]."_".$RAW_TOKEN."_".$bin_size."_spike";
    if(!(-e $out_name) or exists($opt{w}))
    {
        if (exists($opt{d}))
        {
            $command = "recoverChIPlevels -dns $MEGA $chrominfo_file $bam_spike $spike_bins $out_name";
        }
        else
        {
            $command = "recoverChIPlevels -ns $MEGA $chrominfo_file $bam_spike $spike_bins $out_name";
        }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file ", $out_name, " already exist. Skipping the spike Raw Normalization\n");
    }
    #
    $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
    $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    
    if(!(-e $file_avg) or exists($opt{w}))
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file ", $file_avg, " already exist. Skipping creation of the spike avg\n");
    }

    if(!(-e $file_max) or exists($opt{w}))
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file ", $file_max, " already exist. Skipping creation of the spike max\n");
    }
    

    # Sample bins
    $out_name = $NAMES[$i]."_".$RAW_TOKEN."_".$bin_size."_sample";
    if(!(-e $out_name) or exists($opt{w}))
    {
        if (exists($opt{d}))
        {
            $command = "recoverChIPlevels -dns $MEGA $chrominfo_file $bam_sample $sample_bins $out_name";
        }
        else
        {
            $command = "recoverChIPlevels -ns $MEGA $chrominfo_file $bam_sample $sample_bins $out_name";
        }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file ", $out_name, " already exist. Skipping the sample Raw Normalization\n");
    }

    #
    $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
    $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    
    if(!(-e $file_avg) or exists($opt{w}))
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file ", $file_avg, " already exist. Skipping creation of the sample avg\n");
    }

    if(!(-e $file_max) or exists($opt{w}))
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command");
        system($command);
    }else{
        print_mess("\t The file ", $file_max, " already exist. Skipping creation of the sample max\n");
    }
}

sub NormalizationTraditional
{
    my $i = $_[0];
    my ($bam_sample,$bam_spike);
    my $out_name;
    my ($file_all,$file_avg,$file_max);
    my $total_reads;
    
    
    # input bam files
    $bam_sample = $BAM_SAMPLES[$i];
    $bam_spike = $BAM_SPIKES[$i];
    $total_reads = $READS_SAMPLES[$i] + $READS_SPIKES[$i];

    # Spike bins
    $out_name = $NAMES[$i]."_".$TRADITIONAL_TOKEN."_".$bin_size."_spike";

    if (exists($opt{d}))
    {
	$command = "recoverChIPlevels -dns $total_reads $chrominfo_file $bam_spike $spike_bins $out_name";
    }
    else
    {
	$command = "recoverChIPlevels -ns $total_reads $chrominfo_file $bam_spike $spike_bins $out_name";
    }
    print_mess("$command\n");
    system($command);
    #
    $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
    $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
    print_mess("$command\n");
    system($command);
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
    print_mess("$command\n");
    system($command);

    # Sample bins
    $out_name = $NAMES[$i]."_".$TRADITIONAL_TOKEN."_".$bin_size."_sample";

    if (exists($opt{d}))
    {
	$command = "recoverChIPlevels -dns ".$READS_SAMPLES[$i]." $chrominfo_file $bam_sample $sample_bins $out_name";
    }
    else
    {
	$command = "recoverChIPlevels -ns ".$READS_SAMPLES[$i]." $chrominfo_file $bam_sample $sample_bins $out_name";
    }
    print_mess("$command\n");
    system($command);
    #
    $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
    $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
    print_mess("$command\n");
    system($command);
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
    print_mess("$command");
    system($command);
}

sub NormalizationChIPRX
{
    my $i = $_[0];
    my ($bam_sample,$bam_spike);
    my $out_name;
    my ($file_all,$file_avg,$file_max);
    my $total_reads;
    
    
    # input bam files
    $bam_sample = $BAM_SAMPLES[$i];
    $bam_spike = $BAM_SPIKES[$i];
    $total_reads = $READS_SAMPLES[$i] + $READS_SPIKES[$i];

    # Spike bins
    $out_name = $NAMES[$i]."_".$CHIPRX_TOKEN."_".$bin_size."_spike";

    if (exists($opt{d}))
    {
	$command = "recoverChIPlevels -dns ".$READS_SPIKES[$i]." $chrominfo_file $bam_spike $spike_bins $out_name";
    }
    else
    {
	$command = "recoverChIPlevels -ns ".$READS_SPIKES[$i]." $chrominfo_file $bam_spike $spike_bins $out_name";
    }
    print_mess("$command\n");
    system($command);
    #
    $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
    $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
    print_mess("$command\n");
    system($command);
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
    print_mess("$command\n");
    system($command);

    # Sample bins
    $out_name = $NAMES[$i]."_".$CHIPRX_TOKEN."_".$bin_size."_sample";

    if (exists($opt{d}))
    {
	$command = "recoverChIPlevels -dns ".$READS_SPIKES[$i]." $chrominfo_file $bam_sample $sample_bins $out_name";
    }
    else
    {
	$command = "recoverChIPlevels -ns ".$READS_SPIKES[$i]." $chrominfo_file $bam_sample $sample_bins $out_name";
    }
    print_mess("$command\n");
    system($command);
    #
    $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
    $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
    print_mess("$command\n");
    system($command);
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
    print_mess("$command");
    system($command);
}

sub NormalizationTagRemoval
{
    my $i = $_[0];
    my ($bam_sample,$bam_spike);
    my $out_name;
    my ($file_all,$file_avg,$file_max);

    
    # input bam files
    $bam_sample = $BAM_SAMPLES_DOWN[$i];
    $bam_spike = $BAM_SPIKES_DOWN[$i];

    # Spike bins
    $out_name = $NAMES[$i]."_".$TAGREMOVAL_TOKEN."_".$bin_size."_spike";

    if (exists($opt{d}))
    {
	$command = "recoverChIPlevels -dns $MEGA $chrominfo_file $bam_spike $spike_bins $out_name";
    }
    else
    {
	$command = "recoverChIPlevels -ns $MEGA $chrominfo_file $bam_spike $spike_bins $out_name";
    }
    print_mess("$command\n");
    system($command);
    #
    $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
    $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
    print_mess("$command\n");
    system($command);
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
    print_mess("$command\n");
    system($command);

    # Sample bins
    $out_name = $NAMES[$i]."_".$TAGREMOVAL_TOKEN."_".$bin_size."_sample";

    if (exists($opt{d}))
    {
	$command = "recoverChIPlevels -dns $MEGA $chrominfo_file $bam_sample $sample_bins $out_name";
    }
    else
    {
	$command = "recoverChIPlevels -ns $MEGA $chrominfo_file $bam_sample $sample_bins $out_name";
    }
    print_mess("$command\n");
    system($command);
    #
    $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
    $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
    print_mess("$command\n");
    system($command);
    $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
    print_mess("$command");
    system($command);
}

sub JoinNormValues
{
    my $token = $_[0];
    my $i;
    my $n_experiments;
    my ($out_name0,$file_avg0,$file_max0);
    my ($out_name1,$file_avg1,$file_max1);
    my ($out_name,$file_avg,$file_max);
    my ($commandAVG,$commandMAX);
    my ($final_avg,$final_max);
    my $folder;
    

    # join the values of all experiments
    $n_experiments = scalar(@NAMES);

    # First, spike values
    # at least, two experiments
    $out_name0 = $NAMES[0]."_".$token."_".$bin_size."_spike";
    $file_avg0 = "$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_avg.bed";
    $file_max0 = "$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_max.bed";
    $out_name1 = $NAMES[1]."_".$token."_".$bin_size."_spike";
    $file_avg1 = "$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_avg.bed";
    $file_max1 = "$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_max.bed";
    $commandAVG = "JoinNFiles.pl -v $file_avg0 $file_avg1 ";
    $commandMAX = "JoinNFiles.pl -v $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
	$out_name = $NAMES[$i]."_".$token."_".$bin_size."_spike";
	$file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
	$file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
	$commandAVG = $commandAVG." $file_avg ";
	$commandMAX = $commandMAX." $file_max ";
    }

    $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_spike_avg.txt";
    $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_spike_max.txt";
    #
    $commandAVG = $commandAVG." | sort -k 1,1 > $final_avg";
    print_mess("$commandAVG\n");
    system($commandAVG);
    $commandMAX = $commandMAX." | sort -k 1,1 > $final_max";
    print_mess("$commandMAX\n");
    system($commandMAX);  
    #
    SaveFile($final_avg);
    SaveFile($final_max);

    # Second, sample values
    # at least, two experiments
    $out_name0 = $NAMES[0]."_".$token."_".$bin_size."_sample";
    $file_avg0 = "$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_avg.bed";
    $file_max0 = "$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_max.bed";
    $out_name1 = $NAMES[1]."_".$token."_".$bin_size."_sample";
    $file_avg1 = "$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_avg.bed";
    $file_max1 = "$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_max.bed";
    $commandAVG = "JoinNFiles.pl -v $file_avg0 $file_avg1 ";
    $commandMAX = "JoinNFiles.pl -v $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
	$out_name = $NAMES[$i]."_".$token."_".$bin_size."_sample";
	$file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
	$file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
	$commandAVG = $commandAVG." $file_avg ";
	$commandMAX = $commandMAX." $file_max ";
    }

    $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_sample_avg.txt";
    $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_sample_max.txt";
    #
    $commandAVG = $commandAVG." | sort -k 1,1 > $final_avg";
    print_mess("$commandAVG\n");
    system($commandAVG);
    $commandMAX = $commandMAX." | sort -k 1,1 > $final_max";
    print_mess("$commandMAX\n");
    system($commandMAX);
    #
    SaveFile($final_avg);
    SaveFile($final_max);
    
    # remove the seqcode folders (spike and sample)
    print_mess("Removing SeqCode folders\n");
    for($i=0; $i<$n_experiments; $i++)
    {
	$folder = $NAMES[$i]."_".$token."_".$bin_size."_spike_recoverChIPlevels/";
	$command = "rm -rf $folder";
	system($command);
	$folder = $NAMES[$i]."_".$token."_".$bin_size."_sample_recoverChIPlevels/";
	$command = "rm -rf $folder";
	system($command);
    }
}

sub PreparespikChIPValues
{
    my ($spike_values,$sample_values);
    my ($output_file1,$output_file2,$output_file3,$output_file4);
    my $command;
    my ($i,$fields);
    

    # (A) normalization based in avg values
    $spike_values = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$bin_size."_spike_avg.txt";
    $sample_values = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$bin_size."_sample_avg.txt";

    # output files
    $output_file1 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_avg.txt";
    $output_file2 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_avg_values.txt";
    $output_file3 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_avg_names.txt";
    $output_file4 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_avg_names2.txt";
    #
    CleanFile($output_file1);
    CleanFile($output_file2);
    CleanFile($output_file3);
    CleanFile($output_file4);
    
    $command = "cat $spike_values $sample_values > $output_file1";
    print_mess("$command\n");
    system($command);
    
    # extract all the columns with values (one per experiment)
    $fields = "";
    for($i=0; $i<$n_experiments-1; $i++)
    {
	$fields = $fields." \$".($i+2).",";
    }
    $fields = $fields." \$".($i+2);
    
    $command = "gawk 'BEGIN{OFS=\"\\t\"}{print $fields}' $output_file1 > $output_file2";
    print_mess("$command\n");
    system($command);

    $command = "gawk 'BEGIN{OFS=\"\\t\"}{print NR}' $output_file1 > $output_file3";
    print_mess("$command\n");
    system($command);

    $command = "gawk 'BEGIN{OFS=\"\\t\"}{print NR,\$1}' $output_file1 > $output_file4";
    print_mess("$command\n");
    system($command);

    # (B) normalization based in max values
    $spike_values = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$bin_size."_spike_max.txt";
    $sample_values = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$bin_size."_sample_max.txt";

    # output files
    $output_file1 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_max.txt";
    $output_file2 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_max_values.txt";
    $output_file3 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_max_names.txt";
    $output_file4 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_max_names2.txt";
    #
    CleanFile($output_file1);
    CleanFile($output_file2);
    CleanFile($output_file3);
    CleanFile($output_file4);
    
    $command = "cat $spike_values $sample_values > $output_file1";
    print_mess("$command\n");
    system($command);

    $command = "gawk 'BEGIN{OFS=\"\\t\"}{print $fields}' $output_file1 > $output_file2";
    print_mess("$command\n");
    system($command);

    $command = "gawk 'BEGIN{OFS=\"\\t\"}{print NR}' $output_file1 > $output_file3";
    print_mess("$command\n");
    system($command);

    $command = "gawk 'BEGIN{OFS=\"\\t\"}{print NR,\$1}' $output_file1 > $output_file4";
    print_mess("$command\n");
    system($command);
}

sub RunspikChIPValues
{
    my $Rfile;
    my ($output_file1,$output_file2,$output_file3,$output_file4);
    my ($final_avg,$final_avg_spike,$final_avg_sample);
    my ($final_max,$final_max_spike,$final_max_sample);
    my $n_total_bins;
    my $Routput_file;
    my $line;
    my ($i,$fields);
    my $binname;
    

    # spikChIP on avg values
    print_mess("Performing the analysis on average values\n");
    #
    $Rfile = $RSCRIPTS.join("-",@NAMES)."_".$bin_size."_avg.R";
    $output_file1 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_avg.txt";
    $output_file2 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_avg_values.txt";
    $output_file3 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_avg_names.txt";
    $output_file4 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_avg_names2.txt";
    #
    $final_avg = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_avg_normalized.txt";
    $final_avg_spike = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_avg_normalized_spike.txt";
    $final_avg_sample = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_avg_normalized_sample.txt";
    #
    CleanFile($final_avg);
    SaveFile($final_avg_spike);
    SaveFile($final_avg_sample);
    #
    $n_total_bins = $n_spike_bins + $n_sample_bins;
    #
    # R code to perform the loess regression on the whole set of bins in all conditions
    (open(RFILE,'>',$Rfile)) or print_error("R SCRIPT: FILE $Rfile file can not be opened to write");
    print RFILE "library(affy)\n";
    print RFILE "library(MASS)\n";
    print RFILE "c <- scan(\"$output_file2\",sep=\"\\t\")\n";
    print RFILE "c <- c + $PSEUDOCOUNT\n";
    print RFILE "m <-matrix(c,$n_total_bins,$n_experiments,byrow=TRUE)\n";
    print RFILE "d <- read.table(\"$output_file3\")\n";
    print RFILE "rownames(m) <- d[,1]\n";
    print RFILE "s <- seq(1,$n_spike_bins)\n";
    print RFILE "mn <- normalize.loess(m,subset=s)\n";
    print RFILE "write.table(mn,file=\"$final_avg\",sep=\"\\t\",row.names=TRUE,col.names=FALSE,quote=FALSE)\n";
    close(RFILE);
    #
    # execute R script
    $command = "R CMD BATCH $Rfile";
    print_mess("$command\n");
    system($command);
    #
    # error check in Rout file
    $Routput_file = join("-",@NAMES)."_".$bin_size."_avg.Rout";
    (open(ROUT,$Routput_file)) or print_error("R SCRIPTS (avg): FILE $Routput_file can not be opened");
    while($line=<ROUT>)
    {
	if ($line=~/Error/)
	{
	    print_error("R running: $line\n");
	}
    }
    close(ROUT);
    $command = "rm -f $Routput_file";
    print_mess("$command\n");
    system($command);
    #
    # extract all the columns with values (one per experiment)
    $fields = "";
    for($i=0; $i<$n_experiments-1; $i++)
    {
	$fields = $fields." \$".($i+2).",";
    }
    $fields = $fields." \$".($i+2);
    $binname = " \$".($i+3);
    #
    # distinguish spike from sample bins
    $command = "join $final_avg $output_file4 | grep FLY | gawk 'BEGIN{OFS=\"\\t\"}{print $binname,$fields}'> $final_avg_spike";
    print_mess("$command\n");
    system($command);
    $command = "join $final_avg $output_file4 | grep -v FLY | gawk 'BEGIN{OFS=\"\\t\"}{print $binname,$fields}'> $final_avg_sample";
    print_mess("$command\n");
    system($command);

    # spikChIP on max values
    print_mess("Performing the analysis on average values\n");
    #
    $Rfile = $RSCRIPTS.join("-",@NAMES)."_".$bin_size."_max.R";
    $output_file1 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_max.txt";
    $output_file2 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_max_values.txt";
    $output_file3 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_max_names.txt";
    $output_file4 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_spike-sample_max_names2.txt";
    #
    $final_max = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_max_normalized.txt";
    $final_max_spike = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_max_normalized_spike.txt";
    $final_max_sample = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_max_normalized_sample.txt";
    #
    CleanFile($final_max);
    SaveFile($final_max_spike);
    SaveFile($final_max_sample);
    #
    $n_total_bins = $n_spike_bins + $n_sample_bins;
    #
    # R code to perform the loess regression on the whole set of bins in all conditions
    (open(RFILE,'>',$Rfile)) or print_error("R SCRIPT: FILE $Rfile file can not be opened to write");
    print RFILE "library(affy)\n";
    print RFILE "library(MASS)\n";
    print RFILE "c <- scan(\"$output_file2\",sep=\"\\t\")\n";
    print RFILE "c <- c + $PSEUDOCOUNT\n";
    print RFILE "m <-matrix(c,$n_total_bins,$n_experiments,byrow=TRUE)\n";
    print RFILE "d <- read.table(\"$output_file3\")\n";
    print RFILE "rownames(m) <- d[,1]\n";
    print RFILE "s <- seq(1,$n_spike_bins)\n";
    print RFILE "mn <- normalize.loess(m,subset=s)\n";
    print RFILE "write.table(mn,file=\"$final_max\",sep=\"\\t\",row.names=TRUE,col.names=FALSE,quote=FALSE)\n";
    close(RFILE);
    #
    # execute R script
    $command = "R CMD BATCH $Rfile";
    print_mess("$command\n");
    system($command);
    #
    # error check in Rout file
    $Routput_file = join("-",@NAMES)."_".$bin_size."_max.Rout";
    (open(ROUT,$Routput_file)) or print_error("R SCRIPTS (max): FILE $Routput_file can not be opened");
    while($line=<ROUT>)
    {
	if ($line=~/Error/)
	{
	    print_error("R running: $line\n");
	}
    }
    close(ROUT);
    $command = "rm -f $Routput_file";
    print_mess("$command\n");
    system($command);
    #
    # distinguish spike from sample bins
    $command = "join $final_max $output_file4 | grep FLY | gawk 'BEGIN{OFS=\"\\t\"}{print $binname,$fields}'> $final_max_spike";
    print_mess("$command\n");
    system($command);
    $command = "join $final_max $output_file4 | grep -v FLY | gawk 'BEGIN{OFS=\"\\t\"}{print $binname,$fields}'> $final_max_sample";
    print_mess("$command");
    system($command);
}

sub ClassifyBins
{
    my $i;
    my $command;
    my ($folder,$out_name,$out_common,$out_peaks_file);
    
    
    # for each experiment, compare spike/sample genome bins to peaks in the corresponding spike/genome
    for ($i=0; $i<scalar(@NAMES); $i++)
    {
	print_mess("Working with sample $NAMES[$i]: peaks Vs. bins of spike\n");
	
	# use SeqCode to identify the spike bins overlapping with spike peaks
	$out_name = $NAMES[$i]."_".$bin_size."_spike_bins";
	$folder = $PEAKS_TOKEN."_".$out_name."_matchpeaks";
	$out_common = $folder."/common_".$PEAKS_TOKEN."_".$out_name.".bed";
	$out_peaks_file = $RESULTS.$NAMES[$i]."_".$bin_size."_bins_spike_peaks.bed";
	$BINS_PEAKS_SPIKE[$i] = $out_peaks_file;
	#
	CleanFile($out_peaks_file);
	#
	$command = "matchpeaks -v $PEAKS_SPIKES[$i] $spike_bins $PEAKS_TOKEN $out_name";
	print_mess("$command\n");
	system($command);
	$command = "grep spike_bins $out_common > $out_peaks_file";
	print_mess("$command\n");
	system($command);
	#
	print_mess("Removing SeqCode folder\n");
	$command = "rm -rf $folder";
	system($command);
	print_mess("\n");
	#
	print_mess("Working with sample $NAMES[$i]: peaks Vs. bins of sample\n");
	
	# use SeqCode to identify the sample bins overlapping with sample peaks
	$out_name = $NAMES[$i]."_".$bin_size."_sample_bins";
	$folder = $PEAKS_TOKEN."_".$out_name."_matchpeaks";
	$out_common = $folder."/common_".$PEAKS_TOKEN."_".$out_name.".bed";
	$out_peaks_file = $RESULTS.$NAMES[$i]."_".$bin_size."_bins_sample_peaks.bed";
	$BINS_PEAKS_SAMPLE[$i] = $out_peaks_file;
	#
	CleanFile($out_peaks_file);
	#
	$command = "matchpeaks -v $PEAKS_SAMPLES[$i] $sample_bins $PEAKS_TOKEN $out_name";
	print_mess("$command\n");
	system($command);
	$command = "grep sample_bins $out_common > $out_peaks_file";
	print_mess("$command\n");
	system($command);
	#
	print_mess("Removing SeqCode folder\n");
	$command = "rm -rf $folder";
	system($command);
	#
    }
}

sub ClassifyNormalizationValues
{
    my $token = $_[0];
    my $n_experiments;
    my $name;
    my $command;
    my ($input_file,$output_file);
    my ($n_field,$field_notnull);
    my ($file_avg0,$file_avg1,$file_max0,$file_max1,$file_avg,$file_max);
    my ($final_avg,$final_max);
    my ($commandAVG,$commandMAX);
    

    # (A) generate the files of values for bins of peaks and bins of bg for each experiment
    $n_experiments = scalar(@NAMES);
    for($i=0; $i<$n_experiments; $i++)
    {
	$name = $NAMES[$i];
	# select the corresponding column in the results file for the current experiment
	$n_field = "\$".($i + 2);
	#
	# extract values for spike with avg
	if ($token eq $SPIKCHIP_TOKEN)
	{
	    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_avg_normalized_spike.txt";
	}
	else
	{
	    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_spike_avg.txt";
	}
	$output_file = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_spike_avg_peaks.txt";
	$command = "grep -v track ".$BINS_PEAKS_SPIKE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
	print_mess("$command\n");
	system($command);
	#
	CleanFile($output_file);
	#
	$output_file = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_spike_avg_bg.txt";
	$command = "grep -v track ".$BINS_PEAKS_SPIKE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join -v 2 - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
	print_mess("$command\n");
	system($command);
	#
	CleanFile($output_file);
	#
	# extract values for sample with avg
	if ($token eq $SPIKCHIP_TOKEN)
	{
	    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_avg_normalized_sample.txt";
	}
	else
	{
	    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_sample_avg.txt";
	}
	$output_file = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_sample_avg_peaks.txt";
	$command = "grep -v track ".$BINS_PEAKS_SAMPLE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
	print_mess("$command\n");
	system($command);
	#
	CleanFile($output_file);
	#
	$output_file = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_sample_avg_bg.txt";
	$command = "grep -v track ".$BINS_PEAKS_SAMPLE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join -v 2 - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
	print_mess("$command\n");
	system($command);
	#
	CleanFile($output_file);
	#
	# extract values for spike with max
	if ($token eq $SPIKCHIP_TOKEN)
	{
	    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_max_normalized_spike.txt";
	}
	else
	{
	    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_spike_max.txt";
	}
	$output_file = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_spike_max_peaks.txt";
	$command = "grep -v track ".$BINS_PEAKS_SPIKE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
	print_mess("$command\n");
	system($command);
	#
	CleanFile($output_file);
	#
	$output_file = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_spike_max_bg.txt";
	$command = "grep -v track ".$BINS_PEAKS_SPIKE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join -v 2 - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
	print_mess("$command\n");
	system($command);
	#
	CleanFile($output_file);
	#
	# extract values for sample with max
	if ($token eq $SPIKCHIP_TOKEN)
	{
	    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_max_normalized_sample.txt";
	}
	else
	{
	    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_sample_max.txt";
	}
	$output_file = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_sample_max_peaks.txt";
	$command = "grep -v track ".$BINS_PEAKS_SAMPLE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
	print_mess("$command\n");
	system($command);
	#
	CleanFile($output_file);
	#
	$output_file = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_sample_max_bg.txt";
	$command = "grep -v track ".$BINS_PEAKS_SAMPLE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join -v 2 - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
	print_mess("$command\n");
	system($command);
	#
	CleanFile($output_file);
    }

    # (B) join the files of values for bins of peaks and bins of bg of each experiment into a single final experiment
    # generate the not null comparison string for the final filter
    $field_notnull = "gawk '{if (";
    for($i=0; $i<$n_experiments-1; $i++)
    {
	$n_field = "\$".($i + 2);
	$field_notnull = $field_notnull."$n_field !=0 && ";
    }
    $n_field = "\$".($i + 2);
    $field_notnull = $field_notnull."$n_field !=0) print \$0}'";

    # First, spike values for peaks
    # at least, two experiments
    $file_avg0 = $RESULTS.$NAMES[0]."_".$token."_".$bin_size."_spike_avg_peaks.txt";
    $file_max0 = $RESULTS.$NAMES[0]."_".$token."_".$bin_size."_spike_max_peaks.txt";
    $file_avg1 = $RESULTS.$NAMES[1]."_".$token."_".$bin_size."_spike_avg_peaks.txt";
    $file_max1 = $RESULTS.$NAMES[1]."_".$token."_".$bin_size."_spike_max_peaks.txt";
    $commandAVG = "join $file_avg0 $file_avg1 ";
    $commandMAX = "join $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
	$file_avg = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_spike_avg_peaks.txt";
	$file_max = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_spike_max_peaks.txt";
	$commandAVG = $commandAVG." | join - $file_avg ";
	$commandMAX = $commandMAX." | join - $file_max ";
    }

    $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_spike_avg_peaks.txt";
    $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_spike_max_peaks.txt";
    #
    $commandAVG = $commandAVG." | $field_notnull > $final_avg";
    print_mess("$commandAVG\n");
    system($commandAVG);
    $commandMAX = $commandMAX." | $field_notnull > $final_max";
    print_mess("$commandMAX\n");
    system($commandMAX);
    #
    SaveFile($final_avg);
    SaveFile($final_max);

    # Second, sample values for peaks
    # at least, two experiments
    $file_avg0 = $RESULTS.$NAMES[0]."_".$token."_".$bin_size."_sample_avg_peaks.txt";
    $file_max0 = $RESULTS.$NAMES[0]."_".$token."_".$bin_size."_sample_max_peaks.txt";
    $file_avg1 = $RESULTS.$NAMES[1]."_".$token."_".$bin_size."_sample_avg_peaks.txt";
    $file_max1 = $RESULTS.$NAMES[1]."_".$token."_".$bin_size."_sample_max_peaks.txt";
    $commandAVG = "join $file_avg0 $file_avg1 ";
    $commandMAX = "join $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
	$file_avg = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_sample_avg_peaks.txt";
	$file_max = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_sample_max_peaks.txt";
	$commandAVG = $commandAVG." | join - $file_avg ";
	$commandMAX = $commandMAX." | join - $file_max ";
    }

    $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_sample_avg_peaks.txt";
    $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_sample_max_peaks.txt";
    #
    $commandAVG = $commandAVG." | $field_notnull > $final_avg";
    print_mess("$commandAVG\n");
    system($commandAVG);
    $commandMAX = $commandMAX." | $field_notnull > $final_max";
    print_mess("$commandMAX\n");
    system($commandMAX);
    #
    SaveFile($final_avg);
    SaveFile($final_max);

    # Third, spike values for bg
    # at least, two experiments
    $file_avg0 = $RESULTS.$NAMES[0]."_".$token."_".$bin_size."_spike_avg_bg.txt";
    $file_max0 = $RESULTS.$NAMES[0]."_".$token."_".$bin_size."_spike_max_bg.txt";
    $file_avg1 = $RESULTS.$NAMES[1]."_".$token."_".$bin_size."_spike_avg_bg.txt";
    $file_max1 = $RESULTS.$NAMES[1]."_".$token."_".$bin_size."_spike_max_bg.txt";
    $commandAVG = "join $file_avg0 $file_avg1 ";
    $commandMAX = "join $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
	$file_avg = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_spike_avg_bg.txt";
	$file_max = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_spike_max_bg.txt";
	$commandAVG = $commandAVG." | join - $file_avg ";
	$commandMAX = $commandMAX." | join - $file_max ";
    }

    $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_spike_avg_bg.txt";
    $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_spike_max_bg.txt";
    #
    $commandAVG = $commandAVG." | $field_notnull > $final_avg";
    print_mess("$commandAVG\n");
    system($commandAVG);
    $commandMAX = $commandMAX." | $field_notnull > $final_max";
    print_mess("$commandMAX\n");
    system($commandMAX);
    #
    SaveFile($final_avg);
    SaveFile($final_max);
    
    # Forth, sample values for bg
    # at least, two experiments
    $file_avg0 = $RESULTS.$NAMES[0]."_".$token."_".$bin_size."_sample_avg_bg.txt";
    $file_max0 = $RESULTS.$NAMES[0]."_".$token."_".$bin_size."_sample_max_bg.txt";
    $file_avg1 = $RESULTS.$NAMES[1]."_".$token."_".$bin_size."_sample_avg_bg.txt";
    $file_max1 = $RESULTS.$NAMES[1]."_".$token."_".$bin_size."_sample_max_bg.txt";
    $commandAVG = "join $file_avg0 $file_avg1 ";
    $commandMAX = "join $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
	$file_avg = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_sample_avg_bg.txt";
	$file_max = $RESULTS.$NAMES[$i]."_".$token."_".$bin_size."_sample_max_bg.txt";
	$commandAVG = $commandAVG." | join - $file_avg ";
	$commandMAX = $commandMAX." | join - $file_max ";
    }

    $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_sample_avg_bg.txt";
    $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$bin_size."_sample_max_bg.txt";
    #
    $commandAVG = $commandAVG." | $field_notnull > $final_avg";
    print_mess("$commandAVG\n");
    system($commandAVG);
    $commandMAX = $commandMAX." | $field_notnull > $final_max";
    print_mess("$commandMAX");
    system($commandMAX);
    #
    SaveFile($final_avg);
    SaveFile($final_max);
}

sub GenerateBoxplot
{
    my $experiment = $_[0];
    my $value = $_[1];
    my ($Rfile,$Routput_file);
    my $PDF_file;
    my $input_file;
    my $expression;
    my ($i,$j);
    my $midpoint;
    my @CLASSES;

    
    # (A) Boxplot with labelling
    # R code to generate the corresponding final boxplot
    $Rfile = $RSCRIPTS.join("-",@NAMES)."_".$bin_size."_".$experiment."_".$value."_boxplot.R";
    $PDF_file = $PLOTS.join("-",@NAMES)."_".$bin_size."_".$experiment."_".$value.".pdf";
    
    (open(RFILE,'>',$Rfile)) or print_error("R SCRIPT (boxplots): FILE $Rfile file can not be opened to write");
    print RFILE "pdf(\"$PDF_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$RAW_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c1<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c2<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$CHIPRX_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c3<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TAGREMOVAL_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c4<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c5<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$RAW_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c6<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c7<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$CHIPRX_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c8<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TAGREMOVAL_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c9<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c10<-read.table(\"$input_file\")\n";

    print RFILE "p = colorRampPalette($palette)\n";
    
    $expression = "boxplot(";
    for($i=0; $i<10; $i++)
    {
	for($j=0; $j<$n_experiments; $j++)
	{
	    $expression = $expression."log(c".($i+1)."[,".($j+2)."],base=2),";
	}
	$expression = $expression."\n";
    }

    print RFILE "$expression";
    print RFILE "xaxt=\"n\",xlab =\"\",lwd=1,\n";
    print RFILE "ylab=\"log2 ($value)\",outline=FALSE,\n";
    print RFILE "col=rev(p($n_experiments)),\n";
    print RFILE "notch=TRUE,\n";
    print RFILE "main=\"($experiment,$value) ".join(",",@NAMES)."\")\n";

    $expression = "c(";
    for($j=0; $j<$n_experiments-1; $j++)
    {
	$expression = $expression."\"".$NAMES[$j]."\",";
    }
    $expression = $expression."\"".$NAMES[$j]."\")";

    print RFILE "labels <- $expression\n";
    print RFILE "legend(\"topright\",labels,fill=rev(p($n_experiments)))\n";
    
    print RFILE "axis(1, labels = FALSE,tick=FALSE)\n";
    $midpoint = ($n_experiments*$N_STRATEGIES) + 0.5;
    print RFILE "lines(c($midpoint,$midpoint),c(10,-10),lwd=1,col=\"darkblue\")\n";
    #
    $midpoint = ($n_experiments*$N_STRATEGIES)/2;
    $expression = $midpoint + 0.5;
    print RFILE "text(x=$expression,y=0,\"PEAKS\",cex=3,col=\"blue\")\n";
    #
    $expression = ($n_experiments*$N_STRATEGIES) + $midpoint + 0.5;
    print RFILE "text(x=$expression,y=0,\"BG\",cex=3,col=\"blue\")\n";   
    #
    $midpoint = $n_experiments/2;
    @CLASSES=($RAW_TOKEN,$TRADITIONAL_TOKEN,$CHIPRX_TOKEN,$TAGREMOVAL_TOKEN,$SPIKCHIP_TOKEN);
    for($j=0; $j<$n_experiments*$N_STRATEGIES; $j++)
    {
	$expression = ($j*$n_experiments) + $midpoint + 0.5;
	print RFILE "text(x=$expression,y=2.5,\"".$CLASSES[($j % $N_STRATEGIES)]."\",srt=90,cex=0.5,col=\"blue\")\n";   
    }
    print RFILE "dev.off()\n";
    close(RFILE);
    #
    # execute R script
    $command = "R CMD BATCH $Rfile";
    print_mess("$command\n");
    system($command);
    #
    # error check in Rout file
    $Routput_file = join("-",@NAMES)."_".$bin_size."_".$experiment."_".$value."_boxplot.Rout";
    (open(ROUT,$Routput_file)) or print_error("R SCRIPTS (avg): FILE $Routput_file can not be opened");
    while($line=<ROUT>)
    {
	if ($line=~/Error/)
	{
	    print_error("R running: $line\n");
	}
    }
    close(ROUT);
    $command = "rm -f $Routput_file";
    print_mess("$command\n");
    system($command);

    # (B) Boxplot without labelling
    # R code to generate the corresponding final boxplot
    $Rfile = $RSCRIPTS.join("-",@NAMES)."_".$bin_size."_".$experiment."_".$value."_boxplot_".$NOLABELS_TOKEN.".R";
    $PDF_file = $PLOTS.join("-",@NAMES)."_".$bin_size."_".$experiment."_".$value."_".$NOLABELS_TOKEN.".pdf";
    
    (open(RFILE,'>',$Rfile)) or print_error("R SCRIPT (boxplots no labels): FILE $Rfile file can not be opened to write");
    print RFILE "pdf(\"$PDF_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$RAW_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c1<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c2<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$CHIPRX_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c3<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TAGREMOVAL_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c4<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_".$experiment."_".$value."_peaks.txt";
    print RFILE "c5<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$RAW_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c6<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c7<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$CHIPRX_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c8<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TAGREMOVAL_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c9<-read.table(\"$input_file\")\n";
    $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$bin_size."_".$experiment."_".$value."_bg.txt";
    print RFILE "c10<-read.table(\"$input_file\")\n";

    print RFILE "p = colorRampPalette($palette)\n";
    
    $expression = "boxplot(";
    for($i=0; $i<10; $i++)
    {
	for($j=0; $j<$n_experiments; $j++)
	{
	    $expression = $expression."log(c".($i+1)."[,".($j+2)."],base=2),";
	}
	$expression = $expression."\n";
    }

    print RFILE "$expression";
    print RFILE "xaxt=\"n\",yaxt=\"n\",xlab =\"\",lwd=1,\n";
    print RFILE "ylab=\"\",outline=FALSE,\n";
    print RFILE "col=rev(p($n_experiments)),\n";
    print RFILE "notch=TRUE)\n";

    print RFILE "axis(1, labels = FALSE,tick=FALSE)\n";
    print RFILE "axis(2, labels = FALSE,tick=TRUE)\n";
    
    print RFILE "dev.off()\n";
    close(RFILE);
    #
    # execute R script
    $command = "R CMD BATCH $Rfile";
    print_mess("$command\n");
    system($command);
    #
    # error check in Rout file
    $Routput_file = join("-",@NAMES)."_".$bin_size."_".$experiment."_".$value."_boxplot_".$NOLABELS_TOKEN.".Rout";
    (open(ROUT,$Routput_file)) or print_error("R SCRIPTS (avg): FILE $Routput_file can not be opened");
    while($line=<ROUT>)
    {
	if ($line=~/Error/)
	{
	    print_error("R running: $line\n");
	}
    }
    close(ROUT);
    $command = "rm -f $Routput_file";
    print_mess("$command");
    system($command);
}

sub CleanFile
{
    my $file = $_[0];
    my $command;
    
    # info for cleaning intermediate files (option -c)
    $command = "rm -f $file";
    push(@CLEAN_PROCEDURE,$command);
}

sub CleaningFiles
{
    my $i;
    my $command;
    
    
    if (exists($opt{c}))
    {
	for ($i=0; $i<scalar(@CLEAN_PROCEDURE); $i++)
	{
	    $command = $CLEAN_PROCEDURE[$i];
	    system($command);
	}
    }
}

sub SaveFile
{
    my $file = $_[0];
    my $command;
    
    # info for compressing final resulting files (after boxplot generation)
    $command = "gzip -f $file";
    push(@SAVE_PROCEDURE,$command);
}

sub SavingFiles
{
    my $i;
    my $command;
    
    
    for ($i=0; $i<scalar(@SAVE_PROCEDURE); $i++)
    {
	$command = $SAVE_PROCEDURE[$i];
	print_mess("$command\n");
	system($command);
    }
}
