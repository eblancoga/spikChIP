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
use Getopt::Long;
use Term::ANSIColor;

Getopt::Long::Configure("gnu_getopt", "auto_abbrev", "ignore_case_always");

#DEFINEs
my $PROGRAM = "spikChIP";
my $TRUE = 1;
my $FALSE = 0;
my $SUCCESS = 0;
my $DEFAULT_BIN_SIZE = 10000;
my $MIN_BIN_SIZE = 50;
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
my $DEFAULT_RESULTS = "results/";
my $PLOTS = "plots/";
my $RSCRIPTS = "Rscripts/";
my $NOLABELS_TOKEN = "nolabels";
my $N_STRATEGIES = 5;
my $DEFAULT_CHROMKEY = "FLY";

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
my $newdir;
my @CLEAN_PROCEDURE;
my @SAVE_PROCEDURE;

# --clean|-c:  remove intermediate files to reduce the size of the output folder
# --lessMillion|-l: allow the process of BAM files of < 1 Million reads
# --binsize|-b: bin size (default: 10000 bps)
# --chromkey | -k: Key word used to indicate the spike-in species in the chromInfo file (default: FLY)
# --palette|-p: palette (1 for reds, 2 for greens, 3 for blues and 0 for B&W)
# --help|-h: short help
# --verbose|-v: verbose option
# --overwrite|-w: overwrite existing result files
# --raw|-r: Perform the raw normalization
# --traditional|-t: Perform the traditional normalization
# --chiprx|-x: Perform the ChIPRx normalization
# --tagremoval|-g: Perform the tag removal normalization
# --spikchip|s: Perform the spikchip normalization with loess
# --outputfolder|-o: Path to the result folder (default: results/)

my $CLEAN = 0;
my $LESSMILLION = 0;
my $BIN_SIZE = $DEFAULT_BIN_SIZE;
my $PALETTE = $DEFAULT_PALETTE;
my $RESULTS = $DEFAULT_RESULTS;
my $HELP = 0;
my $VERBOSE = 0;
my $OVERWRITE = 0;
my $RAW = 0;
my $TRADITIONAL = 0;
my $CHIPRX = 0;
my $TAGREMOVAL = 0;
my $SPIKCHIP = 0;
my $CHROMKEY = $DEFAULT_CHROMKEY;

# 0.1 Acquire options
Getopt::Long::GetOptions(
    
    'clean|c' => \$CLEAN,
    'lessMillion|l' => \$LESSMILLION,
    'binsize|b=i' => \$BIN_SIZE,
    'chromkey|k=s' => \$CHROMKEY,
    'palette|p=i' => \$PALETTE,
    'help|h' => \$HELP,
    'verbose|v' => \$VERBOSE,
    'overwrite|w' => \$OVERWRITE,
    'raw|r' => \$RAW,
    'traditional|t' => \$TRADITIONAL,
    'chiprx|x' => \$CHIPRX,
    'tagremoval|g' => \$TAGREMOVAL,
    'spikchip|s' => \$SPIKCHIP,
    'outputfolder|o=s' => \$RESULTS,
);

print_mess("$PROGRAM.pl by Enrique Blanco @ CRG (2021)");
print_help();
print_mess("\n");

# 0.2 Save the starting time
$start = time();
$date = localtime();

if(!$RAW and !$TRADITIONAL and !$CHIPRX and !$TAGREMOVAL and !$SPIKCHIP){
        print_error("\t No normalization method was selected. Run $PROGRAM -h to see available options.")
}

# 0.3 Acquire options and arguments
print_mess("[$date] Stage 0.  Configuration of the pipeline\n");
print_mess("Reading arguments of the command line\n");
$n_files = $#ARGV+1;
($n_files == 2) or print_error("ERROR: Two arguments are required but $n_files are provided!\n$PROGRAM -vhcdk <bin_size_kbp> -p <0|1|2|3> <configuration_file> <chrominfo_file>\nPlease, type spikChIP -h for further assistance or spikChIP -v for verbose mode\n");
($configuration_filename,$chrominfo_file) = @ARGV;
print_mess("Syntax is correct and options are acquired");
print_ok();
print_mess("\n");

# 0.4 Checking binaries of SeqCode and samtools are available
print_mess("Checking for external software availability");
CheckExternalSoftware();
print_mess("Closing external software check...");
print_ok();
print_mess("\n");

# 0.5 Prepare the output folders
print_mess("Generating output folders (if necessary)");
# create the results/ folder
if($RESULTS ne "results/"){
    $newdir = $RESULTS."/results";
}else{
    $newdir = $RESULTS;
}
print_mess("Trying the $RESULTS directory");
mkdir($newdir) or print_mess("It is already existing");
# create the Rscripts/ folder
if($RESULTS ne "results/"){
    $newdir = $RESULTS."/".$RSCRIPTS;
}else{
    $newdir = $RSCRIPTS;
}
print_mess("Trying the $RSCRIPTS directory");
mkdir($newdir) or print_mess("It is already existing");
# create the plots/ folder
if($RESULTS ne "results/"){
    $newdir = $RESULTS."/".$PLOTS;
}else{
    $newdir = $PLOTS;
}
print_mess("Trying the $PLOTS directory");
mkdir($newdir) or print_mess("It is already existing");
print_ok();
print_mess("\n");

# 0.6 Check the chrominfo file 
print_mess("Checking information on chromosome sizes");
CheckChromInfoFile($chrominfo_file);
print_mess("Closing ChromInfo file check...");
print_ok();
print_mess("\n");

# 0.7 Checking bin size option (if selected)
print_mess("Establishing the bin size\n");
$BIN_SIZE = SettingBinsize();
print_mess("Effective bin size: $BIN_SIZE");
print_ok();
print_mess("\n");

# 0.8 Checking palette (if selected)
print_mess("Establishing the palette\n");
$PALETTE = SettingPalette();
print_mess("Effective palette: $PALETTE");
print_ok();
print_mess("\n");

# 0.9 Loading filenames from configuration file
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
print_mess("[$date] Stage 1.  Producing the segmentation of both genomes in bins ($BIN_SIZE bps)\n");

# spike segmentation
if($RESULTS ne "results/"){
    $spike_bins = $RESULTS."/results/".join("-",@NAMES)."_"."spike_".$BIN_SIZE.".bed";
}else{
    $spike_bins = join("-",@NAMES)."_"."spike_".$BIN_SIZE.".bed";
}

if(!(-e $spike_bins) or $OVERWRITE)
{
    $command = "grep $CHROMKEY $chrominfo_file | gawk 'BEGIN{OFS=\"\\t\";offset=$BIN_SIZE;}{for(i=1;i<\$2-offset;i=i+offset) print \$1,i,i+offset;}' > $spike_bins";
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
if($RESULTS ne "results/"){
    $sample_bins = $RESULTS."/results/".join("-",@NAMES)."_"."sample_".$BIN_SIZE.".bed";
}else{
    $sample_bins = $RESULTS.join("-",@NAMES)."_"."sample_".$BIN_SIZE.".bed";
}

if(!(-e $sample_bins) or $OVERWRITE)
{
    $command = "grep -v $CHROMKEY $chrominfo_file | gawk 'BEGIN{OFS=\"\\t\";offset=$BIN_SIZE;}{for(i=1;i<\$2-offset;i=i+offset) print \$1,i,i+offset;}' > $sample_bins";
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
    if($RAW)
    {
        print_mess("Starting RAW normalization...\n");
        NormalizationRaw($i);
        print_ok();
        print_mess("\n");
    }
    #
    if($TRADITIONAL){
        print_mess("Starting TRADITIONAL normalization...\n");
        NormalizationTraditional($i);
        print_ok();
        print_mess("\n");
    }
    #
    if($CHIPRX){
        print_mess("Starting CHIP-RX normalization...\n");
        NormalizationChIPRX($i);
        print_ok();
        print_mess("\n");
    }
    #
    if($TAGREMOVAL){
        print_mess("Starting TAG_REMOVAL normalization...\n");
        NormalizationTagRemoval($i);
        print_ok();
        print_mess("\n");
    }

    if(!$RAW and !$TRADITIONAL and !$CHIPRX and !$TAGREMOVAL){
        print_mess("\t None of the following normalization were selected: raw, traditional, chiprx, or tag removal.")
    }
}

# join the values of all experiments in a single file per normalization class
if($RAW)
{
    print_mess("Join raw data values from all experiments\n");
    JoinNormValues($RAW_TOKEN);
    print_mess("\n");
}
if($TRADITIONAL){
    print_mess("Join traditional data values from all experiments\n");
    JoinNormValues($TRADITIONAL_TOKEN);
    print_mess("\n");
}
if($CHIPRX){
    print_mess("Join ChIP-RX data values from all experiments\n");
    JoinNormValues($CHIPRX_TOKEN);
    print_mess("\n");
}
if($TAGREMOVAL)
{
    print_mess("Join tag removal data values from all experiments\n");
    JoinNormValues($TAGREMOVAL_TOKEN);
    print_mess("\n");
}

if($SPIKCHIP)
{
    print_mess("Starting $PROGRAM local normalization...\n");
    print_mess("Prepare $PROGRAM data values for local normalization\n");
    PreparespikChIPValues();
    print_mess("Run local regression adjusting sample values with spike values\n");
    RunspikChIPValues();
    print_ok();
}

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
if($RAW){
    print_mess("Working with $RAW_TOKEN values\n");
    ClassifyNormalizationValues($RAW_TOKEN);
    print_ok();
    print_mess("\n");
}
if($TRADITIONAL){
    print_mess("Working with $TRADITIONAL_TOKEN values\n");
    ClassifyNormalizationValues($TRADITIONAL_TOKEN);
    print_ok();
    print_mess("\n");
}
if($CHIPRX){
    print_mess("Working with $CHIPRX_TOKEN values\n");
    ClassifyNormalizationValues($CHIPRX_TOKEN);
    print_ok();
    print_mess("\n");
}
if($TAGREMOVAL){
    print_mess("Working with $TAGREMOVAL_TOKEN values\n");
    ClassifyNormalizationValues($TAGREMOVAL_TOKEN);
    print_ok();
    print_mess("\n");
}
if($SPIKCHIP){
    print_mess("Working with $SPIKCHIP_TOKEN values\n");
    ClassifyNormalizationValues($SPIKCHIP_TOKEN);
    print_ok();
}
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
    if ($HELP)
    {
        print STDERR color("bold blue"),"spikChIP_v1.0\t\t\t\tUser commands\t\t\t\tspikChIP\n\n";
        print STDERR color("bold blue"),"NAME\n\tspikChIP, a tool to normalize ChIP-seq experiments by spike-in correction\n\n";
        print STDERR color("bold blue"),"SYNOPSIS:\n\t$PROGRAM -rtxgs -vclw -b <bin_size_kbp> -p <0|1|2|3> <configuration_file> <chrominfo_file>\n\n";
        print STDERR color("bold blue"),"OPTIONS:\n";
        print STDERR color("bold blue"),"\t--clean|-c:  remove intermediate files to reduce the size of the output folder\n";
        print STDERR color("bold blue"),"\t--lessMillion|-l: allow the process of BAM files of < 1 Million reads\n";
        print STDERR color("bold blue"),"\t--binsize|-b: bin size (default: 10000 bps)\n";
        print STDERR color("bold blue"),"\t--palette|-p: palette (1 for reds, 2 for greens, 3 for blues and 0 for B&W)\n";
        print STDERR color("bold blue"),"\t--help|-h: short help\n";
        print STDERR color("bold blue"),"\t--verbose|-v: verbose option\n";
        print STDERR color("bold blue"),"\t--overwrite|-w: overwrite existing result files\n";
        print STDERR color("bold blue"),"\t--raw|-r: Perform the raw normalization\n";
        print STDERR color("bold blue"),"\t--traditional|-t: Perform the traditional normalization\n";
        print STDERR color("bold blue"),"\t--chiprx|-x: Perform the ChIPRx normalization\n";
        print STDERR color("bold blue"),"\t--tagremoval|-g: Perform the tag removal normalization\n";
        print STDERR color("bold blue"),"\t--spikchip|s: Perform the spikchip normalization with loess\n";
        print STDERR color("bold blue"),"\t--outputfolder|-o: Path to the result folder (default: results/)\n";
	exit(0);
    }
}

sub print_mess
{
        my @mess = @_;
        print STDERR color("bold green"),"%%%% @mess\n" if $VERBOSE;
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
    if($RESULTS ne "results/"){
        $out_file = $RESULTS."/results/".$name."_flagstat.txt";
    }else{
        $out_file = $RESULTS.$name."_flagstat.txt";
    }
    
    if(!(-e $out_file) or $OVERWRITE)
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

	    if ($n_reads < $MEGA && !$LESSMILLION)
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
	if(!(-e $output_file) or $OVERWRITE)
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
	if(!(-e $output_file) or $OVERWRITE)
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
    # range of valid bin sizes
    if (($BIN_SIZE < $MIN_BIN_SIZE) || ($BIN_SIZE > $MAX_BIN_SIZE))
    {
	print_error("Please, set the bin size in the valid range (current $BIN_SIZE): $MIN_BIN_SIZE - $MAX_BIN_SIZE");
    }

    return($BIN_SIZE);
}

sub SettingPalette
{

    # range of valid paletter
    if (($PALETTE < $MIN_PALETTE) || ($PALETTE > $MAX_PALETTE))
    {
	print_error("Please, set the palette in the valid range (current $PALETTE): $MIN_PALETTE - $MAX_PALETTE");
    }

    if ($PALETTE == 1)
    {
	$PALETTE = "c(\"pink\",\"red\",\"brown\")";
    }
    else
    {
	if ($PALETTE == 2)
	{
	    $PALETTE = "c(\"darkolivegreen1\",\"green4\",\"darkgreen\")";
	}
	else
	{
	    if ($PALETTE == 3)
	    {
		$PALETTE = "c(\"lightskyblue\",\"blue\",\"darkblue\")";
	    }
	    else
	    {
		$PALETTE = "c(\"white\")";
	    }
	}
    }
    
    return($PALETTE);
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
    my $prefix;
    my ($file_all,$file_avg,$file_max);

    # input bam files
    $bam_sample = $BAM_SAMPLES[$i];
    $bam_spike = $BAM_SPIKES[$i];

    # Spike bins
    $out_name = $NAMES[$i]."_".$RAW_TOKEN."_".$BIN_SIZE."_spike";
    if($RESULTS ne "results/"){
        $file_all = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }else{
        $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }
    
    if(!(-e $file_all) or $OVERWRITE)
    {
        if ($LESSMILLION)
        {
            if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -dns $MEGA -x $prefix $chrominfo_file $bam_spike $spike_bins $out_name";
            }else{
                $command = "recoverChIPlevels -dns $MEGA $chrominfo_file $bam_spike $spike_bins $out_name";
            }
        }
        else
        {
            if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -ns $MEGA -x $prefix $chrominfo_file $bam_spike $spike_bins $out_name";
            }else{
                $command = "recoverChIPlevels -ns $MEGA $chrominfo_file $bam_spike $spike_bins $out_name";
            }
        }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file", $file_all, "already exist. Skipping the spike Raw Normalization\n");
    }
    #
    if($RESULTS ne "results/"){
        $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }else{
        $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }
    
    
    if(!(-e $file_avg) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_avg, "already exist. Skipping creation of the spike avg\n");
    }

    if(!(-e $file_max) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_max, "already exist. Skipping creation of the spike max\n");
    }

    # Sample bins
    $out_name = $NAMES[$i]."_".$RAW_TOKEN."_".$BIN_SIZE."_sample";
    if($RESULTS ne "results/"){
        $file_all = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }else{
        $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }
    if(!(-e $file_all) or $OVERWRITE)
    {
        if ($LESSMILLION)
        {
            if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -dns $MEGA -x $prefix $chrominfo_file $bam_sample $sample_bins $out_name";
            }else{
                $command = "recoverChIPlevels -dns $MEGA $chrominfo_file $bam_sample $sample_bins $out_name";
            }
        }
        else
        {
            if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -ns $MEGA -x $prefix $chrominfo_file $bam_sample $sample_bins $out_name";
            }else{
                $command = "recoverChIPlevels -ns $MEGA $chrominfo_file $bam_sample $sample_bins $out_name";
            }
        }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file", $out_name, "already exist. Skipping the sample Raw Normalization\n");
    }

    #
    if($RESULTS ne "results/"){
        $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }else{
        $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }
    
    if(!(-e $file_avg) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_avg, "already exist. Skipping creation of the sample avg\n");
    }

    if(!(-e $file_max) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command");
        system($command);
    }else{
        print_mess("\t The file", $file_max, "already exist. Skipping creation of the sample max\n");
    }
}

sub NormalizationTraditional
{
    my $i = $_[0];
    my ($bam_sample,$bam_spike);
    my $prefix;
    my $out_name;
    my ($file_all,$file_avg,$file_max);
    my $total_reads;
    
    
    # input bam files
    $bam_sample = $BAM_SAMPLES[$i];
    $bam_spike = $BAM_SPIKES[$i];
    $total_reads = $READS_SAMPLES[$i] + $READS_SPIKES[$i];

    # Spike bins
    $out_name = $NAMES[$i]."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_spike";
    if($RESULTS ne "results/"){
        $file_all = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }else{
        $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }
    if(!(-e $file_all) or $OVERWRITE)
    {
        if ($LESSMILLION)
        {
            if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -dns $total_reads -x $prefix $chrominfo_file $bam_spike $spike_bins $out_name";
            }else{
                $command = "recoverChIPlevels -dns $total_reads $chrominfo_file $bam_spike $spike_bins $out_name";
            }
        }
        else
        {
            if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -ns $total_reads -x $prefix $chrominfo_file $bam_spike $spike_bins $out_name";
            }else{
                $command = "recoverChIPlevels -ns $total_reads $chrominfo_file $bam_spike $spike_bins $out_name";
            }
        }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file", $file_all, "already exist. Skipping the spike Traditional Normalization\n");
    }
    
    #
    if($RESULTS ne "results/"){
        $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }else{
        $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }
    
    if(!(-e $file_avg) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_avg, "already exist. Skipping creation of the spike avg\n");
    }

    if(!(-e $file_max) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_max, "already exist. Skipping creation of the spike max\n");
    }

    # Sample bins
    $out_name = $NAMES[$i]."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_sample";
    if($RESULTS ne "results/"){
        $file_all = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }else{
        $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }
    if(!(-e $file_all) or $OVERWRITE)
    {
      if ($LESSMILLION)
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -dns ".$READS_SAMPLES[$i]." -x $prefix $chrominfo_file $bam_sample $sample_bins $out_name";
            }else{
                $command = "recoverChIPlevels -dns ".$READS_SAMPLES[$i]." $chrominfo_file $bam_sample $sample_bins $out_name";
            }
      }
      else
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -ns ".$READS_SAMPLES[$i]." -x $prefix $chrominfo_file $bam_sample $sample_bins $out_name";
            }else{
                $command = "recoverChIPlevels -ns ".$READS_SAMPLES[$i]." $chrominfo_file $bam_sample $sample_bins $out_name";
            }
      }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file", $file_all, "already exist. Skipping the sample Traditional Normalization\n");
    }
    
    #
    if($RESULTS ne "results/"){
        $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }else{
        $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }
    
    if(!(-e $file_avg) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_avg, "already exist. Skipping creation of the sample avg\n");
    }

    if(!(-e $file_max) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command");
        system($command);
    }else{
        print_mess("\t The file", $file_max, "already exist. Skipping creation of the sample max\n");
    }
}

sub NormalizationChIPRX
{
    my $i = $_[0];
    my ($bam_sample,$bam_spike);
    my $out_name;
    my $prefix;
    my ($file_all,$file_avg,$file_max);
    my $total_reads;
    
    
    # input bam files
    $bam_sample = $BAM_SAMPLES[$i];
    $bam_spike = $BAM_SPIKES[$i];
    $total_reads = $READS_SAMPLES[$i] + $READS_SPIKES[$i];

    # Spike bins
    $out_name = $NAMES[$i]."_".$CHIPRX_TOKEN."_".$BIN_SIZE."_spike";
    if($RESULTS ne "results/"){
        $file_all = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }else{
        $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }
    if(!(-e $file_all) or $OVERWRITE)
    {
      if ($LESSMILLION)
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -dns ".$READS_SPIKES[$i]." -x $prefix $chrominfo_file $bam_spike $spike_bins $out_name";
            }else{
                $command = "recoverChIPlevels -dns ".$READS_SPIKES[$i]." $chrominfo_file $bam_spike $spike_bins $out_name";
            }
      }
      else
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -ns ".$READS_SPIKES[$i]." -x $prefix $chrominfo_file $bam_spike $spike_bins $out_name";
            }else{
                $command = "recoverChIPlevels -ns ".$READS_SPIKES[$i]." $chrominfo_file $bam_spike $spike_bins $out_name";
            }
      }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file", $file_all, "already exist. Skipping the spike ChIPRX Normalization\n");
    }
    #
    if($RESULTS ne "results/"){
        $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }else{
        $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }
    if(!(-e $file_avg) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_avg, "already exist. Skipping creation of the spike avg\n");
    }
    if(!(-e $file_max) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_max, "already exist. Skipping creation of the spike max\n");
    }

    # Sample bins
    $out_name = $NAMES[$i]."_".$CHIPRX_TOKEN."_".$BIN_SIZE."_sample";
    if($RESULTS ne "results/"){
        $file_all = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }else{
        $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }
    if(!(-e $file_all) or $OVERWRITE)
    {
      if ($LESSMILLION)
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -dns ".$READS_SPIKES[$i]." -x $prefix $chrominfo_file $bam_sample $sample_bins $out_name";
            }else{
                $command = "recoverChIPlevels -dns ".$READS_SPIKES[$i]." $chrominfo_file $bam_sample $sample_bins $out_name";
            }
      }
      else
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -ns ".$READS_SPIKES[$i]." -x $prefix $chrominfo_file $bam_sample $sample_bins $out_name";
            }else{
                $command = "recoverChIPlevels -ns ".$READS_SPIKES[$i]." $chrominfo_file $bam_sample $sample_bins $out_name";
            }
      }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file", $file_all, "already exist. Skipping the sample ChIPRX Normalization\n");
    }
    
    #
    if($RESULTS ne "results/"){
        $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }else{
        $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }
    if(!(-e $file_avg) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_avg, "already exist. Skipping creation of the sample avg\n");
    }
    if(!(-e $file_max) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_max, "already exist. Skipping creation of the sample max\n");
    }
}

sub NormalizationTagRemoval
{
    my $i = $_[0];
    my ($bam_sample,$bam_spike);
    my $out_name;
    my $prefix;
    my ($file_all,$file_avg,$file_max);

    
    # input bam files
    $bam_sample = $BAM_SAMPLES_DOWN[$i];
    $bam_spike = $BAM_SPIKES_DOWN[$i];

    # Spike bins
    $out_name = $NAMES[$i]."_".$TAGREMOVAL_TOKEN."_".$BIN_SIZE."_spike";
    if($RESULTS ne "results/"){
        $file_all = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }else{
        $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }
    if(!(-e $file_all) or $OVERWRITE)
    {
      if ($LESSMILLION)
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -dns $MEGA -x $prefix $chrominfo_file $bam_spike $spike_bins $out_name";
            }else{
                $command = "recoverChIPlevels -dns $MEGA $chrominfo_file $bam_spike $spike_bins $out_name";
            }
      }
      else
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -ns $MEGA -x $prefix $chrominfo_file $bam_spike $spike_bins $out_name";
            }else{
                $command = "recoverChIPlevels -ns $MEGA $chrominfo_file $bam_spike $spike_bins $out_name";
            }
      }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file", $file_all, "already exist. Skipping the spike TagRemoval Normalization\n");
    }
    #
    if($RESULTS ne "results/"){
        $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }else{
        $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }
    if(!(-e $file_avg) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_avg, "already exist. Skipping creation of the spike avg\n");
    }
    if(!(-e $file_max) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_max, "already exist. Skipping creation of the spike max\n");
    }
    
    # Sample bins
    $out_name = $NAMES[$i]."_".$TAGREMOVAL_TOKEN."_".$BIN_SIZE."_sample";
    if($RESULTS ne "results/"){
        $file_all = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }else{
        $file_all = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name".".bed";
    }
    if(!(-e $file_all) or $OVERWRITE)
    {
      if ($LESSMILLION)
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -dns $MEGA -x $prefix $chrominfo_file $bam_sample $sample_bins $out_name";
            }else{
                $command = "recoverChIPlevels -dns $MEGA $chrominfo_file $bam_sample $sample_bins $out_name";
            }
      }
      else
      {
          if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "recoverChIPlevels -ns $MEGA -x $prefix $chrominfo_file $bam_sample $sample_bins $out_name";
            }else{
                $command = "recoverChIPlevels -ns $MEGA $chrominfo_file $bam_sample $sample_bins $out_name";
            }
      }
    print_mess("$command\n");
    system($command);
    }else{
        print_mess("\t The file", $file_all, "already exist. Skipping the sample TagRemoval Normalization\n");
    }
    #
    if($RESULTS ne "results/"){
        $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }else{
        $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
        $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
    }
    if(!(-e $file_avg) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$5}' $file_all | sort > $file_avg";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_avg, "already exist. Skipping creation of the sample avg\n");
    }
    if(!(-e $file_max) or $OVERWRITE)
    {
        $command = "gawk '{print \$1\"*\"\$2\"*\"\$3,\$6}' $file_all | sort > $file_max";
        print_mess("$command\n");
        system($command);
    }else{
        print_mess("\t The file", $file_max, "already exist. Skipping creation of the sample max\n");
    }
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
    $out_name0 = $NAMES[0]."_".$token."_".$BIN_SIZE."_spike";
    if($RESULTS ne "results/"){
        $file_avg0 = $RESULTS."/"."$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_avg.bed";
        $file_max0 = $RESULTS."/"."$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_max.bed";
    }else{
        $file_avg0 = "$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_avg.bed";
        $file_max0 = "$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_max.bed";
    }
    $out_name1 = $NAMES[1]."_".$token."_".$BIN_SIZE."_spike";
    if($RESULTS ne "results/"){
        $file_avg1 = $RESULTS."/"."$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_avg.bed";
        $file_max1 = $RESULTS."/"."$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_max.bed";
    }else{
        $file_avg1 = "$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_avg.bed";
        $file_max1 = "$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_max.bed";
    }
    $commandAVG = "JoinNFiles.pl -v $file_avg0 $file_avg1 ";
    $commandMAX = "JoinNFiles.pl -v $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
        $out_name = $NAMES[$i]."_".$token."_".$BIN_SIZE."_spike";
        if($RESULTS ne "results/"){
            $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
            $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
        }else{
            $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
            $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
        }
        $commandAVG = $commandAVG." $file_avg ";
        $commandMAX = $commandMAX." $file_max ";
    }
    
    if($RESULTS ne "results/"){
        $final_avg = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_avg.txt";
        $final_max = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_max.txt";
    }else{
        $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_avg.txt";
        $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_max.txt";
    }
    #
    if(!(-e $final_avg) or $OVERWRITE)
    {
        $commandAVG = $commandAVG." | sort -k 1,1 > $final_avg";
        print_mess("$commandAVG\n");
        system($commandAVG);
    }else{
        print_mess("\t The file", $final_avg, "already exist. Skipping joining the spike files for avg\n");
    }
    if(!(-e $final_max) or $OVERWRITE)
    {
        $commandMAX = $commandMAX." | sort -k 1,1 > $final_max";
        print_mess("$commandMAX\n");
        system($commandMAX);
    }else{
        print_mess("\t The file", $final_max, "already exist. Skipping joining the spike files for max\n");
    }
    #
    SaveFile($final_avg);
    SaveFile($final_max);


    # Second, sample values
    # at least, two experiments
    $out_name0 = $NAMES[0]."_".$token."_".$BIN_SIZE."_sample";
    if($RESULTS ne "results/"){
        $file_avg0 = $RESULTS."/"."$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_avg.bed";
        $file_max0 = $RESULTS."/"."$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_max.bed";
    }else{
        $file_avg0 = "$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_avg.bed";
        $file_max0 = "$out_name0"."_recoverChIPlevels/PEAKsignal_"."$out_name0"."_max.bed";
    }
    $out_name1 = $NAMES[1]."_".$token."_".$BIN_SIZE."_sample";
    if($RESULTS ne "results/"){
        $file_avg1 = $RESULTS."/"."$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_avg.bed";
        $file_max1 = $RESULTS."/"."$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_max.bed";
    }else{
        $file_avg1 = "$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_avg.bed";
        $file_max1 = "$out_name1"."_recoverChIPlevels/PEAKsignal_"."$out_name1"."_max.bed";
    }
    $commandAVG = "JoinNFiles.pl -v $file_avg0 $file_avg1 ";
    $commandMAX = "JoinNFiles.pl -v $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
        $out_name = $NAMES[$i]."_".$token."_".$BIN_SIZE."_sample";
        if($RESULTS ne "results/"){
            $file_avg = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
            $file_max = $RESULTS."/"."$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
        }else{
            $file_avg = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_avg.bed";
            $file_max = "$out_name"."_recoverChIPlevels/PEAKsignal_"."$out_name"."_max.bed";
        }
        $commandAVG = $commandAVG." $file_avg ";
        $commandMAX = $commandMAX." $file_max ";
    }
    if($RESULTS ne "results/"){
        $final_avg = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_avg.txt";
        $final_max = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_max.txt";
    }else{
        $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_avg.txt";
        $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_max.txt";
    }
    #
    if(!(-e $final_avg) or $OVERWRITE)
    {
        $commandAVG = $commandAVG." | sort -k 1,1 > $final_avg";
        print_mess("$commandAVG\n");
        system($commandAVG);
    }else{
        print_mess("\t The file", $final_avg, "already exist. Skipping joining the sample files for avg\n");
    }
    if(!(-e $final_max) or $OVERWRITE)
    {
        $commandMAX = $commandMAX." | sort -k 1,1 > $final_max";
        print_mess("$commandMAX\n");
        system($commandMAX);
    }else{
        print_mess("\t The file", $final_max, "already exist. Skipping joining the sample files for max\n");
    }
    #
    SaveFile($final_avg);
    SaveFile($final_max);
    
    # remove the seqcode folders (spike and sample)
    for($i=0; $i<$n_experiments; $i++)
    {
        if($RESULTS ne "results/"){
            $folder = $RESULTS."/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_recoverChIPlevels/";
            CleanFolder($folder);
            $folder = $RESULTS."/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_recoverChIPlevels/";
            CleanFolder($folder);
        }else{
            $folder = $NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_recoverChIPlevels/";
            CleanFolder($folder);
            $folder = $NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_recoverChIPlevels/";
            CleanFolder($folder);
        }
    }
}


sub PreparespikChIPValues
{
    my ($spike_values,$sample_values);
    my ($output_file1,$output_file2,$output_file3,$output_file4);
    my $command;
    my ($i,$fields);
    

    # (A) normalization based in avg values
    if($RESULTS ne "results/"){
        $spike_values = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_spike_avg.txt";
        $sample_values = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_sample_avg.txt";
        # output files
        $output_file1 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg.txt";
        $output_file2 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_values.txt";
        $output_file3 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_names.txt";
        $output_file4 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_names2.txt";
    }else{
        $spike_values = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_spike_avg.txt";
        $sample_values = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_sample_avg.txt";
        # output files
        $output_file1 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg.txt";
        $output_file2 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_values.txt";
        $output_file3 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_names.txt";
        $output_file4 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_names2.txt";
    }
    #
    CleanFile($output_file1);
    CleanFile($output_file2);
    CleanFile($output_file3);
    CleanFile($output_file4);
    
    if(-e $output_file1 and -e $output_file2 and -e $output_file3 and -e $output_file4){
        print_mess("\t The files are already prepared for spike\n")
    }
    
    if(!(-e $output_file1) or $OVERWRITE)
    {
        $command = "cat $spike_values $sample_values > $output_file1";
        print_mess("$command\n");
        system($command);
    }
    
    # extract all the columns with values (one per experiment)
    $fields = "";
    for($i=0; $i<$n_experiments-1; $i++)
    {
	$fields = $fields." \$".($i+2).",";
    }
    $fields = $fields." \$".($i+2);
    
    if(!(-e $output_file2) or $OVERWRITE)
    {
        $command = "gawk 'BEGIN{OFS=\"\\t\"}{print $fields}' $output_file1 > $output_file2";
        print_mess("$command\n");
        system($command);
    }

    if(!(-e $output_file3) or $OVERWRITE)
    {
        $command = "gawk 'BEGIN{OFS=\"\\t\"}{print NR}' $output_file1 > $output_file3";
        print_mess("$command\n");
        system($command);
    }

    if(!(-e $output_file4) or $OVERWRITE)
    {
        $command = "gawk 'BEGIN{OFS=\"\\t\"}{print NR,\$1}' $output_file1 > $output_file4";
        print_mess("$command\n");
        system($command);
    }

    # (B) normalization based in max values
    if($RESULTS ne "results/"){
        $spike_values = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_spike_max.txt";
        $sample_values = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_sample_max.txt";
        # output files
        $output_file1 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max.txt";
        $output_file2 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_values.txt";
        $output_file3 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_names.txt";
        $output_file4 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_names2.txt";
    }else{
        $spike_values = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_spike_max.txt";
        $sample_values = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$TRADITIONAL_TOKEN."_".$BIN_SIZE."_sample_max.txt";
        # output files
        $output_file1 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max.txt";
        $output_file2 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_values.txt";
        $output_file3 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_names.txt";
        $output_file4 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_names2.txt";
    }
    #
    CleanFile($output_file1);
    CleanFile($output_file2);
    CleanFile($output_file3);
    CleanFile($output_file4);
    
    if(-e $output_file1 and -e $output_file2 and -e $output_file3 and -e $output_file4){
        print_mess("\t The files are already prepared for sample\n")
    }

    if(!(-e $output_file1) or $OVERWRITE)
    {
        $command = "cat $spike_values $sample_values > $output_file1";
        print_mess("$command\n");
        system($command);
    }

    if(!(-e $output_file2) or $OVERWRITE)
    {
        $command = "gawk 'BEGIN{OFS=\"\\t\"}{print $fields}' $output_file1 > $output_file2";
        print_mess("$command\n");
        system($command);
    }

    if(!(-e $output_file3) or $OVERWRITE)
    {
        $command = "gawk 'BEGIN{OFS=\"\\t\"}{print NR}' $output_file1 > $output_file3";
        print_mess("$command\n");
        system($command);
    }

    if(!(-e $output_file4) or $OVERWRITE)
    {
        $command = "gawk 'BEGIN{OFS=\"\\t\"}{print NR,\$1}' $output_file1 > $output_file4";
        print_mess("$command\n");
        system($command);
    }
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
    if($RESULTS ne "results/"){
        $Rfile = $RESULTS."/".$RSCRIPTS.join("-",@NAMES)."_".$BIN_SIZE."_avg.R";
        $output_file1 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg.txt";
        $output_file2 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_values.txt";
        $output_file3 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_names.txt";
        $output_file4 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_names2.txt";
        #
        $final_avg = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_avg_normalized.txt";
        $final_avg_spike = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_avg_normalized_spike.txt";
        $final_avg_sample = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_avg_normalized_sample.txt";
    }else{
        $Rfile = $RSCRIPTS.join("-",@NAMES)."_".$BIN_SIZE."_avg.R";
        $output_file1 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg.txt";
        $output_file2 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_values.txt";
        $output_file3 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_names.txt";
        $output_file4 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_avg_names2.txt";
        #
        $final_avg = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_avg_normalized.txt";
        $final_avg_spike = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_avg_normalized_spike.txt";
        $final_avg_sample = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_avg_normalized_sample.txt";
    }
    CleanFile($final_avg);
    SaveFile($final_avg_spike);
    SaveFile($final_avg_sample);
    # R code to perform the loess regression on the whole set of bins in all conditions
    if(!(-e $Rfile) or $OVERWRITE)
    {
        #
        #
        $n_total_bins = $n_spike_bins + $n_sample_bins;
        #
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
        if($RESULTS ne "results/"){
            $Routput_file = $RESULTS."/".join("-",@NAMES)."_".$BIN_SIZE."_avg.Rout";
        }else{
            $Routput_file = join("-",@NAMES)."_".$BIN_SIZE."_avg.Rout";
        }
        
        # execute R script
        $command = "R CMD BATCH $Rfile";
        print_mess("$command\n");
        system($command);
        # error check in Rout file
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
        if(!(-e $final_avg_spike) or $OVERWRITE)
        {
            $command = "join $final_avg $output_file4 | grep $CHROMKEY | gawk 'BEGIN{OFS=\"\\t\"}{print $binname,$fields}'> $final_avg_spike";
            print_mess("$command\n");
            system($command);
        }
        if(!(-e $final_avg_sample) or $OVERWRITE)
        {
            $command = "join $final_avg $output_file4 | grep -v $CHROMKEY | gawk 'BEGIN{OFS=\"\\t\"}{print $binname,$fields}'> $final_avg_sample";
            print_mess("$command\n");
            system($command);
        }
    }else{
        print_mess("The Rscript", $Rfile, "already exists");
    }

    # spikChIP on max values
    print_mess("Performing the analysis on max values\n");
    if($RESULTS ne "results/"){
        $Rfile = $RESULTS."/".$RSCRIPTS.join("-",@NAMES)."_".$BIN_SIZE."_max.R";
        $output_file1 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max.txt";
        $output_file2 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_values.txt";
        $output_file3 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_names.txt";
        $output_file4 = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_names2.txt";
        #
        $final_max = $RESULTS."/results/".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_max_normalized.txt";
        $final_max_spike = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_max_normalized_spike.txt";
        $final_max_sample = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_max_normalized_sample.txt";
    }else{
        $Rfile = $RSCRIPTS.join("-",@NAMES)."_".$BIN_SIZE."_max.R";
        $output_file1 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max.txt";
        $output_file2 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_values.txt";
        $output_file3 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_names.txt";
        $output_file4 = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_spike-sample_max_names2.txt";
        #
        $final_max = $RESULTS.join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_max_normalized.txt";
        $final_max_spike = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_max_normalized_spike.txt";
        $final_max_sample = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$SPIKCHIP_TOKEN."_".$BIN_SIZE."_max_normalized_sample.txt";
    }
    #
    CleanFile($final_max);
    SaveFile($final_max_spike);
    SaveFile($final_max_sample);
    if(!(-e $Rfile) or $OVERWRITE)
    {
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
        if($RESULTS ne "results/"){
            $Routput_file = $RESULTS."/".join("-",@NAMES)."_".$BIN_SIZE."_max.Rout";
        }else{
            $Routput_file = join("-",@NAMES)."_".$BIN_SIZE."_max.Rout";
        }
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
        if(!(-e $final_max_spike) or $OVERWRITE)
        {
            # distinguish spike from sample bins
            $command = "join $final_max $output_file4 | grep $CHROMKEY | gawk 'BEGIN{OFS=\"\\t\"}{print $binname,$fields}'> $final_max_spike";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $final_max_spike already exists");
        }
        if(!(-e $final_max_sample) or $OVERWRITE)
        {
            $command = "join $final_max $output_file4 | grep -v $CHROMKEY | gawk 'BEGIN{OFS=\"\\t\"}{print $binname,$fields}'> $final_max_sample";
            print_mess("$command");
            system($command);
        }else{
            print_mess("\t $final_max_sample already exists");
        }
    }else{
        print_mess("The Rscript", $Rfile, "already exists");
    }
}


sub ClassifyBins
{
    my $i;
    my $command;
    my $prefix;
    my ($folder,$out_name,$out_common,$out_peaks_file);
    
    
    # for each experiment, compare spike/sample genome bins to peaks in the corresponding spike/genome
    for ($i=0; $i<scalar(@NAMES); $i++)
    {
    	print_mess("Working with sample $NAMES[$i]: peaks Vs. bins of spike\n");
    	
    	# use SeqCode to identify the spike bins overlapping with spike peaks
    	if($RESULTS ne "results/"){
            $out_name = $NAMES[$i]."_".$BIN_SIZE."_spike_bins";
            $folder = $PEAKS_TOKEN."_".$out_name."_matchpeaks";
            $out_common = $RESULTS."/".$folder."/common_".$PEAKS_TOKEN."_".$out_name.".bed";
            $out_peaks_file = $RESULTS."/results/".$NAMES[$i]."_".$BIN_SIZE."_bins_spike_peaks.bed";
            $BINS_PEAKS_SPIKE[$i] = $out_peaks_file;
        }else{
            $out_name = $NAMES[$i]."_".$BIN_SIZE."_spike_bins";
            $folder = $PEAKS_TOKEN."_".$out_name."_matchpeaks";
            $out_common = $folder."/common_".$PEAKS_TOKEN."_".$out_name.".bed";
            $out_peaks_file = $RESULTS.$NAMES[$i]."_".$BIN_SIZE."_bins_spike_peaks.bed";
            $BINS_PEAKS_SPIKE[$i] = $out_peaks_file;
        }
    	#
    	CleanFile($out_peaks_file);

    	if(!(-e $out_peaks_file) or $OVERWRITE)
        {
            if($RESULTS ne "results/"){
                $prefix = $RESULTS."/";
                $command = "matchpeaks -v -x $prefix $PEAKS_SPIKES[$i] $spike_bins $PEAKS_TOKEN $out_name";
            }else{
                $command = "matchpeaks -v $PEAKS_SPIKES[$i] $spike_bins $PEAKS_TOKEN $out_name";
            }
        	print_mess("$command\n");
        	system($command);
        	$command = "grep spike_bins $out_common > $out_peaks_file";
        	print_mess("$command\n");
        	system($command);
        	#
    	}
        CleanFolder($folder);
        print_mess("\n");
    	#
    	print_mess("Working with sample $NAMES[$i]: peaks Vs. bins of sample\n");
    	
    	# use SeqCode to identify the sample bins overlapping with sample peaks
    	if($RESULTS ne "results/"){
            $out_name = $NAMES[$i]."_".$BIN_SIZE."_sample_bins";
            $folder = $PEAKS_TOKEN."_".$out_name."_matchpeaks";
            $out_common = $RESULTS."/".$folder."/common_".$PEAKS_TOKEN."_".$out_name.".bed";
            $out_peaks_file = $RESULTS."/results/".$NAMES[$i]."_".$BIN_SIZE."_bins_sample_peaks.bed";
            $BINS_PEAKS_SAMPLE[$i] = $out_peaks_file;
        }else{
            $out_name = $NAMES[$i]."_".$BIN_SIZE."_sample_bins";
            $folder = $PEAKS_TOKEN."_".$out_name."_matchpeaks";
            $out_common = $folder."/common_".$PEAKS_TOKEN."_".$out_name.".bed";
            $out_peaks_file = $RESULTS.$NAMES[$i]."_".$BIN_SIZE."_bins_sample_peaks.bed";
            $BINS_PEAKS_SAMPLE[$i] = $out_peaks_file;
        }
    	#
    	CleanFile($out_peaks_file);
    	if(!(-e $out_peaks_file) or $OVERWRITE)
        {
        	#
        	$command = "matchpeaks -v $PEAKS_SAMPLES[$i] $sample_bins $PEAKS_TOKEN $out_name";
        	print_mess("$command\n");
        	system($command);
        	$command = "grep sample_bins $out_common > $out_peaks_file";
        	print_mess("$command\n");
        	system($command);
        	#
    	}
        CleanFolder($folder);
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
            if($RESULTS ne "results/"){
                $input_file = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_avg_normalized_spike.txt";
            }else{
                $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_avg_normalized_spike.txt";
            }
        }
        else
        {
            if($RESULTS ne "results/"){
                $input_file = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_avg.txt";
            }else{
                $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_avg.txt";
            }
        }
        if($RESULTS ne "results/"){
                $output_file = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
        }else{
                $output_file = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
        }
        if(!(-e $output_file) or $OVERWRITE)
        {
            $command = "grep -v track ".$BINS_PEAKS_SPIKE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $output_file already exists");
        }
        #
        CleanFile($output_file);
        #
        if($RESULTS ne "results/"){
                $output_file = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
        }else{
                $output_file = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
        }
        if(!(-e $output_file) or $OVERWRITE)
        {
            $command = "grep -v track ".$BINS_PEAKS_SPIKE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join -v 2 - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $output_file already exists");
        }
        #
        CleanFile($output_file);
        #
        # extract values for sample with avg
        if ($token eq $SPIKCHIP_TOKEN)
        {
            if($RESULTS ne "results/"){
                $input_file = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_avg_normalized_sample.txt";
            }else{
                $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_avg_normalized_sample.txt";
            }
        }
        else
        {
            if($RESULTS ne "results/"){
                $input_file = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_avg.txt";
                $output_file = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
            }else{
                $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_avg.txt";
                $output_file = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
            }
        }
        if(!(-e $output_file) or $OVERWRITE)
        {
            $command = "grep -v track ".$BINS_PEAKS_SAMPLE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $output_file already exists");
        }
        #
        CleanFile($output_file);
        #
        if($RESULTS ne "results/"){
                $output_file = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
        }else{
            $output_file = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
        }
        if(!(-e $output_file) or $OVERWRITE)
        {
            $command = "grep -v track ".$BINS_PEAKS_SAMPLE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join -v 2 - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $output_file already exists");
        }
        #
        CleanFile($output_file);
        #
        # extract values for spike with max
        if ($token eq $SPIKCHIP_TOKEN)
        {
            if($RESULTS ne "results/"){
                $input_file = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_max_normalized_spike.txt";
            }else{
                $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_max_normalized_spike.txt";
            }
        }
        else
        {
            if($RESULTS ne "results/"){
                $input_file = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_max.txt";
            }else{
                $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_max.txt";
            }
        }
        if($RESULTS ne "results/"){
            $output_file = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
        }else{
            $output_file = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
        }
        if(!(-e $output_file) or $OVERWRITE)
        {
            $command = "grep -v track ".$BINS_PEAKS_SPIKE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $output_file already exists");
        }
        #
        CleanFile($output_file);
        #
        if($RESULTS ne "results/"){
            $output_file = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
        }else{
            $output_file = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
        }
        if(!(-e $output_file) or $OVERWRITE)
        {
            $command = "grep -v track ".$BINS_PEAKS_SPIKE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join -v 2 - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $output_file already exists");
        }
        #
        CleanFile($output_file);
        #
        # extract values for sample with max
        if ($token eq $SPIKCHIP_TOKEN)
        {
            if($RESULTS ne "results/"){
                $input_file = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_max_normalized_sample.txt";
            }else{
                $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_max_normalized_sample.txt";
            }
        }
        else
        {
            if($RESULTS ne "results/"){
                $input_file = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_max.txt";
            }else{
                $input_file = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_max.txt";
            }
        }
        if($RESULTS ne "results/"){
            $output_file = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
        }else{
            $output_file = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
        }
        if(!(-e $output_file) or $OVERWRITE)
        {
            $command = "grep -v track ".$BINS_PEAKS_SAMPLE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $output_file already exists");
        }
        #
        CleanFile($output_file);
        #
        if($RESULTS ne "results/"){
            $output_file = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
        }else{
            $output_file = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
        }
        if(!(-e $output_file) or $OVERWRITE)
        {
            $command = "grep -v track ".$BINS_PEAKS_SAMPLE[$i]." | gawk '{print \$1\"*\"\$2\"*\"\$3}' | sort | join -v 2 - $input_file | gawk '{print \$1,$n_field;}' > $output_file";
            print_mess("$command\n");
            system($command);
        }else{
            print_mess("\t $output_file already exists");
        }
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
    if($RESULTS ne "results/"){
        $file_avg0 = $RESULTS."/results/".$NAMES[0]."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
        $file_max0 = $RESULTS."/results/".$NAMES[0]."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
        $file_avg1 = $RESULTS."/results/".$NAMES[1]."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
        $file_max1 = $RESULTS."/results/".$NAMES[1]."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
    }else{
        $file_avg0 = $RESULTS.$NAMES[0]."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
        $file_max0 = $RESULTS.$NAMES[0]."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
        $file_avg1 = $RESULTS.$NAMES[1]."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
        $file_max1 = $RESULTS.$NAMES[1]."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
    }
    $commandAVG = "join $file_avg0 $file_avg1 ";
    $commandMAX = "join $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
        if($RESULTS ne "results/"){
            $file_avg = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
            $file_max = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
        }else{
            $file_avg = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
            $file_max = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
        }
        $commandAVG = $commandAVG." | join - $file_avg ";
        $commandMAX = $commandMAX." | join - $file_max ";
    }

    if($RESULTS ne "results/"){
        $final_avg = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
        $final_max = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
    }else{
        $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_avg_peaks.txt";
        $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_max_peaks.txt";
    }
    #
    if(!(-e $final_avg) or $OVERWRITE)
    {
        $commandAVG = $commandAVG." | $field_notnull > $final_avg";
        print_mess("$commandAVG\n");
        system($commandAVG);
    }else{
        print_mess("\t $final_avg already exists");
    }
    if(!(-e $final_max) or $OVERWRITE)
    {
        $commandMAX = $commandMAX." | $field_notnull > $final_max";
        print_mess("$commandMAX\n");
        system($commandMAX);
    }else{
        print_mess("\t $final_max already exists");
    }
    #
    SaveFile($final_avg);
    SaveFile($final_max);

    # Second, sample values for peaks
    # at least, two experiments
    if($RESULTS ne "results/"){
        $file_avg0 = $RESULTS."/results/".$NAMES[0]."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
        $file_max0 = $RESULTS."/results/".$NAMES[0]."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
        $file_avg1 = $RESULTS."/results/".$NAMES[1]."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
        $file_max1 = $RESULTS."/results/".$NAMES[1]."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
    }else{
        $file_avg0 = $RESULTS.$NAMES[0]."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
        $file_max0 = $RESULTS.$NAMES[0]."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
        $file_avg1 = $RESULTS.$NAMES[1]."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
        $file_max1 = $RESULTS.$NAMES[1]."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
    }
    $commandAVG = "join $file_avg0 $file_avg1 ";
    $commandMAX = "join $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
        if($RESULTS ne "results/"){
            $file_avg = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
            $file_max = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
        }else{
            $file_avg = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
            $file_max = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
        }
        $commandAVG = $commandAVG." | join - $file_avg ";
        $commandMAX = $commandMAX." | join - $file_max ";
    }

    if($RESULTS ne "results/"){
        $final_avg = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
        $final_max = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
    }else{
        $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_avg_peaks.txt";
        $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_max_peaks.txt";
    }
    #
    if(!(-e $final_avg) or $OVERWRITE)
    {
        $commandAVG = $commandAVG." | $field_notnull > $final_avg";
        print_mess("$commandAVG\n");
        system($commandAVG);
    }else{
        print_mess("\t $final_avg already exists");
    }
    if(!(-e $final_max) or $OVERWRITE)
    {
        $commandMAX = $commandMAX." | $field_notnull > $final_max";
        print_mess("$commandMAX\n");
        system($commandMAX);
    }else{
        print_mess("\t $final_max already exists");
    }
    #
    SaveFile($final_avg);
    SaveFile($final_max);

    # Third, spike values for bg
    # at least, two experiments
    if($RESULTS ne "results/"){
        $file_avg0 = $RESULTS."/results/".$NAMES[0]."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
        $file_max0 = $RESULTS."/results/".$NAMES[0]."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
        $file_avg1 = $RESULTS."/results/".$NAMES[1]."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
        $file_max1 = $RESULTS."/results/".$NAMES[1]."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
    }else{
        $file_avg0 = $RESULTS.$NAMES[0]."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
        $file_max0 = $RESULTS.$NAMES[0]."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
        $file_avg1 = $RESULTS.$NAMES[1]."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
        $file_max1 = $RESULTS.$NAMES[1]."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
    }
    $commandAVG = "join $file_avg0 $file_avg1 ";
    $commandMAX = "join $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
        if($RESULTS ne "results/"){
            $file_avg = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
            $file_max = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
        }else{
            $file_avg = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
            $file_max = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
        }
        $commandAVG = $commandAVG." | join - $file_avg ";
        $commandMAX = $commandMAX." | join - $file_max ";
    }

    if($RESULTS ne "results/"){
        $final_avg = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
        $final_max = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
    }else{
        $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_avg_bg.txt";
        $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_spike_max_bg.txt";
    }
    #
    if(!(-e $final_avg) or $OVERWRITE)
    {
        $commandAVG = $commandAVG." | $field_notnull > $final_avg";
        print_mess("$commandAVG\n");
        system($commandAVG);
    }else{
        print_mess("\t $final_avg already exists");
    }
    if(!(-e $final_max) or $OVERWRITE)
    {
        $commandMAX = $commandMAX." | $field_notnull > $final_max";
        print_mess("$commandMAX\n");
        system($commandMAX);
    }else{
        print_mess("\t $final_max already exists");
    }
    #
    SaveFile($final_avg);
    SaveFile($final_max);
    
    # Forth, sample values for bg
    # at least, two experiments
    if($RESULTS ne "results/"){
        $file_avg0 = $RESULTS."/results/".$NAMES[0]."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
        $file_max0 = $RESULTS."/results/".$NAMES[0]."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
        $file_avg1 = $RESULTS."/results/".$NAMES[1]."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
        $file_max1 = $RESULTS."/results/".$NAMES[1]."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
    }else{
        $file_avg0 = $RESULTS.$NAMES[0]."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
        $file_max0 = $RESULTS.$NAMES[0]."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
        $file_avg1 = $RESULTS.$NAMES[1]."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
        $file_max1 = $RESULTS.$NAMES[1]."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
    }
    $commandAVG = "join $file_avg0 $file_avg1 ";
    $commandMAX = "join $file_max0 $file_max1 ";

    # the rest of experiments (if any)
    for($i=2; $i<$n_experiments; $i++)
    {
        if($RESULTS ne "results/"){
            $file_avg = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
            $file_max = $RESULTS."/results/".$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
        }else{
            $file_avg = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
            $file_max = $RESULTS.$NAMES[$i]."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
        }
        $commandAVG = $commandAVG." | join - $file_avg ";
        $commandMAX = $commandMAX." | join - $file_max ";
    }

    if($RESULTS ne "results/"){
        $final_avg = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
        $final_max = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
    }else{
        $final_avg = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_avg_bg.txt";
        $final_max = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$token."_".$BIN_SIZE."_sample_max_bg.txt";
    }
    #
    if(!(-e $final_avg) or $OVERWRITE)
    {
        $commandAVG = $commandAVG." | $field_notnull > $final_avg";
        print_mess("$commandAVG\n");
        system($commandAVG);
    }else{
        print_mess("\t $final_avg already exists");
    }
    if(!(-e $final_max) or $OVERWRITE)
    {
        $commandMAX = $commandMAX." | $field_notnull > $final_max";
        print_mess("$commandMAX");
        system($commandMAX);
    }else{
        print_mess("\t $final_max already exists");
    }
    #
    SaveFile($final_avg);
    SaveFile($final_max);
}



##################
# GenerateBoxplot
##################


sub NumberLines
{
    my $path = $_[0];
    my $count = `wc -l < $path`;
    die "wc failed: $?" if $?;
    chomp($count);
    if($count == 0){
        print_mess("\t\t The file is empty:", $path);
    }
    return $count;
}

# addPath($norm_token, $exp, $val, $extension, @files_array)
sub addPath
{
    my $norm_token = $_[0];
    my $exp = $_[1];
    my $val = $_[2];
    my $extension = $_[3];
    my @files_array = @{$_[4]};
    my $input;
    my $nb_line;
    
    if($RESULTS ne "results/"){
        $input = $RESULTS."/results/".$FINAL_TOKEN."_".join("-",@NAMES)."_".$norm_token."_".$BIN_SIZE."_".$exp."_".$val.$extension;
    }else{
        $input = $RESULTS.$FINAL_TOKEN."_".join("-",@NAMES)."_".$norm_token."_".$BIN_SIZE."_".$exp."_".$val.$extension;
    }
    $nb_line = NumberLines($input);
    if($nb_line != 0){push(@files_array, $input);}
    return @files_array;
}

# peaksBgTXT($exp, $val, $extension, \@files_array)
sub peaksBgTXT
{
    my $exp = $_[0];
    my $val = $_[1];
    my $extension = $_[2];
    my @files_array = @{$_[3]};
    
    if($RAW){@files_array = addPath($RAW_TOKEN, $exp, $val, $extension, \@files_array);}
    if($TRADITIONAL){@files_array = addPath($TRADITIONAL_TOKEN, $exp, $val, $extension, \@files_array);}
    if($CHIPRX){@files_array = addPath($CHIPRX_TOKEN, $exp, $val, $extension, \@files_array);}
    if($TAGREMOVAL){@files_array = addPath($TAGREMOVAL_TOKEN, $exp, $val, $extension, \@files_array);}
    if($SPIKCHIP){@files_array = addPath($SPIKCHIP_TOKEN, $exp, $val, $extension, \@files_array);}
    
    return @files_array;
}


sub pathVector
{
    my @files_array = @{$_[0]};
    my $last_element;
    my $pathvec;
    my $file;
    
    if((scalar @files_array) > 1){
        $last_element = pop @files_array;
        $pathvec = "paths \<\- c(";
        for my $el (@files_array) {
                $pathvec = $pathvec."\"".$el."\"\,";
        }
        $pathvec = $pathvec."\"".$last_element."\"\)\n\n\n\n";
    }else{
        $pathvec = "paths \<\- \"".$files_array[0]."\"\n\n";
    }
    return $pathvec;
}

sub retrieveNorm
{
    my @norm_array;
    
    if($RAW){push(@norm_array, $RAW_TOKEN);}
    if($TRADITIONAL){push(@norm_array, $TRADITIONAL_TOKEN);}
    if($CHIPRX){push(@norm_array, $CHIPRX_TOKEN);}
    if($TAGREMOVAL){push(@norm_array, $TAGREMOVAL_TOKEN);}
    if($SPIKCHIP){push(@norm_array, $SPIKCHIP_TOKEN);}
    
    return @norm_array;
}

sub GenerateBoxplot
{
    my $experiment = $_[0];
    my $value = $_[1];
    my ($Rfile,$Routput_file);
    my $PDF_file;
    my @input_files;
    my @Norm_array;
    my $input_file;
    my $expression;
    my ($i,$j);
    my $midpoint;
    my @CLASSES;
    my $length;
    my $last_element;
    my $nb_line;
    my $towrite;

    
    # (A) Boxplot with labelling
    # R code to generate the corresponding final boxplot
    if($RESULTS ne "results/"){
        $Rfile = $RESULTS."/".$RSCRIPTS.join("-",@NAMES)."_".$BIN_SIZE."_".$experiment."_".$value."_boxplot.R";
        $PDF_file = $RESULTS."/".$PLOTS."/".join("-",@NAMES)."_".$BIN_SIZE."_".$experiment."_".$value.".pdf";
    }else{
        $Rfile = $RSCRIPTS.join("-",@NAMES)."_".$BIN_SIZE."_".$experiment."_".$value."_boxplot.R";
        $PDF_file = $PLOTS.join("-",@NAMES)."_".$BIN_SIZE."_".$experiment."_".$value.".pdf";
    }
    print_mess("Writing ", $Rfile);
    if(!(-e $Rfile) or $OVERWRITE)
    {
        (open(RFILE,'>',$Rfile)) or print_error("R SCRIPT (boxplots): FILE $Rfile file can not be opened to write");
        @input_files = peaksBgTXT($experiment, $value, "_peaks.txt", \@input_files);
        @input_files = peaksBgTXT($experiment, $value, "_bg.txt", \@input_files);
        @Norm_array = retrieveNorm();
        
        # Creating the vector of files path and the palette
        $towrite = pathVector(\@input_files);
        $towrite = $towrite."p \<\- colorRampPalette($PALETTE)\n\n";
        
        # Preparing values and labels
        $towrite = $towrite."filist <- lapply(paths,function(x) read.table(x, row.names=1))\n";
        $towrite = $towrite."filist <- unlist(lapply(filist, function(x) as.list(data.frame(x))), recursive=FALSE)\n";
        $towrite = $towrite."filist <- lapply(filist, log2)\n\n"; 
        $towrite = $towrite."labels <- \"".join("-",@NAMES)."\"\n";
        $towrite = $towrite."labels <- unlist(strsplit(labels,\"-\"))\n\n";
        
        #Determine position of ticks
        $towrite = $towrite."normNames <- \"".join("-", @Norm_array)."\"\n";
        $towrite = $towrite."normNames <- unlist(strsplit(normNames,\"-\"))\n\n";
        $towrite = $towrite."normNb <- length(unique(normNames))\n";
        $towrite = $towrite."if($n_experiments%%2 == 0){\n";
        $towrite = $towrite."\tpositions <- seq(($n_experiments/2)+1, normNb*$n_experiments, by =$n_experiments)\n";
        $towrite = $towrite."}else\n";
        $towrite = $towrite."\tpositions <- seq(floor($n_experiments/2)+1, normNb*$n_experiments, by = $n_experiments)\n\n";
        
        #Determine label of each boxplot
        $towrite = $towrite."peaksbgvec <- strsplit(paths, \"_\")\n";
        $towrite = $towrite."idx <- unique(sapply(peaksbgvec, function(x) grep(\"txt\",x)))\n";
        $towrite = $towrite."peaksbgvec <- gsub(\".txt\", \"\", sapply(peaksbgvec, \"[\", idx))\n\n";
        
        # Generating code for boxplot
        $towrite = $towrite."pdf(\"$PDF_file\")\n";
        $towrite = $towrite."boxplot(filist, xaxt=\"n\",xlab =\"\",lwd=1, ylab=\"log2 (avg)\",outline=FALSE, col=rev(p(4)), notch=TRUE, main=paste0(\"$experiment,$value \", paste(labels, collapse=\",\")))\n";
        $towrite = $towrite."legend(\"topright\",labels,fill=rev(p(4)))\n";
        $towrite = $towrite."axis(1, at = positions, labels = paste(normNames, peaksbgvec, sep=\"-\"), las=2, cex.axis=0.5)\n";
        $towrite = $towrite."dev.off()\n";
        print RFILE $towrite;
        close(RFILE);
    }else{
        print_mess("\t The R script ",$Rfile, "already exists.") 
    }
    #
    if(!(-e $PDF_file) or $OVERWRITE)
    {
        if($RESULTS ne "results/"){
            $Routput_file = $RESULTS."/".$RSCRIPTS.join("-",@NAMES)."_".$BIN_SIZE."_".$experiment."_".$value."_boxplot.Rout";
        }else{
            $Routput_file = join("-",@NAMES)."_".$BIN_SIZE."_".$experiment."_".$value."_boxplot.Rout";
        }
        # execute R script
        $command = "R CMD BATCH $Rfile $Routput_file";
        print_mess("$command\n");
        system($command);
        #
        # error check in Rout file
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
    }else{
        print_mess("\t The plot ",$PDF_file, "already exists.") 
    }
}

sub CleanFile
{
    my $file = $_[0];
    my $command;
    
    # info for cleaning intermediate files (option -c)
    $command = "rm -f $file";
    push(@CLEAN_PROCEDURE,$command);
}

sub CleanFolder
{
    my $folder = $_[0];
    my $command;
    # info for cleaning intermediate folders (option -c)
    $command = "rm -rf $folder";
    push(@CLEAN_PROCEDURE,$command);
}

sub CleaningFiles
{
    my $i;
    my $command;
    
    
    if ($CLEAN)
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
