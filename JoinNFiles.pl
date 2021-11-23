#!/usr/bin/perl -w                                
#                                                 #
# Script JoinNFiles.pl                            #
#                                                 #
#                                                 #
# Input : (key1,value1),...,(keyN,valueN)         # 
# Output: (key,value1,...,valueN)                 #
#                                                 #
#                                                 #
# by Enrique Blanco (2021)                        #
###################################################

use strict;
use Getopt::Std;
use Term::ANSIColor;

#DEFINEs
my $TRUE = 1;
my $FALSE = 0;
my $PSEUDOCOUNT = 0.001;

## Step 1. Reading arguments
my %opt;

(getopts('v',\%opt)) or print_error("parser: Problems reading options\n");

print_mess("JoinNFiles by Enrique Blanco (CRG 2021)\n\n");
print_mess("Stage 1.  Reading options");

my (@files);
my $n_files;

$n_files = $#ARGV+1;
($n_files >= 2) or print_error("USAGE: Two filenames are at least required but $n_files are provided!");
@files = @ARGV;

print_ok();
##

## Step 2. Loop to load the information from each file in @files
my $i;
my $ifile;
my $n_lines;
my @record;
my $line;
my ($key,$value);
my @RESULTS;
my %KEYS;


print_mess("Stage 2.  Looping along the list of $n_files\n");
for($i=0; $i<scalar(@files); $i++)
{
    $ifile = $files[$i];
    print_mess("Stage 2 - Processing file ($i, $ifile)\n");

    (open(FILE,$ifile)) or print_error("FILE $ifile can not be opened");
    $n_lines = 0;
    while($line=<FILE>)
    {
	chomp($line);

	@record = split(/\s+/,$line);
	$key = $record[0];
	$value = $record[1];

	# register the value for this key and experiment
	$RESULTS[$i]{$key} = $value;

	# save this key in the global catalog of keys
	if (!defined($KEYS{$key}))
	{
	    $KEYS{$key}=1;
	}
	
	$n_lines++;
    }
    close(FILE);

    print_mess("$n_lines elements from $ifile have been acquired");
    print_ok();
}


## Step 3. Reporting the final table of results
my $string;
my $total_keys;


print_mess("Stage 3. Generating the final table of results\n");

$string="";
foreach $key (sort keys(%KEYS))
{
    $string = $key."\t";

    for ($i=0; $i<scalar(@files); $i++)
    {
	if (defined($RESULTS[$i]{$key}))
	{
	    $value = $RESULTS[$i]{$key};
	}
	else
	{
	    $value = $PSEUDOCOUNT;
	}
	$string = $string."\t".$value;
    }

    print $string,"\n";
}

$total_keys = scalar(keys(%KEYS)); 
print_mess("$total_keys elements identified in total among all samples");
print_ok();

## Step X. Finishing successful program execution

print_mess("Successful termination:");
print_ok();
exit(0);
##


############ Subroutines

sub print_mess
{
        my @mess = @_;

        print STDERR color("bold green"),"%%%% @mess" if (exists($opt{v}));
	print STDERR color("reset");
}

sub print_error
{
        my @mess = @_;

        print STDERR color("bold green"),"%%%% @mess\n";
	print STDERR color("reset");
	exit();
}

sub print_ok
{
    if (exists($opt{v}))
    {
	print STDERR color("bold green"), "\t\t[OK]\n\n";
	print STDERR color("reset");
    }
}

