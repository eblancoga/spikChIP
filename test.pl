#!/usr/bin/perl -w                                #

=pod
use Getopt::Long;

Getopt::Long::Configure("gnu_getopt", "auto_abbrev", "ignore_case_always");

my $update='';
my ($configuration_filename,$chrominfo_file);
my $n_files;

   Getopt::Long::GetOptions(
      'normalize|n' => \my $normalize,
      'exclude|e=s' => \my @exclude,
      'help|h'      => \my $help,
      'include|i=s' => \my @include,
      'recurse|r'   => \my $recurse,
      'update|u'  => \$update,
      '2update|2'   => \my $update2,
      'copy|c'      => \my $copy,
      'move|m'      => \my $move,
   );

if($update){
    print "$update\n";
}

$n_files = $#ARGV+1;
($n_files == 2) or print "ERROR\n";
($configuration_filename,$chrominfo_file) = @ARGV;

print "configuration_filename:$configuration_filename\n";
=cut

my @test = (1,2,3);
my $a = 1;
mytest(\@test, $a);

sub mytest
{
    my @all = @{$_[0]};
    my $b = $_[1];
    
    print "The array is: @all \n";
    print "The argument is: $b \n";
    
}