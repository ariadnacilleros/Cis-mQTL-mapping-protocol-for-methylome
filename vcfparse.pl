#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

# vcfparse.pl
# program to extract the first 8 columns from a vcf for the IC.pl program 
# and optionally gzip the output
# Created W.Rayner 2016-2017 wrayner@well.ox.ac.uk
#
# version 1.0 created
# version 1.1 added option for gzip on output and command line options
#
#
#
#

my $folder = '.';
my $gzip = 0;
my $stem;

GetOptions
 (
 "d|directory=s" => \$folder,   # top level directory for the imputed data files
 "g|gzip"        => \$gzip,
 "o|output=s"    => \$stem      # Set path to an output directory, will be created if it doesn't exist
 );

print "Options:\n";
print "Input Folder:    $folder\n";
print "Output Folder:   $stem\n";

if (-e $stem)
 {
 print "Output folder $stem already exists\n";
 }
else
 {
 mkdir $stem;
 print "Creating Output folder $stem\n";
 }

#find vcfs
my @filelisting = find_all_files($folder);
my $total = @filelisting;
print "Found $total files in $folder\n";

foreach my $file (@filelisting)
 {
 my $toprocess;
 if ($file =~ /(.*\.vcf)\.gz$/)
  {
  $toprocess = $1;
  }
 elsif ($file =~ /(.*\.vcf)$/) 
  {
  $toprocess = $file;
  }
 else
  {
  }
  
 if ($toprocess)
  {
  my $cmd;
  my($filename, $dirs, $suffix) = fileparse($toprocess);
  my $outfile = $stem.'/'.$filename.'.cut';
  if (-e $outfile)
   {
   print "ERROR: $outfile Exists, skipping\n";
   }
  else
   {
   print "Processing $file\n";
   if (!$gzip)
    {
    $cmd = "zgrep -v '##'  $file \| cut -f 1-8 \> $outfile"; 
    }
   else
    {
    $outfile = $outfile.'.gz';
    $cmd = "zgrep -v '##'  $file \| cut -f 1-8 \| gzip \> $outfile"; #gzip output
    }
   system $cmd;
   print "Output in $outfile\n";
   }
  }
 }

sub find_all_files
 {
 my $path = $_[0];
 my @filelist;
 
 # Open the directory.
 if (opendir (DIR, $path)) #or die "Unable to open $path: $!";
  {
  # Read files
  my @files = grep { !/^\.{1,2}$/ } readdir (DIR);
  # Close the directory.
  closedir (DIR);
 
  @files = map { $path . '/' . $_ } @files;

  foreach my $file (@files)
   {
   if (-d $file) # If it is a directory
    {
    push @filelist, find_all_files($file);
    }
   else # If it isn't a directory 
    {
    push @filelist, $file;
    }
   } 
  }
 return @filelist; 
 }
 






