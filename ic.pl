#!/usr/bin/perl

# If any file header changes then search for calls to get_header to update the names, change in future to read from a file
# requires file with the mapping of directory to study/cohort name
#
#
# W.Rayner 2016 
# wrayner@well.ox.ac.uk
#
# Requires GD to be installed on the system for the plotting
# - needs libgd, for ubuntu: sudo apt-get -y install libgd2-xpm-dev build-essential
# - needs GD                 cpan GD::Graph
# -
# -


use strict;
use warnings;
use GD::Graph::points;
use Getopt::Long;
use GD::Graph::hbars;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Time::localtime;
use File::stat;
use File::Basename;
use File::Path qw(make_path);


$| = 1;
my %altaf;
my %directorylist;
my @gwsummary;
my $reffile;
my $hrcflag = 0;
my $kgflag = 0;
my $population = 'ALL';
my $directory; # Top level directory containing files or sub-folders with files
my $infile;
my %directorymapping;
my $outputpath = 'outputDir';
my @allGWinfofiles;
my $columns = get_width(); 
my $mid = int($columns/2+0.5);

print "\n\n";
printf("%*s", $mid+12, "IC: a program to check imputed data files\n");
printf("%*s", $mid+5, "William Rayner 2016\n");
printf("%*s", $mid+6, "wrayner\@well.ox.ac.uk\n");
print "\n";
printf("%*s", $mid+5, "Version 1.0.2\n\n\n");

GetOptions
 (
 "f|file=s"      => \$infile,      # file containing folder to name mappings
 "d|directory=s" => \$directory,   # top level directory for the imputed data files
 "h|hrc"         => \$hrcflag,     # flag to set HRC reference panel is being used
 "r|ref=s"       => \$reffile,     # input reference file name (1000G or HRC) 
 "g|1000g"       => \$kgflag,      # flag to set 1000G reference panel is being used
 "p|pop=s"       => \$population,  # 1000G population frequency to use in the checks
 "o|output=s"    => \$outputpath   # Set path to an output directory, will be created if it doesn't exist
 );

if (!$reffile)
 {
 print "No reference file specified\n";
 usage();
 exit;
 }

if (!$directory)
 {
 print "No directory containing imputed summary files specified\n";
 usage();
 exit;
 }

if (!$hrcflag and !$kgflag)
 {
 $hrcflag = 1;
 print "No reference panel type specified, assuming HRC\n";
 }

# get run date
my $dt = get_time();

# open Logfile
my $logfile = 'ic-'.$dt.'.log';

make_directory($outputpath);
$logfile = $outputpath.'/'.$logfile;

open L, ">$logfile" or die $!;

print "Options Specified:\n";
print "Directory path:    $directory\n";
if ($hrcflag)
 {
 print "Reference:          HRC\n";
 print "Reference path:     $reffile\n";
 }
elsif ($kgflag) 
 {
 print "Reference:          1000G\n";
 print "Reference path:     $reffile\n";
 print "Population:         $population\n";
 }
 
print "Output Saved in:    $outputpath\n";
 
if ($infile) 
 {
 print "Mapping file:       $infile\n";
 } 
print "\nLog file:           $logfile\n\n"; 

# find all files
my @filestoprocess = find_all_files($directory);

if ($#filestoprocess > 0) #if there are no files to process, don't continue
 {
 
 # list mapping directory names to study names
 if ($infile)
  {
  print "Directory mapping file: $infile\n";
  %directorymapping = get_listing($infile);
  }
 else #if no mapping file for directory is given take top most directory as name
  {
  print "No directory mapping file, creating a directory mapping\n";
  %directorymapping = get_mapping(\@filestoprocess);
  my @k = keys %directorymapping;
  print "Using: (name -> path)\n";
  foreach my $s (@k)
   {
   print "$directorymapping{$s} -> $s\n";
   }
  print "\n";
  }
 
 if ($hrcflag)
  {
  my($filename, $dirs, $suffix) = fileparse($reffile);
  print "Reading HRC file $filename from $dirs\n";
  %altaf = read_hrc($reffile);
  }
 elsif ($kgflag)
  {
  my($filename, $dirs, $suffix) = fileparse($reffile);
  print "Reading 1000G file $filename from $dirs\nUsing $population population\n";
  %altaf = read_kg($reffile, $population);
  }

 # extract the data from each per chromosome file and invoke the plotting
 %directorylist = process_data(\@filestoprocess, $outputpath);

 my @directorylisting = keys %directorylist;
 foreach my $subdir (@directorylisting)
  {
  print "Final processing: $subdir\n";
  my $infofile = create_genome_wide_summary($subdir); #create summaries and return info summary file name for across cohort summary 
  push (@allGWinfofiles, $infofile); 
  study_date_upload($subdir, \@filestoprocess, \%directorymapping); #adds a file with the date of upload of the cohort based on date modified time
  }

 #plot genome wide info for across cohort summary
 my $outfile = $outputpath.'/Summary-info-test-temp.txt';
 my $allinfofile = $outputpath.'/Summary-info-test.txt';
  
 #my @files = find_info_files($outputpath, $outfile);
 
 open AO, ">$allinfofile" or die $!;
 foreach my $file (@allGWinfofiles)
  {
  my $line = get_gw_counts($file);
  print AO "$line";
  }
 close AO;
  
 plot_summary($allinfofile, $outputpath); 
 
 my $summaryfolder = $outputpath.'/summaryOutput';
 my $pathtoscript = $0;
 my($scriptfile, $scriptdirs, $scriptsuffix) = fileparse($pathtoscript); 

 my $cmd = "java -Xmx4g -classpath $scriptdirs/Junk.jar uk.ac.ox.well.t2d.reports.imputation.Main $outputpath $summaryfolder";
 print "$cmd\n";
 system ($cmd);
 }
else
 {
 print "No files found to process, exiting\n";
 }
 
sub get_mapping
 { #derive from a directory name
 my @files = @{$_[0]};
 my %mapping;
 my %map;
 my $k = 0;
 my @all;
 my $max = 0;
 my @maps;
 
 foreach my $path (@files)
  {
  my($filename, $dirs, $suffix) = fileparse($path); 
  $dirs =~ s/^\///;
  #print "$dirs\n";
  $map{$dirs} = 1;
  }

 @maps = keys(%map);
 
 foreach my $m (@maps)
  {
  my @directories = split(/\//, $m);
  for (my $i = 0; $i <= $#directories; $i++)
   {
   $all[$i][$k] = $directories[$i];
   if ($i > $max)
    {
    $max = $i;
    }
   }
  $k++; 
  }
 
 for (my $j = 0; $j <= $max; $j++)
  {
  my %temp;
  for (my $i = 0; $i < $k; $i++)
   {
   if ($all[$j][$i])
    {
    $temp{$all[$j][$i]} = 1;
    }
   }
  
  my @temparray = keys %temp;
  my $count = @temparray;
  if  ($count == $k)
   {
   foreach my $m (@maps)
    {
    my @directories = split(/\//, $m);
    $mapping{$m} = $directories[$j];
    }
   }
  }
 return %mapping; 
 } 
 
sub study_date_upload
 {
 my $directory = $_[0];
 my @files = @{$_[1]};
 my %mapping = %{$_[2]};
 
 my @directories = split(/\//, $directory);
 my $finaldir = $directories[$#directories];
 $finaldir =~ /(.*?)\.(.*?)/;
 #my $study = $1; 
 
 #$directory =~ /(.*?)\.(.*?)/;
 my $studyfromdirname = $1;
 
 foreach my $file (@files)
  {
  my $studyfromfile = get_study($file, \%mapping);
  if ($studyfromfile eq $studyfromdirname)
   {
   my $month = (localtime(stat($file)->mtime)->mon) + 1;
   my $day = localtime(stat($file)->mtime)->mday;
   my $year = (localtime(stat($file)->mtime)->year) + 1900;
   my $timestamp = "$day.$month.$year";
   my $studyfile = $directory.'/Meta-'.$studyfromdirname.'.txt';
   open O, ">>$studyfile" or die $!;
   print O "UPLOADED\t$timestamp\n";
   close O;
   }
  }
 close IN;
 }
 
sub create_genome_wide_summary
 {
 my $directory = $_[0];
  
 opendir(DIR, "$directory") or die "Not a directory $directory";
 my @all_files =  grep !/^\.\.?\z/, readdir DIR;
  
 my %info;
 my %af;
 my %counts;

 foreach my $file (@all_files)
  {
  my $filepath = $directory.'/'.$file;
  if ($file =~ /^(.*?)-(.*?)\.chr(.*?)\.txt$/)
   {
   #print "$file\n";
   my $type = $1;
   my $study = $2;
   my $chr = $3;
   my $flag = 0;
   
   open IN, "$filepath" or die $!;
   while (<IN>)
    {
    chomp;
    my @temp = split/\t/;
     
    if ($type eq 'AF')
     {
     $af{$temp[0]} += $temp[1];
     }
    elsif ($type eq 'INFO')
     {
     $info{$temp[0]} += $temp[1];    
     }
    elsif ($type eq 'COUNTS')
     {
     if ($flag)
      {
      $counts{'NOMATCH'} += $temp[0];
      }
     else
      {
      $counts{'MATCH'} += $temp[0];
      $flag = 1;
      }
     }
    else 
     {
     # unknown file type
     print "WARNING: $file unknown type\n";
     }
    }
   close IN; 
   }
  }
 #print All summaries
 
 my @directories = split(/\//, $directory);
 my $finaldir = $directories[$#directories];
 $finaldir =~ /(.*?)\.(.*?)/;
 my $study = $1;
 my $affile = $directory.'/AF-'.$study.'.GW.txt';
 my $infofile = $directory.'/INFO-'.$study.'.GW.txt';
 my $countsfile = $directory.'/COUNTS-'.$study.'.GW.txt';
 
 print_summary($affile, \%af);
 print_summary($infofile, \%info);
 print_summary($countsfile, \%counts);
 return $infofile;
 }
 
sub print_summary
 {
 my $outfile = $_[0];
 my %data = %{$_[1]};
 my $total = 0;
 
 my @items = keys %data;
 my @sorted = sort {$a cmp $b} @items;
  
 if ($data{'Total'})
  {
  $total = $data{'Total'};
  }
 else
  {
  for (my $i = 0; $i <= $#sorted; $i++)
   {
   $total += $data{$sorted[$i]};
   }
  }
 
 open OUT, ">$outfile" or die $!;
 for (my $i = 0; $i <= $#sorted; $i++)
  {
  my $pct = sprintf("%0.3f", $data{$sorted[$i]}/$total * 100);
  print OUT "$sorted[$i]\t$data{$sorted[$i]}\t$pct\n";
  }
 close OUT;
 }

  
 
sub read_hrc
 {
 my %altaf;
 my $file = $_[0];
 open IN, "$file" or die $!;
 while (<IN>)
  {
  chomp;
  if (!/\#.*/)
   {
   if ($. % 100000 == 0)
    {
    print L " $.\n";
    }
   my @temp = split/\s+/;
   my $refaltchrpos = $temp[0].':'.$temp[1].'_'.$temp[3].'_'.$temp[4];
   $altaf{$refaltchrpos} = $temp[7];
   }
  }
 close IN;
 return %altaf;
 }

sub process_data #add new calls to subroutines here for further functionality
 {
 my @filelisting = @{$_[0]};
 my $outputdir = $_[1];
 my %genomewidelist;
 
 for (my $i = 0; $i <= $#filelisting; $i++)
  {
  my $format = detect_format($filelisting[$i]); #format can be OXFORD, MINIMAC, SANGER, MINIMAC2 or UNKNOWN. 
  #print "$filelisting[$i]\t$format\n";
  if ($format ne 'UNKNOWN' and $format ne 'MINIMAC2') #MINIMAC2 used when processing MINIMAC so ignore here
   {
   my $file;
   my $study = get_study($filelisting[$i], \%directorymapping);
  
   my $newdir = $outputdir.'/'.$study.'.'.$dt;
   make_directory($newdir);
   my $datadir = $newdir.'/RAW';
   make_directory($datadir);
   my $formatfile = $newdir.'/Meta-'.$study.'.txt';
   open O, ">>$formatfile" or die $!;
   print O "FORMAT\t$format\n";
   close O;
   $genomewidelist{$newdir} = 1;
   
   if ($format eq 'OXFORD')
    {
    print L "$filelisting[$i]\t$format\t$study\n";
    $file = process_oxford($filelisting[$i], $study, $datadir, $newdir);
    }
   elsif ($format eq 'SANGER')
    {
    print L "$filelisting[$i]\t$format\t$study\n";
    $file = process_sanger($filelisting[$i], $study, $datadir, $newdir);
    }
   elsif ($format eq 'MINIMAC')
    {
    print L "$filelisting[$i]\t$format\t$study\n";
    $file = process_minimac($filelisting[$i], $study, $datadir, $newdir);
    }
   print "$filelisting[$i]\t$format\t$study\n";
   
   #start summarising file
   my $afplotfile = plot_af($file, $study, $newdir);
   my $manhfile = plot_manh_info($file, $study, $newdir); 
   my $afinfofile = plot_af_info($file, $study, $newdir);
      
   my $summaryinfofile = get_counts($file, $study, $newdir, 1, 'INFO');
   my $summaryinfoimage = plot_counts($summaryinfofile, $study, $newdir, 'INFO Score');
   my $summaryaffile = get_counts($file, $study, $newdir, 3, 'AF');
   my $summaryafimage = plot_counts($summaryaffile, $study, $newdir, 'Alt Allele Frequency');
   
   my $positionfile = plot_position_row($file, $study, $newdir);
   
   print L "$file\t$study\t$newdir\n";
   
   print L "Manhattan Info plot:           $manhfile\n";
   print L "AF plot:                       $afplotfile\n";
   print L "MAF vs INFO plot:              $afinfofile\n";
   print L "Info summary:                  $summaryinfofile\n";
   print L "Allele Frequency summary:      $summaryaffile\n";
   print L "Position check file:           $positionfile\n";
   print  "Manhattan Info plot:            $manhfile\n";
   print  "AF plot:                        $afplotfile\n";
   print  "MAF vs INFO plot:               $afinfofile\n";
   print  "Info summary:                   $summaryinfofile\n";
   print  "Allele Frequency summary:       $summaryaffile\n";
   print  "Position check file:            $positionfile\n";
   
   }
  }
 return %genomewidelist;
 }
 
sub get_counts
 {
 my $file = $_[0];
 my $study = $_[1];
 my $outdir = $_[2];
 my $col = $_[3];
 my $type = $_[4];
 my %counts;
 #my $chr = get_chr($file, '2');
 my $total = 0;
 
 $file =~ /.*\/(.*?)\.txt$/;
 
 my $outfile = $outdir.'/'.$type.'-'.$1.'.txt';
 
 open IN, "$file" or die $!;
 while (<IN>)
  {
  chomp;
  my @temp = split/\t/;
  for (my $i = 10; $i >= 0; $i--)
   {
   if (int($temp[$col]*10) == $i)
    {
    $counts{$i}++;
    $total++;
    $i = 0;
    }
   }
  }
 open OUT, ">$outfile" or die $!;
 my @c = keys %counts;
 my @cs = sort {$a <=> $b} @c;
 for (my $i = 0; $i <= $#cs; $i++)
  {
  my $band = sprintf("%0.1f",$cs[$i]/10);
  my $pct = sprintf("%0.4f", $counts{$cs[$i]}/$total*100);
  print OUT "$band\t$counts{$cs[$i]}\t$pct\n";
  }
 print OUT "Total\t$total\n"; 
 close IN;
 close OUT;
 return $outfile;
 }
 
 
sub make_directory
 {
 my $dir = $_[0];
 if (-e $dir)
  {
  #print L "Output directory $dir exists, will not recreate\n";
  }
 else
  {
  make_path($dir);
  #print L "Creating output directory $dir\n";
  }
 }
 
sub get_listing
 {
 my $file = $_[0];
 my %directories;
 open IN, "$file" or die $!; 
 while (<IN>)
  {
  chomp;
  my @temp = split/\t/;
  $temp[1] =~ s/\s+//g; #remove all spaces from name as these will cause problems
  $directories{$temp[0]} = $temp[1];
  }
 return %directories; 
 }
 
sub get_study
 {
 my $filename = $_[0];
 my %directories= %{$_[1]};
 my $study;
 my @path = split(/\//, $filename);
 
 for(my $i = 0; $i <= $#path; $i++)
  {
  if ($directories{$path[$i]})
   {
   $study = $directories{$path[$i]};
   }
  }

 if (!$study)
  {
  my($file, $dirs, $suffix) = fileparse($filename); 
  $dirs =~ s/^\///;
  #print "$dirs\n";
  $study = $directories{$dirs};
  }

 if (!$study) # all attempts at lookup have not worked, assign a name
  {
  $study = 'STUDY';
  }

 return $study; 
 }
 
sub get_time
 {
 my $y = localtime->year() + 1900;
 my $d = localtime->mday();
 my $m = localtime->mon() + 1;
 my $h = localtime->hour();
 my $min = localtime->min();
 my $dt = $d.'-'.$m.'-'.$y; #.'.'.$h.'h'.$min.'m';
 return $dt;
 }
 
sub detect_format
 {
 my $file = $_[0];
 my $gzipped = checkgz($file);
 my $format; 
 my $z;
 my $dose = 0;
 my $oxf = 0;
 my $chrom = 0;
 
 if ($gzipped eq -1)
  {
  $format = 'ZIPPED';
  return $format;
  }
 elsif ($gzipped)
  {
  $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
 else
  {
  open $z, "$file" or die $!;
  }
 
 my $header = <$z>;
 chomp $header;
 my @temp = split(/\t/, $header); 
 
 for (my $i = 0; $i <= $#temp; $i++)
  {
  if ($temp[$i] eq 'Dose1')
   {
   $dose = 1;
   }
  if ($temp[$i] eq 'SNPID')
   {
   $oxf = 1;
   }
  if ($temp[$i] eq '#CHROM')
   {
   $chrom = 1;
   }
  }
 
 #if ($temp[0] eq 'SNPID')#changed to allow column order to be random
 if ($oxf)
 {
  $format = 'OXFORD';
  }
 #elsif ($temp[0] eq '#CHROM')
 elsif ($chrom)
  {
  my $dataline = <$z>;
  chomp $dataline;
  my @entries = split(/\t/, $dataline);
  my @info = split(/;/, $entries[7]);
  
  for (my $i = 0; $i <= $#info; $i++)
   {
   if ($info[$i] =~ /^R2\=.*/)
    {
    $format = 'MINIMAC';
    }
   elsif ($info[$i] =~ /^INFO\=.*/)
    {
    $format = 'SANGER';
    }
   }
  if (!$format)
   {
   $format = 'UNKNOWN';
   }
  }
 elsif ($temp[0] eq 'SNP' and $dose)
  {
  $format = 'MINIMAC2';
  }
 else
  {
  $format = 'UNKNOWN';
  }
 close $z; 
 return $format; 
 }
 
sub checkgz
 {
 my $file = $_[0];
 my @filecomponents = split(/\./, $file);
 my $zipped = 0;
 
 if ($filecomponents[$#filecomponents] eq 'gz')
  {
  $zipped = 1;
  }
 elsif ($filecomponents[$#filecomponents] eq 'zip')
  {
  print L "WARNING: .zip format is not supported, skipping $file\n";
  $zipped = -1;
  } 
 else
  {
  $zipped = 0;
  }
 return $zipped; 
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
 
sub get_header
 {
 my $line = $_[0];
 my $header = $_[1];
 my $column;
 
 my @allheaders = split(/\s+/, $line);
 for (my $i = 0; $i <= $#allheaders; $i++)
  {
  if ($header eq $allheaders[$i])
   {
   $column = $i;
   }
  }
 return $column; 
 }
 
sub get_chr
 {
 my $file = $_[0];
 my $column = $_[1];
 my $chr;
 my $z;
 my $idcol = 0; #default is first column
 my $gzipped = checkgz($file);
 my $flag = 0;
 
 if ($gzipped)
  {
  $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
 else
  {
  open $z, "$file" or die $!;
  }

 my $titleline = <$z>;
 chomp $titleline;
 my @titles = split(/\t/, $titleline);
 
 if ($column > $#titles) #if column index passed is too large, try finding the id column, based on format, all are currently column 0
  {
  $flag = 1;
  $column = 0;
  my $format = detect_format($file);
  if ($format eq 'OXFORD' or $format eq 'MINIMAC2' or $format eq 'SANGER' or $format eq 'MINIMAC')
   {
   $idcol = 0;
   }
  }
 
 my $line = <$z>;
 chomp $line;
 my @temp = split(/\t/, $line);
   
 if ($flag or $temp[$column] eq 'NA') #using ID as chr identifier
  {
  $temp[$idcol] =~ /^(.*?)\:.*$/;
  $chr = $1;    
  if (!$chr) #no match as SNP is directly typed (so no position added to id)
   {
   $chr = $temp[$idcol];
   }
  }
  
 else
  {
  $chr = $temp[$column];
  }
  
 close $z;
 $chr =~ s/^0//;
 print L "$file is Chromosome $chr\n";
 return $chr;
 }

sub process_minimac
 {
 my $file = $_[0];
 my $study = $_[1];
 my $datadir = $_[2];
 my $outdir = $_[3];
 #my @allfiles = @{$_[4]};
 my $match = 0;
 my $nomatch = 0;
 
 my $afflag = 'AF';
 my $infoflag = 'R2';
 my $af;
 my $info;
 my $an;
 
 my $gzipped = checkgz($file);
 my $z;
  
 if ($gzipped)
  {
  $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
 else
  {
  open $z, "$file" or die $!;
  }
 my $titleline = <$z>;
 chomp $titleline;
 
 my $chrcol = get_header($titleline, '#CHROM');
 my $poscol = get_header($titleline, 'POS');
 my $infocol = get_header($titleline, 'INFO');
 my $refcol =  get_header($titleline, 'REF');
 my $altcol = get_header($titleline, 'ALT');
 my $idcol = get_header($titleline, 'ID');
 
 my $chr = get_chr($file, $chrcol);
 
 my $newfile = $datadir.'/'.$study.'.chr'.$chr.'.txt';
 my $outfile = set_outfilename($newfile);
 open OUT, ">$outfile" or die $!;
 
 while (<$z>)
  {
  chomp;
  my @temp = split/\s+/;
  my @info = split(/\;/, $temp[$infocol]);
 
  for (my $i = 0; $i <= $#info; $i++)
   {
   if ($info[$i] =~ /$infoflag\=(.*)/)
    {
    $info = $1;
    }
   if ($info[$i] =~ /^$afflag\=(.*)/)
    {
    $af = $1;
    } 
   }

  if ($info < 0) #required as some info scores are negative
   {
   $info = 0;
   }
  #check id is present in HRC (%altaf) and add Alternate AF to the data file
  if (!$af)
   {
   print L "WARNING: allele ";  
   }
  my $id = $temp[$chrcol].':'.$temp[$poscol].'_'.$temp[$refcol].'_'.$temp[$altcol];
  print OUT "$temp[$poscol]\t$info\t$temp[$chrcol]\t$af\t$temp[$refcol]\t$temp[$altcol]"; 
  
  if ($altaf{$id})
   {
   print OUT "\t$altaf{$id}\n";
   $match++;
   }
  else
   {
   print OUT "\t\n";
   $nomatch++;
   }
  }
 $outfile =~ /.*\/(.*?)\.txt$/;
 my $countsfile = $outdir.'/COUNTS-'.$1.'.txt';
 open O, ">$countsfile" or die;
 print O "$match\n$nomatch\n";
 close O; 
 close $z;
 close OUT;
 return $outfile;
 }  
 
sub process_sanger
 {
 my $file = $_[0];
 my $study = $_[1];
 my $datadir = $_[2];
 my $outdir = $_[3];
 my $acflag = 'AC';
 my $allelenumber = 'AN';
 my $infoflag = 'INFO';
 my $ac;
 my $info;
 my $an;
 my $z;
 my $match = 0;
 my $nomatch = 0;
 
 my $gzipped = checkgz($file);
 my $af;
 
 if ($gzipped)
  {
  $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
 else
  {
  open $z, "$file" or die $!;
  }
 my $titleline = <$z>;
 chomp $titleline;
 
 my $chrcol = get_header($titleline, '#CHROM');
 my $poscol = get_header($titleline, 'POS');
 my $infocol = get_header($titleline, 'INFO');
 my $refcol =  get_header($titleline, 'REF');
 my $altcol = get_header($titleline, 'ALT');
 my $chr = get_chr($file, $chrcol);

 my $newfile = $datadir.'/'.$study.'.chr'.$chr.'.txt';
 my $outfile = set_outfilename($newfile);
 open OUT, ">$outfile" or die $!;
 
 while (<$z>)
  {
  chomp;
  my @temp = split/\s+/;
  my @info = split(/\;/, $temp[$infocol]);
  for (my $i = 0; $i <= $#info; $i++)
   {
   if ($info[$i] =~ /$infoflag\=(.*)/)
    {
    $info = $1;
    }
   if ($info[$i] =~ /$acflag\=(.*)/)
    {
    $ac = $1;
    }
   if ($info[$i] =~ /$allelenumber\=(.*)/)
    {
    $an = $1;
    } 
   }
  
  if ($an)
   {
   #print "$allelenumber\n";
   $af = sprintf("%0.7g", $ac/$an);
   }
  else
   {
   $af = 0;
   }
  
  if ($info < 0)
   {
   $info = 0;
   }
   
  #check id is present in HRC (%altaf) and add Alternate AF to the data file
  my $id = $temp[$chrcol].':'.$temp[$poscol].'_'.$temp[$refcol].'_'.$temp[$altcol];
  print OUT "$temp[$poscol]\t$info\t$temp[$chrcol]\t$af\t$temp[$refcol]\t$temp[$altcol]";
  
  if ($altaf{$id})
   {
   print OUT "\t$altaf{$id}\n";
   $match++;
   }
  else
   {
   print OUT "\t\n";
   $nomatch++;
   }
  }
 $outfile =~ /.*\/(.*?)\.txt$/;
 my $countsfile = $outdir.'/COUNTS-'.$1.'.txt';
 open O, ">$countsfile" or die;
 print O "$match\n$nomatch\n";
 close O; 
 close $z;
 close OUT;
 return $outfile;
 } 
 
sub set_outfilename
 {
 my $dirfile = $_[0];
 
 if (-e $dirfile) #file already exists
  {
  print L "$dirfile already exists\n";
  if ($dirfile =~ /^(.*)\/(\d+?)\.(.*?)\.chr(.*?)\.txt$/)
   {
   my $counter = $2+1;
   $dirfile = $1.'/'.$counter.'.'.$3.'.chr'.$4.'.txt';
   $dirfile = set_outfilename($dirfile);
   }
  elsif ($dirfile =~ /^(.*)\/(.*?)\.chr(.*?)\.txt$/)
   {
   $dirfile = $1.'/1.'.$2.'.chr'.$3.'.txt';
   $dirfile = set_outfilename($dirfile);
   }
  }
 return $dirfile; 
 }
 
 
sub process_oxford
 {
 my $file = $_[0];
 my $study = $_[1];
 my $datadir = $_[2];
 my $outdir = $_[3];
 my $z;
 my $gzipped = checkgz($file);
 my $af = 0;
 my $match = 0;
 my $nomatch = 0;
 
 if ($gzipped)
  {
  $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
 else
  {
  open $z, "$file" or die $!;
  }
 my $titleline = <$z>;
 chomp $titleline;
 
 
 my $chrcol = get_header($titleline, 'chromosome');
 my $poscol = get_header($titleline, 'position');
 my $infocol = get_header($titleline, 'information');
 my $refcol =  get_header($titleline, 'A_allele');
 my $altcol = get_header($titleline, 'B_allele');
 my $afcol = get_header($titleline, 'MAF');
 my $idcol = get_header($titleline, 'SNPID');
 my $minorallelecol = get_header($titleline, 'minor_allele');
 
 my $chr = get_chr($file, $chrcol);
  
 my $newfile = $datadir.'/'.$study.'.chr'.$chr.'.txt';
 my $outfile = set_outfilename($newfile);
 open OUT, ">$outfile" or die $!;
 
 while (<$z>)
  {
  chomp;
  my @temp = split/\t/;
  
  if ($temp[$minorallelecol] eq $temp[$altcol])
   {
   $af = $temp[$afcol];
   }
  else
   {
   $af = 1 - $temp[$afcol];
   }
  
  if ($temp[$infocol] < 0) #required as there are some negative info scores
   {
   $temp[$infocol] = 0;
   }
   
  print OUT "$temp[$poscol]\t$temp[$infocol]\t$chr\t$af\t$temp[$refcol]\t$temp[$altcol]"; 
  
  #check id is present in HRC (%altaf) and add Alternate AF to the data file
  #my $id = $temp[$idcol]; DO not use SNPID column as only contains chromosome when SNP is directly typed
  my $id = $chr.':'.$temp[$poscol].'_'.$temp[$refcol].'_'.$temp[$altcol];
  
 
  if ($altaf{$id})
   {
   print OUT "\t$altaf{$id}\n";
   $match++;
   }
  else
   {
   print OUT "\t\n";
   $nomatch++;
   }
  }
 $outfile =~ /.*\/(.*?)\.txt$/;
 my $countsfile = $outdir.'/COUNTS-'.$1.'.txt';
 open O, ">$countsfile" or die;
 print O "$match\n$nomatch\n";
 close O;
 
 close $z; 
 close OUT; 
 return $outfile;
 }
 
 
sub plot_manh_info
 {
 my $file = $_[0];
 my $study = $_[1];
 my $dir = $_[2];
 my $z;
 my @data;
 my @ids;
 my $xmax = 0;
 my $ymax = 0;
 my $i = 0;
 my $chr;
 my $old = 0;
 my $clr = 'dblue';
 my $counter = 0;
 
 my $gzipped = checkgz($file);

 if ($gzipped eq -1)
  {
  return -1;
  }
 elsif ($gzipped)
  {
  $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
 else
  {
  open $z, "$file" or die $!;
  }
  
 #read data file
 while (<$z>)
  {
  my $j = 0;
  chomp;
  my @temp = split/\t/;
  
  if ($old)
   {
   if ($temp[0] - $old > 950000)
    {
    $counter++
    }
   $old = $temp[0]; 
   }

  if ($i == 0)  #sets last position to check for gaps in the chromosome, should find at most 1 (centromere)
   {
   $old = $temp[0];
   }
   
  $data[0][$i] = $temp[0];
  $data[1][$i] = $temp[1];
  if ($chr)
   {
   if ($chr ne $temp[2])
    {
    print L "WARNING: More than one chromsome in $file\n";
    }
   }
  else
   {
   $chr = $temp[2];
   }
  $i++;
  }

 print L "Read $i Lines $file\n";
 $ymax = 1.1;
 $xmax = $data[0][$i-1];
 print L "Max values Read\n x: $xmax\n y: $ymax\n";
 
 if ($counter > 1)
  {
  print L "WARNING: $file seems to have missing chunks\n";
  $clr = 'red';
  }
 
 #derive the new file name, accounting for any duplicate names
 $file =~ /.*\/(.*?)\.txt$/;
 my $imagename = $dir.'/Manhattan.Info-'.$1.'.png';
 
 my $skip = int($xmax/10);
 $i = 5;
  
 #print "plotting\n";
 my $graph = GD::Graph::points->new(2000, 1000);
 $graph->set_title_font('14459.ttf', 18);
 $graph->set_x_label_font('14459.ttf', 18);
 $graph->set_y_label_font('14459.ttf', 18);
 $graph->set_x_axis_font('14459.ttf', 14);
 $graph->set_y_axis_font('14459.ttf', 14);
 
 $graph->set
  ( 
  x_label           => 'Chromosome '.$chr,
  y_label           => "Info Score",
  y_max_value       => $ymax,
  y_min_value       => 0,
  x_tick_number     => $i,
  transparent       => 0,
  x_max_value       => $xmax,
  x_min_value       => 0,
  marker_size       => 0.01,
  t_margin          => 20,
  b_margin          => 40,
  l_margin          => 10,
  r_margin          => 40,
  x_label_skip      => $skip,
  x_tick_offset     => $skip,
  dclrs             => [ $clr ],
  labelclr          => 'black',
  axislabelclr      => 'black',
  legendclr         => 'black',
  textclr           => 'black',
  fgclr             => 'black'
  )
 or die $graph->error;
 #print "$imagename\n";
 my $gd = $graph->plot(\@data) or die $graph->error;
 open(IMG, ">$imagename") or die $!;
 binmode IMG;
 print IMG $gd->png;
 close IMG;
 return $imagename;
 } 
 
 
sub plot_af
 {
 my $file = $_[0];
 my $study = $_[1];
 my $dir = $_[2];
 
 my $tempfile = $dir.'/temp.txt';
 my @data;
 my $i = 0;
 my $xmax = 1;
 my $ymax = 1;
 my $chr;
 
 #sort the file, required as plotting doesn't handle true X/Y plots 
 my $args = "sort -g -k 7 $file > $tempfile";
 system ($args);
 
 #read data file
 open PL, "$tempfile" or die $!;
 while (<PL>)
  {
  my $j = 0;
  chomp;
  my @temp = split/\t/;
  if ($temp[6]) #if no entry for HRC then don't plot
   {
   $data[0][$i] = $temp[6];
   $data[1][$i] = $temp[3];
   if ($chr)
    {
    if ($chr ne $temp[2])
     {
     print L "WARNING: More than one chromsome in $file\n";
     }
    }
   else
    {
    $chr = $temp[2];
    }
   $i++;
  }
 } 
 close PL;
 
 #derive the new file name, accounting for any duplicate names
 $file =~ /.*\/(.*?)\.txt$/;
 my $imagename = $dir.'/AF-'.$1.'.png';
 
 my $skip = 0.1;
  
 my $graph = GD::Graph::points->new(1200, 1200);
 $graph->set_title_font('14459.ttf', 18);
 $graph->set_x_label_font('14459.ttf', 18);
 $graph->set_y_label_font('14459.ttf', 18);
 $graph->set_x_axis_font('14459.ttf', 14);
 $graph->set_y_axis_font('14459.ttf', 14);

 $graph->set
   ( 
   x_label           => 'HRC Allele Frequency',
   y_label           => "$study Chr $chr Allele Frequency",
   y_max_value       => $ymax,
   y_min_value       => 0,
   x_ticks           => 1,
   x_tick_number     => 10,
   transparent       => 0,
   x_max_value       => $xmax,
   x_min_value       => 0,
   marker_size       => 0.1,
   markers           => [4, 3, 8, 6],
   t_margin          => 20,
   b_margin          => 40,
   l_margin          => 10,
   r_margin          => 40,
   x_label_skip      => $skip,
   x_tick_offset     => $skip,
   dclrs             => [ qw(dpurple) ],
   labelclr          => 'black',
   axislabelclr      => 'black',
   legendclr         => 'black',
   textclr           => 'black',
   fgclr             => 'black'
   )
 or die $graph->error;
 
 my $gd = $graph->plot(\@data) or die $graph->error;
 open(IMG, ">$imagename") or die $!;
 binmode IMG;
 print IMG $gd->png;
 close IMG; 
 unlink $tempfile or print L "Couldn't remove $tempfile\n";
 return $imagename;
 }


sub plot_position_row
 {
 my $file = $_[0];
 my $study = $_[1];
 my $dir = $_[2];
  
 my @data;
 my $i = 0;
 my $xmax = 0;
 my $ymax = 0; 
 my $xmin = 0;
 
 #read data file
 open PL, "$file" or die $!;
 while (<PL>)
  {
  my $j = 0;
  chomp;
  my @temp = split/\t/;
  $data[0][$i] = $temp[0];
  $data[1][$i] = $i;
  if ($temp[0] > $xmax)
   {
   $xmax = $temp[0];
   }
  if ($i == 0)
   {
   $xmin = $temp[0]; 
   }
  $i++;
  }
 $ymax = $i;
 close PL;
  
 #derive the new file name, accounting for any duplicate names
 $file =~ /.*\/(.*?)\.txt$/;
 my $imagename = $dir.'/POSITION.ROW-'.$1.'.png';
  
 my $skip = 0.1;
   
 my $graph = GD::Graph::points->new(1800, 1200);
 $graph->set_title_font('14459.ttf', 18);
 $graph->set_x_label_font('14459.ttf', 18);
 $graph->set_y_label_font('14459.ttf', 18);
 $graph->set_x_axis_font('14459.ttf', 14);
 $graph->set_y_axis_font('14459.ttf', 14);
 
 $graph->set
  ( 
  x_label           => 'Position on chromosome',
  y_label           => "Line number",
  y_max_value       => $ymax,
  x_max_value       => $xmax,
  y_min_value       => 0,
  x_min_value       => $xmin,
  x_ticks           => 1,
  x_tick_number     => 10,
  transparent       => 0,
  x_min_value       => 0,
  marker_size       => 2,
  markers           => [7, 8, 6],
  t_margin          => 20,
  b_margin          => 40,
  l_margin          => 10,
  r_margin          => 40,
  x_label_skip      => $skip,
  x_tick_offset     => $skip,
  dclrs             => [ qw(marine) ],
  labelclr          => 'black',
  axislabelclr      => 'black',
  legendclr         => 'black',
  textclr           => 'black',
  fgclr             => 'black'
  )
 or die $graph->error;
  
 my $gd = $graph->plot(\@data) or die $graph->error;
 open(IMG, ">$imagename") or die $!;
 binmode IMG;
 print IMG $gd->png;
 close IMG; 
 return $imagename;
 }


sub plot_af_info
 {
 my $file = $_[0];
 my $study = $_[1];
 my $dir = $_[2];
 
 my $tempfile = $dir.'/temp.txt';
 my @data;
 my $i = 0;
 my $xmax = 0.5;
 my $ymax = 1;
 my $chr;
 
 #sort the file, required as plotting doesn't handle true X/Y plots 
 my $args = "sort -g -k 4 $file > $tempfile";
 system ($args);
 
 #read data file
 open PL, "$tempfile" or die $!;
 while (<PL>)
  {
  my $j = 0;
  chomp;
  my @temp = split/\t/;
  if ($temp[3] > 0.5)
   {
   $temp[3] = 1 - $temp[3];
   }
  $data[0][$i] = $temp[3];
  $data[1][$i] = $temp[1];
  if ($chr)
   {
   if ($chr ne $temp[2])
    {
    print L "WARNING: More than one chromsome in $file\n";
    }
   }
  else
   {
   $chr = $temp[2];
   }
  $i++;
  } 
 close PL;
 
 #derive the new file name, accounting for any duplicate names
 $file =~ /.*\/(.*?)\.txt$/;
 my $imagename = $dir.'/MAF.INFO-'.$1.'.png';
 
 my $skip = 0.1;
  
 my $graph = GD::Graph::points->new(1800, 1200);
 $graph->set_title_font('14459.ttf', 18);
 $graph->set_x_label_font('14459.ttf', 18);
 $graph->set_y_label_font('14459.ttf', 18);
 $graph->set_x_axis_font('14459.ttf', 14);
 $graph->set_y_axis_font('14459.ttf', 14);

 $graph->set
  ( 
  x_label           => 'Minor Allele Frequency',
  y_label           => "Info Score",
  y_max_value       => $ymax,
  y_min_value       => 0,
  x_ticks           => 1,
  x_tick_number     => 10,
  transparent       => 0,
  x_max_value       => $xmax,
  x_min_value       => 0,
  marker_size       => 2,
  markers           => [7, 8, 6],
  t_margin          => 20,
  b_margin          => 40,
  l_margin          => 10,
  r_margin          => 40,
  x_label_skip      => $skip,
  x_tick_offset     => $skip,
  dclrs             => [ qw(marine) ],
  labelclr          => 'black',
  axislabelclr      => 'black',
  legendclr         => 'black',
  textclr           => 'black',
  fgclr             => 'black'
  )
 or die $graph->error;
 
 my $gd = $graph->plot(\@data) or die $graph->error;
 open(IMG, ">$imagename") or die $!;
 binmode IMG;
 print IMG $gd->png;
 close IMG; 
 unlink $tempfile or print L "Couldn't remove $tempfile\n";
 return $imagename;
 }


sub plot_counts
 {
 my $file = $_[0];
 my $study = $_[1];
 my $dir = $_[2];
 my $type = $_[3];
 
 my @data;
 my $i = 0;
 my $xmax = 100;
 my $ymax = 100;
 my $chr;
 my $cl = 'dgreen';
 
 #print "$file\n"; 
 #read data file
 open PL, "$file" or die $!;
 while (<PL>)
  {
  my $j = 0;
  chomp;
  my @temp = split/\t/;
  if ($temp[2])
   {
   $data[0][$i] = $temp[0];
   $data[1][$i] = $temp[2];
   $i++;
   }
  } 
 close PL;
 
 $file =~ /.*\/.*\.chr(.*?)\.txt$/;
 $chr = $1;
 
 my $graph = GD::Graph::hbars->new(1800, 1200);
 $graph->set_title_font('14459.ttf', 18);
 $graph->set_x_label_font('14459.ttf', 18);
 $graph->set_y_label_font('14459.ttf', 18);
 $graph->set_x_axis_font('14459.ttf', 10);
 $graph->set_y_axis_font('14459.ttf', 10);

 if ($type eq 'INFO Score')
  {
  $cl = 'dgreen';
  }
 elsif ($type eq 'Alt Allele Frequency')
  {
  $cl = 'dred';
  }
 else
  {
  $cl = 'dyellow';
  }
  
 $graph->set
  ( 
  x_label           => $type,
  y_label           => 'Percentage',
  x_ticks           => 1,
  transparent       => 0,
  t_margin          => 20,
  y_label_position  => 0.85,
  x_label_position  => 0.85,
  b_margin          => 40,
  y_max_value       => $ymax,
  #x_max_value       => $xmax,
  x_min_value       => 0,
  l_margin          => 10,
  r_margin          => 40,
  shadow_depth      => -5,
  bar_spacing       => 20,
  dclrs             => [ $cl ],
  labelclr          => 'black',
  axislabelclr      => 'black',
  legendclr         => 'black',
  textclr           => 'black',
  fgclr             => 'black'
  )
 or die $graph->error;
 
 #derive the new file name, accounting for any duplicate names
 $type =~ s/\s/\./g; 
 
 $file =~ /.*\/(.*?)\.txt$/;
 my $imagename = $dir.'/'.$type.'.'.$1.'.png';

 my $gd = $graph->plot(\@data) or die $graph->error;
 open(IMG, ">$imagename") or die $!;
 binmode IMG;
 print IMG $gd->png;
 close IMG; 
 return $imagename;
 }

sub usage
 {
 print "\nUsage:\nic -d <directory> -r <Reference panel> [-h | -g -p <population>] [-f <Mapping file>] [-o <output path>]\n";
 print "\n\n";
 printusage("-d --directory", "Directory", "Top level directory containing either one set of per chromosome files, or multiple directories of per chromosome files");
 printusage("-f --file", "Mapping file", "Mapping file of directory name to cohort name");
 printusage("-r --ref", "Reference panel", "Path to reference Panel summary file, either 1000G or HRC");
 printusage("-h --hrc", "", "Flag to indicate Reference panel file given is HRC, defaults to HRC if not supplied");
 printusage("-g --1000g", "", "Flag to indicate Reference panel file given is 1000G");
 printusage("-p --pop", "Population", "Population to check frequency against, 1000G only. Default ALL, options ALL, EUR, AFR, AMR, SAS, EAS");
 print "\n\n";
 }

sub get_width
 {
 my $output = `stty size`;
 my @rowcols = split(/\s/, $output);
 my $cols = $rowcols[1];
 if (!$cols)
  {
  $cols = 80;
  }
 return $cols;
 } 
 
sub printusage
 {
 my $option = $_[0];
 my $file = $_[1];
 my $text = $_[2];
 
 my $cols = getwidth(); 

 printf("%-18s", $option);
 printf("%-15s", $file);
 my $textwidth = $cols - 36; #33 plus 3 spaces
 my $newtext = substr ($text, 0, $textwidth);
 
 if (length($text) > $textwidth)
  {
  printf("  %*s", $textwidth, $newtext);
  for (my $i = $textwidth; $i < length($text); $i+= $textwidth)
   {
   my $newtext = substr ($text, $i, $textwidth);
   $newtext =~ s/^\s//;
   printf("%33s", " ");
   printf("   %*s", -$textwidth, $newtext);
   }
  }
 else
  {
  printf("  %*s", -$textwidth, $newtext);
  }
 print "\n"; 
 } 

sub read_kg
 {
 my $file = $_[0];
 my $pop = $_[1];
 my %altaf;

 open IN, "$file" or die $!;
 my $header = <IN>;
 chomp $header;
 
 my $freqcol = get_header($header, $pop);
 my $typecol = get_header($header, 'TYPE');
 
 print "Found Frequency for population $pop in column $freqcol\n";
 
 if (!$freqcol)
  {
  print "ERROR: Population specified, $pop not found in $file\n";
  return -1;
  }
 
 while (<IN>)
  {
  chomp;
  my @temp = split/\s+/;
  my $refaltchrpos = $temp[1].':'.$temp[2].'_'.$temp[3].'_'.$temp[4];
  $altaf{$refaltchrpos} = $temp[$freqcol];
  }
 
 close IN;
 return %altaf;
 }
  
sub plot_summary
 {
 my $file = $_[0];
 my $path = $_[1];
 my @data;
 my $i = 0;
 my $ymax = 100;
 my @cl = [ qw(orange gold lorange) ];

 #read data file
 open PL, "$file" or die $!;

 while (<PL>)
  {
  chomp;
  my @temp = split/\t/;
  if ($temp[2])
   {
   $data[0][$i] = $temp[0]; # cohort name
   $data[1][$i] = $temp[1]/$temp[3]*100; # counts if 0.4 info score used
   $data[2][$i] = ($temp[1]-$temp[2])/$temp[3]*100; # extra counts if 0.3 info score used
   $i++;
   }
  } 
 close PL;
 my $graph = GD::Graph::hbars->new(1800, 1200);
 $graph->set_title_font('14459.ttf', 18);
 $graph->set_x_label_font('14459.ttf', 18);
 $graph->set_y_label_font('14459.ttf', 18);
 $graph->set_x_axis_font('14459.ttf', 10);
 $graph->set_y_axis_font('14459.ttf', 10);
   
 $graph->set
   ( 
   x_label           => 'Cohort',
   y_label           => 'Percent of variants with given info score',
   x_ticks           => 1,
   transparent       => 0,
   t_margin          => 20,
   y_label_position  => 0.5,
   x_label_position  => 0.5,
   b_margin          => 40,
   y_max_value       => $ymax,
   x_min_value       => 0,
   l_margin          => 10,
   r_margin          => 40,
   shadow_depth      => -10,
   bar_spacing       => 20,
   dclrs             => @cl,
   labelclr          => 'black',
   axislabelclr      => 'black',
   legendclr         => 'black',
   textclr           => 'black',
   fgclr             => 'black',
   cumulate          => 1
   )
 or die $graph->error;

 my $imagename = $path.'/All-Cohort-Info-Summary.png';

 my $gd = $graph->plot(\@data) or die $graph->error;
 my $white = $gd->colorAllocate(255,255,255);
 my $black = $gd->colorAllocate(0,0,0);
 my $orange = $gd->colorAllocate(255,120,0);
 my $yellow = $gd->colorAllocate(230,230,0);
 
 my $fontpath = '14459.ttf';
 $gd -> rectangle(1450, 90, 1715, 200, $black);
 $gd -> filledRectangle(1490, 100, 1550, 140, $orange);
 $gd -> filledRectangle(1490, 150, 1550, 190, $yellow);
 # my @bounds = $gd->stringFT($black, $fontpath, 20, 0, 1490, 85, 'Legend');
 my @bounds = $gd->stringFT($black, $fontpath, 20, 0, 1560, 125, 'Info > 0.4');
 @bounds = $gd->stringFT($black, $fontpath, 20, 0, 1560, 175, 'Info > 0.3');
 if ($@)
  {
  print "String not written: $@\n";
  }
 open(IMG, ">$imagename") or die $!;
 binmode IMG;
 print IMG $gd->png;
 close IMG; 
 }

sub get_gw_counts
 {
 my $file = $_[0];
 my $cf = 0;
 my $ct = 0;
 my $tt = 0;

 $file =~ /^(.*)\/(.*?)-(.*?)\.GW\.txt$/;
  
 my $directory = $1;
 my $type = $2;
 my $study = $3;
 my $flag = 0;
 
 open IN, "$file" or die $!;
 while (<IN>)
  {
  chomp;
  my @temp = split/\t/;
  
  if ($temp[0] ne 'Total')
   {
   if ($temp[0] >= 0.4)
    {
    $cf += $temp[1];
    }
   if ($temp[0] >= 0.3)
    {
    $ct += $temp[1];
    }
   $tt += $temp[1];
   }
  }
 close IN; 
 my $line = "$study\t$ct\t$cf\t$tt\n";
 return $line;
 }
