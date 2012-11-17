#! /usr/bin/perl -w
#####################################################################################
#
# velvetReadMatcher.pl
# Extract reads, contigs, and tiling from a velvet assembly.
#
# Requires a velvet assembly from a set of Illumina GA short reads. The GA pipeline
# with a Gerald ANALYSIS sequence or ANALYSIS sequence_pairs options, or else the
# equivalent fastq files with the standard Illumina naming conventions must be 
# produced by other means. See the buildExpFileNames subroutine for further details
# about the file name.
# 
# In addition, Velvet must be run with read-tracking and amos file  production options
# turned on. The program also requires connection to a Postgres database with a schema
# as described in the text file. The SQL commands are simplistic, so the connections
# should work with any DBMS -- simply edit the DBI->connect settings. You might also 
# need to do something the 'set standard_conforming_strings = 1' line: Postgres has
# string interpolation which does bad things if your quality strings contain a 
# backslash. The set command turns this off for Postgres. If your database does this
# sort of non-standard interpolation, you will want to find an equivalent command.
#
#####################################################################################
#
# Copyright 2009, Children's Hospital of Cincinnati
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
###################################################################################

use strict;
use Getopt::Long;
use DBI;
use Term::ReadKey;

# get parameters
my $vdir = ".";
my $fdir = ".";
my $rdir = ".";
my $odir = "foobar-dummy-dir";
my $result = GetOptions("velvet:s" => \$vdir, 
			"fasta:s" => \$fdir,
			"rundir:s" => \$rdir,
			"outdir:s" => \$odir, );

print "FILES ". localtime(). "\n";

# get input files
my $log = "$vdir/Log";
open (LOG, "<$log") || die "Can't open $log";
my @vhCommand;
my @fl;
while (<LOG>) {
  chomp;
  if (/velveth/) {
    @vhCommand = split;
  }
}
close LOG;
unless ($#vhCommand > 0) {die "Fatal: No velveth command in $log!\n"};

splice(@vhCommand, 0, 3);
while ($_ = shift(@vhCommand)) {
  unless (/^-/) {
    # these are our files
    push(@fl, $_);
  }
}

# work on collating and storing reads
my %sample;
my %gDir;
foreach my $file (@fl) {
  unless (-e "$fdir/$file") {
    die "Couldn't find $fdir/$file\n";
  }
  $sample{$file} = `head -n 2 $fdir/$file`;
  print "Finding Gerald directory for $file\n";
  $gDir{$file} = drillDown($rdir, $sample{$file});
  unless ($gDir{$file}) {die "No experiment directory for $file rooted at $rdir\n"};
}
my $seqFile = "$vdir/Sequences";
my $tempdir = "$vdir/tmp";
unless (-e $tempdir) {`mkdir $tempdir`};
unless (-e $tempdir) {die "$tempdir cannot be created\n"};
my $fc = 0; ## fc = file counter -- indexes into the file name array to track progress
my ($seq1, $seq2) = buildExpFileNames($gDir{$fl[$fc]}, $sample{$fl[$fc]});
open (SEQ, "<$seqFile") || die "Couldn't open $seqFile\n";
print "opening new fasta file: $fl[$fc]\n";
open (FA, "<$fdir/$fl[$fc]") || die "Couldn't open $fdir/$fl[$fc]\n";
open (SEQ1, "<$seq1") || die "Couldn't open $seq1\n";
open (SEQ2, "<$seq2") || die "Couldn't open $seq2\n";
open (RMAP, ">$tempdir/RDMap") || die "Couldn't open $vdir/RDMap\n";
unless (-e $odir) {
  $odir = $vdir;
}
open (PHD, ">$odir/phd.ball") || die "Couldn't open $odir/phd.ball\n";
$fc++;

print STDERR "READS ". localtime(). "\n";
## our files exist, we can start parsing using the velvet Sequences file to coordinate us
my ($vDef, $vSeq, $iDef, $iSeq, $fDef1, $fSeq, $fQual);
my $AID = 1; #afg read numbering starts at 1 and runs sequentially with the Sequences file (which starts at 0)
local $| = 1;
print "Prepping reads";
while ($vDef = <SEQ>) {
  $vSeq = <SEQ>;
  unless (defined($iDef = <FA>)) {
    ## eof, we need to move to the next set of input files
    print "\nopening new fasta file: $fl[$fc]\n";
    open (FA, "<$fdir/$fl[$fc]") || die "Couldn't open $fdir/$fl[$fc]\n";
    ($seq1, $seq2) = buildExpFileNames($gDir{$fl[$fc]}, $sample{$fl[$fc]});
    $fc++;
    open (SEQ1, "<$seq1") || die "Couldn't open $seq1\n";
    open (SEQ2, "<$seq2") || die "Couldn't open $seq2\n";
    $iDef =<FA>;
    print "Prepping reads";
  }
  $iSeq = <FA>;
  chomp($vDef);
  chomp($vSeq);
  chomp($iDef);
  chomp($iSeq);
  $iSeq =~ s/N/A/g;
  unless ($vSeq =~ /$iSeq/) {die "Files are unsynched: $vSeq ne $iSeq\n"};
  
  if ($iDef =~ /\/2/) {
    # get the qualities for a second read
    do {
      $fDef1 = <SEQ2>;
      chomp($fDef1);
      $fDef1 =~ s/@/>/;
    } until ($fDef1 eq $iDef);
    $fSeq = <SEQ2>;
    $fQual = <SEQ2>;
    $fQual = <SEQ2>;
  } else {
    # get the qualities for a first read or for a single read
    do {
      $fDef1 = <SEQ1>;
      chomp($fDef1);
      $fDef1 =~ s/@/>/;
    } until ($fDef1 eq $iDef);
    $fSeq = <SEQ1>;
    $fQual = <SEQ1>;
    $fQual = <SEQ1>;
  }
  chomp($fSeq);
  chomp($fQual);
  if ($vDef =~ /SEQUENCE_(\d*)_length/) {
    $vDef = $1;
  }
  
  ## phd file
  $iDef =~ s/^>//;
  $iDef =~ s/\//_/;
  my $l = length($fSeq);
  my $chromat = "CHROMAT_FILE: $iDef";
  my $time = "TIME: ".localtime();
  my $chem = "CHEM: Solexa";
  my $outfile =  "$iDef.phd.1";
  print RMAP "$AID\t" . join(' ','RD', $iDef, $l, '0', '0') . "?$fSeq?QA 1 $l 1 $l?DS $chromat PHD_FILE: $outfile $time $chem\n";
  print PHD "BEGIN_SEQUENCE $iDef 1\nBEGIN_COMMENT\n" . "CHROMAT_FILE: $iDef\nCALL_METHOD: Bustard\nQUALITY_LEVELS: 50\n" .
    "$time\nTRACE_ARRAY_MIN_INDEX: 0\nTRACE_ARRAY_MAX_INDEX: " . ($l - 1) . "\n$chem\n" . "END_COMMENT\nBEGIN_DNA\n";
  for (my $i = 0; $i < $l; $i++) {
    my $quality = ord(substr($fQual, $i, 1)) - 64;
    $quality = ($quality<0)?0:$quality;
    print PHD join(' ', substr($fSeq, $i, 1), $quality, $i + 1) . "\n";
  }
  print PHD "END_DNA\nEND_SEQUENCE\n";
  $AID++;
  unless ($AID%10000) {print "."};
}
$| = 0;
print "\n";
close SEQ;
close FA;
close SEQ1;
close SEQ2;
close RMAP;
close PHD;

# Work on asm file
my @block;
my %values;
my $scfIdx = 0;
my $pid = - 1;
my $multi = 0;

print STDERR "INDEX ". localtime(). "\n";
print "indexing reads\n";
open (RDMAP, "<$tempdir/RDMap") || die "Couldn't open RDMap\n";
my @index = (pack('J', 0));
while (<RDMAP>) { 
  push (@index, (pack 'J', tell RDMAP));
}
pop @index;
my $index = join ('', @index);
print 'Size of index: ', length $index, "\n";
close RDMAP;

print STDERR "CONTIGS ". localtime(). "\n";
my ($file) = "$vdir/velvet_asm.afg";
open (IN, "<$file") || die "Couldn't open $file\n";
print "Working on contigs and tiling from $file\n";
my (%ctg, %reads, %offset, $ctgCount, $readCount);
# parse
while (<IN>) {
  chomp;
  
  #start of a block
  if (s/^{//) {
    if ($#block > -1) {
      ## preflush the outer block and set parent id
      my $label = pop(@block);
      if ($label eq 'SCF') {
 	## Save the scaffold info to include with tiling
 	$values{'iid'} = "s_" . $values{'eid'};
      }
      if ($label eq 'CTG') {
  	## flush previous contig
	if (defined($ctg{eid})) {
	  outputContig(\%ctg, \%reads, \%offset, $index);
	  $ctgCount++;
	  undef %ctg;
	  undef %reads;
	  undef %offset;
	}
	$ctg{eid} = $values{eid};
	$ctg{seq} = $values{seq};
      }
      $pid = $values{'iid'};
      undef %values;      
    }
    push(@block, $_);
    next;
  }
  
  #end of a block
  if (/^}/) {
    my $label = pop(@block);
    if ($label) {
      if (($pid =~ /s_/) || ($pid > -1)) {
  	$values{'parentID'} = $pid;
      }
      if ($pid =~ /s_(\d*)/) {
 	# scaffold tile; reset the parent values;
 	$values{'iid'} = $scfIdx++;
 	$values{'eid'} = $1;
 	$label = 'SCF-TLE';
      }
      if ($label eq 'TLE') {
	my ($start, $stop) = split(/,/, $values{'clr'});
	$reads{$values{src}} = join(" ", "AF", $values{src}, ($start)?'C':'U', $values{off} + 1) . "\n";
	$offset{$values{src}} = $values{off};
	$readCount++;
      }
      undef %values;
    } else {
      $pid = -1;
    }
    next;
  }
  
  #multi-line values
  if (/^seq|^qlt/) {
    s/://g;
    $multi = $_;
    next;
  }
  if (/\./) {
    $multi = 0;
    next;
  }
  if ($multi) {
    $values{$multi} .= $_;
    next;
  }
  
  #everything else
  my ($key, $value) = split(/:/);
  $values{$key} = $value;
}
close IN;
if (defined($ctg{eid})) {
  outputContig(\%ctg, \%reads, \%offset, $index);
  $ctgCount++;
}

print STDERR "ACE ". localtime(). "\n";
my $acefile = "$odir/velvet.ace.";
my $suffix = 1;
while (-e "$acefile$suffix") {
  $suffix++;
}
$acefile .= $suffix;
open (ACE, ">$acefile") || die "Couldn't open $acefile\n";
print ACE "AS $ctgCount $readCount\n";
close ACE;
`cat $tempdir/*.tmp >> $acefile`;
`rm -r $tempdir`;
print STDERR "DONE ". localtime(). "\n";

#### Subroutines #####
# find the Gerald directory
sub drillDown {
  my ($dir, $example) = @_;
  my $exp = `ls $dir`;
  my ($mach, $run) = split(/_|:/, $example);
  $mach =~ s/>//;
  if ($exp =~ /(\d*_${mach}_0*${run})/) {
    my $edir = $rdir . "/" . $1 . "/Data";
    $edir = findDir($edir, "Firecrest");
    $edir = findDir($edir, "Bustard");
    $edir = findDir($edir, "GERALD");
    return $edir;
  }
  return 0;
}

sub findDir {
  # matches a directory phrase vs an ls of possibles
  my ($dir, $phrase) = @_;
  my @dirs = split(/\n/, `ls -d $dir/*$phrase*`);
  if ($#dirs < 0) {die "no $phrase in $dir\n"};
  if ($#dirs == 0) {return $dirs[0]};
  return choose(\@dirs);
}

sub choose {
  # present the user with a way to choose a directory from multiple possibilties
  my ($aref) = @_;
  print "Choose one of the following directories:\n";
  my $i = 0;
  foreach my $dir (@$aref) {
    print "\t", $i++, "\t", $dir, "\n";
  }
  print "\t", $i, "\tAbort\n";
  print ">";
  ReadMode 4;
  my $char = ReadKey(0);
  ReadMode 0;
  print " $char\n";
  if ($char == $i) {die "Program aborted by user\n"};
  return($$aref[$char]);
}

# build file names for Illumina fastq files
sub buildExpFileNames {
  my ($dir, $example) = @_;
  my $paired = ($example =~ /\/1/)?1:0;
  my ($mach, $run, $lane, $tile) = split(/_|:/, $example);
  my ($fn1, $fn2);
  $fn1 = "/s_". $lane;
  $fn2 = $fn1;
  if ($paired) {
    $fn1 .= "_1";
    $fn2 .= "_2";
  }
  $fn1 .= "_sequence.txt";
  $fn2 .= "_sequence.txt";
  unless (-e "$dir/$fn1") {
    $dir .= "/Temp";
  }
  return ("$dir/$fn1", "$dir/$fn2");
}

# turn a .afg CTG block into a .ace CO block
sub outputContig {
  my ($ctgRef, $readsRef, $offsetRef, $index) = @_;
  my $file = "$tempdir/ctg_" . $ctgRef->{eid} . ".tmp";
  open (OUT, ">$file") || die "couldn't open $file\n";
  open (RDMAP, "<$tempdir/RDMap") || die "Couldn't open RDMap\n";
  my $length = length($ctgRef->{seq});
  my @keys = sort(keys(%$readsRef));
  print OUT join(" ", "CO", $ctgRef->{eid}, $length, ($#keys + 1), '0', 'U') . "\n";
  while (length($ctgRef->{seq})) {
    print OUT substr($ctgRef->{seq}, 0, 70, '') . "\n";
  }
  print OUT "\nBQ\n";
  while ($length > 0) {
    my $limit = ($length > 70)?70:$length;
    for (my $i = 0; $i < $limit; $i++) {
      print OUT "30 ";
    }
    print OUT "\n";
    $length -= 70;
  }
  print OUT "\n";
  for my $i (@keys) {
    my $afline = $readsRef->{$i};
    my $line = unpack('J', substr($index, $i*4, 4));
    seek (RDMAP, $line, 0);
    my $map = <RDMAP>;
    my (undef, undef, $read) = split(" ", $map);
    $afline =~ s/${i}/${read}/;
    print OUT $afline;
  }
  print OUT "\n";

  for my $i (@keys) {
    $i--;
    my $line = unpack('J', substr($index, $i*4, 4));
    seek (RDMAP, $line, 0);
    my $map = <RDMAP>;
    chomp($map);
    my ($id, $RDStr) = split("\t", $map);
#    unless ($id == $i-1) { die "mismatch between desired index $i and fetched line $id\n"};
    my ($RD, $seq, $QA, $DS) = split(/\?/, $RDStr);
    if ($offsetRef->{$id}) {
      $seq = revComp($seq);
    }
    print OUT "$RD\n$seq\n\n$QA\n$DS\n";
  }
  close RDMAP;
  close OUT;
}

# write sequence to a fasta file
sub fastaPrint {
  my ($ctg, $seq) = @_;
  my $file = $odir . "/Contig_" . $ctg . ".fa";
  unless (open (OUT, ">$file")) {return "NULL"};
  print OUT ">Contig_$ctg " . length($seq) . "\n";
  while (length($seq)) {
    print OUT substr($seq, 0, 70, '') . "\n";
  }
  close OUT;
  return $file;
}

# general sequence data munging routines
sub revComp {
  my ($seq) = @_;
  $seq = reverse($seq);
  $seq =~ tr/ACGT/TGCA/;
  return $seq;
}

sub convertFromSolexaSymbol {
  my ($sym) = @_;
  return (ord($sym) - 64);
}

sub convertToSolexaSymbol {
  my ($num) = @_;
  if ($num >= 62) {
    $num = 62;
  }
  return chr($num + 64);
}
sub convertFromSangerSymbol {
  my ($sym) = @_;
  return (ord($sym) - 33);
}

sub convertToSangerSymbol {
  my ($num) = @_;
  if ($num >= 93) {
    $num = 93;
  }
  return chr($num + 33);
}

sub convertSolexaToPhred {
  my ($score) = @_;
  return 10 * log(1 + 10 ** ($score / 10.0)) / log(10);
}

# perldoc

=pod

=head1 NAME


=head1 SYNOPSIS



=head1 DESCRIPTION



=head1 POSSIBLE TODO LIST


=head1 LICENSE

This software is released under the GNU General Public License (see <http://www.gnu.org/licenses/>).

=head1 AUTHOR

Curtis Jamison



                 Copyright (C) 2009
         Children's Hospital of Cincinnati
                 All rights Reserved

=cut
