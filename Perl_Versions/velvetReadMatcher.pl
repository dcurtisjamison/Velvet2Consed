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
use Term::ReadKey;
use DBI;

# get parameters
my $vdir = ".";
my $fdir = ".";
my $rdir = ".";
my $odir = "contig_dir";
my $gDB = "velvet";
my $dbUser = "postgres";
my $schema = 0;
my $result = GetOptions("velvet:s" => \$vdir, 
			"fasta:s" => \$fdir,
			"rundir:s" => \$rdir,
			"outdir:s" => \$odir,
			"genome:s" => \$gDB,
			"dbusername:s" => \$dbUser,
			"showschema" => \$schema, );
if ($schema) {
  die (showSchemaSQL());
}

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

# figure out where we're putting the fasta & quality files
unless ($odir =~ /^\//) {
  my $rootdir = `pwd`;
  chomp $rootdir;
  $odir = $rootdir . '/' . $odir;
}
unless (-e $odir) {
  `mkdir $odir`;
  unless (-w $odir) {die "Can't create files in $odir\n"};
}

# open database
print "Database password for $dbUser:";
ReadMode 2; ## no echo mode for password entry, in case your cubicle mate is a spy...
my $password = ReadLine(0);
chomp $password;
ReadMode 0;
print "\n";
my $dbh = DBI->connect("dbi:Pg:dbname=$gDB",$dbUser, $password,{AutoCommit => 1});
my $retval = $dbh->do('set standard_conforming_strings = 1');
unless ($retval) {die $retval}; ## psql may not be returning 0/1 for commands?

# work on collating and storing reads
my %sample;
my %gDir;
foreach my $file (@fl) {
  unless (-e "$fdir/$file") {
    die "Couldn't find $fdir/$file\n";
  }
  $sample{$file} = `head -n 2 $fdir/$file`;
  if ($rdir eq ".") {
    # our files are local files
    $gDir{$file} = $rdir;
  } else {
    #try to get files from Gerald directory
    print "Finding Gerald directory for $file\n";
    $gDir{$file} = drillDown($rdir, $sample{$file});
    unless ($gDir{$file}) {die "No experiment directory for $file rooted at $rdir\n"};
  }
}
my $seqFile = "$vdir/Sequences";
my $fc = 0; ## fc = file counter -- indexes into the file name array to track progress
my ($seq1, $seq2) = buildExpFileNames($gDir{$fl[$fc]}, $sample{$fl[$fc]});
open (SEQ, "<$seqFile") || die "Couldn't open $seqFile\n";
open (FA, "<$fdir/$fl[$fc]") || die "Couldn't open $fdir/$fl[$fc]\n";
open (SEQ1, "<$seq1") || die "Couldn't open $seq1\n";
open (SEQ2, "<$seq2") || die "Couldn't open $seq2\n";
$fc++;

## our files exist, we can start parsing using the velvet Sequences file to coordinate us
## this would all be so much simpler if Velvet didn't tromp on the deflines...
my ($vDef, $vSeq, $iDef, $iSeq, $fDef1, $fSeq, $fQual);
my $AID = 1; #afg read numbering starts at 1 and runs sequentially with the Sequences file (which starts at 0)
while ($vDef = <SEQ>) {
  $vSeq = <SEQ>;
  unless (defined($iDef = <FA>)) {
    ## eof, we need to move to the next set of input files
    print "opening new fasta file: $fl[$fc]\n";
    open (FA, "<$fdir/$fl[$fc]") || die "Couldn't open $fdir/$fl[$fc]\n";
    ($seq1, $seq2) = buildExpFileNames($gDir{$fl[$fc]}, $sample{$fl[$fc]});
    $fc++;
    open (SEQ1, "<$seq1") || die "Couldn't open $seq1\n";
    open (SEQ2, "<$seq2") || die "Couldn't open $seq2\n";
    $iDef =<FA>;
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
  $iDef =~ s/^>//;
  store('RED', $dbh, {'AID'=>$AID, 'VID'=>$vDef, 'GAID'=>$iDef, 'Sequence'=>$fSeq, 'Quality'=>$fQual});
  $AID++;
}
close SEQ;
close FA;
close SEQ1;
close SEQ2;
## mate pairing
## TBD

# Work on asm file
my @block;
my %values;
my $scfIdx = 0;
my $pid = - 1;
my $multi = 0;

my ($file) = "$vdir/velvet_asm.afg";
open (IN, "<$file") || die "Couldn't open $file\n";
print "Working on contigs and tiling from $file\n";

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
  	## store reads, tiles and scaffolds later
  	store($label, $dbh, \%values);
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
      unless ($label eq 'RED') {
  	## don't store reads from afg -- we did that earlier
  	store($label, $dbh, \%values);
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

# Quality scoring
print "Computing contig quality scores\n";
# create some indices to make our program work a bit better
$retval = $dbh->do(createIndexSQL('Read', 'aid'));
$retval = $dbh->do(createIndexSQL('Tiling', 'contig'));
$retval = $dbh->do(createIndexSQL('Tiling', 'contig, origin')),
# get contig list from database
my $sth = $dbh->prepare(getCtgSQL());
$retval = $sth->execute();
unless ($retval) {die $retval};
my $cDef = $sth->fetchall_arrayref({});

foreach my $cRef (@$cDef) {
  # read sequence from fasta file
  open (SEQ, "<" . $cRef->{sequence}) || die "Couldn't open " . $cRef->{sequence} ."\n" ;
  my $sequence;
  while (<SEQ>) {
    chomp;
    unless (/>/) {$sequence .= $_};
  }
  close SEQ;

  # get reads from database
  $sth = $dbh->prepare(getReadSQL($cRef->{aid}));
  $retval = $sth->execute();
  unless ($retval) {die "psql returned $retval\n"};
  my ($rowRef, @sum, @count);

  while ($rowRef = $sth->fetchrow_hashref) {
    $rowRef->{quality} =~ s/^\"//; ## deal with a potential wierdness from early verison of database...
    if ($rowRef->{clr_start} > $rowRef->{clr_stop}) {
      $rowRef->{sequence} = revComp($rowRef->{sequence});
      $rowRef->{quality} = reverse($rowRef->{quality});
    }
    if ($rowRef->{origin} < 0) {
      $rowRef->{sequence} = substr($rowRef->{sequence}, 0 - $rowRef->{origin});
      $rowRef->{quality} = substr($rowRef->{quality}, 0 - $rowRef->{origin});
      $rowRef->{origin} = 0;
    }
    for (my $i = $rowRef->{origin}; $i < length($rowRef->{quality}) + $rowRef->{origin}; $i++) {
      $sum[$i] += convertSolexaToPhred(convertFromSolexaSymbol(substr($rowRef->{quality}, $i - $rowRef->{origin}, 1)));
      $count[$i]++;
    }
  }
 
  print "updating quality for contig " . $cRef->{vid} . " from $retval bases\n";
  # send quality to an output file
  my $qFile = $cRef->{sequence};
  $qFile =~ s/\.fa/\.qual/;
  open (OUT, ">$qFile") || die "Couldn't open $qFile\n";
  for (my $i = 0; $i <= length($sequence) - 1; $i++) {
    print OUT $count[$i]?int($sum[$i]/$count[$i]):0;
    if (($i+1)%70) {
      print OUT " ";
    } else {
      print OUT "\n";
    }
  }
  print OUT "\n";
  close OUT;
  my $retval = $dbh->do(updateQualSQL($cRef->{aid}, $qFile));
}

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
  if ($#dirs < 0) {
    if ($phrase eq "Firecrest") {
      return findDir("IPAR");
    } 
    die "no $phrase in $dir\n";
  }
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

# build file names for fastq files
sub buildExpFileNames {
  my ($dir, $example) = @_;
  my ($fn1, $fn2);
  if ($dir eq ".") {
    # local files
    # I assume the root name is the first thing before the colon in the defline, and the fastq
    # are named root_1.fastq and root_2.fastq
    # eg., >root:read:whatever other stuff
    my ($root) = split(":", $example);
    $fn1 = $root . "_1.fastq";
    $fn2 = $root . "_2.fastq";
  } else {
    # Find file names in an Illumina directory
    my $paired = ($example =~ /\/1/)?1:0;
    my ($mach, $run, $lane, $tile) = split(/_|:/, $example);
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
  }
  return ("$dir/$fn1", "$dir/$fn2");
}

#store records in database
sub store {
  my ($label, $dbh, $hashRef) = @_;
  my $sqlStr;
  if ($label eq 'RED') {
    # Reads extracted from the velvet Sequence file
    $sqlStr = "INSERT INTO Read (aid, vid, gaid, sequence, quality) VALUES (" . 
      join(', ', $hashRef->{'AID'},  "\'" . $hashRef->{'VID'} . "\'",  "\'" . $hashRef->{'GAID'} . "\'", "\'" . 
	   $hashRef->{'Sequence'} . "\'",  "\'" . $hashRef->{'Quality'} . "\'") . ")";
  } elsif ($label eq 'CTG') {
    # Contigs extracted from the agf file
    print "Storing Contig " . $hashRef->{'eid'} . "\n";
    my $ofile =  fastaPrint($hashRef->{'eid'}, $hashRef->{'seq'});
    $sqlStr = "INSERT INTO Contig (aid, vid, sequence, length) VALUES (" .
      join(',', $hashRef->{'iid'}, "\'" . $hashRef->{'eid'} . "\'","\'" . $ofile . "\'", length($hashRef->{'seq'})) . ")";
  } elsif ($label eq 'TLE') {
    # tiling information extracted from the afg file
    my $table = 'Tiling';
    my ($start, $stop) = split(/,/, $hashRef->{'clr'});
    $sqlStr = "INSERT INTO $table (read, contig, origin, clr_start, clr_stop) VALUES (" .
      join(',', $hashRef->{'src'}, $hashRef->{'parentID'}, $hashRef->{'off'}, $start, $stop) . ")";
  } elsif ($label eq 'SCF-TLE') {
    my ($start, $stop) = split(/,/, $hashRef->{'clr'});
    $sqlStr = "INSERT INTO Scaffold (aid, vid, contig, origin, clr_start, clr_stop) VALUES (" .
      join(',', $hashRef->{'iid'}, "\'" . $hashRef->{'eid'} . "\'", $hashRef->{'src'}, $hashRef->{'off'}, $start, $stop) . ")";
  }else {
    # skip anything else
    return;
  }
  my $rtv = $dbh->do($sqlStr);
  unless ($rtv) {die "$rtv"};
}

# SQL statement templates

sub showSchemaSQL {

  return "\nCREATE TABLE contig (\n\taid integer primary key, \n\tvid character varying(50), \n\tsequence text, \n\tquality text, \n\tlength int);\n\nCREATE TABLE read (\n\taid integer primary key, \n\tvid character varying(50), \n\tgaid character varying(50), \n\tsequence text, \n\tquality text, \n\tmate integer);\n\nCREATE TABLE scaffold (\n\taid integer primary key, \n\tvid character varying(50), \n\tcontig integer, \n\torigin integer, \n\tclr_start integer, \n\tclr_stop integer);\n\nCREATE TABLE tiling (\n\tread integer, \n\tcontig integer, \n\torigin integer, \n\tclr_start integer, \n\tclr_stop integer);\n\n";

}

sub createContigTableSQL {
  return "CREATE TABLE contig (aid integer primary key, vid character varying(50), sequence text, quality text, length int)";
}

sub createReadTableSQL {
  return "CREATE TABLE read (aid integer primary key, vid character varying(50), gaid character varying(50), sequence text, quality text, mate integer)";
}

sub createScaffoldTableSQL {
  return "CREATE TABLE scaffold (aid integer primary key, vid character varying(50), contig integer, origin integer, clr_start integer, clr_stop integer)";
}

sub createTilingTableSQL {
  return "CREATE TABLE tiling (read integer, contig integer, origin integer, clr_start integer, clr_stop integer)":
}

sub getCtgSQL {
  return "select aid, vid, sequence from Contig";
}

sub getReadSQL {
  my ($ctg) = @_;
  return "select r.gaid, r.sequence, r.quality, t.origin, t.clr_start, t.clr_stop from read as r, tiling as t where t.contig = $ctg and t.read = r.aid order by t.origin"
}

sub updateQualSQL {  
  my ($aid, $qstr) = @_;
  return "update Contig set quality = '$qstr' where aid = $aid";
}

sub createIndexSQL {
  my ($tbl, @col) = @_;
  my $idxName = join('_', @col) . "_idx";
  my $col = join(", ", @col);
  return "create index $idxName on $tbl ($col)";
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

velvetReadMatcher.pl - Create a database for the output from velvet

=head1 SYNOPSIS

velvetReadMatcher.pl [-v <velvet directory>] [-f <input fasta directory>] [-r <run folder root>] [-o <contig output directory>] [-g <assembly database>] [-d <db user name>]


=head1 DESCRIPTION

velvetReadMatcher.pl extract reads, contigs, and tiling from a velvet assembly. The program requires a Velvet assembly based on a set of Illumina GA short reads: Velvet must be run with read-tracking and amos file production options turned on. The GA pipeline must have been with a Gerald ANALYSIS sequence or ANALYSIS sequence_pairs options. The program also requires connection to a Postgres database with a schema as described in the README file. For more implementation details, see the README document and the header in the program file. We also require the strict, Getopt::Long, Term::ReadKey and DBI perl modules (again, see the README file).

All the parameters are optional. There are default values for everything, but you probably won't like them much. The parameters are as follows:
 -s prints out the current schema as SQL create table statements, suitable for cut and paste into your psql (or other) command-line interface
 -v <velvet directory> The path to the directory in which you ran velvet. It will contain the Log, Sequence, and velvet_asm.afg files. The default path is the current directory.
 -f <input fasta directory> The path to the files you used as input into the velveth command. The defaut path is the current directory.
 -r <run folder root> The path to the Illumina run folder directory. This is where the experiment folders are. The default path is the current directory.
 -o <contig output directory> The place where you want the fasta files and quality files to be stored. The default directory is "./contig_dir".
 -g <assembly database> The database to load the current assembly into. The default is 'velvet', you probably want to use your organism name.
 -d <db user name> The user name to connect to the database as. Must have write permission to the database. The default user is 'postgres'.

The goal of this program is to link the velvet output back into the Illumina namespace. It first rummages through the velvet Sequences file and the input files to get the more informative Illumina GA Pipeline names, and then grabs the quality string for the reads from the Gerald folder. These are put into the database. Next, it reads the contig and tiling information from the afg file, storing everything except the sequence in the database. The sequence is written to an individual fasta file, and the file name stored in the database. The contig base qualities are computed from the read qualities, and stored as a quality file, and the file names are stored in the database. Note that the quality values are converted from Q(Solexa) to Q(phred).

=head1 POSSIBLE TODO LIST

There is always too many things to do with a program, most of which become clear only after it has been released upon an unsuspecting world. Here are a few of the things that I can think of:
 1. More hardening to guard against misuse. 
 2. Is an average the best way to compute the quality? I could use Statistics::Descriptive to provide other potential methods.
 3. Store the sequences & quality strings in the DB rather than files: insert on long text works fine, but update truncates the command string after a certain point. fixing this will take some work.
 4. For real masochism, the quality step could be threaded and run parallel. DB connection would have to be record locking. Ick.

Note there is no mention of reworking to include other sequencer formats. First off, if you're using anything other than a GA or Solid system, you probably should be using something other than Velvet. Secondly, this script was written to support a specific sequencing project using our GA, not as a general coding exercise. Check out the README file for more ranting on this topic.

=head1 LICENSE

This software is released under the GNU General Public License (see <http://www.gnu.org/licenses/>).

=head1 AUTHOR

Curtis Jamison



                 Copyright (C) 2009
         Children's Hospital of Cincinnati
                 All rights Reserved

=cut

