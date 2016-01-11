#!/usr/bin/perl

=head1 NAME

seq_unique.pl - removes sequences of organisms with few sequences

=head1 USAGE

head -1 input.fa
>OPRX_RAT source=uniprot;query=OPRX_HUMAN;organism=Rattus norvegicus;taxid=10116

perl seq_min_group.pl -t 100 input.fa > output.fa

=head1 DESCRIPTION

Organisms with not a full GPCR family are filtered out.

The number of sequences for each organism is counted.
And sequences are printed of organism with are above or on the limit

=head1 OPTIONS

=item t

Minimum number of sequences belonging to a organism.
If less than organism is skipped.

Default 1 (1=no filtering).
Sets $min_nr_seqs_per_species.

=item d

Debug. If true then more verbose.
Default is false.

=head1 AUTHOR

Stefan Verhoeven <stefan.verhoeven@organon.com>

=cut

use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;
use Carp;
use Getopt::Long;

sub split_desc {
  my $d=shift;

  my %a;
  my @annots=split(/;/,$d);
  foreach my $annot (@annots) {
    my ($key,$value)=split(/=/,$annot);
    $a{$key}=$value;
  }
  return \%a;
}

sub join_desc {
  my $a=shift;

  my @d;
  while (my ($key,$value)=each %$a) {
    push(
      @d,
      "$key=$value"
    );
  }

  return join(';',@d);
}


my $stream = Bio::SeqIO->newFh(-format => 'Fasta',
                           -fh     => \*ARGV);

my $min_nr_seqs_per_species=1;

my %s; #hash of {species}->{sequence}->{id}=1
my %d;
my $seq_in_nr;
my $debug;

GetOptions(
	't=s'=>\$min_nr_seqs_per_species,
	'd=s'=>\$debug,
);

while (my $seq = <$stream> ) {
  $seq_in_nr++;

  my $id=$seq->display_id();
  my $d=$seq->desc();
  my $a=split_desc($d);
  my $organism=$a->{'organism'};

  push  
    @{$s{$organism}},
    '>'.$id.' '.$d."\n".$seq->seq()."\n"
}

my %o2;
my $seq_out_nr;
while ( my ($organism,$r)=each %s) {
  my $nr_seqs=scalar @{$r};
  if ($nr_seqs>=$min_nr_seqs_per_species) {
    print @{$r};
    $o2{$organism}=$nr_seqs;
    $seq_out_nr+=$nr_seqs;
  }
}

print STDERR "Number of sequences in: $seq_in_nr\n";
print STDERR "Number of sequences out: $seq_out_nr\n";
print STDERR "Number of organisms in: ",scalar keys(%s),"\n";
print STDERR "Number of organisms out: ",scalar keys(%o2),"\n";

print STDERR "Organism distribution :\norg|nr\n";

foreach my $i (sort { $o2{$b}<=>$o2{$a} } keys(%o2)) {
  print STDERR "$i|",$o2{$i},"\n";
}

