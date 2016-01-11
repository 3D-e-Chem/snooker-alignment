#!/usr/bin/perl

=head1 NAME

highly_similar_reducer.pl - Remove highly similar sequences

=head1 USAGE

head -1 in.fa
>ENSOGAP00000005678 source=ensembl;organism=Otolemur garnettii;taxid=30611

perl highly_similar_reducer.pl -m 9 -s skipped.ids in.fa > out.fa


perl highly_similar_reducer.pl -i 'Homo sapiens' -m 9 -s skipped.ids in.fa > out.fa


=head1 DESCRIPTION

Removes equal and highly similar sequences.

Look within the same organism for sequence pairs which have a Levenshtein distance
of less or equal than minumum distance option.
If sequence pair has distance<=minumum distance then one sequence is skipped.

It preferes to remove the sequence which looks like trembl or ensembl.

=head1 OPTIONS

=item m

Minimum distance, if distance between 2 sequences is less or eqaul than this.
Then one of the sequences will be skipped
Optional.
Default is 0.

=item s

File with skipped ids.
Format is skipped_id\tkept_id\t$distance\t$reason\n.

Of a pair of ids one must be kept and one must be kept.
Reasons for skipping can be:
1 =  prefer shorter or alpabeticly first with seqs of same source
2 =  prefer source with lowest score

=item include

Apply reduction only to included organisms.

=back

=head1 AUTHOR

Stefan Verhoeven <stefan.verhoeven@spcorp.com>

=cut

use strict;
use warnings;
use Bio::SeqIO;
use Text::LevenshteinXS qw(distance);
use Data::Dumper;
use Getopt::Long;
use Carp;

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

my $min_dist=0;
my $skip_file;
my $skip_fh;
my @include_taxs;

GetOptions(
	'm=s'=>\$min_dist,
    's:s'=>\$skip_file,
    'include:s'=>\@include_taxs,
);

if ($skip_file) {
  open($skip_fh,'>',$skip_file);
}

my %include_taxs = map { $_ => 1 } @include_taxs;

my $stream = Bio::SeqIO->newFh(-format => 'Fasta',
                           -fh     => \*ARGV);

my %s; # {species}->{seqid}=sequence
my %desc; # {seqid}=desc lookup table of description
my %sources;

my %source_scores = (
    'refseq' => 3,
    'ensembl' => 4,
    'trembl' => 2,
    'swissprot' => 1,
    'uniprot' => 2,
    'unknown' => 9999,
);

print STDERR "Reading file\n";
while (my $seq = <$stream> ) {

  my $id=$seq->display_id();
  my $a=split_desc($seq->desc());
  my $organism=$a->{'organism'};

  if (exists($desc{$id})) {
    my $sequence=$seq->seq();
    unless ($sequence eq $s{$organism}->{$id}) {
      die "Same id '$id' but different sequence found,'$sequence' and '$s{$organism}->{$id}', please remove incorrect sequence\n";
    }

    warn "Same id '$id' and the same sequence, keeping first, ignoring others\n";
    next;
  }

  $desc{$id}=$seq->desc();

  # for uniprot do a swissprot or trembl distinction
  if ($a->{'source'} eq 'uniprot') {
    if ($id=~/^[QABO]\d/) {
        $a->{'source'} = 'trembl';
    } else {
        $a->{'source'} = 'swissprot';
    }
  }
  $sources{$id} = 'unknown' unless ($seq->desc() && exists($a->{'source'}) && exists($source_scores{$a->{'source'}}));

  $sources{$id} = $a->{'source'};

  $s{$organism}->{$id}=$seq->seq();
}

my %d; # {seqid1}->{seqid2}=distance
my %skip_id; #store ids which are thrown away

my %nr_per_organism; #{skipped}++ {kept}++
print STDERR "Calculating distances\n";

while ( my ($organism,$seqs)=each %s) {
  #calc distance of
  #all organisms without include
  #or only organisms which are included
  unless (%include_taxs && !(exists($include_taxs{$organism}))) {
      print STDERR "Current organism = $organism\n";
      foreach my $id1 (keys(%$seqs)) {
        #don't need to calc distance of $id1 against something else
        #if we already know $id1 should be skipped
        next if (exists($skip_id{$id1}));
        foreach my $id2 (keys(%$seqs)) {
          #don't need to calc distance of $id1 against $id2
          #if we already know $id2 should be skipped
          next if (exists($skip_id{$id2}));

          #only do below the diagonal
          next if ($id1 ge $id2);

          my $dist;
          #do not calculate levenstein dist when min_dist==0
          #only use 'eq' comparison
          if ($min_dist==0) {
            if ($seqs->{$id1} eq $seqs->{$id2}) {
              $dist=0;
            } else {
              $dist=$min_dist+1;
            }
          } else {
            # eq is faster then distance() if 1==2
            if ($seqs->{$id1} eq $seqs->{$id2}) {
              $dist=0;
            } else {
              $dist=distance(
                $seqs->{$id1},
                $seqs->{$id2}
              );
            }
          }

          #test dist
          if ($dist<=$min_dist) {
            #id1 should be kept

            my $perform_swap=0;

            #seqs of same source
            if ($sources{$id1} eq $sources{$id2}) {
                # keep seq which is shorter or alphabeticly first
                if ((length($id2) < length($id1)) || ($id1 ge $id2)) {
                    $perform_swap=1;
                }
            } else {
                # keep seq with lowest scoring source
                if (
                    $source_scores{$sources{$id1}} > $source_scores{$sources{$id2}}
                ) {
                    $perform_swap=2;
                } else {
                  if ($source_scores{$sources{$id1}} == $source_scores{$sources{$id2}} && length($id1) < length($id2)) {
                    # when source the same prefer shortest id
                    $perform_swap=2;
                  }
                }
            }

            if ($perform_swap) {
               my $swap=$id1;
              $id1=$id2;
              $id2=$swap;
              print STDERR "Swapping $id1 <=> $id2 : $perform_swap\n";
            }

            #store distance
            $d{$id1}->{$id2}=$dist;
            #store id2 which is highly similar to $id1
            $skip_id{$id2}=1;
            print STDERR "$id2 is highly similar ($dist aa diff) to $id1, skipping $id2 \n";
            if ($skip_file) {
              print $skip_fh $id2,"\t",$id1,"\t",$dist,"\t",$perform_swap,"\n";
            }
            $nr_per_organism{organism}->{skipped_dist}->{$dist}++;
            $nr_per_organism{organism}->{skipped}++;
          }
        }
      }
  }
  foreach my $id1 (keys(%$seqs)) {
    next if (exists($skip_id{$id1}));
#    print STDERR "\n\n\n===== Dumper $id1 $desc{$id1} $seqs->{$id1} ==\n\n" if length($seqs->{$id1})<100;
    print ">$id1 $desc{$id1}\n$seqs->{$id1}\n";
    $nr_per_organism{organism}->{kept}++;
  }
}

print STDERR Dumper \%nr_per_organism;
