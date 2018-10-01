Run in g/ directory:

1. See gpcrdb_gapped_tm_numbering.csv

3 + 4
```bash
python ../scripts/gpcrdb_spreadsheet_filter.py -n gpcrdb_gapped_tm_numbering.csv gpcrdb_human_swissprot_tm_aligment.csv > gpcrdb_human_swissprot.fa
```
5
```bash
export BLASTDB=../ensembl
# swissprot
blastp -db uniprot_sprot -query gpcrdb_human_swissprot.fa -out gpcrlike_swissprot.idlist -evalue 10 -outfmt '6 sseqid' -num_threads 4 -num_alignments 1000
sort -u gpcrlike_swissprot.idlist > gpcrlike_swissprot.uniq.idlist
blastdbcmd -db uniprot_sprot -entry_batch gpcrlike_swissprot.uniq.idlist -outfmt '>%i taxid=%T;organism=%S;source=swissprot;source_version=2015_08\n%s' | grep -v taxid=0 | perl -pe 's/^>.*\|.*\|(.*?) />$1 /;s/\xFF/*/g;s/\\n/\n/' > gpcrlike_swissprot.fa

# trembl
blastp -db uniprot_trembl -query gpcrdb_human_swissprot.fa -out gpcrlike_trembl.idlist -evalue 10 -outfmt '6 sseqid' -num_threads 4 -num_alignments 1000
sort -u gpcrlike_trembl.idlist > gpcrlike_trembl.uniq.idlist
blastdbcmd -db uniprot_trembl -entry_batch gpcrlike_trembl.uniq.idlist -outfmt '>%i taxid=%T;organism=%S;source=trembl;source_version=2015_08\n%s' | grep -v taxid=0 | perl -pe 's/^>.*\|.*\|(.*?) />$1 /;s/\xFF/*/g;s/\\n/\n/' > gpcrlike_trembl.fa

# ensembl
blastp -db ensembl -query gpcrdb_human_swissprot.fa -out gpcrlike_ensembl.idlist -evalue 10 -outfmt '6 sseqid' -num_threads 4 -num_alignments 1000
sort -u gpcrlike_ensembl.idlist > gpcrlike_ensembl.uniq.idlist
# replace xFF with *
blastdbcmd -db ensembl -entry_batch gpcrlike_ensembl.uniq.idlist -outfmt '>%a taxid=%T;organism=%S;source=ensembl;source_version=81\n%s' | perl -pe 's/\xFF/*/g;s/\\n/\n/' > gpcrlike_ensembl.fa

cat gpcrlike_swissprot.fa gpcrlike_trembl.fa gpcrlike_ensembl.fa > source.fa
../scripts/highly_similar_reducer.pl -m 0 source.fa > uniq.fa
../scripts/seq_min_group.pl -t 100 uniq.fa > uniqt100.fa
../scripts/alignment.py -run config-gapped.dat uniqt100.fa
cp files/total_alignment.dat ali.fa
../scripts/highly_similar_reducer.pl -m 9 ali.fa > ali.s9.fa
../scripts/seq_min_group.pl -t 100 ali.s9.fa > ali.s9t100.fa
clustalo -i ali.s9t100.fa -o ali.s9t100.clustalo.fa --guidetree-out=ali.s9t100.ph -v
# correct clustalo -0 distances to positive 0.0, biopython chokes on them
perl -p -e 's/:-0/:0.0/g' ali.s9t100.ph > ali.s9t100.ph-0
../scripts/entropy.py -a ali.s9t100.fa -n gpcrdb_gapped_tm_numbering.csv -t ali.s9t100.ph-0 > ali.s9t100.entropies

# Check if all seed sequences are present in output
grep '>' gpcrdb_human_swissprot2.fa|sort > gpcrdb_human_swissprot2.ids
wc -l gpcrdb_human_swissprot2.ids
310
grep '>' ali.s9t100.fa | awk '{print $1}' |sort > ali.s9t100.ids
diff  gpcrdb_human_swissprot2.ids ali.s9t100.ids | grep '<' |wc -l
87
# So 87 missing in final alignment
grep '>' gpcrlike_swissprot.fa | awk '{print $1}' |sort > gpcrlike_swissprot.fa.ids
diff  gpcrdb_human_swissprot2.ids gpcrlike_swissprot.fa.ids | grep '<' |wc -l
0
# So no missing in swissprot blast
grep '>' uniqt100.fa | awk '{print $1}' |sort > uniqt100.fa.ids
diff  gpcrdb_human_swissprot2.ids uniqt100.fa.ids | grep '<' |wc -l
0
# So no missing in input for alignment
grep '>' ali.fa | awk '{print $1}' |sort > ali.fa.ids
diff  gpcrdb_human_swissprot2.ids ali.fa.ids | grep '<' |wc -l
85
# Alignment looses most human seed seqs

# hypotheses seqs with gaps inside tms where dropped because they have lower scores
```
