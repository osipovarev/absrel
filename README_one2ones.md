# Preparing alignments of one2one orthologs

## 1. Extract alignemnts of one2one orthologs from TOGA

```
transcripts=$(cut -f2 /projects/hillerlab/genome/gbdb-HL/galGal6/TOGA/final.isoformes.galGal6.ncbi_appris.tsv)
refbed=/projects/hillerlab/genome/gbdb-HL/galGal6/TOGA/final.galGal6.ncbi_appris.anno.bed

for db in $(cat  nectar_birds.lst | grep -v galGal); do echo -e "/projects/hillerlab/genome/gbdb-HL/galGal6/TOGA/vs_${db}/"; done > toga_dirs.lst

for t in $transcripts; do echo -e "extract_codon_alignments_from_toga.py toga_dirs.lst $refbed $t --skip_dups > $t.fa"; done > jobs.extract_msa_toga

md batch_extract_toga
shuf jobs.extract_msa_toga | splitFile stdin 100 batch_extract_toga/batch_
ls batch_extract_toga/* | xargs -i echo "bash {}" >> batch.extract
para make batch.extract batch.extract
```
Most failed due to the time limit

```
for t in $transcripts; do if [ ! -f one2one_extracted_toga/$t.fa ]; then echo -e "extract_codon_alignments_from_toga.py toga_dirs.lst $refbed $t --skip_dups > one2one_extracted_toga/$t.fa"; fi; done > jobs2.extract_msa_toga

md batch_extract_toga
shuf jobs.extract_msa_toga | splitFile stdin 10 batch_extract_toga/batch_
ls batch_extract_toga/* | xargs -i echo "bash {}" >> batch.extract
para make batch.extract batch.extract
```

delete empty files
```
find one2one_extracted_toga/ -type f -name "*.fa" -size -1k | xargs -i rm -r {}
```
51245 files



## 2. Clean with HmmCleaner
```
for f in $(ls one2one_extracted_toga/); do echo -e "HmmCleaner.pl -costs -0.25 -0.20 0.15 0.45 one2one_extracted_toga/$f"; done > jobs.hmmcleaner

md batch_hmm
shuf jobs.hmmcleaner | splitFile stdin 30 batch_hmm/batch_
ls batch_hmm/* | xargs -i echo "bash {}" >> batch.hmm
para make batch.hmm batch.hmm
```

### 3. Get one2one orthlogs only
```
for db in $(cat ../nectar_birds.lst); do grep "one2one\|one2zero" $genomePath/gbdb-HL/galGal6/TOGA/vs_${db}/orthology_classification.tsv | cut -f1 | s -u > transcripts_count/$db.one2ones.lst; done
intersect_multiple_files.py  -f transcripts_count/
* > one2ones.inall.lst
```
11688 one2ones.inall.lst


## 4. Run manual filtering
```
hmm_dir=one2one_extracted_toga_hmm/
manual_dir=one2one_extracted_toga_hmm_manual
tanscripts=$(cat trans.one2ones.inall.lst)
for t in $tanscripts; do if [ -s $hmm_dir/${t}_hmm.fasta ]; then echo "manual_filter_msa.py -m  -mc 0.5 -ml 30 -a $hmm_dir/${t}_hmm.fasta > $manual_dir/$t.hmm.manual.fa"; fi; done  > jobs.manual
```

## 5. Check alignments
```
remove vs_ from names; replace REFERENCE with galGal6

for f in $(ls one2one_extracted_toga_hmm_manual/); do echo -e " sed -i 's/vs_//g' one2one_extracted_toga_hmm_manual/$f; sed -i 's/REFERENCE/galGal6/g' one2one_extracted_toga_hmm_manual/$f"; done > jobs.sed

clades="HLcalAnn5,HLfloFus1,HLphaSup1 HLtriMol2 HLlicMelCas1,HLphyNov1,HLlicPen1,HLgraPic1"

for clade in  $clades; do for t in $tanscripts; do if [ -s $manual_dir/$t.hmm.manual.fa ]; then echo "check_ali_file.py -f $manual_dir/$t.hmm.manual.fa -n 20 -r $clade -t $t"; fi; done; done > jobs.check.ali
```


## 6. Prepare directory per gene
```
gene_dir=by_gene_one2one_MSAs
md $gene_dir
for g in $(cat one2ones.inall.lst); do mkdir -p $gene_dir/$g; done

isoformes=/projects/hillerlab/genome/gbdb-HL/galGal6/TOGA/final.isoformes.galGal6.ncbi_appris.tsv
```

### [optional] Translate MSAs to prot
```
md prot_one2one_extracted_toga
for t in $(cat checked.transcripts.lst); do echo -e "translate.py -f one2one_extracted_toga/$t.fa | sed 's/REFERENCE/galGal6/g' | sed 's/vs_//g' > prot_one2one_extracted_toga/prot.$t.fa"; done > jobs.translate
```

next step:

##
