# Enrichment analysis with ClusterProfiler

## 1. Get genes under selection per nectar clade
```
TARGET=nectar
clades=$(echo "hmmbrds honeyeaters nectar_parrots sunbirds")

TARGET=nonnectar
clades=$(echo "swifts falcons lyrebirds passerides")

for i in $clades; \
do \
 grep ^$i under_selection_per_clade_0.01.$TARGET.tsv | cut -f2 > $i.under_selection_per_clade_0.01.txt; \
done
```



## 2. Run enrichment GO analysis 
### NB: Replace gene names with human gene symbols for the enrichcment analysis
```
clades=$(echo hmmbrds honeyeaters nectar_parrots sunbirds falcons lyrebirds swifts passerides)
RENAMEDICT=~/Documents/LabDocs/Chicken_resources/galGal6_gene.hg38_gene_symbol.tsv

md ClusterProfiler/

for i in $clades; \
do \
 echo $i; \

 renameToHLscaffolds.py -c 1 -a $i.under_selection_per_clade_${P}.txt -d <(sed 's/\t/,/' $RENAMEDICT) > hg38.$i.under_selection_per_clade_${P}.txt; \ 
 
 goenrich_genelist.R -w $(pwd) -g hg38.$i.under_selection_per_clade_${P}.txt -o ClusterProfiler/hg38.goenrich.$i.under_selection_per_clade_${P}.tsv; \
done
```



## 3. Run enrichment GO analysis for genes of different ranks
```
TARGET=nectar
TARGET=nonnectar
#TARGET=nectar_norelax
#TARGET=nonnectar_norelax
```

### 3.1. Make gene lists
```
P=0.01

for TARGET in nectar nonnectar; \
do \
	for i in {1..4}; do echo $i; cut -f2 under_selection_per_clade_${P}.$TARGET.tsv | sort | uniq -c | awk -v var=$i '$1>=var{print $2}' > rank${i}.$TARGET.genes_selection.lst; done; \
done
```

### 3.2. Run enrichGO by rank
```
RENAMEDICT=~/Documents/LabDocs/Chicken_resources/galGal6_gene.hg38_gene_symbol.tsv

for f in $(ls rank*genes_selection.lst); \
do \
	echo $f; \
	renameToHLscaffolds.py -c 1 -a $f -d <(sed 's/\t/,/' $RENAMEDICT)  > hg38.$f; \
done

for t in nectar nonnectar; \
do \
	for i in  {1..4}; \
	do \
		echo $i; \
		goenrich_genelist.R -w $(pwd) -g hg38.rank${i}.$t.genes_selection.lst -o ClusterProfiler/hg38.goenrich.rank${i}.$t.tsv; \
	done; \
done
```



## 4. Test convergence: get 3-4 way convergent terms
```
all_nectaf="hg38.goenrich.hmmbrds.under_selection_per_clade_${P}.tsv hg38.goenrich.honeyeaters.under_selection_per_clade_${P}.tsv  hg38.goenrich.nectar_parrots.under_selection_per_clade_${P}.tsv hg38.goenrich.sunbirds.under_selection_per_clade_${P}.tsv"


for g in $(cat $all_nectaf| cut -f1 | grep -v ^ID | s | uniq -c | awk '$1>2{print $2}'); \
do \
	grep $g $all_nectaf | cut -f1,2,8; \
done | s -u > 3_4_way_convergent_terms.nectar.tsv
```


### 4.1. Exlcude children GO terms; exlucde terms with gene count < 3
```
GOOBO=~/Documents/LabDocs/GO_terms_genes/go.obo
f=hg38.goenrich.rank2.nonnectar.tsv
c=1

golist=$(cut -f${c} $f | tail -n +2|  tr '\n' ',')

for g in $(echo $golist | tr ',' '\n'); do echo $g; get_go_children.py -f $GOOBO  -go $g -l $golist ; done| grep "has parents" | awk '{print $1}' > to_exclude_go.lst

filter_annotation_with_list.py -b -c $c -l to_exclude_go.lst -a $f |  awk -F"\t" '$9>=3{print}'> noChildren.gt2genes.$f 
```


## 5. Find representative (ancestral) GO terms
```
for i in {1..4}; \
do \
	for c in nectar nonnectar; \
	do \
		echo $i; \
		find_representative_go.R -w $(pwd) -e hg38.goenrich.rank${i}.$c.tsv -o represent_GO.hg38.goenrich.rank${i}.$c.tsv; \
	done; \
done
```


### 5.1. Make breakdown of selected terms in representative GO terms
```
for i in {1..4}; \
do \
	for g in $(grep "homeostas\|blood\|cardi\|system process\|carbo\|macromolec" represent_GO.rank${i}.nectar.tsv  | cut -f1 | cut -d: -f1,2); \
	do  \
		children=$(grep "^$g\t" represent_GO.rank${i}.nectar.tsv | cut -f5 | sed 's/,/ /g'); \
		for c in $children; \
		do \
			printf "$g\t"; \
			grep $c gene_hg38.goenrich.rank${i}.nectar.tsv; \
		done; \
	done; \
done > breakdown.interesting_rank_terms.nectar.tsv


for g in $(echo $gos | tr ',' ' '); do g $g noChildren.hg38.goenrich.rank2.nectar.tsv; done | cut -f1,2
```



## 6. Test convergence with Fisher-exact test
```
P=0.01

for i in 2 3 4; \
do \
	paste <(for p in $(grep ^nectar${i}way pairs_test_control.tsv | cut -f2); do shared=$(grep "$p" under_selection_per_clade_${P}.*nectar.tsv | cut -f2 | s | uniq -c | awk -v var=$i '$1==var{print}' | wc -l| awk '{print $1}'); universe=$(grep "$p" under_selection_per_clade_${P}.*nectar.tsv | cut -f2 | s | uniq -c | wc -l| awk '{print $1}'); echo -e "$p\t$shared\t$universe"; done) <(for p in $(grep ^nonnectar${i}way pairs_test_control.tsv | cut -f2); do shared=$(grep "$p" under_selection_per_clade_${P}.*nectar.tsv | cut -f2 | s | uniq -c | awk -v var=$i '$1==var{print}' | wc -l| awk '{print $1}'); universe=$(grep "$p" under_selection_per_clade_${P}.*nectar.tsv | cut -f2 | s | uniq -c | wc -l| awk '{print $1}'); echo -e "$p\t$shared\t$universe"; done); \

done | awk '{OFS="\t"; print $1,$2,$3-$2,$4,$5,$6-$5}' > matching_test_control.tsv


paste <(cat header.matching_test_control.tsv | tail +2) <(for line in $(cat matching_test_control.tsv | awk '{print $2","$3","$5","$6}'); do simple_fisher_exact.R $line; done | cut -d: -f2 | sed 's/^ //' | sed 's/ $//' | tr ' ' '\n') > file;

mv file summary.matching_test_control.tsv
```




## 7. Reviews: Running RELAX

### 7.1. get transcripts under selection in nectar or nonnectar
```
 filter_annotation_with_list.py -c 4 -l <(cut -f1 under_selection_ranked_genes_0.5_n*) -a all.genes.pval.table.tsv |  awk '$2<${P}{print}' | cut -f3 | sort -u > transcripts_for_relax.lst
```

## 8. Check semantic similarity between terms in GO enrichments of the default and SRV models
```
## rank2
OUT1=represent_GO.hg38.goenrich.rank2.nectar.tsv
OUT2=../../combined_default_srv_mh_v1/ClusterProfiler/represent_GO.rank2.nectar.tsv

make_all_against_all_table.py -l1 $(cut -f1 $OUT1 | tail -n +2 | tr '\n' ',')  -l2 $(cut -f1 $OUT2 | tail -n +2 |tr '\n' ',') > all_against_all.rank2_versions.txt

awk '$4>0.2{print}' similarity.all_against_all.rank2_versions.txt

## 3-4 way convergence
OUT1=3_4_way_convergent_terms.nectar.tsv
OUT2=../../combined_default_srv_mh_v1/ClusterProfiler/3_4_way_convergent_terms.nectar.tsv

make_all_against_all_table.py -l1 $(cut -f2 $OUT1 | tail -n +2 | tr '\n' ',')  -l2 $(cut -f1 $OUT2 | tail -n +2 |tr '\n' ',') > all_against_all.3_4_way_versions.txt

perl go_comb.pl all_against_all.3_4_way_versions.txt similarity.all_against_all.3_4_way_versions.txt
``` 



