# Running positive selection screen with aBSREL


## 0. Set up for TOGA projections-based run:
```
gene_file=validated.final.isoformes.galGal6.ncbi_appris.csv
isoformes=validated.final.isoformes.galGal6.ncbi_appris.tsv
tree=/projects/hillerlab/genome/gbdb-HL/galGal6/maf/multiz45way_2021/45_bird_tree.philofit_4d.mod
suff=".hmm.manual.fa"
msa_final_dir=/projects/project-eosipova/NectarivoryProject/one2one_MSA_2024/one2one_MSAs/
absrel_dir=/projects/project-eosipova/NectarivoryProject/absrel_2024/run_absrel_toga_by_gene/
results_dir=/projects/project-eosipova/NectarivoryProject/absrel_2024/results_absrel_toga/
qual_trans=/projects/project-eosipova/NectarivoryProject/one2one_MSA_2024/good.transcripts.lst
```

## 1. Write tree_doctor jobs; prepare output dir structure
```
write_tree_doctor_jobs.sh $gene_file $msa_final_dir $absrel_dir $tree $qual_trans $suff > jobs.tree_doctor

# split jobs
md batch_tree_doc
shuf jobs.tree_doctor | splitFile stdin 30 batch_tree_doc/batch_
ls batch_tree_doc/* | xargs -i echo "bash {}" > batch.tree.doc
para make batch.tree.doc batch.tree.doc
# 1 min
```

## 2. Write aBSREL jobs
```
write_absrel_jobs.sh $gene_file $msa_final_dir $absrel_dir $qual_trans $suff > jobs.absrel

md batch_absrel
shuf jobs.absrel | splitFile stdin 30 batch_absrel/batch_
ls batch_absrel/* | xargs -i echo "bash {}" > batch.absrel
para make batch.absrel batch.absrel
# 20 hours
```

### Run aBSREL for four differnet models independently:
1) -defalut: 	with just "ENV='TOLERATE_NUMERICAL_ERRORS=1;"
2) -SRV: 		"--srv Yes"
3) -MH:  		"--multiple-hits Double+Triple"
3) -SRV+MH:  	"--multiple-hits Double+Triple --srv Yes"



## 3. Parse aBSREL results
```
write_parse_absrel_jobs.sh $gene_file $msa_final_dir $absrel_dir $qual_trans $suff > jobs.parse.absrel

md batch_parse
shuf jobs.parse.absrel | splitFile stdin 30 batch_parse/batch_
ls batch_parse/* | xargs -i echo "bash {}" > batch.parse
para make batch.parse batch.parse
# 1 min
```

## 4. Combine all p-values in one table
```
write_pval_table.sh isoforms.ncbi.galGal6.tsv $absrel_dir > $results_dir/all.genes.pval.table.tsv
```


## 5. Find the best model fit for each transcript

### 5.1. Make table AIC-c v58 default VS srv
```
for t in $(cat ../one2one_MSA_2024/good.transcripts.lst); do \
    g=$(grep $t validated.final.isoformes.galGal6.ncbi_appris.tsv | cut -f1); \
    aic=$(grep "AIC-c" run_absrel_toga_by_gene/$g/$t/out_v58.$t.json | tail -1 | awk -F"[:,]" '{print $2}'); \
    aic_srv=$(grep "AIC-c" run_absrel_toga_by_gene/$g/$t/out_v58_srv.$t.json | tail -1 | awk -F"[:,]" '{print $2}'); \
    echo -e "$g\t$t\t$aic\t$aic_srv"; \
done > AICc.v58_VS_v58_srv.tsv
```

### 5.2. Get AIC-c for v58 srv mh
```
for t in $(cat ../one2one_MSA_2024/good.transcripts.lst); do \
    g=$(grep $t validated.final.isoformes.galGal6.ncbi_appris.tsv | cut -f1); \
    aic_srv_mh=$(grep "AIC-c" run_absrel_toga_by_gene/$g/$t/out_v58_srv_mh.$t.json | tail -1 | awk -F"[:,]" '{print $2}'); \
    echo -e "$t\t$aic_srv_mh"; \
done > AICc.v58_srv_mh.tsv
```

### 5.3. Find the best model for each transcript
```
find_best_values.py -c 2 -i AICc.v58_VS_v58_srv_VS_v58_srv_mh.tsv > best_model.AIC.tsv
awk '$3==1{print $1}' best_model.AIC.tsv > srv_better.transcripts.lst	# 262
awk '$3==2{print $1}' best_model.AIC.tsv > mh_better.transcripts.lst	# 87
awk '$3==3{print $1}' best_model.AIC.tsv > srv_mh_better.transcripts.lst 	# 38
# default = 31321
```

### 5.4. Make combined table
```
filter_annotation_with_list.py -c 3 -a all.genes.pval.table.v58_srv_mh.tsv -l srv_mh_better.transcripts.lst > srv_mh_better.genes.pval.table.v58_srv_mh.tsv
filter_annotation_with_list.py -c 3 -a all.genes.pval.table.v58_tolerate_mh.tsv -l mh_better.transcripts.lst > mh_better.genes.pval.table.v58_tolerate_mh.tsv
filter_annotation_with_list.py -c 3 -a
all.genes.pval.table.v58_srv_mh.tsv -l srv_mh_better.transcripts.lst > srv_mh_better.genes.pval.table.v58_srv_mh.tsv
filter_annotation_with_list.py -b -c 3 -a all.genes.pval.table.v58_tolerate.tsv -l <(cat srv_better.transcripts.lst srv_mh_better.transcripts.lst mh_better.transcripts.lst) > default_better.genes.pval.table.v58_tolerate.tsv

cat default_better.genes.pval.table.v58_tolerate.tsv srv_better.genes.pval.table.v58_srv.tsv mh_better.genes.pval.table.v58_tolerate_mh.tsv srv_mh_better.genes.pval.table.v58_srv_mh.tsv > combined.all.genes.pval.table.v58.tsv
```


## 6. Next: see [absrel_stats.ipynb](https://github.com/osipovarev/absrel/blob/main/absrel_analysis.ipynb)


***

## Label trees with p-values
```
for g in $(awk '{if (NF == 4) {print $2}}' genes.to_analse.summary.all_clades.tsv); do t=$(grep $g results_absrel_annotations/lowest.lowest.all.genes.pval.table.tsv | head -1 | cut -f3); grep -w $t results_absrel_annotations/lowest.lowest.all.genes.pval.table.tsv > Label_trees/$g.anno.lowest_pval.tsv; tree_add_attribute.py -t run_absrel_by_gene/$g/$t/$t.ans.tree -l Label_trees/$g.anno.lowest_pval.tsv > Label_trees/labeled.$g.anno.tree; done
```

## Plot p-value trees for each gene
```
for g in $(awk '{if (NF == 4) {print $2}}' genes.to_analse.summary.all_clades.tsv); do Rscript plot_pvalue_tree.R Label_trees/labeled.$g.anno.tree Label_trees/tree.$g.anno.lowest_pval.pdf; done
```


