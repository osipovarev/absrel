# Running positive selection screen with aBSREL

### for annotation-based run: 
```
gene_file=/projects/project-osipova/NectarivoryProject/absrel/isoforms.ncbi.galGal6.csv
msa_final_dir=/projects/project-osipova/NectarivoryProject/make_one2one_MSA_iso_ncbi/by_gene_orthologs_MSAs_final/
absrel_dir=/projects/project-osipova/NectarivoryProject/absrel/run_absrel_by_gene/
tree=/projects/project-osipova/NectarivoryProject/bird_phylogeny/30_bird_tree.phyloFit.Maggie.mod
qual_trans=/projects/project-osipova/NectarivoryProject/absrel/results_absrel_annotations/quality.transcripts.lst
results_dir=/projects/project-osipova/NectarivoryProject/absrel/results_absrel_annotations/
suff=".30birds.macse.hmm.manual.fa"
```

### for TOGA projections-absed run:
```
md run_absrel_toga_by_gene/
gene_file=/projects/project-osipova/NectarivoryProject/absrel/isoforms.ncbi.galGal6.csv
msa_final_dir=/projects/project-osipova/NectarivoryProject/make_one2one_MSA_iso_ncbi/by_gene_toga_orthologs_MSAs_final/
absrel_dir=/projects/project-osipova/NectarivoryProject/absrel/run_absrel_toga_by_gene/
tree=/projects/project-osipova/NectarivoryProject/bird_phylogeny/30_bird_tree.phyloFit.Maggie.mod
results_dir=/projects/project-osipova/NectarivoryProject/absrel/results_absrel_toga/
qual_trans=/projects/project-osipova/NectarivoryProject/absrel/results_absrel_toga/toga.quality.transcripts.lst
suff=".30birds.macse.hmm.manual.fa"
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

## 2. Write absrel jobs
```
write_absrel_jobs.sh $gene_file $msa_final_dir $absrel_dir $qual_trans $suff > jobs.absrel

md batch_absrel
shuf jobs.absrel | splitFile stdin 30 batch_absrel/batch_
ls batch_absrel/* | xargs -i echo "bash {}" > batch.absrel
para make batch.absrel batch.absrel
# 20 hours
```
1738 absrel jobs failed due to time limit -> resubmitted.
8 hours



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

## 5. Next: see absrel_stats.ipynb)




***

## Label trees with p-values
```
for g in $(awk '{if (NF == 4) {print $2}}' genes.to_analse.summary.all_clades.tsv); do t=$(grep $g results_absrel_annotations/lowest.lowest.all.genes.pval.table.tsv | head -1 | cut -f3); grep -w $t results_absrel_annotations/lowest.lowest.all.genes.pval.table.tsv > Label_trees/$g.anno.lowest_pval.tsv; tree_add_attribute.py -t run_absrel_by_gene/$g/$t/$t.ans.tree -l Label_trees/$g.anno.lowest_pval.tsv > Label_trees/labeled.$g.anno.tree; done
```

## Plot p-value trees for each gene
```
for g in $(awk '{if (NF == 4) {print $2}}' genes.to_analse.summary.all_clades.tsv); do Rscript plot_pvalue_tree.R Label_trees/labeled.$g.anno.tree Label_trees/tree.$g.anno.lowest_pval.pdf; done
```


