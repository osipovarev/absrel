## Get the p-val trees
## Get the MEME anno csv files
```
meme_dir=run_absrel_toga_by_gene

for t in $(cut -f2 selected_iso.gene_list.meme.tsv); do g=$(grep $t selected_iso.gene_list.meme.tsv | cut -f1);  scp eosipova@10.10.50.202:/projects/project-eosipova/NectarivoryProject/absrel_2021/$meme_dir/$g/$t/site.pval.jalview.$t.fg.csv $g/; scp eosipova@10.10.50.202:/projects/project-eosipova/NectarivoryProject/absrel_2021/$meme_dir/$g/$t/pval.$t.ans.tree $g; done
```

## Plot trees with p-values
```
for t in $(cut -f2 selected_iso.gene_list.meme.tsv); do g=$(grep $t selected_iso.gene_list.meme.tsv | cut -f1); plot_pvalue_tree.R $g/pval.$t.ans.tree  $g/pval.$t.tree.pdf; done
```


### Make Figures with trees and alignemtns
```
# tree=~/Documents/LabDocs/Birds_Phylogeny/45_bird_code_tree_topology.nh
tree=pval.HK3_ENSGALT00000046805.2.ans.tree
bird_tab=~/Documents/LabDocs/Birds_Phylogeny/45_birds.master_table.tsv
prot=prot.HK3_ENSGALT00000046805.2.hmm.manual.fa

rename_fasta_from_dict.py -f $prot -d <(awk -F"\t" '{print $1","$5}' $bird_tab) | sed 's/ /_/g' > common.$prot
```

## rename to common name in trees
```
for t in $(cut -f2 selected_iso.gene_list.meme.tsv); do g=$(grep $t selected_iso.gene_list.meme.tsv | cut -f1); substituteNamesTree.py -t $g/rev.pval.$t.ans.tree -n <(cut -f1,4,5 $bird_tab) -c -o 4 | sed 's/ /_/g' > $g/common.rev.pval.$t.ans.tree; done
```

## plot reverse p-value trees with common names
```
for t in $(cut -f2 selected_iso.gene_list.meme.tsv); do g=$(grep $t selected_iso.gene_list.meme.tsv | cut -f1); plot_pvalue_tree.R $g/common.rev.pval.$t.ans.tree $g/common.rev.pval.$t.ans.tree.pdf; done
```

## Sort alignments
```
awk -F"\t" '{print $1","$5}' /Users/osipova/Documents/LabDocs/Birds_Phylogeny/45_birds.master_table.tsv > /Users/osipova/Documents/LabDocs/Birds_Phylogeny/45_birds.dict.csv

dict=/Users/osipova/Documents/LabDocs/Birds_Phylogeny/45_birds.dict.csv
order=/Users/osipova/Documents/LabDocs/Birds_Phylogeny/45_birds.phylo_ordered.lst

for t in $(cut -f2 selected_iso.gene_list.meme.tsv);
do
 g=$(grep $t selected_iso.gene_list.meme.tsv | cut -f1);
 rename_fasta_from_dict.py -f $g/prot.$t.hmm.manual.fa -d $dict > $g/commonNames.prot.$t.fa;
 sort_fasta_by_list.py -f $g/commonNames.prot.$t.fa -l $order -r chicken > $g/ordered.commonNames.prot.$t.fa;
done
```

