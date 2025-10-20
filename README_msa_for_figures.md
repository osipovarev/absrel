#  Describes how to prepare gene MSA for Figures

```
gene=MLXIPL
dic=mammals_birds.dict.csv
order=mammals_birds.phylo_ordered.lst


nectar=HLlepAsp1,HLcinPul1,HLlicMelCas1,HLlicPen1,HLphyNov1,HLfloFus1,HLphaSup1,HLcalAnn5,HLparPun1,HLdicExi1,HLamaAes1,HLaraSol1,HLtriMol2,HLlorGal1

gene=HK3
gene=SLC2A5
```

## Get prot alignments of full gene and BEB.nogaped
```
scp eosipova@10.10.50.202:/projects/project-eosipova/NectarivoryProject/$gene/prot* .
```

## Combine these files together and align
```
#cat prot.nectar.*nogaps.fa <(get_seq_by_name_fasta.py -ab -n $(echo $nectar) -f prot.*.$gene.fa) > combined.prot.$gene.fa
```
## better: add only one nogaps nectar seq as reference and align 
```
mafft combined.prot.$gene.fa > mafft.combined.prot.$gene.fa
```

## Replace assembly codes with common names and sort MSA by phylogeny
```
dict=/Users/osipova/Documents/LabDocs/Birds_Phylogeny/45_birds.dict.csv
order=/Users/osipova/Documents/LabDocs/Birds_Phylogeny/45_birds.phylo_ordered.lst

rename_fasta_from_dict.py -f mafft.combined.prot.$gene.fa -d $dict > commonNames.combined.prot.$gene.fa

sort_fasta_by_list.py -f commonNames.combined.prot.$gene.fa -l $order -r chicken > ordered.commonNames.combined.prot.$gene.fa
```