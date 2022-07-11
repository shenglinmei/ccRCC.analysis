

args <- commandArgs(trailingOnly = T);
fname=args[1]

fano=paste('inferCNA.',fname,'.ano.txt',sep='')
fexp=paste('inferCNA.',fname,'.exp.txt',sep='')

library("infercnv")

ano=read.csv(fano,sep='\t',header=F)

tmp=as.character(unique(ano[,2]))

# Normal reference
control=tmp[grepl('Proximal|podocytes',tmp)]
print(control)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=fexp,
                                    annotations_file=fano,
                                    delim="\t",
                                    gene_order_file="hg19.sort.inferCNV.txt",   #,'Epitheial_Hillock'
                                    ref_group_names=control)





# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=fname, 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             use_zscores = TRUE,
			     mask_nonDE_pval = 1,
                             mask_nonDE_genes = F,
                             include.spike=F  # used for final scaling to fit range (0,2) centered at 1.
)



