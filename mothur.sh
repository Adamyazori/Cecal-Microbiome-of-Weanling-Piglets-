
#make a directory "mothur"
mkdir mothur

#go into the directory

cd mothur

#Load mothur silva trained database

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138.tgz

#unzip- 
tar -zxvf silva.nr_v138.tgz

#copy your .fa file into mothur
cp /work/samodha/sadams23/pig_colonization_history_Study/Analysis/ASVs.fa /work/samodha/sadams23/pig_colonization_history_Study/mothur

#call mothur then you can start aligning (for me I call mothur as such: ~/mothur/mothur)
mothur

pcr.seqs(fasta=silva.nr_v138.align, start=11894, end=25319, keepdots=F, procesors=16) 

#next rename silva
system(mv silva.nr_v138.pcr.align silva.v4.fasta) #you can rename this how you want #should only need to do once

#align
align.seqs(fasta=ASVs.fa, reference=silva.v4.fasta) 

#exit mothur
quit()

#after aligning must make each ASV have a minimum of 10 characters to be recongized (do this outside of mothur... ie need a new terminal)
sed -i -e 's/>/>AAAAAAAAAA/g' ASVs.align
#also doesnt like ..... must take out (do this outside of mothur... ie need a new terminal)
sed -i -e 's/\./-/g' ASVs.align

#back into mothur
mothur

#create distances in mothur
dist.seqs(fasta=ASVs.align, processors=16, cutoff=.10, output=phylip)

#last step. finalize tree in mothur 
clearcut(phylip=ASVs.phylip.dist)
quit()
#change ASV back (do this outside of mothur... ie need a new terminal)
sed -i -e 's/AAAAAAAAAA//g' mothur/ASVs_mockpipeline.phylip.tre
#move finalized tree back to R working directory before proceeding
cp /work/samodha/sadams23/pig_colonization_history_Study/mothur/ASVs.phylip.tre /work/samodha/sadams23/pig_colonization_history_Study/Analysis/



#final phyloseq object
library(phyloseq)
tree= read_tree("Analysis/ASVs.phylip.tre")
library("readxl")
Multistate_pig_Study_meta= readxl::read_xlsx("Mapping_pig_colonization.xlsx")
Multistate_pig_Study_meta= as.data.frame(Multistate_pig_Study_meta)
rownames(Multistate_pig_Study_meta) =Multistate_pig_Study_meta$SampleID
summary(Multistate_pig_Study_meta)

ps= phyloseq(otu_table(asv_tab, taxa_are_rows=T),
                     sample_data(Multistate_pig_Study_meta),
		             tax_table(asv_tax)                
       	         	     phy_tree(tree))
save(ps, file= "ps.RData")
