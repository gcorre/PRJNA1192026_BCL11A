#### CORRE GUILLAUME @ GENETHON
#### Nov 11th, 2024
#### Analysis of ATF4 and GATA1 loci


# libraries 'template*' are fake reads generated to validate the pipeline.

  # template_BCL11A_WT is the WT sequence : NNNNNN site1 NNNN site2 NNNNNNN
  # template_BCL11A_DEL contains a deletion between the 2 editing sites : NNNNNN site1/site2 NNNNNNN
  # template_BCL11A_overDEL contains a deletion between the 2 editing sites and extra deletion  NNNNNN sit../..te2 NNNNNNN
  # template_BCL11A_INV contains an inversion between the 2 editing sites
  # template_BCL11A_microDEL contains small deletions at each loci :   NNNN sit.. NNNN ..te2 NNNN


#### Generate SV analysis

library(tidyverse)
library(ggseqlogo)


enrichment_start = 60490000

CBE_edit_5 = 60495266 - enrichment_start
CBE_edit_3=60498490 - enrichment_start
DBE_edit_5 = 60495265 - enrichment_start
DBE_edit_3 = 60498493 - enrichment_start

# nucleotide composition analysis using Seqlogo

compo <- list.files(path = "basecalling_hac3.3/",pattern = "composition.txt", recursive = T,full.names = T)
names(compo) <- basename(dirname(dirname(compo)))



compo_list <- lapply(compo, read.delim)


#### 5' editing
compo_logo <- lapply(compo_list, function(x) {
  x %>%
    filter(between(pos, 5255,5280)) %>% 
    select(pos,A,T,C,G) %>% 
    pivot_longer(cols = c(A:G),names_to = "base") %>% 
    pivot_wider(names_from = "pos", values_from = "value") %>% 
    column_to_rownames('base') %>% as.matrix
}
)

ggseqlogo(compo_logo[c(1:4)],facet = "wrap",ncol = 1, method = 'prob') +  
  annotate('rect', xmin = 10.5, xmax = 15.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow')+
  labs(title = "gRNA20 & gRNA40") + 
  scale_x_continuous(breaks = 1:26,labels = 5255:5280)

ggsave("basecalling_hac3.3/for_the_paper//GATA1_seqlogo.pdf", device = "pdf",width = 14,height = 6,units = "in")

GATA1 <- lapply(compo_logo[c(1:4)], function(x) {
  x %>% as.table %>% prop.table(margin = 2)  %>% as.data.frame %>%  pivot_wider(names_from = 'Var2',values_from = "Freq")}
  )


writexl::write_xlsx(x = GATA1,path = "basecalling_hac3.3/for_the_paper/GATA1_nucleotide_composition.xlsx",col_names = T)


#### 3' editing
compo_logo <- lapply(compo_list, function(x) {
  x %>%
  filter(between(pos, 8480,8505)) %>% 
  select(pos,A,T,C,G) %>% 
    pivot_longer(cols = c(A:G),names_to = "base") %>% 
    pivot_wider(names_from = "pos", values_from = "value") %>% 
    column_to_rownames('base') %>% as.matrix
  }
)


ggseqlogo(compo_logo[c(1:4)],facet = "wrap",ncol = 1, method = 'prob') +  
  annotate('rect', xmin = 11.5, xmax = 15.5, ymin = -0.05, ymax = 1.05, alpha = .1, col='black', fill='yellow')+ labs(title = "gRNA53 & gRNA66") + 
  scale_x_continuous(breaks = 1:26,labels = 8480:8505)

ggsave("basecalling_hac3.3/for_the_paper/ATF4_seqlogo.pdf", device = "pdf",width = 14,height = 6,units = "in")

ATF4 <- lapply(compo_logo[c(1:4)], function(x) {
  x %>% as.table %>% prop.table(margin = 2)  %>% as.data.frame %>%  pivot_wider(names_from = 'Var2',values_from = "Freq")}
)

writexl::write_xlsx(x = ATF4,path = "basecalling_hac3.3/for_the_paper/ATF4_nucleotide_composition.xlsx",col_names = T)






# INDELS analysis -------------------------------------------------------------------


small_indel = 4

cut_left = 5267
cut_right = 8493

metadata <- readxl::read_xlsx("M:\\DB_Unite de Genomique\\Alexandra\\ONT\\metadata_librairies.xlsx",sheet = 1)

# summary of libraries ---------

#l oad bed with locus aligned reads
bed <- list.files(path = "./basecalling_hac3.3/",pattern = ".locusAligned.sorted.end2end.annotated.bed", recursive = T ,full.names = T)

names(bed) <- basename(str_match(bed,pattern = "./(.+)/locus")[,2])


beds<- lapply(bed, function(x){
  if(file.size(x)>0){
    read.delim(x, header=F,comment.char = "#")
  } else {
    NULL
  }
})

beds_df <- bind_rows(beds,.id='library')

## number of unique reads aligned on the locus (with or without deletion)
beds_df_unique_all <- beds_df %>% 
  distinct(library,V4) 

rm(bed); rm(beds);



# Process reads aligned on the locus sequence ---------------------------------------


## bed files -------------------------------------------------------------------------

# calculate order of fragments for split reads using CIGAR string and orientation
beds_df <- beds_df %>% 
  mutate(left_trim = case_when(V6 == "+" ~ as.numeric(str_match(V7, "(^[0-9]+)S")[,2]),
                               TRUE ~ as.numeric(str_match(V7, "([0-9]+)S$")[,2])),
         overlap_feature = round(V12 / (V10-V9)*100,1)) %>%  
  replace_na(list(left_trim = 0)) %>%  
  arrange(library, V4,left_trim,V10) %>%
  group_by(library,V4) %>% 
  mutate(fragment = as.numeric(fct_inorder(V7)))


# get single read structure from all aligned fragments
structures_bed  <- beds_df %>% 
  mutate(fragment_strand = case_when(first(V6)=="+" ~ fragment,
                   TRUE ~ -fragment)) %>% 
  arrange(fragment_strand) %>% 
  group_by(library,V4) %>%
  summarise(structure = toString(V11),
            n_SplitFragments= n_distinct(fragment),
            strand = toString(V6),
            read_orientation = first(V6),   
            inversion = n_distinct(V6)>1,            # check if one fragment is in a different order
            duplication = str_count(structure,"deletion_part")>1,  # check if there are multiple copies of the fragment between the 2 editing sites
            deletion_between_gRNA = str_count(structure,"deletion_part")<1)




## VCF files -------------------------------------------------------------------------
## read VCF file to analyse SV

files <- list.files(path = "./basecalling_hac3.3", recursive = T, full.names = T, pattern = "locusAligned.sorted.end2end.sniffles.vcf")
files

files <- files[file.size(files)>0]

names(files) <- basename(str_match(files,pattern = "./(.+)/locus")[,2])

vcfs<- lapply(files, function(x){
  read.delim(x, header=F,comment.char = "#",colClasses = c("character","integer","integer",rep("character",7)))
})

vcfs_df <- bind_rows(vcfs,.id='library')

rm(files); rm(vcfs)

## extract metadatas and format vcf file
vcfs_df <- vcfs_df %>% 
  mutate(SVTYPE = str_match(V8,"SVTYPE=([A-Z/]+)")[,2]) %>%
  group_by(library,V1,SVTYPE) %>% 
  arrange(V2,.by_group = T) %>% 
  mutate(SVLEN = as.numeric(str_match(V8,"SVLEN=([-0-9]+)")[,2]),
         RE = as.numeric(str_match(V8,"RE=([-0-9]+)")[,2]),
         ID = row_number(),
         CHR2 = str_match(V8,"CHR2=([[:alnum:]-_]+)")[,2],
         RNAMES = str_match(V8,"RNAMES=([[:alnum:]-_,]+)")[,2]) %>% 
  filter(V7 != "UNRESOLVED")


## Get SV end position
vcfs_df2 <- vcfs_df %>%  mutate(SV_start = V2  ,
       SV_end = case_when(SVTYPE %in% c("DEL","INV","DEL/INV") ~ (V2+ abs(SVLEN)-1),
                          SVTYPE== "INS" ~ V2  + 10))

# anotate SV to editing site
vcfs_df2 <- vcfs_df2 %>%
  rowwise %>% 
  mutate(position = case_when(between(cut_left,SV_start-100,SV_end+100) & between(cut_right,SV_start-100,SV_end+100)~ "GATA1-ATF4",
                              between(cut_left,SV_start-20,SV_end+20) ~ "GATA1",
                              between(cut_right,SV_start-20,SV_end+20) ~ "ATF4",
                              TRUE ~NA))

# calculat erelative position to cut site according to locus annotation
vcfs_df2 <- vcfs_df2 %>% 
  mutate(SV_start_relative = case_when(position == "GATA1" ~  SV_start-cut_left,
                                       position == "ATF4" ~  SV_start-cut_right,
                                       position == "GATA1-ATF4" ~ SV_start-cut_left),
         SV_end_relative = case_when(position == "GATA1" ~  SV_end-cut_left,
                                     position == "ATF4" ~  SV_end-cut_right,
                                     position == "GATA1-ATF4" ~ SV_end-cut_right))

vcfs_df2 <- vcfs_df2 %>% filter(!is.na(position))


# classify SV by size categories

vcfs_df2 <- vcfs_df2 %>% 
  mutate(SV_size = case_when(abs(SVLEN)<=small_indel~ "Unmodified",
                             between(abs(SVLEN),small_indel+1,50)~ "small (<=50bp)",
                             between(abs(SVLEN),51,200)~ "intermediate(<=200bp)",
                             abs(SVLEN)>200~ "Large (>200bp)")) %>% 
  mutate(SVTYPE_size = case_when(SV_size == "Unmodified"~"Unmodified",
                                 TRUE ~ SVTYPE)) %>% 
  separate_rows(RNAMES, sep=",")


# generate missing locus/library combinations
all_combinations <- vcfs_df2 %>%
  filter(position %in% c("GATA1","ATF4")) %>% 
  distinct(library,RNAMES) %>% 
  expand(nesting(library,RNAMES),position=c("GATA1","ATF4")) 


# get all reads in bed file
bed_combinations <- beds_df_unique_all %>% 
  anti_join(vcfs_df2 %>% 
              distinct(library,RNAMES), by = c("library","V4"="RNAMES")) %>% 
  expand(nesting(library,RNAMES=V4),position=c("GATA1","ATF4"))


# generate complete table with all analyzed reads
vcfs_df2 <- vcfs_df2 %>% 
  full_join(all_combinations) %>%
  full_join(bed_combinations) %>% 
  replace_na(list(SVTYPE_size = "Unmodified",
                  SV_size = "Unmodified",
                  SV_start = 0, SV_end = 0,
                  SV_start_relative=0,
                  SV_end_relative = 0))



## keep largest SV at each position
vcfs_df2 <- vcfs_df2 %>% 
  group_by(library,RNAMES,position) %>% 
  slice_max(abs(SVLEN),n = 1) %>% 
  mutate(position = factor(position, levels = c("GATA1", "ATF4","GATA1-ATF4")))

lib_counts <- vcfs_df2 %>% ungroup %>% distinct(library,RNAMES) %>% count(library,name = 'total_reads')


# Make some plots and tables

data_plot <- vcfs_df2 %>% ungroup %>% count(library,position,SVTYPE_size) %>% left_join(lib_counts) %>% mutate(proportion = round(n/total_reads,digits = 2))

ggplot(data_plot, aes(position,proportion,fill=SVTYPE_size)) + 
  geom_col(col = "black")+
  facet_wrap(~library,scales="free_y",ncol=4)+
  ggprism::theme_prism()+
  labs(x= "",y = "Percent of reads")+scale_y_continuous(limits = c(0,1),labels = scales::percent_format())



ggsave("basecalling_hac3.3/for_the_paper/SVtype_per_site.pdf", device = "pdf",width = 20,height = 8,units = "in")

write.table(data_plot, "basecalling_hac3.3/for_the_paper/SVtype_per_site.csv", sep=";", row.names = F, quote = F)


data_plot <- vcfs_df2 %>% ungroup %>% count(library,position,SV_size) %>% left_join(lib_counts) %>% mutate(proportion = round(n/total_reads,digits = 2))

ggplot(data_plot, aes(position,proportion,fill=SV_size)) + 
  geom_col(col = "black")+
  facet_wrap(~library,scales="free_y",ncol = 4)+
  ggprism::theme_prism()+
  labs(x= "",y = "Percent of reads")+scale_y_continuous(limits = c(0,1),labels = scales::percent_format())


ggsave("basecalling_hac3.3/for_the_paper/SVsize_per_site.pdf", device = "pdf",width = 20,height = 8,units = "in")
write.table(data_plot, "basecalling_hac3.3/for_the_paper/SVsize_per_site.csv", sep=";", row.names = F, quote = F)



# reads with deletion

table_data <- vcfs_df2 %>% ungroup %>%
  filter(position %in% c("GATA1","GATA1-ATF4")) %>% 
  count(library,
        locus=position, 
        SV_type=SVTYPE_size, 
        SV_size,name = "Read_Count") %>% left_join(lib_counts)

write.table(table_data, "basecalling_hac3.3/for_the_paper/table_for_graph.csv", sep=";", row.names = F, quote = F)



ggplot(vcfs_df2 %>% 
         #filter(str_detect(SVTYPE_size, "DEL|INV")) %>% 
         filter(SVTYPE_size!="Unmodified") %>% 
         group_by(library) %>%  
         arrange(desc(SV_start_relative-SV_end_relative),.by_group = T) %>% 
         mutate(rank = row_number()),
       aes(x=SV_start_relative,xend=SV_end_relative,y=rank,yend=rank, col = position)) +
  geom_vline(xintercept = 0, lty =2, col = "grey")+
  geom_segment(lwd =1)+
  facet_wrap(~library, scales = "free", ncol = 3,as.table = F,dir = "v")+
  ggprism::theme_prism()+
  labs("Relative position to cutting site", y = "Reads", title = "Deletion size around cutting sites")

ggsave("basecalling_hac3.3/for_the_paper/SVsize_size_around_cut_site.pdf", device = "pdf",width = 20,height = 12,units = "in")



ggplot(vcfs_df2 %>% 
         #filter(str_detect(SVTYPE_size, "DEL|INV")) %>% 
         filter(SVTYPE_size!="Unmodified") %>% 
         group_by(library) %>%  
         arrange(desc(abs(SVLEN)),.by_group = T) %>% 
         mutate(rank = row_number()),
       aes(x=SV_start,xend=SV_end,y=rank,yend=rank, col = SVTYPE_size)) +
  geom_segment()+
  geom_vline(xintercept = c(cut_left,cut_right), lty =2, col = "grey")+
  facet_wrap(~library, scales = "free_y", ncol = 3, dir = "v",as.table = F)+
  ggprism::theme_prism()+
  labs("Relative position to cutting site", y = "Reads", title = "Deletion size ")+
  scale_x_continuous(limits = c(5000,10000))

ggsave("basecalling_hac3.3/for_the_paper/DeletionSize_at_cut_sites.pdf", device = "pdf",width = 20,height = 12,units = "in")


write.table(vcfs_df2, "basecalling_hac3.3/for_the_paper/complete_SV_table.csv", sep = ";", quote = F, row.names = F)

###### END #######