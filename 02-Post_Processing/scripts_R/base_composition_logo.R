library(tidyverse)
library(ggseqlogo)

enrichment_start = 60490000

CBE_edit_5 = 60495266 - enrichment_start
CBE_edit_3=60498490 - enrichment_start
DBE_edit_5 = 60495265 - enrichment_start
DBE_edit_3 = 60498493 - enrichment_start

# nucleotide composition

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


# BCL11A locus --------------------------------------------------------------------------
library(tidyverse)
small_indel = 4

cut_left = 5267
cut_right = 8493

metadata <- readxl::read_xlsx("M:\\DB_Unite de Genomique\\Alexandra\\ONT\\metadata_librairies.xlsx",sheet = 1)

# summary of libraries ---------

#load bed with locus aligned reads
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

beds_df <- beds_df %>% 
  mutate(left_trim = case_when(V6 == "+" ~ as.numeric(str_match(V7, "(^[0-9]+)S")[,2]),
                               TRUE ~ as.numeric(str_match(V7, "([0-9]+)S$")[,2])),
         overlap_feature = round(V12 / (V10-V9)*100,1)) %>%  
  replace_na(list(left_trim = 0)) %>%  
  arrange(library, V4,left_trim,V10) %>%
  group_by(library,V4) %>% 
  mutate(fragment = as.numeric(fct_inorder(V7)))

structures_bed  <- beds_df %>% 
  mutate(fragment_strand = case_when(first(V6)=="+" ~ fragment,
                   TRUE ~ -fragment)) %>% 
  arrange(fragment_strand) %>% 
  group_by(library,V4) %>%
  summarise(structure = toString(V11),
            n_SplitFragments= n_distinct(fragment),
            strand = toString(V6),
            read_orientation = first(V6),
            inversion = n_distinct(V6)>1,
            duplication = str_count(structure,"deletion_part")>1,
            deletion_between_gRNA = str_count(structure,"deletion_part")<1)







## VCF files -------------------------------------------------------------------------


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





vcfs_df2 <- vcfs_df %>%  mutate(SV_start = V2  ,
       SV_end = case_when(SVTYPE %in% c("DEL","INV","DEL/INV") ~ (V2+ abs(SVLEN)-1),
                          SVTYPE== "INS" ~ V2  + 10))

vcfs_df2 <- vcfs_df2 %>% rowwise %>% mutate(position = case_when(between(cut_left,SV_start-100,SV_end+100) & between(cut_right,SV_start-100,SV_end+100)~ "GATA1-ATF4",
                                         between(cut_left,SV_start-20,SV_end+20) ~ "GATA1",
                                         between(cut_right,SV_start-20,SV_end+20) ~ "ATF4",
                    TRUE ~NA))

vcfs_df2 <- vcfs_df2 %>% mutate(SV_start_relative = case_when(position == "GATA1" ~  SV_start-cut_left,
                                                  position == "ATF4" ~  SV_start-cut_right,
                                                  position == "GATA1-ATF4" ~ SV_start-cut_left),
                    SV_end_relative = case_when(position == "GATA1" ~  SV_end-cut_left,
                                                position == "ATF4" ~  SV_end-cut_right,
                                                position == "GATA1-ATF4" ~ SV_end-cut_right))

vcfs_df2 <- vcfs_df2 %>% filter(!is.na(position))


vcfs_df2 <- vcfs_df2 %>% mutate(SV_size = case_when(abs(SVLEN)<=small_indel~ "Unmodified",
                                                    between(abs(SVLEN),small_indel+1,50)~ "small (<=50bp)",
                                                    between(abs(SVLEN),51,200)~ "intermediate(<=200bp)",
                                                    abs(SVLEN)>200~ "Large (>200bp)")) %>% 
  mutate(SVTYPE_size = case_when(SV_size == "Unmodified"~"Unmodified",
                                 TRUE ~ SVTYPE)) %>% 
  separate_rows(RNAMES, sep=",")


# reads in vcf
all_combinations <- vcfs_df2 %>%
  filter(position %in% c("GATA1","ATF4")) %>% 
  distinct(library,RNAMES) %>% 
  expand(nesting(library,RNAMES),position=c("GATA1","ATF4")) 


#reads in bed
bed_combinations <- beds_df_unique_all %>% 
  anti_join(vcfs_df2 %>% 
              distinct(library,RNAMES), by = c("library","V4"="RNAMES")) %>% 
  expand(nesting(library,RNAMES=V4),position=c("GATA1","ATF4"))



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


# KEEP BELOW BUT DO NOT USE ---------------------------------------------------------



























## 04 - Process reads without vector  -------------------------------------------------

vcfs_df_HBA21_novec <- vcfs_df %>% 
  filter(V1 == "BCL11A_deleted", SVTYPE %in% c("DEL","INS")) %>% 
  separate_rows(RNAMES, sep = ",") 


vcfs_df_HBA21_novec <- vcfs_df_HBA21_novec %>% 
  mutate(SV_start = V2  ,
         SV_end = case_when(SVTYPE == "DEL" ~ (V2 ) + abs(SVLEN)-1,
                            SVTYPE== "INS" ~ V2  + 10),
         SV_start_relative = SV_start - cut_hba2,
         SV_end_relative =  SV_end - cut_hba2) %>% 
  filter((SV_start_relative -2 < 0) & (SV_end_relative + 2) > 0) %>% 
  distinct(library,RNAMES, SVTYPE,SV_start,SVLEN,SV_end,SV_start_relative,SV_end_relative) %>% 
  mutate(         locus = "BCL11A",
                  site = "BCL11A2/1",
                  vector = "no",
                  HBA15Del = "yes")





## 05 - Rescue reads with other SV ---------------------------------------------------------

vcfs_df_HBA21_other <- vcfs_df %>% 
  ungroup %>%
  filter(V1 == "BCL11A_deleted",
         SVTYPE %in% c("INV","DUP")) %>% 
  mutate(SV_start = (V2) ,
         SV_end = (V2 ) + abs(SVLEN),
         SV_start_relative = SV_start - cut_hba2,
         SV_end_relative =  SV_end - cut_hba2) %>% 
  filter((SV_start_relative-2)<0 & (SV_end_relative+2) >0) %>% ## remove INDEL far from cutting site
  separate_rows(RNAMES, sep = ",") %>% 
  group_by(library,RNAMES) %>% 
  slice_head(n = 1) %>%  ### for INV, there are somtimes 2 almost identical SVs, keep the first
  mutate(         locus = "BCL11A",
                  site = "BCL11A2/1",
                  vector = "no",
                  HBA15Del = "yes")


## 06 - Aggregate HBA2/1 Reads  ------------------------------------------------------------

# reads with SV : with vector, without vector, with other SV
# then add reads with HBA deletion but no other SV

vcfs_df_HBA21 <- vcfs_df_HBA21_novec %>%
              ungroup %>% 
              distinct(library,RNAMES,SV_start,SV_start_relative,SV_end,SV_end_relative,locus,site, HBA15Del, vector,SVLEN,SVTYPE) %>% 
  bind_rows(vcfs_df_HBA21_other %>%
              ungroup %>% 
              distinct(library,RNAMES,SV_start,SV_start_relative,SV_end,SV_end_relative,locus,site, HBA15Del, vector,SVLEN,SVTYPE))  %>% 
  full_join(beds_df_unique) %>% ## this is the complete list of reads with HBA deletion
  replace_na(list(SV_start = 0,
                  SV_end = 1,
                  locus = "BCL11A",
                  site = "BCL11A2/1",
                  vector = "no",
                  HBA15Del = "yes",
                  #type = "noINDELs",
                  SVLEN = 0,
                  SVTYPE = "none"))



# Reads without HBA deletion ---------------------------------------------

## 01 - Get all aligned reads without HBA deletion -----------------------------------------

bed <- list.files(path = "./",pattern = "locusAligned.sorted.end2end.withoutDEL.bed", recursive = T ,full.names = T)

names(bed) <- str_match(bed,pattern = "./(.+)/locus")[,2]

beds<- lapply(bed, function(x){
  if(file.size(x)>0){
    read.delim(x, header=F,comment.char = "#")
  } else {
    NULL
  }
})


beds_df <- bind_rows(beds,.id='library')



## Get reads with vector
vector_containing_reads_HBA21 <- beds_df %>% 
  filter(V1 != "BCL11A") %>%
  distinct(library, RNAMES=V4) %>%
  mutate(vector_bed = "yes")


## get unique reads per library with HBA deletion
beds_df_unique <- beds_df %>% 
  distinct(library,RNAMES=V4) %>% 
  left_join(vector_containing_reads_HBA21) %>% 
  replace_na(list(vector_bed="no"))

rm(vector_containing_reads_HBA21)
rm(beds)
rm(beds_df)




## 02 - Load variant calling results -----------------------------------

files <- list.files(path = "./", recursive = T, full.names = T, pattern = "withoutDEL.sniffles.vcf")
files

files <- files[file.size(files)>0]

names(files) <- str_match(files,pattern = "./(.+)/locus")[,2]

vcfs<- lapply(files, function(x){
  read.delim(x, header=F,comment.char = "#",colClasses = c("character","integer","integer",rep("character",7)))
})

vcfs_df <- bind_rows(vcfs,.id='library')

rm(files); rm(vcfs)

vcfs_df <- vcfs_df %>% 
  #filter(str_detect(V8, "SVTYPE=DEL|SVTYPE=INS|SVTYPE=BND")) %>%
  mutate(SVTYPE = str_match(V8,"SVTYPE=([A-Z/]+)")[,2]) %>%
  group_by(library,V1,SVTYPE) %>% 
  arrange(V2,.by_group = T) %>% 
  mutate(SVLEN = as.numeric(str_match(V8,"SVLEN=([-0-9]+)")[,2]),
         RE = as.numeric(str_match(V8,"RE=([-0-9]+)")[,2]),
         ID = row_number(),
         CHR2 = str_match(V8,"CHR2=([[:alnum:]-_]+)")[,2],
         RNAMES = str_match(V8,"RNAMES=([[:alnum:]-_,]+)")[,2]) %>% 
  filter(V7 != "UNRESOLVED")





## 04 - @HBA reads without vector  -------------------------------------------------
vcfs_df_HBA2_novec <- vcfs_df %>% 
  filter(V1 == "BCL11A",SVTYPE %in% c("DEL","INS")) %>% 
  separate_rows(RNAMES, sep = ",") 


vcfs_df_HBA2_novec <- vcfs_df_HBA2_novec %>% 
  mutate(SV_start = V2 ,
         SV_end = case_when(SVTYPE == "DEL" ~ (V2 ) + abs(SVLEN)-1,
                            SVTYPE== "INS" ~ V2  + 10)) %>% 
  mutate(site = case_when(((SV_end + SV_start)/2 - (cut_hba1 + cut_hba2)/2) > 0 ~ "BCL11A1",
                          ((SV_end + SV_start)/2 - (cut_hba1 + cut_hba2)/2) < 0 ~ "BCL11A2",
                          TRUE ~ "unknown_site")) %>% 
  mutate(SV_start_relative = case_when(site == "BCL11A2" ~ SV_start - cut_hba2,
                                       site == "BCL11A1" ~ SV_start - cut_hba1,
                                       TRUE ~ NA),
         SV_end_relative = case_when(site == "BCL11A2" ~ SV_end - cut_hba2,
                                     site == "BCL11A1" ~ SV_end - cut_hba1,
                                     TRUE ~ NA)) %>% 
  filter( (SV_start_relative-2) < 0 & (SV_end_relative+2) > 0) %>% 
  distinct(library,RNAMES, SVTYPE,SV_start,SVLEN,SV_end,SV_start_relative,SV_end_relative,site) %>% 
  mutate(locus = "BCL11A",
         vector = "no",
         HBA15Del = "no")


## 05 - @HBA reads with other SV ---------------------------------------------------
vcfs_df_HBA2_other <- vcfs_df %>% 
  ungroup %>%
  filter(V1 == "BCL11A",
         SVTYPE %in% c("INV","DUP")) %>% 
  mutate(SV_start = V2  ,
         SV_end = V2 + abs(SVLEN)) %>% 
  separate_rows(RNAMES, sep = ",") %>% 
  group_by(library,RNAMES) %>% 
 # slice_head(n = 1) %>%  ### for INV, there are somtimes 2 almost identical SVs, keep the first
  mutate(site = case_when(SVTYPE=="INV" & (SV_start-10 < cut_hba2) & (SV_end +10 > cut_hba1) ~ "BCL11A2",
                          ((SV_end + SV_start)/2 - (cut_hba1 + cut_hba2)/2) > 100 ~ "BCL11A1",
                          ((SV_end + SV_start)/2 - (cut_hba1 + cut_hba2)/2) < 100 ~ "v2",
                          TRUE ~ "unknown_site")) %>% 
  mutate(SV_start_relative = case_when(site == "BCL11A2" ~ SV_start - cut_hba2,
                                       site == "BCL11A1" ~ SV_start - cut_hba1,
                                       site == "BCL11A2" ~ SV_start-cut_hba2,
                                       TRUE ~ NA),
         SV_end_relative = case_when(site == "BCL11A2" ~ SV_end - cut_hba2,
                                     site == "BCL11A1" ~ SV_end - cut_hba1,
                                     site == "BCL11A2" ~ SV_end-cut_hba1,
                                     TRUE ~ NA)) %>% 
  filter( (SV_start_relative-2) < 20 & (SV_end_relative+2) > 0) %>% 
  mutate(         locus = "BCL11A",
                  vector = "no",
                  HBA15Del = "no")


## 06 - @HBA reads aggregated ---------------------------------------------------

vcfs_df_HBA <- vcfs_df_HBA2_novec %>%
              ungroup %>% 
              distinct(library,RNAMES,SV_start,SV_end,locus,site, HBA15Del, vector,SVLEN,SV_start_relative,SV_end_relative,SVTYPE) %>% 
  bind_rows(vcfs_df_HBA2_other %>%
              ungroup %>% 
              distinct(library,RNAMES,SV_start,SV_end,locus,site, HBA15Del, vector,SVLEN, SV_start_relative,SV_end_relative,SVTYPE)) 

filing_table <- beds_df_unique %>% 
  tidyr::expand(nesting(library,RNAMES),site=c("BCL11A2","BCL11A1")) %>%
  anti_join(vcfs_df_HBA %>% filter(site!="BCL11A2"& site!="BCL11A1") %>%  select(library,RNAMES))

vcfs_df_HBA = vcfs_df_HBA %>%
  full_join(filing_table) %>% 
  full_join(beds_df_unique ) %>% 
  replace_na(list(locus = "BCL11A",
                  vector = "no",
                  HBA15Del = "no",
                  #type = "noINDELs",
                  SVLEN = 0,
                  SVTYPE = "none",
                  # SVTYPE_paper = "none",
                  #SVTYPE_category = "Unmodified",
                  SV_start = 0,
                  SV_end = 0, SV_start_relative = 0,
                  SV_end_relative=0) ) 





# Aggregate reads with and without HBA del ------------
vcfs_df_HBA <- vcfs_df_HBA %>%
  bind_rows(vcfs_df_HBA21)%>% 
  left_join(metadata %>% 
              select(library, pooled_libraries, Cells,Target,donor_template )) 


vcfs_df_HBA %>% ungroup %>% distinct(library,RNAMES) %>% count()


# Vector integration ---------------------------------------------------------------
##01 - select reads with vector ----------------------------------

## reads with the HBA15 deletion
bed_del <- list.files(path = "./",pattern = "withDEL.locusDelAligned.annotated.bed", recursive = T ,full.names = T)

names(bed_del) <- str_match(bed_del,pattern = "./(.+)/locus")[,2]

beds_del<- lapply(bed_del, function(x){
  if(file.size(x)>0){
    read.delim(x, header=F,comment.char = "#")
  } else {
    NULL
  }
})

beds_del_df <- bind_rows(beds_del,.id='library')
rm(bed_del)
rm(beds_del)



## all reads
bed <- list.files(path = "./",pattern = "locusAligned.sorted.end2end.annotated.bed", recursive = T ,full.names = T)

names(bed) <- str_match(bed,pattern = "./(.+)/locus")[,2]


beds<- lapply(bed, function(x){
  if(file.size(x)>0){
    read.delim(x, header=F,comment.char = "#")
  } else {
    NULL
  }
})

beds_nodel_df <- bind_rows(beds,.id='library')

rm(beds);rm(bed)

## reads without deletion
beds_nodel_df <- beds_nodel_df %>% # remove reads with HBA deletion
  anti_join(beds_del_df %>%
              distinct(library,V4))


# merge reads with the deletion (aligned on the deleted locus) and reads without the HBA15 deletion (aligned on the locus)
beds_df <-  bind_rows("noHBA15del"= beds_nodel_df, "HBA15del"=beds_del_df, .id="type")
rm(beds_nodel_df)



## 04 - Detect HBA15 deletion/ duplication  ---------------------
hba15_del <- beds_df %>%  #count number of HBA15 fragments (identify deletions or duplications)
  group_by(library,V4) %>% 
  dplyr::summarise(hba15_fragment = length(which(str_detect(V12, "deletion")))) 







## 07 - Add missing combination of HBA sites (with no vector) --------------------------------------------------------

missing_reads <- beds_df %>% distinct(type,library,V4) %>% 
  mutate(site = case_when(type == "HBA15del" ~ "BCL11A2/1",
                          TRUE ~ "BCL11A2,BCL11A1")) %>% 
  separate_rows(site,sep = ",")


# Annotate INDELS --------------------------------------------------------------------


vcfs_df_HBA <- vcfs_df_HBA %>% mutate(SVTYPE= case_when(SVTYPE=="BND" & SVLEN > small_indel ~ "INS",
                                                        SVTYPE=="BND" & SVLEN < -small_indel ~ "DEL",
                                                        TRUE ~ SVTYPE),
                                      
                                      SVTYPE_paper = case_when(abs(SVLEN)<=  small_indel ~ "none",
                                                               SVLEN < -small_indel ~ "DEL",
                                                               SVLEN > small_indel & SVTYPE !="DUP" & SVTYPE != "INV" &SVTYPE != "INVDUP" ~ "INS",
                                                               TRUE ~ SVTYPE),
                                      
                                      SVTYPE_category = case_when(between(abs(SVLEN), small_indel+1,50) ~ "small",
                                                                  between(abs(SVLEN), 51, 200) ~ "intermediate",
                                                                  abs(SVLEN)> 200 ~ "large",
                                                                  TRUE ~ "Unmodified"), 
                                      
                                      typeIndels = case_when(abs(SVLEN)<= small_indel ~ "noINDELs",
                                                             TRUE ~ "INDELs"))


# Collapse everything ---------------------------------------------------------------


yyy <- vcfs_df_HBA %>% 
  full_join(missing_reads, by = c("RNAMES" = 'V4', "library","site")) %>% 
  left_join(hba15_del, by = c("RNAMES"="V4","library"))

yyy  <- yyy %>%  
  mutate(SVTYPE_paper = case_when(hba15_fragment>1 & site == "BCL11A2" ~ "DUP", TRUE ~ SVTYPE_paper))



save(list=c("yyy","metadata","vcfs_df_HBA"),file = "BCL11A.rdata")


write.table(yyy , "BCL11A_complete_table.csv", sep = ";", quote = F, row.names  = F)



ggplot(yyy, aes(site, fill=SVTYPE_paper)) + 
  geom_bar(position = "stack", col = "black") + 
  facet_wrap(~library) + 
  theme_bw()


