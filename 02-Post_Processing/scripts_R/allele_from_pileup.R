#### CORRE GUILLAUME @ GENETHON
#### Nov 11th, 2024
#### Analysis of ATF4 and GATA1 locus on the same reads

# We here work on files generated with samtools mpileup program. 
# For each library, we have the nucleotide composition at 5 nucleotides around the editing site.
# for each position and each read, we have the corresponding base, incuding INDELS

# libraries 'template****' are fake reads generated to validate the pipeline.

  # template_BCL11A_WT is the WT sequence : NNNNNN site1 NNNN site2 NNNNNNN
  # template_BCL11A_DEL contains a deletion between the 2 editing sites : NNNNNN site1/site2 NNNNNNN
  # template_BCL11A_overDEL contains a deletion between the 2 editing sites and extra deletion  NNNNNN sit../..te2 NNNNNNN
  # template_BCL11A_INV contains an inversion between the 2 editing sites
  # template_BCL11A_microDEL contains small deletions at each loci :   NNNN sit.. NNNN ..te2 NNNN




# ================================================================
# Mpileup file generation for the 2 loci from the bam file
# ================================================================

# for lib in *; do 
#   samtools mpileup $lib/locus/$lib.locusAligned.sorted.end2end.bam \
#   -r BCL11A:8490-8495 \
#   --output-QNAME  \
#   --no-output-ins \
#   --no-output-del \
#   --no-output-ins-mods \
#   --no-output-ends \
#   --min-BQ 0 \
#   -aa \
#   -o $lib"_ATF4.mpileup" \ 
#   echo $lib;\
# done


# for lib in *; do \
#   samtools mpileup $lib/locus/$lib.locusAligned.sorted.end2end.bam \
#   -r BCL11A:5264-5270 
#   --output-QNAME  \
#   --no-output-ins \
#   --no-output-del \
#   --no-output-ins-mods \
#   --no-output-ends \
#   --min-BQ 0 \
#   -aa \
#   -o $lib"_GATA1.mpileup"; \
#   echo $lib;\
# done


library(tidyverse);
library(ggpubr);
library(patchwork);
library(writexl)
library(ggprism)


## Load mpileup files 

files <- list.files(path = "basecalling_hac3.3/",pattern = "mpileup",full.names = T)

#files <- grep(files2,pattern = "pooled", invert = T,value=T)
#files <- c(files,grep(files2,pattern = "WT", value=T) )

files_list <- lapply(files, function(x) {
  read.delim(x, sep="", header = F)
  }
  )
names(files_list) <- str_remove(basename(files),".mpileup")

files_df <- bind_rows(files_list,.id = "library")



## extract and format data for ATF4 locus

## We end up with a dataframe with 3 columns, including the library name, read name and 5nt sequence
# library                         V7                                   ATF4_seq
# FAV39284_pass_dca5b723_41ce3bd7 34113f84-25fe-4033-b163-51cd4e96ce84 CATCC   
# FAV39284_pass_dca5b723_41ce3bd7 5b46cacd-5272-4d46-ac3b-e9c4ffc93a42 CATCC   
# FAV39284_pass_dca5b723_41ce3bd7 7783bf13-2ba0-454b-a3c8-9bc8979138b8 CATCC   
# FAV39284_pass_dca5b723_41ce3bd7 998a6299-ff43-4759-81bd-b6bf9049492f *****   
# FAV39284_pass_dca5b723_41ce3bd7 d6faf348-bd59-48c8-9505-aa3875818530 CATCC 
# ...


seq_ATF4 <- files_df %>% 
  filter(str_detect(library,"ATF4")) %>%  
  mutate(V5 = str_replace_all(V5,pattern = "[-\\+0-9]+",replacement = "")) %>% 
  rowwise() %>% 
  mutate(V5 = paste(unlist(strsplit(V5, "", fixed = TRUE)),collapse = ",")) %>% 
  separate_rows(V5,V7,sep = ",") %>% 
  pivot_wider(names_from = V2,values_from = V5,id_cols = c(library,V7),values_fill = "*") %>% 
  select(library,V7,`8491`:`8495`) %>% 
  unite(col = "ATF4_seq",3:last_col(),sep = "") %>% 
  mutate(ATF4_seq = toupper(ATF4_seq)) %>%
  mutate(library = str_remove(library,"_ATF4"))


## extract and format data for GATA1 locus
seq_GATA1 <- files_df %>%
  filter(str_detect(library,"GATA1")) %>%  
  mutate(V5 = str_replace_all(V5,pattern = "[-\\+0-9]+",replacement = "")) %>% 
  rowwise() %>% 
  mutate(V5 = paste(unlist(strsplit(V5, "", fixed = TRUE)),collapse = ",")) %>% 
  separate_rows(V5,V7,sep = ",") %>% 
  pivot_wider(names_from = V2,values_from = V5,id_cols = c(library,V7),values_fill = "*") %>% 
  select(library,V7,`5265`:`5269`) %>% 
  unite(col = "GATA1_seq",3:last_col(),sep = "") %>% 
  mutate(GATA1_seq = toupper(GATA1_seq))%>% 
  mutate(library = str_remove(library,"_GATA1"))


## Combine data for the 2 loci

## Check if sequence is different from WT.
## Classify each locus as unmodified or modified
## Classify pair of loci as unmodified, edited at 1 or 2 loci

seq <- seq_GATA1 %>% left_join(seq_ATF4)%>% 
  drop_na() %>% 
  mutate(GATA1_type =case_when(GATA1_seq == "GTGAT" ~ "unmodified",
                               TRUE ~ "modified"),
         ATF4_type = case_when(ATF4_seq == "CATCC" ~ "unmodified",
                               TRUE ~ "modified")) %>% 
  filter(str_detect(library,"WT") | str_starts(library, "FA")) %>% 
  mutate(GATA1_ATF4 = case_when(GATA1_type == "modified" & ATF4_type == "modified" ~ "Dual edition",
                               GATA1_type == "modified" | ATF4_type == "modified" ~ "Mono edition",
                               TRUE ~ "Unmodified"))



# summarise : count nucleotide sequence per locus in 5 bp window
seq_stat <- seq %>%
  group_by(library) %>%
  count(GATA1_seq,ATF4_seq) %>%
  arrange(desc(n)) %>%
  group_by(library) %>%
  mutate(proportion = n/sum(n) *100) 


write.table(seq_stat, "basecalling_hac3.3/for_the_paper/GATA1_ATF4_pattern_5nt.tsv", sep="\t",quote = F, row.names = F)


# summarise : count edited/unedited reads per locus in 5bp window
seq_stat_bin_class <- seq %>%
  group_by(library) %>%
  count(GATA1_type,ATF4_type) %>%
  arrange(desc(n)) %>%
  group_by(library) %>%
  mutate(proportion = n/sum(n) *100) 


write.table(seq_stat_bin_class, "basecalling_hac3.3/for_the_paper/GATA1_ATF4_editing.tsv", sep="\t",quote = F, row.names = F)


# Summarise : count co-edited reads @ GATA1 and/or ATF4

seq_stat_bin_coediting <- seq_stat_bin_class %>%
  mutate(GATA1_ATF4 = case_when(GATA1_type == "modified" & ATF4_type == "modified" ~ "Dual edition",
                                GATA1_type == "modified" | ATF4_type == "modified" ~ "Mono edition",
                                TRUE ~ "Unmodified")) %>% 
  group_by(library,GATA1_ATF4) %>%
  summarise(count = sum(n)) %>%
  group_by(library) %>%
  mutate(proportion =count/sum(count) *100) 


write.table(seq_stat_bin_coediting, "basecalling_hac3.3/for_the_paper/GATA1_ATF4_coediting.tsv", sep="\t",quote = F, row.names = F)


## make graph of co-editing

p_gata <- ggplot(seq, aes(library,fill=interaction(GATA1_type))) +
  geom_bar(col = "black", position = position_fill()) +
  ggprism::theme_prism() + 
  coord_flip() + 
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::percent_format())+
  labs(y = "% reads", x = NULL, title = 'GATA1 position (BCL11A:5265-5269)')+
  theme(legend.position = "none")



p_atf4 <- ggplot(seq, aes(library,fill=interaction(ATF4_type))) +
  geom_bar(col = "black", position = position_fill()) + 
  ggprism::theme_prism() + 
  coord_flip() + 
  scale_fill_viridis_d()+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y = "% reads", x = NULL, title = 'ATF4 position (BCL11A:8491-8495)')+
  theme(axis.text.y = element_blank())


p_combined <- ggplot(seq, aes(library,fill=GATA1_ATF4)) +
  geom_bar(col = "black", position = position_fill()) + 
  ggprism::theme_prism() + 
  coord_flip() + 
  scale_fill_viridis_d()+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y = "% reads", x = NULL , title = 'GATA1.ATF4 position')


library(patchwork)
x11()
(p_gata | p_atf4 )/p_combined 





## single graphs at each position of the 5bp window 


ATF4_pattern = "CATCC"
GATA1_pattern = "GTGAT"

# Position on the enrichment region
ATF4_position <- 8491:8495
GATA1_position <- 5265:5269


GATA1_plot_seq <- list()
GATA1_plot_type <- list()

ATF4_plot_seq <- list()
ATF4_plot_type <- list()

plot_data <- list()

for(i in 1:length(GATA1_position)){
  
  
  
  seq_GATA1 <- files_df %>% filter(str_detect(library,"GATA1")) %>%  
    mutate(V5 = str_replace_all(V5,pattern = "[-\\+0-9]+",replacement = "")) %>% 
    filter(str_detect(library,"WT") | str_starts(library, "FA")) %>% 
    rowwise() %>% 
    mutate(V5 = paste(unlist(strsplit(V5, "", fixed = TRUE)),collapse = ",")) %>% 
    separate_rows(V5,V7,sep = ",") %>% 
    pivot_wider(names_from = V2,values_from = V5,id_cols = c(library,V7),values_fill = "*") %>% 
    select(library,V7,paste(GATA1_position[i])) %>% 
    unite(col = "read_sequence",3:last_col(),sep = "") %>% 
    mutate(read_sequence = toupper(read_sequence))%>% 
    mutate(locus = "GATA1",
           library = str_remove(library,"_GATA1"),
           ) %>% 
    mutate(type =case_when(read_sequence == str_sub(GATA1_pattern,i,i) ~ "unmodified",
                                 TRUE ~ "modified"),
           position = GATA1_position[i],
           reference_sequence = str_sub(GATA1_pattern,i,i))
  
  
  GATA1_plot_type[[i]] <- ggplot(seq_GATA1, aes(library,fill=interaction(type))) +
    geom_bar(col = "black", position = position_fill()) +
    ggprism::theme_prism(base_size = 9) +
    coord_flip() +
    scale_fill_manual(values = c("red3","green4")) +
    scale_y_continuous(labels = scales::percent_format())+
    labs(y = "% reads", x = NULL, title = paste('GATA1 (BCL11A:',GATA1_position[i],")",sep = ""))

  
  
  GATA1_plot_seq[[i]] <- ggplot(seq_GATA1, aes(library,fill=interaction(read_sequence))) +
    geom_bar(col = "black", position = position_fill()) +
    ggprism::theme_prism(base_size = 9) + 
    coord_flip() + 
    scale_fill_manual(values = c("grey","green4","blue","orange","red3")) +
    scale_y_continuous(labels = scales::percent_format())+
    labs(y = "% reads", x = NULL, title = paste('GATA1 (BCL11A:',GATA1_position[i],")",sep = ""))
  

  
  seq_ATF4 <- files_df %>% filter(str_detect(library,"ATF4")) %>%  
    mutate(V5 = str_replace_all(V5,pattern = "[-\\+0-9]+",replacement = "")) %>% 
    filter(str_detect(library,"WT") | str_starts(library, "FA")) %>% 
    rowwise() %>% 
    mutate(V5 = paste(unlist(strsplit(V5, "", fixed = TRUE)),collapse = ",")) %>% 
    separate_rows(V5,V7,sep = ",") %>% 
    pivot_wider(names_from = V2,values_from = V5,id_cols = c(library,V7),values_fill = "*") %>% 
    select(library,V7,paste(ATF4_position[i])) %>% 
    unite(col = "read_sequence",3:last_col(),sep = "") %>% 
    mutate(read_sequence = toupper(read_sequence))%>% 
    mutate(locus = "ATF4",
           library = str_remove(library,"_ATF4"),
    ) %>% 
    mutate(type =case_when(read_sequence == str_sub(ATF4_pattern,i,i) ~ "unmodified",
                           TRUE ~ "modified"),
           position = ATF4_position[i],
           reference_sequence = str_sub(ATF4_pattern,i,i))
  
  
  ATF4_plot_type[[i]] <- ggplot(seq_ATF4, aes(library,fill=interaction(type))) +
    geom_bar(col = "black", position = position_fill()) +
    ggprism::theme_prism(base_size = 9) +
    coord_flip() +
    scale_fill_manual(values = c("red3","green4")) +
    scale_y_continuous(labels = scales::percent_format())+
    labs(y = "% reads", x = NULL, title = paste('ATF4 (BCL11A:',ATF4_position[i],")",sep = ""))
  
  
  
  ATF4_plot_seq[[i]] <- ggplot(seq_ATF4, aes(library,fill=interaction(read_sequence))) +
    geom_bar(col = "black", position = position_fill()) +
    ggprism::theme_prism(base_size = 9) + 
    coord_flip() + 
    scale_fill_manual(values = c("grey","green4","blue","orange","red3")) +
    scale_y_continuous(labels = scales::percent_format())+
    labs(y = "% reads", x = NULL, title = paste('ATF4 (BCL11A:',ATF4_position[i],")",sep = ""))
  
  
  GATA1_data <- seq_GATA1 %>% count(library,locus,position,reference_sequence,read_sequence,type) %>% 
    group_by(library) %>% 
    mutate(proportion = n / sum(n) * 100)
  ATF4_data <- seq_ATF4 %>% count(library,locus,position,reference_sequence,read_sequence,type) %>% 
    group_by(library) %>% 
    mutate(proportion = n / sum(n) * 100)
  
  plot_data[[paste("GATA1",GATA1_position[i],sep = "_")]] <- GATA1_data
  plot_data[[paste("ATF4",ATF4_position[i],sep = "_")]] <- GATA1_data
  
}

writexl::write_xlsx(plot_data,col_names = T,format_headers = T,path = "basecalling_hac3.3/for_the_paper/Single_position_plot_data.xlsx")



x11()
p1 <- ggarrange(plotlist = GATA1_plot_seq, ncol = 1,nrow=5,common.legend = T)
p2 <- ggarrange(plotlist = GATA1_plot_type, ncol = 1,nrow=5,common.legend = T)

p1+p2

p3 <- ggarrange(plotlist = ATF4_plot_seq, ncol = 1,nrow=5,common.legend = T)
p4 <- ggarrange(plotlist = ATF4_plot_type, ncol = 1,nrow=5,common.legend = T)

x11()
p3+p4




###### END #######