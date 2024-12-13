## set file path to current directory

setwd(system("pwd", intern = T) )


library(dplyr)
library(tibble)
library(R.utils)


args <- commandArgs(TRUE)

#Control_name = 

Sample_name = args[1]


Sample_1_base = as.numeric(args[2])
Sample_2_base = as.numeric(args[3])

Control_1_base = as.numeric(args[4])
Control_2_base = as.numeric(args[5])

Control_name = args[6]

# open files

Sample_1 <- read.delim(file =(paste("kmer_centromere/ChIP-seq/kmer_dbs/51/",Sample_name,"_r1_normkcounts.txt", 
                                     sep="")), header = F, sep = "")
Sample_2 <- read.delim(file =(paste("kmer_centromere/ChIP-seq/kmer_dbs/51/",Sample_name,"_r2_normkcounts.txt", 
                                     sep="")), header = F, sep = "")


Control_1 <- read.delim(file =(paste("kmer_centromere/ChIP-seq/kmer_dbs/51/control/",Control_name,"_r1_normkcounts.txt",
                                     sep="")), header = F, sep = "")
Control_2 <- read.delim(file =(paste("kmer_centromere/ChIP-seq/kmer_dbs/51/control/",Control_name,"_r2_normkcounts.txt",
                                     sep="")), header = F, sep = "")


# set column names

colnames(Control_1) <- c('K-mer', 'Control_1_count', 'Control_1_normalised')
colnames(Control_2) <- c('K-mer', 'Control_2_count', 'Control_2_normalised')

colnames(Sample_1) <- c('K-mer','Sample_1_count', 
                        'Sample_1_normalised')
colnames(Sample_2) <- c('K-mer','Sample_2_count', 
                        'Sample_2_normalised')

# full_join #set up script to join K-mer enrichment from different files

# full join replicates together for control 


Control_average <- full_join(Control_1, Control_2, by = "K-mer")

# full_join vs controls

Sample_1_vs_Control <- full_join(Sample_1, Control_average, by = "K-mer")
Sample_2_vs_Control <- full_join(Sample_2, Control_average, by = "K-mer")

# remove normalised columns

Sample_1_vs_Control <- select(Sample_1_vs_Control,!contains("normalised"))
Sample_2_vs_Control <- select(Sample_2_vs_Control,!contains("normalised"))

# I replaced NAs with 1 instead -as possibility of K-mer count of 1 being excluded (underestimating fold change if at all)


Sample_1_vs_Control <- mutate_all(Sample_1_vs_Control, 
                                  ~replace(., is.na(.), 1))
Sample_2_vs_Control <- mutate_all(Sample_2_vs_Control, 
                                  ~replace(., is.na(.), 1))

# total base count values - extracted from all reads (not just filtered for centromere) - needed for next step
#loaded above for sample and control with args


# normalisation step


Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Sample_1_normalised = 
                                 as.numeric(Sample_1_count/Sample_1_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Sample_2_normalised = 
                                 as.numeric(Sample_2_count/Sample_2_base))

Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_1_normalised = 
                                 as.numeric(Control_1_count/Control_1_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_2_normalised = 
                                 as.numeric(Control_2_count/Control_2_base))


Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_1_normalised = 
                                 as.numeric(Control_1_count/Control_1_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_2_normalised = 
                                 as.numeric(Control_2_count/Control_2_base))



### Filter so the sample has more than a normalised count of 5e-09 (row numbers containing normalised count for Sample only)
# In this example there are 2 controls and 1 sample - so "K-mer" column, 3 count columns, and then the Sample normalised count column is number 5

Sample_1_vs_Control <- Sample_1_vs_Control %>% filter(.[[5]]>5e-09)
Sample_2_vs_Control <- Sample_2_vs_Control %>% filter(.[[5]]>5e-09)

### Take average of normalised values for control (in this example, there are 3 controls, but if more divide by the number of control samples)

Sample_1_vs_Control <- mutate(Sample_1_vs_Control, average_Control = as.numeric(Control_1_normalised + 
Control_2_normalised)/2)

Sample_2_vs_Control <- mutate(Sample_2_vs_Control, average_Control = as.numeric(Control_1_normalised + 
Control_2_normalised)/2)


### calculate fold change

Sample_1_vs_Control <- mutate(Sample_1_vs_Control, Sample_1_v_Control = as.numeric(Sample_1_normalised/average_Control))
Sample_2_vs_Control <- mutate(Sample_2_vs_Control, Sample_2_v_Control = as.numeric(Sample_2_normalised/average_Control))


### Filter by fold change - can choose whatever you like here


Sample_1_vs_Control_fc2 <- filter(Sample_1_vs_Control, Sample_1_v_Control >2)
Sample_2_vs_Control_fc2 <- filter(Sample_2_vs_Control, Sample_2_v_Control >2)


write.table(Sample_1_vs_Control_fc2, file=paste("kmer_centromere/ChIP-seq/enriched_kmers/",Sample_name,"_1_vs_",Control_name,"_fc2.txt", sep=""),sep =" ", row.names = F, 
col.names = T, quote = F)
write.table(Sample_2_vs_Control_fc2, file=paste("kmer_centromere/ChIP-seq/enriched_kmers/",Sample_name,"_2_vs_",Control_name,"_fc2.txt", sep=""),sep =" ", row.names = F, 
col.names = T, quote = F)
