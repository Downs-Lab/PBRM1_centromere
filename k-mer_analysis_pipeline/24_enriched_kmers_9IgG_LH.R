## set file path to current directory

setwd(system("pwd", intern = T) )


library(dplyr)
library(tibble)
library(R.utils)


args <- commandArgs(TRUE)

#Control_name = 

Sample_name = args[1]


Sample_L1_base = as.numeric(args[2])
Sample_H1_base = as.numeric(args[3])
Sample_L2_base = as.numeric(args[4])
Sample_H2_base = as.numeric(args[5])
Sample_L3_base = as.numeric(args[6])
Sample_H3_base = as.numeric(args[7])



# open files

Sample_L1 <- read.delim(file =(paste("kmer_centromere/kmer_dbs/51/",Sample_name,"_L_r1_normkcounts.txt", 
                                     sep="")), header = F, sep = "")
Sample_H1 <- read.delim(file =(paste("kmer_centromere/kmer_dbs/51/",Sample_name,"_H_r1_normkcounts.txt", 
                                     sep="")), header = F, sep = "")
Sample_L2 <- read.delim(file =(paste("kmer_centromere/kmer_dbs/51/",Sample_name,"_L_r2_normkcounts.txt", 
                                     sep="")), header = F, sep = "")
Sample_H2 <- read.delim(file =(paste("kmer_centromere/kmer_dbs/51/",Sample_name,"_H_r2_normkcounts.txt", 
                                     sep="")), header = F, sep = "")
Sample_L3 <- read.delim(file =(paste("kmer_centromere/kmer_dbs/51/",Sample_name,"_L_r3_normkcounts.txt", 
                                     sep="")), header = F, sep = "")
Sample_H3 <- read.delim(file =(paste("kmer_centromere/kmer_dbs/51/",Sample_name,"_H_r3_normkcounts.txt", 
                                     sep="")), header = F, sep = "")



Control_L1 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_L_r1_normkcounts.txt", 
                         header = F, sep = "")
Control_H1 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_H_r1_normkcounts.txt", 
                         header = F, sep = "")
Control_L2 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_L_r2_normkcounts.txt", 
                         header = F, sep = "")
Control_H2 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_H_r2_normkcounts.txt", 
                         header = F, sep = "")
Control_L3 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_L_r3_normkcounts.txt", 
                         header = F, sep = "")
Control_H3 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_H_r3_normkcounts.txt", 
                         header = F, sep = "")
Control_L4 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_L_r4_normkcounts.txt", 
                         header = F, sep = "")
Control_H4 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_H_r4_normkcounts.txt", 
                         header = F, sep = "")
Control_L5 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_L_r5_normkcounts.txt", 
                         header = F, sep = "")
Control_H5 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_H_r5_normkcounts.txt", 
                         header = F, sep = "")
Control_L6 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_L_r6_normkcounts.txt", 
                         header = F, sep = "")
Control_H6 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/WT_IgG_H_r6_normkcounts.txt", 
                         header = F, sep = "")
Control_L7 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/PBRM1KO_IgG_L_r1_normkcounts.txt", 
                         header = F, sep = "")
Control_H7 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/PBRM1KO_IgG_H_r1_normkcounts.txt", 
                         header = F, sep = "")
Control_L8 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/PBRM1KO_IgG_L_r2_normkcounts.txt", 
                         header = F, sep = "")
Control_H8 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/PBRM1KO_IgG_H_r2_normkcounts.txt", 
                         header = F, sep = "")
Control_L9 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/PBRM1KO_IgG_L_r3_normkcounts.txt", 
                         header = F, sep = "")
Control_H9 <- read.delim(file ="kmer_centromere/kmer_dbs/51/control/PBRM1KO_IgG_H_r3_normkcounts.txt", 
                         header = F, sep = "")

# set column names

colnames(Control_L1) <- c('K-mer', 'Control_L1_count', 'Control_L1_normalised')
colnames(Control_H1) <- c('K-mer', 'Control_H1_count', 'Control_H1_normalised')
colnames(Control_L2) <- c('K-mer', 'Control_L2_count', 'Control_L2_normalised')
colnames(Control_H2) <- c('K-mer', 'Control_H2_count', 'Control_H2_normalised')
colnames(Control_L3) <- c('K-mer', 'Control_L3_count', 'Control_L3_normalised')
colnames(Control_H3) <- c('K-mer', 'Control_H3_count', 'Control_H3_normalised')
colnames(Control_L4) <- c('K-mer', 'Control_L4_count', 'Control_L4_normalised')
colnames(Control_H4) <- c('K-mer', 'Control_H4_count', 'Control_H4_normalised')
colnames(Control_L5) <- c('K-mer', 'Control_L5_count', 'Control_L5_normalised')
colnames(Control_H5) <- c('K-mer', 'Control_H5_count', 'Control_H5_normalised')
colnames(Control_L6) <- c('K-mer', 'Control_L6_count', 'Control_L6_normalised')
colnames(Control_H6) <- c('K-mer', 'Control_H6_count', 'Control_H6_normalised')
colnames(Control_L7) <- c('K-mer', 'Control_L7_count', 'Control_L7_normalised')
colnames(Control_H7) <- c('K-mer', 'Control_H7_count', 'Control_H7_normalised')
colnames(Control_L8) <- c('K-mer', 'Control_L8_count', 'Control_L8_normalised')
colnames(Control_H8) <- c('K-mer', 'Control_H8_count', 'Control_H8_normalised')
colnames(Control_L9) <- c('K-mer', 'Control_L9_count', 'Control_L9_normalised')
colnames(Control_H9) <- c('K-mer', 'Control_H9_count', 'Control_H9_normalised')


colnames(Sample_L1) <- c('K-mer','Sample_L1_count', 
                        'Sample_L1_normalised')
colnames(Sample_H1) <- c('K-mer','Sample_H1_count', 
                        'Sample_H1_normalised')
colnames(Sample_L2) <- c('K-mer','Sample_L2_count', 
                        'Sample_L2_normalised')
colnames(Sample_H2) <- c('K-mer','Sample_H2_count', 
                        'Sample_H2_normalised')
colnames(Sample_L3) <- c('K-mer','Sample_L3_count', 
                        'Sample_L3_normalised')
colnames(Sample_H3) <- c('K-mer','Sample_H3_count', 
                        'Sample_H3_normalised')

# full_join #set up script to join K-mer enrichment from different files

# full join replicates together for control 


Control_1 <- full_join(Control_L1, Control_H1, by = "K-mer")
Control_2 <- full_join(Control_L2, Control_H2, by = "K-mer")
Control_3 <- full_join(Control_L3, Control_H3, by = "K-mer")
Control_4 <- full_join(Control_L4, Control_H4, by = "K-mer")
Control_5 <- full_join(Control_L5, Control_H5, by = "K-mer")
Control_6 <- full_join(Control_L6, Control_H6, by = "K-mer")
Control_7 <- full_join(Control_L7, Control_H7, by = "K-mer")
Control_8 <- full_join(Control_L8, Control_H8, by = "K-mer")
Control_9 <- full_join(Control_L9, Control_H9, by = "K-mer")

Control_12 <- full_join(Control_1, Control_2, by = "K-mer")
Control_123 <- full_join(Control_12, Control_3, by = "K-mer")
Control_1234 <- full_join(Control_123, Control_4, by = "K-mer")
Control_12345 <- full_join(Control_1234, Control_5, by = "K-mer")
Control_123456 <- full_join(Control_12345, Control_6, by = "K-mer")
Control_1234567 <- full_join(Control_123456, Control_7, by ="K-mer")
Control_12345678 <- full_join(Control_1234567, Control_8, by ="K-mer")
Control_average <- full_join(Control_12345678, Control_9, by = "K-mer")

# full join high and low salt together for sample

Sample_1 <- full_join(Sample_L1, Sample_H1, by = "K-mer")
Sample_2 <- full_join(Sample_L2, Sample_H2, by = "K-mer")
Sample_3 <- full_join(Sample_L3, Sample_H3, by = "K-mer")

# full_join vs controls

Sample_1_vs_Control <- full_join(Sample_1, Control_average, by = "K-mer")
Sample_2_vs_Control <- full_join(Sample_2, Control_average, by = "K-mer")
Sample_3_vs_Control <- full_join(Sample_3, Control_average, by = "K-mer")

# remove normalised columns

Sample_1_vs_Control <- select(Sample_1_vs_Control,!contains("normalised"))
Sample_2_vs_Control <- select(Sample_2_vs_Control,!contains("normalised"))
Sample_3_vs_Control <- select(Sample_3_vs_Control,!contains("normalised"))

# I replaced NAs with 1 instead -as possibility of K-mer count of 1 being excluded (underestimating fold change if at all)


Sample_1_vs_Control <- mutate_all(Sample_1_vs_Control, 
                                  ~replace(., is.na(.), 1))
Sample_2_vs_Control <- mutate_all(Sample_2_vs_Control, 
                                  ~replace(., is.na(.), 1))
Sample_3_vs_Control <- mutate_all(Sample_3_vs_Control, 
                                  ~replace(., is.na(.), 1))

# total base count values - extracted from all reads (not just filtered for centromere) - needed for next step
#loaded above for sample with args and below for the 9 IgG controls


Control_L1_base <- 2869929771
Control_H1_base <- 4087594056
Control_L2_base <- 1575217926
Control_H2_base <- 1254249900
Control_L3_base <- 930748261
Control_H3_base <- 1407113100
Control_L4_base <- 3250859650
Control_H4_base <- 1826984541
Control_L5_base <- 2577798743
Control_H5_base <- 2515203189
Control_L6_base <- 2509780898
Control_H6_base <- 3292206355
Control_L7_base <- 2146562152
Control_H7_base <- 2115313008
Control_L8_base <- 1873976526
Control_H8_base <- 2283733579
Control_L9_base <- 2266606176
Control_H9_base <- 2354414293



# normalisation step


Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Sample_L1_normalised = 
                                 as.numeric(Sample_L1_count/Sample_L1_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Sample_H1_normalised = 
                                 as.numeric(Sample_H1_count/Sample_H1_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Sample_L2_normalised = 
                                 as.numeric(Sample_L2_count/Sample_L2_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Sample_H2_normalised = 
                                 as.numeric(Sample_H2_count/Sample_H2_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Sample_L3_normalised = 
                                 as.numeric(Sample_L3_count/Sample_L3_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Sample_H3_normalised = 
                                 as.numeric(Sample_H3_count/Sample_H3_base))

Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L1_normalised = 
                                 as.numeric(Control_L1_count/Control_L1_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H1_normalised = 
                                 as.numeric(Control_H1_count/Control_H1_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L2_normalised = 
                                 as.numeric(Control_L2_count/Control_L2_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H2_normalised = 
                                 as.numeric(Control_H2_count/Control_H2_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L3_normalised = 
                                 as.numeric(Control_L3_count/Control_L3_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H3_normalised = 
                                 as.numeric(Control_H3_count/Control_H3_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L4_normalised = 
                                 as.numeric(Control_L4_count/Control_L4_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H4_normalised = 
                                 as.numeric(Control_H4_count/Control_H4_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L5_normalised = 
                                 as.numeric(Control_L5_count/Control_L5_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H5_normalised = 
                                 as.numeric(Control_H5_count/Control_H5_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L6_normalised = 
                                 as.numeric(Control_L6_count/Control_L6_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H6_normalised = 
                                 as.numeric(Control_H6_count/Control_H6_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L7_normalised = 
                                 as.numeric(Control_L7_count/Control_L7_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H7_normalised = 
                                 as.numeric(Control_H7_count/Control_H7_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L8_normalised = 
                                 as.numeric(Control_L8_count/Control_L8_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H8_normalised = 
                                 as.numeric(Control_H8_count/Control_H8_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_L9_normalised = 
                                 as.numeric(Control_L9_count/Control_L9_base))
Sample_1_vs_Control <-  mutate(Sample_1_vs_Control, Control_H9_normalised = 
                                 as.numeric(Control_H9_count/Control_H9_base))

Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L1_normalised = 
                                 as.numeric(Control_L1_count/Control_L1_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H1_normalised = 
                                 as.numeric(Control_H1_count/Control_H1_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L2_normalised = 
                                 as.numeric(Control_L2_count/Control_L2_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H2_normalised = 
                                 as.numeric(Control_H2_count/Control_H2_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L3_normalised = 
                                 as.numeric(Control_L3_count/Control_L3_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H3_normalised = 
                                 as.numeric(Control_H3_count/Control_H3_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L4_normalised = 
                                 as.numeric(Control_L4_count/Control_L4_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H4_normalised = 
                                 as.numeric(Control_H4_count/Control_H4_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L5_normalised = 
                                 as.numeric(Control_L5_count/Control_L5_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H5_normalised = 
                                 as.numeric(Control_H5_count/Control_H5_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L6_normalised = 
                                 as.numeric(Control_L6_count/Control_L6_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H6_normalised = 
                                 as.numeric(Control_H6_count/Control_H6_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L7_normalised = 
                                 as.numeric(Control_L7_count/Control_L7_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H7_normalised = 
                                 as.numeric(Control_H7_count/Control_H7_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L8_normalised = 
                                 as.numeric(Control_L8_count/Control_L8_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H8_normalised = 
                                 as.numeric(Control_H8_count/Control_H8_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_L9_normalised = 
                                 as.numeric(Control_L9_count/Control_L9_base))
Sample_2_vs_Control <-  mutate(Sample_2_vs_Control, Control_H9_normalised = 
                                 as.numeric(Control_H9_count/Control_H9_base))

Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L1_normalised = 
                                 as.numeric(Control_L1_count/Control_L1_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H1_normalised = 
                                 as.numeric(Control_H1_count/Control_H1_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L2_normalised = 
                                 as.numeric(Control_L2_count/Control_L2_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H2_normalised = 
                                 as.numeric(Control_H2_count/Control_H2_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L3_normalised = 
                                 as.numeric(Control_L3_count/Control_L3_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H3_normalised = 
                                 as.numeric(Control_H3_count/Control_H3_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L4_normalised = 
                                 as.numeric(Control_L4_count/Control_L4_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H4_normalised = 
                                 as.numeric(Control_H4_count/Control_H4_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L5_normalised = 
                                 as.numeric(Control_L5_count/Control_L5_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H5_normalised = 
                                 as.numeric(Control_H5_count/Control_H5_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L6_normalised = 
                                 as.numeric(Control_L6_count/Control_L6_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H6_normalised = 
                                 as.numeric(Control_H6_count/Control_H6_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L7_normalised = 
                                 as.numeric(Control_L7_count/Control_L7_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H7_normalised = 
                                 as.numeric(Control_H7_count/Control_H7_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L8_normalised = 
                                 as.numeric(Control_L8_count/Control_L8_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H8_normalised = 
                                 as.numeric(Control_H8_count/Control_H8_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_L9_normalised = 
                                 as.numeric(Control_L9_count/Control_L9_base))
Sample_3_vs_Control <-  mutate(Sample_3_vs_Control, Control_H9_normalised = 
                                 as.numeric(Control_H9_count/Control_H9_base))

# combine low and high salt normalised values

Sample_1_vs_Control$Sample_1_normalised <- as.numeric(Sample_1_vs_Control$Sample_L1_normalised + 
                                                        Sample_1_vs_Control$Sample_H1_normalised)
Sample_2_vs_Control$Sample_2_normalised <- as.numeric(Sample_2_vs_Control$Sample_L2_normalised + 
                                                        Sample_2_vs_Control$Sample_H2_normalised)
Sample_3_vs_Control$Sample_3_normalised <- as.numeric(Sample_3_vs_Control$Sample_L3_normalised + 
                                                        Sample_3_vs_Control$Sample_H3_normalised)

Sample_1_vs_Control$Control_1_normalised <- as.numeric(Sample_1_vs_Control$Control_L1_normalised + 
                                                         Sample_1_vs_Control$Control_H1_normalised)
Sample_1_vs_Control$Control_2_normalised <- as.numeric(Sample_1_vs_Control$Control_L2_normalised + 
                                                         Sample_1_vs_Control$Control_H2_normalised)
Sample_1_vs_Control$Control_3_normalised <- as.numeric(Sample_1_vs_Control$Control_L3_normalised + 
                                                         Sample_1_vs_Control$Control_H3_normalised)
Sample_1_vs_Control$Control_4_normalised <- as.numeric(Sample_1_vs_Control$Control_L4_normalised + 
                                                         Sample_1_vs_Control$Control_H4_normalised)
Sample_1_vs_Control$Control_5_normalised <- as.numeric(Sample_1_vs_Control$Control_L5_normalised + 
                                                         Sample_1_vs_Control$Control_H5_normalised)
Sample_1_vs_Control$Control_6_normalised <- as.numeric(Sample_1_vs_Control$Control_L6_normalised + 
                                                         Sample_1_vs_Control$Control_H6_normalised)
Sample_1_vs_Control$Control_7_normalised <- as.numeric(Sample_1_vs_Control$Control_L7_normalised + 
                                                         Sample_1_vs_Control$Control_H7_normalised)
Sample_1_vs_Control$Control_8_normalised <- as.numeric(Sample_1_vs_Control$Control_L8_normalised + 
                                                         Sample_1_vs_Control$Control_H8_normalised)
Sample_1_vs_Control$Control_9_normalised <- as.numeric(Sample_1_vs_Control$Control_L9_normalised + 
                                                         Sample_1_vs_Control$Control_H9_normalised)

Sample_2_vs_Control$Control_1_normalised <- as.numeric(Sample_2_vs_Control$Control_L1_normalised + 
                                                         Sample_2_vs_Control$Control_H1_normalised)
Sample_2_vs_Control$Control_2_normalised <- as.numeric(Sample_2_vs_Control$Control_L2_normalised + 
                                                         Sample_2_vs_Control$Control_H2_normalised)
Sample_2_vs_Control$Control_3_normalised <- as.numeric(Sample_2_vs_Control$Control_L3_normalised + 
                                                         Sample_2_vs_Control$Control_H3_normalised)
Sample_2_vs_Control$Control_4_normalised <- as.numeric(Sample_2_vs_Control$Control_L4_normalised + 
                                                         Sample_2_vs_Control$Control_H4_normalised)
Sample_2_vs_Control$Control_5_normalised <- as.numeric(Sample_2_vs_Control$Control_L5_normalised + 
                                                         Sample_2_vs_Control$Control_H5_normalised)
Sample_2_vs_Control$Control_6_normalised <- as.numeric(Sample_2_vs_Control$Control_L6_normalised + 
                                                         Sample_2_vs_Control$Control_H6_normalised)
Sample_2_vs_Control$Control_7_normalised <- as.numeric(Sample_2_vs_Control$Control_L7_normalised + 
                                                         Sample_2_vs_Control$Control_H7_normalised)
Sample_2_vs_Control$Control_8_normalised <- as.numeric(Sample_2_vs_Control$Control_L8_normalised + 
                                                         Sample_2_vs_Control$Control_H8_normalised)
Sample_2_vs_Control$Control_9_normalised <- as.numeric(Sample_2_vs_Control$Control_L9_normalised + 
                                                         Sample_2_vs_Control$Control_H9_normalised)

Sample_3_vs_Control$Control_1_normalised <- as.numeric(Sample_3_vs_Control$Control_L1_normalised + 
                                                         Sample_3_vs_Control$Control_H1_normalised)
Sample_3_vs_Control$Control_2_normalised <- as.numeric(Sample_3_vs_Control$Control_L2_normalised + 
                                                         Sample_3_vs_Control$Control_H2_normalised)
Sample_3_vs_Control$Control_3_normalised <- as.numeric(Sample_3_vs_Control$Control_L3_normalised + 
                                                         Sample_3_vs_Control$Control_H3_normalised)
Sample_3_vs_Control$Control_4_normalised <- as.numeric(Sample_3_vs_Control$Control_L4_normalised + 
                                                         Sample_3_vs_Control$Control_H4_normalised)
Sample_3_vs_Control$Control_5_normalised <- as.numeric(Sample_3_vs_Control$Control_L5_normalised + 
                                                         Sample_3_vs_Control$Control_H5_normalised)
Sample_3_vs_Control$Control_6_normalised <- as.numeric(Sample_3_vs_Control$Control_L6_normalised + 
                                                         Sample_3_vs_Control$Control_H6_normalised)
Sample_3_vs_Control$Control_7_normalised <- as.numeric(Sample_3_vs_Control$Control_L7_normalised + 
                                                         Sample_3_vs_Control$Control_H7_normalised)
Sample_3_vs_Control$Control_8_normalised <- as.numeric(Sample_3_vs_Control$Control_L8_normalised + 
                                                         Sample_3_vs_Control$Control_H8_normalised)
Sample_3_vs_Control$Control_9_normalised <- as.numeric(Sample_3_vs_Control$Control_L9_normalised + 
                                                         Sample_3_vs_Control$Control_H9_normalised)

# remove columns with separate low and high salt normalisation, just keeping combined normalised columns

Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H1_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H2_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H3_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H4_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H5_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H6_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H7_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H8_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("H9_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L1_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L2_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L3_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L4_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L5_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L6_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L7_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L8_normalised"))
Sample_1_vs_Control <- Sample_1_vs_Control %>% select(-ends_with("L9_normalised"))

Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H1_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H2_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H3_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H4_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H5_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H6_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H7_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H8_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("H9_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L1_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L2_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L3_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L4_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L5_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L6_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L7_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L8_normalised"))
Sample_2_vs_Control <- Sample_2_vs_Control %>% select(-ends_with("L9_normalised"))

Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H1_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H2_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H3_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H4_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H5_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H6_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H7_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H8_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("H9_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L1_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L2_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L3_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L4_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L5_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L6_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L7_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L8_normalised"))
Sample_3_vs_Control <- Sample_3_vs_Control %>% select(-ends_with("L9_normalised"))


### Filter so the sample has more than a normalised count of 5e-09 (row numbers containing normalised count for Sample only)
# In this example there are 9 controls (H and L salt) and 1 sample (H+L) - so "K-mer" column, 20 count columns, and then the Sample normalised count column is number 22

Sample_1_vs_Control <- Sample_1_vs_Control %>% filter(.[[22]]>5e-09)
Sample_2_vs_Control <- Sample_2_vs_Control %>% filter(.[[22]]>5e-09)
Sample_3_vs_Control <- Sample_3_vs_Control %>% filter(.[[22]]>5e-09)

### Take average of normalised values for control (in this example, there are 3 controls, but if more divide by the number of control samples)

Sample_1_vs_Control <- mutate(Sample_1_vs_Control, average_Control = as.numeric(Control_1_normalised + 
Control_2_normalised + Control_3_normalised + Control_4_normalised + Control_5_normalised + Control_6_normalised + 
Control_7_normalised + Control_8_normalised + Control_9_normalised)/9)

Sample_2_vs_Control <- mutate(Sample_2_vs_Control, average_Control = as.numeric(Control_1_normalised + 
Control_2_normalised + Control_3_normalised + Control_4_normalised + Control_5_normalised + Control_6_normalised + 
Control_7_normalised + Control_8_normalised + Control_9_normalised)/9)

Sample_3_vs_Control <- mutate(Sample_3_vs_Control, average_Control = as.numeric(Control_1_normalised + 
Control_2_normalised + Control_3_normalised + Control_4_normalised + Control_5_normalised + Control_6_normalised + 
Control_7_normalised + Control_8_normalised + Control_9_normalised)/9)

### calculate fold change

Sample_1_vs_Control <- mutate(Sample_1_vs_Control, Sample_1_v_Control = as.numeric(Sample_1_normalised/average_Control))
Sample_2_vs_Control <- mutate(Sample_2_vs_Control, Sample_2_v_Control = as.numeric(Sample_2_normalised/average_Control))
Sample_3_vs_Control <- mutate(Sample_3_vs_Control, Sample_3_v_Control = as.numeric(Sample_3_normalised/average_Control))


### Filter by fold change - can choose whatever you like here


Sample_1_vs_Control_fc2 <- filter(Sample_1_vs_Control, Sample_1_v_Control >2)
Sample_2_vs_Control_fc2 <- filter(Sample_2_vs_Control, Sample_2_v_Control >2)
Sample_3_vs_Control_fc2 <- filter(Sample_3_vs_Control, Sample_3_v_Control >2)



write.table(Sample_1_vs_Control_fc2, file=paste("kmer_centromere/enriched_kmers/",Sample_name,"_1_vs_Control_fc2.txt", sep=""),sep =" ", row.names = F, 
col.names = T, quote = F)
write.table(Sample_2_vs_Control_fc2, file=paste("kmer_centromere/enriched_kmers/",Sample_name,"_2_vs_Control_fc2.txt", sep=""),sep =" ", row.names = F, 
col.names = T, quote = F)
write.table(Sample_3_vs_Control_fc2, file=paste("kmer_centromere/enriched_kmers/",Sample_name,"_3_vs_Control_fc2.txt", sep=""),sep =" ", row.names = F, 
col.names = T, quote = F)
