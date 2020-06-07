library(dplyr)
library(tidyr)
library(ggplot2)
library(biomaRt)

clinical_fname = "/Volumes/Elements/MYCNAmplicon/Data/TCGA_CopyNumberData/TCGA-UCS/clinical.cart.2020-03-14/clinical.tsv"
sample_sheet_fname = "/Volumes/Elements/MYCNAmplicon/Data/TCGA_CopyNumberData/TCGA-UCS/gdc_sample_sheet.2020-03-14.tsv"
data_dir = "/Volumes/Elements/MYCNAmplicon/Data/TCGA_CopyNumberData/TCGA-UCS/gdc_download_20200314_152554.192163/"

clinical = 
  read.table(clinical_fname, 
             header=T,
             sep="\t") %>% 
  as_tibble() %>%
  dplyr::select(case_id, submitter_id, days_to_death, days_to_last_follow_up, vital_status, primary_diagnosis, ajcc_pathologic_stage) %>%
  mutate(days_to_death = as.numeric(ifelse(as.character(days_to_death) == "--", "NA", as.character(days_to_death))),
         days_to_last_follow_up = as.numeric(ifelse(as.character(days_to_last_follow_up) == "--", "NA", as.character(days_to_last_follow_up))),
         isDead = ifelse(vital_status == "Dead", T, ifelse(vital_status == "Alive", F, NA))) %>%
  mutate(days_to_last_follow_up = ifelse(is.na(days_to_last_follow_up), days_to_death, days_to_last_follow_up),
         days_to_last_follow_up = ifelse(!is.na(days_to_death) & days_to_death>days_to_last_follow_up, days_to_death, days_to_last_follow_up)) %>% 
  distinct() %>%
  dplyr::rename(TCGACase = submitter_id)

# clinical %>% View
table(clinical$vital_status)
nrow(clinical) # 371 for HCC, 515 for LUAD
length(unique(clinical$case_id)) # 371 for HCC, 515 for LUAD
length(unique(clinical$TCGACase)) # 371 for HCC, 515 for LUAD
table(clinical$primary_diagnosis)
table(clinical$ajcc_pathologic_stage)

sample_sheet = 
  read.table(sample_sheet_fname,
             header=T,
             sep="\t") %>%
  as_tibble() %>% 
  dplyr::rename(TCGACase=Case.ID, TCGASample=Sample.ID) %>% 
  filter(Sample.Type == "Primary Tumor")

# Some cases have several primary tumor samples; randomly select one sample
# to get one sample per case 
set.seed(42)
sample_sheet = 
  sample_sheet %>% 
  group_by(TCGACase) %>%
  dplyr::sample_n(1) %>%
  ungroup()

#sample_sheet %>% View
nrow(sample_sheet)  
length(unique(sample_sheet$TCGACase)) 
length(unique(sample_sheet$TCGASample)) 

sample_sheet = sample_sheet %>% 
  inner_join(clinical, by="TCGACase")

nrow(sample_sheet) 

sample_sheet = 
  sample_sheet %>%
  filter(!is.na(File.ID), !is.na(File.Name)) %>% 
  mutate(File = paste0(data_dir, 
                       File.ID,
                       "/",
                       File.Name))

nrow(sample_sheet) 
#sample_sheet = sample_sheet[sample(nrow(sample_sheet),60,replace = F),] # subsampling for debugging
sampleTable = data.frame(sampleName = sample_sheet$TCGACase,
                         fileName = sample_sheet$File)

for (i in 1:nrows(sampleTable)){
  
}


