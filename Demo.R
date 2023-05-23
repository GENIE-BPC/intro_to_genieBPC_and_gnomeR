#devtools::install_github("MSKCC-Epi-Bio/gnomeR") 

library(genieBPC)
library(tidyverse)
library(gnomeR)
library(gtsummary)

###############################################################################
#                         Demo with NSCLC 2.0-public
###############################################################################
# log into Synapse
genieBPC::set_synapse_credentials()

# pull data from Synapse into R
nsclc_synapse_data <- pull_data_synapse(cohort = "NSCLC",
                                  version = "v2.0-public")

# simple example of create_analytic_cohort
# first index cancer, all stages, institutions and histologies
nsclc_cohort <- create_analytic_cohort(data_synapse = nsclc_synapse_data$NSCLC_v2.0,
                                       return_summary = TRUE)

# review list returned
names(nsclc_cohort)
# summary(nsclc_cohort)

# review summary tables returned
nsclc_cohort$tbl_overall_summary
nsclc_cohort$tbl_cohort
nsclc_cohort$tbl_drugs
nsclc_cohort$tbl_ngs

# create analytic cohort for case study
nsclc_cohort <- create_analytic_cohort(
  data_synapse = nsclc_synapse_data$NSCLC_v2.0,
  stage_dx = c("Stage IV"),
  histology = "Adenocarcinoma",
  regimen_drugs = c("Carboplatin, Pemetrexed Disodium",
                    "Cisplatin, Pemetrexed Disodium",
                    "Bevacizumab, Carboplatin, Pemetrexed Disodium",
                    "Bevacizumab, Cisplatin, Pemetrexed Disodium"),
  regimen_type = "Exact",
  regimen_order = 1,
  regimen_order_type = "within cancer",
  return_summary = TRUE
)

# review summary tables
nsclc_cohort$tbl_overall_summary
nsclc_cohort$tbl_cohort
nsclc_cohort$tbl_drugs
nsclc_cohort$tbl_ngs

# create sunburst plot
nsclc_sunburst <- drug_regimen_sunburst(data_synapse = nsclc_synapse_data$NSCLC_v2.0,
                                        data_cohort = nsclc_cohort)

nsclc_sunburst$sunburst_plot


# Genomic Analysis --------------------------------------------------------


# select one sample per patient
nrow(nsclc_cohort$cohort_ngs)

nsclc_samp <- select_unique_ngs(
  data_cohort = nsclc_cohort$cohort_ngs,
  oncotree_code = "LUAD",
  sample_type = "Metastasis",
  min_max_time = "max"
)

nrow(nsclc_samp) 

# get mutation, cna and fusions data
mutations <- nsclc_synapse_data$NSCLC_v2.0$mutations_extended
cna <- nsclc_synapse_data$NSCLC_v2.0$cna
fusions <- nsclc_synapse_data$NSCLC_v2.0$fusions

# Reformat CNA and fusions -----

# reformat fusions and cna to get in standard format
reformat_fusions <- gnomeR::reformat_fusion(fusions)
nrow(reformat_fusions)

reformat_cna <- gnomeR::pivot_cna_longer(cna)
nrow(reformat_cna)


nsclc_panels <- data.frame(
  sample_id = nsclc_samp$cpt_genie_sample_id,
  panel_id = nsclc_samp$cpt_seq_assay_id) %>% 
  mutate(panel_id = ifelse(!is.na(panel_id),
                           panel_id, "no"))

# Create your gene binary data -----
gene_binary <- gnomeR::create_gene_binary(
  mutation = mutations,
  cna = reformat_cna,
  fusion = reformat_fusions,
  samples = nsclc_samp$cpt_genie_sample_id,
  specify_panel = nsclc_panels, 
  recode_aliases = "impact")

# Join clinical data to mutation data -----
# get patient IDs and sample IDs
patient_index <-  nsclc_cohort$cohort_ngs %>%
  select(record_id, cpt_genie_sample_id)

# Join sex data to patient ID index
select_clinical <- nsclc_cohort$cohort_pt_char %>%
  select(record_id, naaccr_sex_code) %>%
  left_join(patient_index) 

# Join all to gene binary data
gene_binary <- gene_binary %>% 
  left_join(select_clinical,
            by = c("sample_id"= "cpt_genie_sample_id")) %>%
  select(-record_id)

gene_binary <- gene_binary %>% 
  select(sample_id, naaccr_sex_code, everything())

gene_binary %>%
  select(naaccr_sex_code) %>%
  tbl_summary() 

# Subset by gene frequency -----
ncol(gene_binary)

nscl_subset <- gene_binary %>%
  subset_by_frequency(t = .4, other_vars = naaccr_sex_code)

ncol(nscl_subset)

# General summary table  -----
nscl_subset %>% 
  select(-naaccr_sex_code) %>%
  tbl_genomic() %>%
  bold_labels() 


# Pathways -----
path_df <- gene_binary %>%
  select(-naaccr_sex_code) %>%
  add_pathways() 

path_df %>%
  select(contains("pathway")) %>%
  tbl_summary() %>% 
  bold_labels() 


# Summarize Frequencies By Clinical Variable -----
tbl_gene <- gene_binary %>%
  subset_by_frequency(
    t = .4,
    other_vars = naaccr_sex_code) %>%
  tbl_genomic(by = naaccr_sex_code) %>%
  bold_labels() %>%
  add_p() %>%
  add_q()

tbl_gene


###############################################################################
#                         Demo with NSCLC test data
###############################################################################
# simple example of create_analytic_cohort
nsclc_cohort <- create_analytic_cohort(data_synapse = nsclc_test_data,
                                       return_summary = TRUE)

# review summary tables
nsclc_cohort$tbl_overall_summary
nsclc_cohort$tbl_cohort
nsclc_cohort$tbl_drugs
nsclc_cohort$tbl_ngs

# create analytic cohort for case study
nsclc_cohort <- create_analytic_cohort(
  data_synapse = genieBPC::nsclc_test_data,
  stage_dx = c("Stage IV"),
  histology = "Adenocarcinoma", 
  regimen_drugs = c("Carboplatin, Pemetrexed Disodium",
                    "Cisplatin, Pemetrexed Disodium",
                    "Bevacizumab, Carboplatin, Pemetrexed Disodium",
                    "Bevacizumab, Cisplatin, Pemetrexed Disodium"),
  regimen_type = "Exact",
  regimen_order = 1,
  regimen_order_type = "within cancer",
  return_summary = TRUE
)

# review summary tables
nsclc_cohort$tbl_overall_summary
nsclc_cohort$tbl_cohort
nsclc_cohort$tbl_drugs
nsclc_cohort$tbl_ngs

# create sunburst plot
nsclc_sunburst <- drug_regimen_sunburst(data_synapse = nsclc_test_data,
                                        data_cohort = nsclc_cohort)

nsclc_sunburst$sunburst_plot
