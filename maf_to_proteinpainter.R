#!/usr/bin/env Rscript


# Libraries ---------------------------------------------------------------

library(data.table)
library(feather)
library(tidyverse)
library(stringr)


# Argument parsing --------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Inputs
file_paths <- list()
file_paths$maf    <- args[1]
file_paths$genes  <- args[2]
file_paths$output <- args[3]

disease <- args[4]
if (is.null(disease)) disease <- "Lymphoma"



# Loading -----------------------------------------------------------------

maf <- fread(file_paths$maf)
genes <- read_lines(file_paths$genes)

nonsyn <- c("Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del",
            "Frame_Shift_Ins", "Nonstop_Mutation", "Translation_Start_Site",
            "In_Frame_Ins", "In_Frame_Del", "Missense_Mutation")


# Tidying -----------------------------------------------------------------

pp <-
  maf %>%
  .[Variant_Classification %in% nonsyn] %>%
  .[Hugo_Symbol %in% genes] %>%
  .[, .(
    disease = disease,
    sampletype = "diagnosis",
    annovar_type = ifelse(Variant_Type == "SNP", "snv", "complex"),
    annovar_gene = Hugo_Symbol,
    annovar_class = recode(Variant_Classification,
                           Nonsense_Mutation = "nonsense",
                           Missense_Mutation = "missense",
                           Splice_Site = "splice",
                           Frame_Shift_Del = "frameshift",
                           Frame_Shift_Ins = "frameshift",
                           In_Frame_Del = "proteinDel",
                           In_Frame_Ins = "proteinIns",
                           Translation_Start_Site = "missense"),
    annovar_aachange = sub("^p[.]", "", HGVSp_Short),
    annovar_isoform = sub("[.][0-9]+$", "", str_split_fixed(RefSeq, ",", n = Inf)[,1]),
    TARGET_CASE_ID = Tumor_Sample_Barcode,
    Trio = "",
    Hugo_Symbol = Hugo_Symbol,
    Variant_Classification = Variant_Classification,
    VariantType = Variant_Type,
    dbSNP_RS = dbSNP_RS,
    Mutation_Status = Mutation_Status,
    PFAM_DOMAIN = "",
    Somatic_Score = "",
    Somatic_Rank = "",
    Somatic_quality = "",
    Tumor_ReadCount_Alt = 0,
    Tumor_ReadCount_Ref = 0,
    Tumor_ReadCount_Total = 0,
    Normal_ReadCount_Alt = 0,
    Normal_ReadCount_Ref = 0,
    Normal_ReadCount_Total = 0,
    Cosmic = "",
    Cosmic_Gene = "",
    Reference_Allele = Reference_Allele,
    TumorSeq_Allele1 = Tumor_Seq_Allele1,
    TumorSeq_Allele2 = Tumor_Seq_Allele2,
    Match_Norm_Seq_Allele1 = Match_Norm_Seq_Allele1,
    Match_Norm_Seq_Allele2 = Match_Norm_Seq_Allele2,
    Tumor_Sample_Barcode = Tumor_Sample_Barcode,
    Match_Normal_Sample_Barcode = Matched_Norm_Sample_Barcode,
    Entrez_Gene_Id = Entrez_Gene_Id,
    Chromosome = Chromosome,
    Start_position = Start_Position,
    End_position = End_Position,
    miRNA = "",
    Verification_Status = "",
    Verification_Method = "",
    FET_Score = "",
    TumorRefCount_VS = "",
    TumorVarCount_VS = "",
    TumorTotalCount_VS = "",
    NormalRefCount_VS = "",
    NormalVarCount_VS = "",
    NormalTotalCount_VS = ""
  )]

# Manual tweaking
pp[Hugo_Symbol == "HIST1H2BC", annovar_isoform := "NM_003526"]
pp[Hugo_Symbol == "TP53" & annovar_class == "splice" &
     (is.na(annovar_aachange) | annovar_aachange == "") &
     Start_position > 7685855,
   annovar_aachange := "X0_splice"]
pp[Hugo_Symbol == "FOXO1" & VariantType == "DEL" &
     (is.na(annovar_aachange) | annovar_aachange == "") &
     Start_position > 40666200,
   annovar_aachange := "M1*"]

# Add alternative columns names
pp[, `:=`(
  gene = annovar_gene,
  refseq = annovar_isoform,
  chromosome = Chromosome,
  start = Start_position,
  aachange = annovar_aachange,
  class = annovar_class
)]


# Outputting --------------------------------------------------------------

fwrite(pp, file_paths$output, sep = "\t")
