#!/usr/local/bin/Rscript

release <- '2025.12'

otp_assoc_data <- list()
otp_assoc_data[['release']] <- release
otp_assoc_data[['overall']] <- data.frame()
otp_assoc_data[['datasource']] <- data.frame()
otp_assoc_data[['datatype']] <- data.frame()
otp_assoc_data[['disease']] <- data.frame()
otp_assoc_data[['target']] <- readRDS(
  file = file.path(
    "output",
    paste0("opentargets_target_",release,".rds")))


####---- ASSOCIATIONS - OVERALL ----####

basepath <- file.path(here::here(), "data", 
                      release, "association_overall_direct")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))
m <- 0
for(parquet_fname in parquet_files){
  cat(paste0('Chunk - association_overall - ',m))
  cat('\n')
  
  # Read the parquet file and iterate through each row
  association_data <- as.data.frame(
    arrow::read_parquet(parquet_fname)) |>
    dplyr::rename(
      disease_id = diseaseId,
      target_id = targetId,
      evidence_count = evidenceCount
    )
  otp_assoc_data[['overall']] <- 
    dplyr::bind_rows(
      otp_assoc_data[['overall']], association_data
    )
  m <- m + 1
}

####---- ASSOCIATIONS - SOURCE ----####
basepath <- file.path(here::here(), "data", 
                      release, "association_by_datasource_direct")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))

m <- 0
OT_datasource_assocs <- data.frame()
for(parquet_fname in parquet_files){

  cat(paste0('Chunk - association_by_datasource - ',m))
  cat('\n')
  
  datasource_assoc_df <- as.data.frame(
    arrow::read_parquet(parquet_fname)) |>
    dplyr::rename(
      datatype_id = datatypeId,
      datasource_id = datasourceId,
      disease_id = diseaseId,
      target_id = targetId,
      evidence_count = evidenceCount
    )
  
  otp_assoc_data[['datasource']] <- 
    dplyr::bind_rows(
      otp_assoc_data[['datasource']], datasource_assoc_df
    )
  m <- m + 1
}

####---- ASSOCIATIONS - TYPE ----####
basepath <- file.path(here::here(), "data", 
                      release, "association_by_datatype_direct")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))

m <- 0
for(parquet_fname in parquet_files){

  cat(paste0('Chunk - association_by_datatype - ',m))
  cat('\n')
  
  assoc_dtype_df <- as.data.frame(
    arrow::read_parquet(parquet_fname)) |>
    dplyr::rename(
      datatype_id = datatypeId,
      disease_id = diseaseId,
      target_id = targetId,
      evidence_count = evidenceCount
    )
  
  otp_assoc_data[['datatype']] <- 
    dplyr::bind_rows(
      otp_assoc_data[['datatype']], assoc_dtype_df
    )
  m <- m + 1
}



####---- DISEASE ----####

basepath <- file.path(here::here(), "data", 
                      release, "disease")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))
m <- 0
for(parquet_fname in parquet_files){
  cat(paste0('Chunk - disease - ',m))
  cat('\n')
  
  disease_assoc_df <- as.data.frame(
    arrow::read_parquet(parquet_fname)) |>
    dplyr::rename(
      disease_id = id) |>
    dplyr::select(
      disease_id, name, description
    )
  
  otp_assoc_data[['disease']] <- 
    dplyr::bind_rows(
      otp_assoc_data[['disease']], disease_assoc_df
    )
  
  m <- m + 1
}

## Aggregate datatype and datasource support per disease-target pair

otp_assoc_data[['datatype']] <- as.data.frame(
  otp_assoc_data[['datatype']] |>
    dplyr::mutate(score = round(score, digits = 12)) |>
    dplyr::mutate(datatype_support = paste0(
      datatype_id,"|",evidence_count,"|",score)) |>
    dplyr::group_by(disease_id, target_id) |>
    dplyr::summarise(datatype_items = paste(datatype_support, collapse=","),
                     .groups = "drop")
)

otp_assoc_data[['datasource']] <- as.data.frame(
  otp_assoc_data[['datasource']] |>
    dplyr::mutate(score = round(score, digits = 12)) |>
    dplyr::mutate(datasource_support = paste0(
      datasource_id,"|",evidence_count,"|",score)) |>
    dplyr::group_by(disease_id, target_id) |>
    dplyr::summarise(
      datasource_items = paste(datasource_support, collapse=","),
      .groups = "drop")
)

OT_association_all <- as.data.frame(
  otp_assoc_data[['overall']] |>
    dplyr::mutate(score = round(score, digits = 12)) |>
    dplyr::left_join(
      otp_assoc_data[['datatype']], 
      by = c("disease_id","target_id"),
      relationship = "many-to-many") |>
    dplyr::left_join(
      otp_assoc_data[['datasource']], 
      by = c("disease_id","target_id"),
      relationship = "many-to-many") |>
    dplyr::rename(
      target_ensembl_gene_id = target_id) |>
    dplyr::left_join(
      dplyr::select(
        otp_assoc_data[['target']], 
        target_symbol, 
        target_ensembl_gene_id),
      by = c("target_ensembl_gene_id"),
      relationship = "many-to-many") |>
    dplyr::left_join(
      dplyr::select(
        otp_assoc_data[['disease']], 
        disease_id, name),
      by = "disease_id", 
      relationship = "many-to-many") |>
    dplyr::rename(disease_label = name)
)

saveRDS(OT_association_all, 
        file = file.path(
          "output",
          paste0("opentargets_association_direct_",
                 release,".rds"))
)
        
OT_associations_multiple_types <- OT_association_all |>
  dplyr::filter(!stringr::str_detect(
    tolower(disease_label),"^(neoplasm|cancer)$")) |>
  dplyr::filter(stringr::str_detect(datatype_items,","))

OT_associations_single_types <- OT_association_all |>
  dplyr::filter(!stringr::str_detect(
    tolower(disease_label),"^(neoplasm|cancer)$")) |>
  dplyr::filter(!stringr::str_detect(datatype_items,",")) |>
  tidyr::separate(datatype_items, 
                  c("type","type_count","type_score"),sep = "\\|", remove = F) |>
  dplyr::filter(type != "rna_expression") |>
  dplyr::filter(type != "known_drug") |>
  dplyr::filter(type != "literature" |
                  (type == "literature" & evidence_count > 4)) |>
  dplyr::filter(type != "somatic_mutation" |
                  (type == "somatic_mutation" & 
                     ("cancer_gene_census" %in% datasource_items |
                     stringr::str_detect(datasource_items,",")))) |>
  dplyr::select(-c("type_count","type_score","type"))
  
OT_association_hc <- OT_associations_single_types |>
  dplyr::bind_rows(OT_associations_multiple_types) |>
  dplyr::arrange(dplyr::desc(score))

saveRDS(OT_association_hc, 
        file = file.path(
          "output",
          paste0("opentargets_association_direct_HC_",
                 release,".rds"))
)

