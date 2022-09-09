release <- '2022.06'

####---- ASSOCIATIONS - OVERALL ----####

basepath <- file.path(here::here(), "data", 
                      release, "associationByOverallDirect")
json_files <- sort(
  list.files(basepath, pattern = ".json", 
             all.files = T, full.names = T))

OT_association <- data.frame()

m <- 0
for(json_chunk_lines in json_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  con <- file(json_chunk_lines, "r")
  association_data <- 
    jsonlite::stream_in(con, flatten = F, simplifyVector = F)
  i <- 1
  
  while(i <= length(association_data)){
    association_item <- association_data[[i]]
    
    df <- data.frame(
      'disease_id' =  association_item$diseaseId,
      'target_id' = association_item$targetId,
      'score' = association_item$score,
      'evidence_count' = association_item$evidenceCount,
      stringsAsFactors = F)

    i <- i + 1
    OT_association <- OT_association |>
      dplyr::bind_rows(df)
  }
  
  close.connection(con)
  saveRDS(OT_association, 
          file = paste0("output/association/OT_overall_association_",
                        release,"_",m,".rds"))
  OT_association <- data.frame()
  m <- m + 1
}

####---- ASSOCIATIONS - SOURCE ----####

basepath <- file.path(here::here(), "data", 
                      release, "associationByDatasourceDirect")
json_files <- sort(
  list.files(basepath, pattern = ".json", 
             all.files = T, full.names = T))

OT_datasource_association <- data.frame()

m <- 0
for(json_chunk_lines in json_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  con <- file(json_chunk_lines, "r")
  association_datasource <- 
    jsonlite::stream_in(con, flatten = F, simplifyVector = F)
  i <- 1
  
  while(i <= length(association_datasource)){
    association_datasource_item <- association_datasource[[i]]
    
    df <- data.frame(
      'datatype_id' =  association_datasource_item$datatypeId,
      'datasource_id' = association_datasource_item$datasourceId,
      'disease_id' = association_datasource_item$diseaseId,
      'target_id' = association_datasource_item$targetId,
      'score' = association_datasource_item$score,
      'evidence_count' = association_datasource_item$evidenceCount,
      stringsAsFactors = F)
    
    i <- i + 1
    
    OT_datasource_association <- OT_datasource_association |>
      dplyr::bind_rows(df)
    
  }
  close.connection(con)
  saveRDS(OT_datasource_association, 
          file = paste0("output/association/OT_datasource_association_",
                        release,"_",m,".rds"))
  OT_datasource_association <- data.frame()
  m <- m + 1
}

####---- ASSOCIATIONS - TYPE ----####

basepath <- file.path(here::here(), "data", 
                      release, "associationByDatatypeDirect")
json_files <- sort(
  list.files(basepath, pattern = ".json", 
             all.files = T, full.names = T))

OT_datatype_association <- data.frame()

m <- 0
for(json_chunk_lines in json_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  con <- file(json_chunk_lines, "r")
  association_datatype <- 
    jsonlite::stream_in(con, flatten = F, simplifyVector = F)
  i <- 1
  
  while(i <= length(association_datatype)){
    association_datatype_item <- association_datatype[[i]]
    
    df <- data.frame(
      'datatype_id' =  association_datatype_item$datatypeId,
      'disease_id' = association_datatype_item$diseaseId,
      'target_id' = association_datatype_item$targetId,
      'score' = association_datatype_item$score,
      'evidence_count' = association_datatype_item$evidenceCount,
      stringsAsFactors = F)
    
    i <- i + 1
    
    OT_datatype_association <- OT_datatype_association |>
      dplyr::bind_rows(df)
    
  }
  close.connection(con)
  saveRDS(OT_datatype_association, 
          file = paste0("output/association/OT_datatype_association_",
                        release,"_",m,".rds"))
  OT_datatype_association <- data.frame()
  m <- m + 1
}



####---- DISEASE ----####

basepath <- file.path(here::here(), "data", 
                      release, "diseases")
json_files <- sort(
  list.files(basepath, pattern = ".json", 
             all.files = T, full.names = T))

OT_diseases <- data.frame()

m <- 0
for(json_chunk_lines in json_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  con <- file(json_chunk_lines, "r")
  disease_assoc <- 
    jsonlite::stream_in(con, flatten = F, simplifyVector = F)
  i <- 1
  
  while(i <= length(disease_assoc)){
    disease_item <- disease_assoc[[i]]
    
    df <- data.frame(
      'disease_id' =  disease_item$id,
      'name' = disease_item$name,
      stringsAsFactors = F)
    
    if(!is.null(disease_item$description)){
      df$description <- disease_item$description
    }else{
      df$description <- NA
    }
    
    OT_diseases <- OT_diseases |>
      dplyr::bind_rows(df)
    
    i <- i + 1
  
  }
  close.connection(con)
  m <- m + 1
}



opentargets_target <- 
  readRDS(file=paste0("output/opentargets_target_",release,".rds"))

i <- 0

OT_association_all <- data.frame()

all_datatypes <- data.frame()

while(i <= 199){
  
  overallAssoc <- 
    readRDS(file=paste0("output/association/OT_overall_association_",release,"_",i,".rds"))
  datatypeAssoc <- 
    readRDS(file=paste0("output/association/OT_datatype_association_",release,"_",i,".rds"))
  datasourceAssoc <- 
    readRDS(file=paste0("output/association/OT_datasource_association_",release,"_",i,".rds"))
  
  datatype_df <- as.data.frame(
    datatypeAssoc |>
      dplyr::mutate(score = round(score, digits = 12),
                    overall_datatype_evidence_count = evidence_count) |>
      dplyr::mutate(datatype_support = paste0(datatype_id,"|",evidence_count,"|",score)) |>
      dplyr::group_by(disease_id, target_id) |>
      dplyr::summarise(datatype_items = paste(datatype_support, collapse=","),
                       .groups = "drop")
  )
  
  all_datatypes <- all_datatypes |>
    dplyr::bind_rows(datatype_df)
  
  datasource_df <- as.data.frame(
    datasourceAssoc |>
      dplyr::mutate(score = round(score, digits = 12),
                    overall_datasource_evidence_count = evidence_count) |>
      dplyr::mutate(datasource_support = paste0(datasource_id,"|",evidence_count,"|",score)) |>
      dplyr::group_by(disease_id, target_id) |>
      dplyr::summarise(datasource_items = paste(datasource_support, collapse=","),
                       .groups = "drop")
  )

  dataoverall_df <- as.data.frame(
    overallAssoc |>
      dplyr::mutate(score = round(score, digits = 12)) |>
      dplyr::mutate(json_chunk = i) |>
      dplyr::left_join(datatype_df, by = c("disease_id","target_id")) |>
      dplyr::left_join(datasource_df, by = c("disease_id","target_id")) |>
      dplyr::rename(target_ensembl_gene_id = target_id) |>
      dplyr::left_join(
        dplyr::select(opentargets_target, target_symbol, target_ensembl_gene_id),
        by = c("target_ensembl_gene_id")) |>
      dplyr::left_join(
        dplyr::select(OT_diseases, disease_id, name),
        by = "disease_id") |>
      dplyr::rename(disease_label = name)
  )

  OT_association_all <- OT_association_all |>
    dplyr::bind_rows(dataoverall_df)
  cat(i,'\n')
  i <- i + 1
}


saveRDS(OT_association_all, 
        file=paste0("output/opentargets_association_direct_",
                    release,".rds"))
        
OT_association_hc <- OT_association_all |>
  dplyr::filter(!stringr::str_detect(
    tolower(disease_label),"^(neoplasm|cancer)$| neoplasm$")) |>
  dplyr::filter(stringr::str_detect(datatype_items,","))

saveRDS(OT_association_hc, 
        file=paste0("output/opentargets_association_direct_HC_",
                    release,".rds"))

