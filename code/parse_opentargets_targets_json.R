release <- '2023.12'

####---- TARGETS ----####

# wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.04/output/etl/json/targets

basepath <- file.path(here::here(), "data", 
                      release, "targets")
json_files <- sort(
  list.files(basepath, pattern = ".json", 
             all.files = T, full.names = T))

OT_target <- data.frame()

m <- 1
for(json_chunk_lines in json_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  con <- file(json_chunk_lines, "r")
  target_data <- 
    jsonlite::stream_in(con, flatten = F, simplifyVector = F)
  i <- 1
  
  while(i <= length(target_data)){
    target_item <- target_data[[i]]
    
    df <- data.frame(
      'target_ensembl_gene_id' =  target_item$id,
      stringsAsFactors = F)
    for(e in c('approvedSymbol',
               'biotype')){
      if(!is.null(target_item[[e]])){
        df[,e] <- target_item[[e]]
      }else{
        df[,e] <- NA
      }
    }
    
    hgnc_id <- NA
    if(is.list(target_item$dbXrefs)){
      if(length(target_item$dbXrefs) > 0){
        for(k in 1:length(target_item$dbXrefs)){
          if(target_item$dbXrefs[[k]]$source == "HGNC"){
            hgnc_id <- target_item$dbXrefs[[k]]$id
          }
        }
      }
    }
    df$hgncId <- hgnc_id
    
    
    df <- df |>
      dplyr::rename(
        target_symbol = approvedSymbol,
        target_biotype = biotype,
        hgnc_id = hgncId)
    
    if(!is.null(target_item$functionDescriptions)){
      df$function_description <- 
        paste(target_item$functionDescriptions, 
              collapse = "|")
      
      df$function_description <- stringr::str_replace(
        stringr::str_replace_all(
          df$function_description,
          "(ECO:[0-9]{1,}\\|(UniProtKB:[A-Z0-9]{1,}|PubMed:[0-9]{1,}))(, )?",
          ""
        ), "\\{\\}\\.$","")
      
    }else{
      df$function_description <- NA
    }
    
    
    
    tractability_small_molecule <- ""
    tractability_antibody <- ""
    if(!is.null(target_item$tractability)){
      for(j in 1:length(target_item$tractability)){
        tractability_item <- target_item$tractability[[j]]
        if(tractability_item$modality == "SM" &
           tractability_item$value == T){
          tractability_small_molecule <- paste(
            tractability_small_molecule, tractability_item$id,
            sep = " <b>|</b> "
          )
        }
        if(tractability_item$modality == "AB" &
           tractability_item$value == T){
          tractability_antibody <- paste(
            tractability_antibody, tractability_item$id,
            sep = " <b>|</b> "
          )
        }
      }
    }
    if(tractability_antibody == ""){
      df$AB_tractability_support <- ""
      df$AB_tractability_category <- "Unknown"
    }else{
      df$AB_tractability_support <- stringr::str_replace(
        tractability_antibody,"^ <b>\\|</b> ","")
      
      df$AB_tractability_category <- "Unknown"
      if(stringr::str_detect(
        tractability_antibody,"Phase 1 Clinical|Approved Drug|Advanced Clinical")){
        df$AB_tractability_category <- "Clinical_Precedence"
      }
      else if(stringr::str_detect(
        tractability_antibody,"UniProt loc high conf|GO CC high conf")){
        df$AB_tractability_category <- "Predicted_Tractable_High_confidence"
      }
      else if(stringr::str_detect(
        tractability_antibody,"UniProt loc med conf|UniProt SigP or TMHMM|Human Protein Atlas loc|GO CC med conf")){
        df$AB_tractability_category <- "Predicted_Tractable_Medium_to_low_confidence"
      }
      
    }
    if(tractability_small_molecule == ""){
      df$SM_tractability_support <- ""
      df$SM_tractability_category <- "Unknown"
    }else{
      df$SM_tractability_support <- stringr::str_replace(
        tractability_small_molecule,"^ <b>\\|</b> ","")
      
      df$SM_tractability_category <- "Unknown"
      if(stringr::str_detect(
        tractability_small_molecule,"Phase 1 Clinical|Approved Drug|Advanced Clinical")){
        df$SM_tractability_category <- "Clinical_Precedence"
      }
      else if(stringr::str_detect(
        tractability_small_molecule,"Structure with Ligand|High-Quality Ligand|High-Quality Pocker")){
        df$SM_tractability_category <- "Discovery_Precedence"
      }
      else if(stringr::str_detect(
        tractability_small_molecule,"Med-Quality Pocket|Druggable Family")){
        df$SM_tractability_category <- "Predicted_Tractable"
      }
      
    }
    
    if(length(target_item$hallmarks) > 0){
      if(!is.null(target_item$hallmarks$cancerHallmarks)){
        
        all_labels <- list()
        for(n in 1:length(target_item$hallmarks$cancerHallmarks)){
          hallmark_eitem <-
            target_item$hallmarks$cancerHallmarks[[n]]
          
          # DEBUG
          #cat(names(hallmark_eitem), paste(df$target_symbol,m,i,sep=" - "), "\n")
          
          hallmark_eitem$promote <- FALSE
          hallmark_eitem$suppress <- FALSE
          
          if(!is.null(hallmark_eitem$impact)){
            if(hallmark_eitem$impact == "promotes"){
              hallmark_eitem$promote <- TRUE
            }
            if(hallmark_eitem$impact == "suppresses"){
              hallmark_eitem$suppress <- TRUE
            }
          }else{
            hallmark_eitem$suppress <- TRUE
          }
          label <- paste(
            hallmark_eitem$label,
            paste0("PROMOTE:", hallmark_eitem$promote),
            paste0("SUPPRESS:", hallmark_eitem$suppress),
            sep = "|")
          if(is.null(all_labels[[label]])){
            all_labels[[label]] <-
              paste(hallmark_eitem$pmid,
                    hallmark_eitem$description,
                    sep="%%%")
          }else{
            all_labels[[label]] <- paste(
              all_labels[[label]],
              paste(hallmark_eitem$pmid,
                    hallmark_eitem$description,
                    sep="%%%"),
              sep="&"
            )
          }
        }
        
        all_hallmark_entries <- c()
        if(length(all_labels) > 0){
          for(u in names(all_labels)){
            all_hallmark_entries <- c(
              all_hallmark_entries,
              paste0(u, "|",all_labels[[u]])
            )
            
          }
          df$cancer_hallmark <-
            paste(all_hallmark_entries,
                  collapse="@@@")
        }
      }
      if(!is.null(target_item$hallmarks$attributes)){
        
        summary_found <- 0
        for(v in 1:length(target_item$hallmarks$attributes)){
          
          if(target_item$hallmarks$attributes[[v]]$attribute_name == 
             "function summary"){
            summary_found <- 1
            df$cancer_hallmark_function_summary <- 
              target_item$hallmarks$attributes[[v]]$description
          }
        }
        if(summary_found == 0){
          df$cancer_hallmark_function_summary <- NA
        }
      }else{
        df$cancer_hallmark_function_summary <- NA
      }
    }else{
      df$cancer_hallmark_function_summary <- NA
      df$cancer_hallmark <- NA
    }
    
    OT_target <- OT_target |>
      dplyr::bind_rows(df)
    
    i <- i + 1
  }
  close.connection(con)
  m <- m + 1
}

duplicate_symbols <- as.data.frame(
  OT_target |> 
  dplyr::group_by(target_symbol) |> 
  dplyr::summarise(n = dplyr::n()) |> 
  dplyr::filter(n > 1)
)

entries_for_removal <- duplicate_symbols |>
  dplyr::inner_join(OT_target) |>
  dplyr::filter(is.na(hgnc_id)) |>
  dplyr::select(target_ensembl_gene_id, target_symbol)


OT_target_clean <- OT_target |>
  dplyr::anti_join(entries_for_removal) |>
  dplyr::mutate(
    AB_tractability_category =
      factor(AB_tractability_category,
             levels = c("Clinical_Precedence",
                        "Predicted_Tractable_High_confidence",
                        "Predicted_Tractable_Medium_to_low_confidence",
                        "Unknown"))
  ) |>
  dplyr::mutate(
    SM_tractability_category =
      factor(SM_tractability_category,
             levels = c("Clinical_Precedence",
                        "Discovery_Precedence",
                        "Predicted_Tractable",
                        "Unknown"))
  )


saveRDS(OT_target_clean, 
        file = file.path(
          "output", 
          paste0("opentargets_target_",
                 release,".rds")))

rm(OT_target)
rm(target_data)

