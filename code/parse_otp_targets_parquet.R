release <- '2026.03'

####---- GENE CROSSREF----####
gene_oncox <- list()
gene_oncox[['basic']] <- geneOncoX::get_basic(cache_dir = file.path(
  here::here(), "data"
))
gene_oncox[['gencode']] <- geneOncoX::get_gencode(cache_dir = file.path(
  here::here(), "data"
))

gene_xref <- dplyr::bind_rows(
  dplyr::inner_join(
    dplyr::select(gene_oncox[['basic']]$records, entrezgene, 
                  symbol, name, gene_biotype),
    dplyr::select(gene_oncox[['gencode']]$records$grch38, 
                  entrezgene, ensembl_gene_id),
    by = c("entrezgene")),
  dplyr::inner_join(
    dplyr::select(gene_oncox[['basic']]$records, 
                  entrezgene, symbol, name, gene_biotype),
    dplyr::select(gene_oncox[['gencode']]$records$grch37, 
                  entrezgene, ensembl_gene_id),
    by = c("entrezgene"))) |>
  dplyr::filter(!is.na(ensembl_gene_id)) |>
  dplyr::rename(target_genename = name,
                target_symbol2 = symbol,
                target_entrezgene = entrezgene,
                target_ensembl_gene_id = ensembl_gene_id) |>
  dplyr::select(target_genename, 
                target_symbol2, 
                target_entrezgene,
                target_ensembl_gene_id) |>
  dplyr::distinct()


####---- TARGETS ----####

basepath <- file.path(here::here(), "data", 
                      release, "target")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))

OT_target <- data.frame()

m <- 1
for(parquet_fname in parquet_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  target_data <- 
    arrow::read_parquet(parquet_fname)
  i <- 1
  
  while(i <= nrow(target_data)){
    target_item <- target_data[i,]
    
    df <- data.frame(
      'target_ensembl_gene_id' =  target_item$id,
      'target_name' = target_item$approvedName,
      'fname' = basename(parquet_fname),
      'idx' = i,
      stringsAsFactors = F)
    for(e in c('approvedSymbol',
               'biotype')){
      if(!is.null(target_item[[e]])){
        df[,e] <- target_item[[e]]
      }else{
        df[,e] <- NA
      }
    }
    
    df <- df |>
      dplyr::left_join(
        dplyr::select(
          gene_xref,
          target_ensembl_gene_id,
          target_entrezgene),
        by = c("target_ensembl_gene_id" = "target_ensembl_gene_id")
      )
    
    hgnc_id <- NA
    
    if (is.list(target_item$dbXrefs) && length(target_item$dbXrefs) > 0) {
      # Keep only non-empty data.frames
      dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, target_item$dbXrefs)
      
      if (length(dfs) > 0) {
        # Bind them together into one data.frame
        all_refs <- do.call(rbind, dfs)
        
        # Extract HGNC id if present
        if ("source" %in% names(all_refs) && "id" %in% names(all_refs)) {
          idx <- which(all_refs$source == "HGNC")
          if (length(idx) > 0) {
            hgnc_id <- all_refs$id[idx[1]]  # first match
          }
        }
      }
    }
    
    
    df$hgnc_id <- hgnc_id
    
    df <- df |>
      dplyr::rename(
        target_symbol = approvedSymbol,
        target_biotype = biotype)
    
    if(!is.null(target_item$functionDescriptions) &
       is.list(target_item$functionDescriptions)){
      df$function_description <- 
        paste(target_item$functionDescriptions[[1]], 
              collapse = "|")
      
      df$function_description <- stringr::str_replace(
        stringr::str_replace_all(
          df$function_description,
          "(ECO:[0-9]{1,}\\|(UniProtKB:[A-Z0-9]{1,}|PubMed:[0-9]{1,}))(, )?|ECO:[0-9]{1,}",
          ""
        ), "\\{\\}\\.$","")
      
    }else{
      df$function_description <- NA
    }
    
    tractability <- stats::setNames(as.list(rep("", 3)), c("SM", "AB", "PR"))
    
    if (is.list(target_item$tractability) && length(target_item$tractability) > 0) {
      # Keep only non-empty data.frames
      dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, target_item$tractability)
      
      if (length(dfs) > 0) {
        # Combine into one data.frame
        all_tract <- do.call(rbind, dfs)
        
        # Process each modality (SM, AB, PR)
        for (t in names(tractability)) {
          ids <- all_tract$id[all_tract$modality == t & all_tract$value == TRUE]
          if (length(ids) > 0) {
            tractability[[t]] <- paste(ids, collapse = " <b>|</b> ")
          }
        }
      }
    }
    
    for(t in c('SM','AB','PR')){
      if(tractability[[t]] == ""){
        df[paste0(t,'_tractability_support')] <- ""
        df[paste0(t,'_tractability_category')] <- "Unknown"
      }else{
        df[paste0(t,'_tractability_support')] <- stringr::str_replace(
          tractability[[t]],"^ <b>\\|</b> ","")
        df[paste0(t,'_tractability_category')] <- "Unknown"
        if(stringr::str_detect(
          tractability[[t]],
          "Phase 1 Clinical|Approved Drug|Advanced Clinical")){
          df[paste0(t,'_tractability_category')] <- "Clinical_Precedence"
        }
        else if(stringr::str_detect(
          tractability[[t]],
          "UniProt loc high conf|GO CC high conf")){
          df[paste0(t,'_tractability_category')] <- 
            "Predicted_Tractable_High_confidence"
        }
        else if(stringr::str_detect(
          tractability[[t]],
          "Structure with Ligand|High-Quality Ligand|High-Quality Pocker")){
          df[paste0(t,'_tractability_category')] <- 
            "Discovery_Precedence"
        }
        else if(stringr::str_detect(
          tractability[[t]],
          "Med-Quality Pocket|Druggable Family")){
          df[paste0(t,'_tractability_category')] <- "Predicted_Tractable"
        }
        else if(stringr::str_detect(
          tractability[[t]],
          "UniProt loc med conf|UniProt SigP or TMHMM|Human Protein Atlas loc|GO CC med conf")){
          df[paste0(t,'_tractability_category')] <- 
            "Predicted_Tractable_Medium_to_low_confidence"
        }
        
      }
    }
    
    df$cancer_hallmark_function_summary <- NA
    df$cancer_hallmark_role <- NA
    df$cancer_hallmark <- NA
    if(length(target_item$hallmarks) > 0){
      if(!is.null(target_item$hallmarks$cancerHallmarks) &
         is.data.frame(target_item$hallmarks$cancerHallmarks[[1]])){
        
        all_labels <- list()
        if(NROW(target_item$hallmarks$cancerHallmarks[[1]]) > 0){
          for(n in 1:nrow(target_item$hallmarks$cancerHallmarks[[1]])){
            hallmark_eitem <-
              target_item$hallmarks$cancerHallmarks[[1]][n,]
            
            # DEBUG
            #cat(names(hallmark_eitem), paste(df$target_symbol,m,i,sep=" - "), "\n")
            
            hallmark_eitem$promote <- FALSE
            hallmark_eitem$suppress <- FALSE
            
            if(!is.null(hallmark_eitem$impact) &
               !is.na(hallmark_eitem$impact)){
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
      }
      if(!is.null(target_item$hallmarks$attributes)){
        
        summary_found <- 0
        role_in_cancer <- c()
        function_summary <- c()
        if(NROW(target_item$hallmarks$attributes[[1]]) > 0){
          for(v in 1:nrow(target_item$hallmarks$attributes[[1]])){
            if(isTRUE(target_item$hallmarks$attributes[[1]][v,]$name ==
                      "function summary")){
              function_summary <- c(
                function_summary,
                stringr::str_trim(
                  target_item$hallmarks$attributes[[1]][v,]$description))
            }
            if(isTRUE(target_item$hallmarks$attributes[[1]][v,]$name ==
                      "role in cancer")){
              if(nchar(stringr::str_trim(
                target_item$hallmarks$attributes[[1]][v,]$description)) < 30){
                role_in_cancer <- c(
                  role_in_cancer,
                  stringr::str_trim(
                    target_item$hallmarks$attributes[[1]][v,]$description))
              }
            }
          }
        }
        if(length(function_summary) > 0){
          df$cancer_hallmark_function_summary <- 
            paste(function_summary, collapse = "|")
        }
        if(length(role_in_cancer) > 0){
          df$cancer_hallmark_role <- 
            paste(unique(sort(role_in_cancer)), collapse = "|")
          df$cancer_hallmark_role <- stringr::str_replace_all(
            df$cancer_hallmark_role,
            c("Oncogene" = "oncogene",
              "TSG" = "tsg")
          )
        }
      }
    }
    
    OT_target <- OT_target |>
      dplyr::bind_rows(df)
    
    i <- i + 1
  }
  m <- m + 1
}


# Remove duplicates
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

