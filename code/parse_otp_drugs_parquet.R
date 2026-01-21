#!/usr/local/bin/Rscript

release <- '2025.12'

OT_datasets <- list()
OT_datasets[['known_drug']] <- data.frame()
OT_datasets[['molecule']] <- data.frame()
OT_datasets[['drug_moa']] <- data.frame()
OT_datasets[['drug_indication']] <- data.frame()
OT_datasets[['drug_target']] <- readRDS(
  file = file.path(
    here::here(), "output", 
    paste0("opentargets_target_", release,".rds")
  )
)

####-- MECHANISM_OF_ACTION --####

basepath <- file.path(here::here(), "data", 
                      release, "drug_mechanism_of_action")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))

m <- 1
for(parquet_fname in parquet_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  moa_data <- arrow::read_parquet(parquet_fname)
  
  i <- 1
  while(i <= NROW(moa_data)){
    moa_item <- moa_data[i,]
    action_type <- NA
    if(!is.null(moa_item$actionType)){
      action_type <- moa_item$actionType
    }
    df <- data.frame(
      'parquet_chunk' = stringr::str_replace(
        basename(parquet_fname), ".snappy.parquet", ""),
      'drug_action_type' =  action_type,
      stringsAsFactors = F)
    for(e in c('mechanismOfAction',
               'targetType')){
      if(!is.null(moa_item[[e]])){
        df[,e] <- moa_item[[e]]
      }else{
        df[,e] <- NA
      }
    }
    
    df <- df |>
      dplyr::rename(drug_moa = mechanismOfAction,
                    target_type = targetType) |>
      dplyr::mutate(
        target_type = stringr::str_replace_all(
          target_type, " |-", "_")
      )
    
    df$drug_moa_references <- NA
    if(is.list(moa_item$references)){
      if(length(moa_item$references) > 0){
        if("source" %in% colnames(moa_item$references[[1]]) &
           "ids" %in% colnames(moa_item$references[[1]])){
          if(NROW(moa_item$references[[1]]) >= 1){
            j <- 1
            all_sources <- c()
            while(j <= NROW(moa_item$references[[1]])){
              source <- moa_item$references[[1]][j,]$source
              source_ids <- NA
              if(!is.null(moa_item$references[[1]][j,]$ids)){
                source_ids <- paste(
                  unique(moa_item$references[[1]][j,]$ids[[1]]),
                  collapse="&"
                )
              }
              source_element <- paste(
                source, source_ids, sep="|"
              )
              all_sources <- 
                c(all_sources,
                  source_element)
              j <- j + 1
            }
            df$drug_moa_references <- 
              paste(sort(all_sources), collapse=",")
          }
        }
      }
    }
    
    
    for(n in 1:length(moa_item$chemblIds[[1]])){
      molecule_chembl_id <- moa_item$chemblIds[[1]][n]
      moa_entry <- df
      moa_entry$molecule_chembl_id <- molecule_chembl_id
      if(length(moa_item$targets[[1]]) > 0){
        for(k in 1:length(moa_item$targets[[1]])){
          ensembl_gene_id <- moa_item$targets[[1]][k]
          moa_entry2 <- moa_entry
          moa_entry2$target_ensembl_gene_id <- ensembl_gene_id
          OT_datasets[['drug_moa']] <- OT_datasets[['drug_moa']] |>
            dplyr::bind_rows(moa_entry2)
        }
      }else{
        moa_entry$target_ensembl_gene_id <- NA
        OT_datasets[['drug_moa']] <- OT_datasets[['drug_moa']] |>
          dplyr::bind_rows(moa_entry)
      }
      
    }
    i <- i + 1
  }
  m <- m + 1
}



####-- KNOWN DRUG --####

basepath <- file.path(here::here(), "data", 
                      release, "known_drug")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))

m <- 1
for(parquet_fname in sort(parquet_files)){
  cat(paste0('Chunk - ',m))
  cat('\n')
  drug_data <- 
    arrow::read_parquet(parquet_fname)
  
  drug_data$tradeNames <- NULL
  #drug_data$synonyms <- NULL
  
  df <- data.frame(
    'parquet_chunk' = stringr::str_replace(
      basename(parquet_fname), ".snappy.parquet", ""),
    'molecule_chembl_id' =  drug_data$drugId,
    'target_ensembl_gene_id' = drug_data$targetId,
    'drug_name' = drug_data$prefName,
    'target_symbol' = drug_data$approvedSymbol,
    'drug_type' = drug_data$drugType)
  
  df$target_class <- 
    sapply(
      drug_data$targetClass, function(x) paste(x, collapse = "|"))
  
  df$drug_synonyms <- 
    sapply(
      drug_data$synonyms, function(x) paste(x, collapse = "|"))
  
  df <- df |>
    dplyr::select(
      molecule_chembl_id,
      drug_name,
      drug_type,
      drug_synonyms,
      target_ensembl_gene_id,
      target_symbol,
      target_class,
      dplyr::everything())
  
  OT_datasets[['known_drug']] <- OT_datasets[['known_drug']] |>
    dplyr::bind_rows(df)

  m <- m + 1
}

    
####-- MOLECULES --####

basepath <- file.path(here::here(), "data", 
                      release, "drug_molecule")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))

m <- 1
for(parquet_fname in sort(parquet_files)){
  cat(paste0('Chunk - ',m))
  cat('\n')
  molecule_data <- 
    arrow::read_parquet(parquet_fname)
  
  i <- 1
  while(i <= NROW(molecule_data)){
    molecule_item <- molecule_data[i,]    
    df <- data.frame(
      #'parquet_chunk' = stringr::str_replace(
      #  basename(parquet_fname), ".snappy.parquet", ""),
      'molecule_chembl_id' =  molecule_item$id)
    for(e in c('name',
               'isApproved',
               'drugType',
               'canonicalSmiles',
               'inchiKey',
               'blackBoxWarning',
               'maximumClinicalTrialPhase',
               'hasBeenWithdrawn',
               'description')){
      if(!is.null(molecule_item[[e]])){
        df[,e] <- molecule_item[[e]]
      }else{
        df[,e] <- NA
      }
    }

    df <- df |>
      dplyr::rename(
        drug_name = name,
        drug_inchi_key = inchiKey,
        drug_is_approved = isApproved,
        drug_canonical_smiles = canonicalSmiles,
        drug_withdrawn = hasBeenWithdrawn,
        drug_max_ct_phase = maximumClinicalTrialPhase,
        drug_blackbox_warning = blackBoxWarning,
        drug_type = drugType,
        drug_description = description) |>
      dplyr::mutate(
        drug_is_approved = dplyr::if_else(
          is.na(drug_is_approved), FALSE, 
          as.logical(drug_is_approved)),
        )

    df$drug_synonyms <- NA
    df$drug_tradenames <- NA
    for(e in c('synonyms','tradeNames')){
      if(!is.null(molecule_item[[e]])){
        if(is.character(molecule_item[[e]][[1]]) &
           length(molecule_item[[e]][[1]]) > 0){
          df[,paste0('drug_',tolower(e))] <-
            paste(unique(sort(molecule_item[[e]][[1]])), collapse = "|")
        }
      }
    }

    df$drug_year_first_approval <- NA
    if(!is.null(molecule_item[['yearOfFirstApproval']])){
      df$drug_year_first_approval <-
        molecule_item[['yearOfFirstApproval']]
    }

    df$drug_linked_diseases <- NA
    if(!is.null(molecule_item[['linkedDiseases']])){
      if(is.data.frame(molecule_item[['linkedDiseases']]) &
         NROW(molecule_item[['linkedDiseases']]) > 0){
        if(length(molecule_item[['linkedDiseases']]$rows[[1]]) > 0){
          df$drug_linked_diseases <-
            paste(unique(molecule_item[['linkedDiseases']]$rows[[1]]),
                  collapse = "|")
        }
      }
    }

    df$drug_description <- NA
    if(!is.null(molecule_item[['description']])){
      df$drug_description <- 
        molecule_item[['description']]
    }

    df$parent_molecule_chembl_id <- NA
    if(!is.null(molecule_item[['parentId']])){
      df$parent_molecule_chembl_id <- 
        molecule_item[['parentId']]
    }

    df$drug_pubchem_xref <- NA
    df$drug_drugbank_xref <- NA
    df$drug_chebi_xref <- NA
    df$drug_wikipedia_xref <- NA
    if(!is.null(molecule_item[['crossReferences']])){
      if(is.data.frame(molecule_item[['crossReferences']][[1]])){
        drug_df <- molecule_item[['crossReferences']][[1]]
        if(NROW(drug_df) > 0){
          for(j in 1:nrow(drug_df)){
            for(db in c('PubChem','drugbank','chEBI','Wikipedia')){
              if(drug_df[j,"source"] == db){
                df[,paste0("drug_",tolower(db),"_xref")] <- 
                  paste(drug_df[j,"ids"][[1]][[1]], collapse="|")
              }
            }
          }
        }
      }
    }

    all_drug_child_entries <- data.frame()
    if(!is.null(molecule_item[['childChemblIds']])){
      for(n in 1:length(molecule_item[['childChemblIds']][[1]])){
        drug_child_entry <- df
        drug_child_entry$child_molecule_chembl_id <-
          molecule_item[['childChemblIds']][[1]][n]

        all_drug_child_entries <- all_drug_child_entries |>
          dplyr::bind_rows(drug_child_entry)
      }
    }

    if(nrow(all_drug_child_entries) == 0){
      all_drug_child_entries <- df
      all_drug_child_entries$child_molecule_chembl_id <- NA
    }


    if(!is.null(molecule_item[['linkedTargets']])){
      if(molecule_item[['linkedTargets']]$count > 0){
          j <- 1
          while(j <= nrow(all_drug_child_entries)){
            drug_target_entry <- all_drug_child_entries[j,]
            k <- 1
            while(k <= length(molecule_item[['linkedTargets']]$rows[[1]])){
              drug_target_entry$target_ensembl_gene_id <-
                molecule_item[['linkedTargets']]$rows[[1]][k]
              OT_datasets[['molecule']] <- OT_datasets[['molecule']] |>
                dplyr::bind_rows(drug_target_entry)
              k <- k + 1
            }
            j <- j + 1
          }
      }else{
        df$target_ensembl_gene_id <- NA
        OT_datasets[['molecule']] <- OT_datasets[['molecule']] |>
          dplyr::bind_rows(df)

      }
    }else{
      df$target_ensembl_gene_id <- NA
      OT_datasets[['molecule']] <- OT_datasets[['molecule']] |>
        dplyr::bind_rows(df) |>
        dplyr::distinct()

    }

    i <- i + 1
  }

  m <- m + 1
}

OT_datasets[['molecule']] <- OT_datasets[['molecule']] |>
  dplyr::distinct() |>
  dplyr::group_by(
    dplyr::across(-c("child_molecule_chembl_id"))) |> 
  dplyr::reframe(
    child_molecule_chembl_id = paste(
      child_molecule_chembl_id, collapse="&")) |> 
  dplyr::distinct()


####-- INDICATIONS --####

basepath <- file.path(here::here(), "data", 
                      release, "drug_indication")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))
m <- 1
for(parquet_fname in parquet_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  indication_data <- 
    arrow::read_parquet(parquet_fname)
  i <- 1
  
  while(i <= NROW(indication_data)){
    indication_item <- indication_data[i,]
    molecule_chembl_id <- NA
    if(!is.null(indication_item$id)){
      molecule_chembl_id <- indication_item$id
    }
    
    df <- data.frame(
      #'parquet_chunk' = stringr::str_replace(
      #  basename(parquet_fname), ".snappy.parquet", ""),
      'molecule_chembl_id' =  molecule_chembl_id,
      stringsAsFactors = F)
    
    
    approved_indications <- data.frame()
    if(!is.null(indication_item$approvedIndications[[1]])){
      if(length(indication_item$approvedIndications[[1]]) > 0){
        for(n in 1:length(indication_item$approvedIndications[[1]])){
          df2 <- 
            data.frame('disease_efo_id' = 
                         indication_item$approvedIndications[[1]][n],
                       'drug_approved_indication' = TRUE,
                       stringsAsFactors = F)
          approved_indications <- dplyr::bind_rows(
            approved_indications, df2)
        }
      }
    }
    
    OT_indication_entries <- data.frame()
    if(NROW(indication_item$indications[[1]]) > 0){
      k <- 1
      while(k <= NROW(indication_item$indications[[1]])){
        indication <- indication_item$indications[[1]][k,]
        indication_entry <- df
        for(e in c('disease',
                   'efoName',
                   'maxPhaseForIndication')){
          if(!is.null(indication[[e]])){
            indication_entry[,e] <- indication[[e]]
          }else{
            indication_entry[,e] <- NA
          }
        }
        
        indication_entry <- indication_entry |>
          dplyr::rename(
            disease_efo_id = disease,
            disease_efo_label = efoName,
            drug_max_phase_indication = maxPhaseForIndication)
        
        
        all_references <- c()
        if(NROW(indication$references[[1]]) > 0){
          for(n in 1:NROW(indication$references[[1]])){
            ref <- indication$references[[1]][n,]
            
            for(u in 1:length(ref$ids[[1]])){
              e <- indication_entry
              e$drug_clinical_source <- ref$source
              e$drug_clinical_id <- ref$ids[[1]][u]
              
              e <- e |>
                dplyr::mutate(
                  drug_clinical_source = dplyr::if_else(
                    drug_clinical_source == "ClinicalTrials",
                    as.character("clinicaltrials.gov"),
                    as.character(drug_clinical_source)
                  ))
              
              OT_indication_entries <- OT_indication_entries |>
                dplyr::bind_rows(e)
            }
          }
        }
        k <- k + 1
      }
    }
    
    if(nrow(approved_indications) > 0){
      OT_indication_entries <- 
        OT_indication_entries |> 
        dplyr::left_join(
          approved_indications, by = "disease_efo_id") |>
        dplyr::mutate(
          drug_approved_indication =
            dplyr::if_else(
              is.na(drug_approved_indication),
              as.logical(FALSE),
              as.logical(drug_approved_indication)
            ))
    }else{
      OT_indication_entries$drug_approved_indication <- FALSE
    }
    
    OT_datasets[['drug_indication']] <- OT_datasets[['drug_indication']] |>
      dplyr::bind_rows(OT_indication_entries)
    
    #cat(i,'\n')
    i <- i + 1
  }
  m <- m + 1
}

## Get all drug synonyms
synonyms_known_drugs <- OT_datasets$known_drug |> 
  dplyr::rename(drug_synonyms2 = drug_synonyms) |> 
  dplyr::select(molecule_chembl_id, drug_synonyms2)

drug_synonyms <- 
  OT_datasets$molecule |> 
  dplyr::select(molecule_chembl_id, drug_synonyms) |> 
  dplyr::left_join(synonyms_known_drugs, relationship = "many-to-many") |> 
  dplyr::mutate(drug_synonyms = dplyr::if_else(
    !is.na(drug_synonyms) & !is.na(drug_synonyms2), 
    paste(drug_synonyms, drug_synonyms2, sep="|"), 
    as.character(drug_synonyms))) |> 
  dplyr::distinct() |> 
  dplyr::select(-c("drug_synonyms2")) |> 
  tidyr::separate_rows(drug_synonyms, sep="\\|") |> 
  dplyr::distinct() |> 
  dplyr::group_by(molecule_chembl_id) |> 
  dplyr::summarise(drug_synonyms = paste(
    drug_synonyms, collapse="|")) |> 
  dplyr::mutate(drug_synonyms = dplyr::if_else(
    drug_synonyms == "NA", "", 
    as.character(drug_synonyms))) |> 
  dplyr::arrange(molecule_chembl_id)

## Combine molecules, mechanism of action data, and drug indication data

OT_drugs_complete <- OT_datasets$molecule |> 
  dplyr::select(-c("drug_synonyms")) |>
  dplyr::left_join(
    drug_synonyms, 
    by = "molecule_chembl_id") |>
  dplyr::left_join(
    OT_datasets$drug_moa, 
    by = c("molecule_chembl_id", 
           "target_ensembl_gene_id")) |> 
  dplyr::mutate(child_molecule_chembl_id = dplyr::if_else(
    child_molecule_chembl_id == "NA",
    NA_character_,
    as.character(child_molecule_chembl_id)
  )) |>
  dplyr::distinct() |> 
  dplyr::left_join(
    OT_datasets$drug_indication, 
    by = "molecule_chembl_id",
    relationship = "many-to-many") |> 
  dplyr::distinct() |>
  dplyr::left_join(
    dplyr::select(
      OT_datasets$drug_target, 
      c("target_ensembl_gene_id", 
        "target_symbol", 
        "target_name",
        "target_entrezgene")),
    by = "target_ensembl_gene_id") |>
  dplyr::rename(target_genename = target_name) |>
  dplyr::distinct() |>
  dplyr::select(-c("parquet_chunk")) |>
  dplyr::distinct()

saveRDS(OT_drugs_complete, 
        file = file.path(
          "output",
          paste0(
            "opentargets_drugs_",
            release,".rds"))
)


        