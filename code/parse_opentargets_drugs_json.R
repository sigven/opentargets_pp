release <- '2023.12'

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
                target_symbol = symbol,
                target_entrezgene = entrezgene,
                target_ensembl_gene_id = ensembl_gene_id) |>
  dplyr::select(target_genename, 
                target_symbol, 
                target_entrezgene,
                target_ensembl_gene_id) |>
  dplyr::distinct()
  


####-- MOLECULES --####

basepath <- file.path(here::here(), "data", 
                      release, "molecule")
json_files <- sort(
  list.files(basepath, pattern = ".json", 
             all.files = T, full.names = T))

OT_drugs <- data.frame()

m <- 1
for(json_chunk_lines in sort(json_files)){
  cat(paste0('Chunk - ',m))
  cat('\n')
  con <- file(json_chunk_lines, "r")
  drug_data <- 
    jsonlite::stream_in(con, flatten = F, simplifyVector = F)
  
  i <- 1
  while(i <= length(drug_data)){
    drug_item <- drug_data[[i]]
    df <- data.frame(
      'molecule_chembl_id' =  drug_item$id,
      stringsAsFactors = F)
    for(e in c('name',
               'isApproved',
               'drugType',
               'canonicalSmiles',
               'inchiKey',
               'blackBoxWarning',
               'maximumClinicalTrialPhase',
               'hasBeenWithdrawn',
               'description')){
      if(!is.null(drug_item[[e]])){
        df[,e] <- drug_item[[e]]
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
        drug_description = description)
    
    df$drug_synonyms <- NA
    if(!is.null(drug_item[['synonyms']])){
      df$drug_synonyms <- 
        paste(drug_item[['synonyms']], collapse = "|")
    }
    
    df$drug_tradenames <- NA
    if(length(drug_item[['tradeNames']]) > 0){
      df$drug_tradenames <- 
        paste(drug_item[['tradeNames']], collapse = "|")
    }
    
    df$drug_year_first_approval <- NA
    if(!is.null(drug_item[['yearOfFirstApproval']])){
      df$drug_year_first_approval <- 
        drug_item[['yearOfFirstApproval']]
    }
    
    df$drug_linked_diseases <- NA
    if(!is.null(drug_item[['linkedDiseases']])){
      df$drug_linked_diseases <- 
        paste(drug_item[['linkedDiseases']]$rows, 
              collapse = "|")
    }
    
    df$drug_description <- NA
    if(!is.null(drug_item[['description']])){
      df$drug_description <- drug_item[['description']]
    }
    
    df$parent_molecule_chembl_id <- NA
    if(!is.null(drug_item[['parentId']])){
      df$parent_molecule_chembl_id <- drug_item[['parentId']]
    }
    
    df$PubChem <- NA
    df$drugbank <- NA
    df$chEBI <- NA
    df$Wikipedia <- NA
    if(!is.null(drug_item[['crossReferences']])){
      for(db in c('PubChem','drugbank','chEBI','Wikipedia')){
        if(!is.null(drug_item[['crossReferences']][[db]])){
          df[,db] <- 
            paste(drug_item[['crossReferences']][[db]], 
                  collapse = "|")
        }else{
          df[,db] <- NA
        }
      }
    }
    
    df <- df |>
      dplyr::rename(drug_pubchem_xref = PubChem,
                    drug_drugbank_xref = drugbank,
                    drug_chebi_xref = chEBI,
                    drug_wikipedia_xref = Wikipedia)
    
    all_drug_child_entries <- data.frame()
    if(!is.null(drug_item[['childChemblIds']])){
      for(n in 1:length(drug_item[['childChemblIds']])){
        drug_child_entry <- df
        drug_child_entry$child_molecule_chembl_id <-
          drug_item[['childChemblIds']][[n]]
        
        all_drug_child_entries <- all_drug_child_entries |>
          dplyr::bind_rows(drug_child_entry)
      }
    }
    
    if(nrow(all_drug_child_entries) == 0){
      all_drug_child_entries <- df
      all_drug_child_entries$child_molecule_chembl_id <- NA
    }
    
    
    if(!is.null(drug_item[['linkedTargets']])){
      if(drug_item[['linkedTargets']]$count > 0){
          j <- 1
          while(j <= nrow(all_drug_child_entries)){
            drug_target_entry <- all_drug_child_entries[j,]
            k <- 1
            while(k <= length(drug_item[['linkedTargets']]$rows)){
              drug_target_entry$target_ensembl_gene_id <- 
                drug_item[['linkedTargets']]$rows[[k]]
              OT_drugs <- OT_drugs |>
                dplyr::bind_rows(drug_target_entry)
              k <- k + 1
            }
            j <- j + 1
          }
      }else{
        df$target_ensembl_gene_id <- NA
        OT_drugs <- OT_drugs |>
          dplyr::bind_rows(df)
        
      }
    }else{
      df$target_ensembl_gene_id <- NA
      OT_drugs <- OT_drugs |>
        dplyr::bind_rows(df)
      
    }
    
    i <- i + 1
  }
  
  m <- m + 1
  close.connection(con)
}

OT_drugs <- OT_drugs |>
  dplyr::left_join(gene_xref) |>
  dplyr::select(
    target_genename, 
    target_symbol,
    target_ensembl_gene_id,
    drug_name,
    drug_synonyms,
    drug_tradenames,
    drug_type, 
    drug_description,
    molecule_chembl_id,
    child_molecule_chembl_id, 
    parent_molecule_chembl_id,
    dplyr::everything()
  ) |>
  dplyr::arrange(target_symbol)


####-- MECHANISM_OF_ACTION --####

basepath <- file.path(here::here(), "data", 
                      release, "mechanismOfAction")
json_files <- sort(
  list.files(basepath, pattern = ".json", 
             all.files = T, full.names = T))

OT_moa <- data.frame()
m <- 1
for(json_chunk_lines in json_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  con <- file(json_chunk_lines, "r")
  moa_data <- 
    jsonlite::stream_in(con, flatten = F, simplifyVector = F)
  
  i <- 1
  while(i <= length(moa_data)){
    moa_item <- moa_data[[i]]
    action_type <- NA
    if(!is.null(moa_item$actionType)){
      action_type <- moa_item$actionType
    }
    df <- data.frame(
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
    
    if(length(moa_item$references) > 0){
      all_sources <- c()
      for(n in 1:length(moa_item$references)){
        source <- moa_item$references[[n]]$source
        source_ids <- NA
        if(!is.null(moa_item$references[[n]]$ids)){
          source_ids <- paste(
            unique(moa_item$references[[n]]$ids),
            collapse="&"
          )
        }
        source_element <- paste(
          source, source_ids, sep="|"
        )
        all_sources <- 
          c(all_sources,
            source_element)
      }
      
      df$drug_moa_references <- paste(all_sources, collapse=",")
    }else{
      df$drug_moa_references <- NA
    }
      
    
    for(n in 1:length(moa_item$chemblIds)){
      molecule_chembl_id <- moa_item$chemblIds[[n]]
      moa_entry <- df
      moa_entry$molecule_chembl_id <- molecule_chembl_id
      if(length(moa_item$targets) > 0){
        for(k in 1:length(moa_item$targets)){
          ensembl_gene_id <- moa_item$targets[[k]]
          moa_entry2 <- moa_entry
          moa_entry2$target_ensembl_gene_id <- ensembl_gene_id
          OT_moa <- OT_moa |>
            dplyr::bind_rows(moa_entry2)
        }
      }else{
        moa_entry$target_ensembl_gene_id <- NA
        OT_moa <- OT_moa |>
          dplyr::bind_rows(moa_entry)
      }
     
    }
    i <- i + 1
  }
  close.connection(con)
  m <- m + 1
}



## Combine drugs and mechanism of action data

OT_moa_ensembl <- 
  OT_moa |> dplyr::filter(!is.na(target_ensembl_gene_id))
OT_moa_no_ensembl <- 
  OT_moa |> dplyr::filter(is.na(target_ensembl_gene_id))
OT_drugs_ensembl <- 
  OT_drugs |> dplyr::filter(!is.na(target_ensembl_gene_id))
OT_drugs_no_ensembl <- 
  OT_drugs |> dplyr::filter(is.na(target_ensembl_gene_id))
OT_drugs_no_ensembl$target_ensembl_gene_id <- NULL
OT_drugs_no_ensembl$target_symbol <- NULL
OT_drugs_no_ensembl$target_genename <- NULL
OT_drugs_no_ensembl$target_entrezgene <- NULL
OT_moa_no_ensembl$target_ensembl_gene_id <- NULL

set1 <- OT_drugs_ensembl |> 
  dplyr::full_join(OT_moa_ensembl) |>
  dplyr::filter(!is.na(target_symbol))
set2 <- OT_drugs_no_ensembl |> 
  dplyr::full_join(OT_moa_no_ensembl) |>
  dplyr::filter(!is.na(drug_name)) |>
  dplyr::arrange(drug_name)

OT_drugs_moa_all <- set1 |>
  dplyr::bind_rows(set2) 

OT_drugs_with_children <- OT_drugs_moa_all |>
  dplyr::filter(!is.na(child_molecule_chembl_id)) |>
  dplyr::group_by(across(-starts_with("child"))) |>
  dplyr::summarise(
    child_molecule_chembl_id = 
      paste(child_molecule_chembl_id,
            collapse="&"),
    .groups = "drop")

OT_drugs_without_children <- OT_drugs_moa_all |>
  dplyr::filter(is.na(child_molecule_chembl_id))

OT_drugs_moa_final <- OT_drugs_with_children |>
  dplyr::bind_rows(OT_drugs_without_children)

rm(set1)
rm(set2)
rm(OT_drugs_moa_all)
rm(OT_drugs_without_children)
rm(OT_drugs_with_children)
rm(OT_moa_ensembl)
rm(OT_moa_no_ensembl)
rm(OT_drugs_ensembl)
rm(OT_drugs_no_ensembl)


####-- INDICATIONS --####

basepath <- file.path(here::here(), "data", 
                      release, "indication")
json_files <- sort(
  list.files(basepath, pattern = ".json", 
             all.files = T, full.names = T))

OT_indication <- data.frame()

m <- 1
for(json_chunk_lines in json_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  con <- file(json_chunk_lines, "r")
  indication_data <- 
    jsonlite::stream_in(con, flatten = F, simplifyVector = F)
  i <- 1
  
  while(i <= length(indication_data)){
    indication_item <- indication_data[[i]]
    molecule_chembl_id <- NA
    if(!is.null(indication_item$id)){
      molecule_chembl_id <- indication_item$id
    }
    
    # if(molecule_chembl_id == "CHEMBL2165191"){
    #   cat(m,i,sep=" - ")
    # }
    
    df <- data.frame(
      'molecule_chembl_id' =  molecule_chembl_id,
      stringsAsFactors = F)
    
    
    approved_indications <- data.frame()
    if(length(indication_item$approvedIndications) > 0){
      for(n in 1:length(indication_item$approvedIndications)){
        df2 <- 
          data.frame('disease_efo_id' = 
                       indication_item$approvedIndications[[n]],
                     'drug_approved_indication' = TRUE,
                     stringsAsFactors = F)
        approved_indications <- dplyr::bind_rows(
          approved_indications, df2)
      }
    }
    
    OT_indication_entries <- data.frame()
    if(length(indication_item$indications) > 0){
      k <- 1
      while(k <= length(indication_item$indications)){
        indication <- indication_item$indications[[k]]
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
        if(length(indication$references) > 0){
          for(n in 1:length(indication$references)){
            ref <- indication$references[[n]]
            
            for(u in 1:length(ref$ids)){
              e <- indication_entry
              e$drug_clinical_source <- ref$source
              e$drug_clinical_id <- ref$ids[[u]]
              
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
        OT_indication_entries |> dplyr::left_join(
        approved_indications, by = "disease_efo_id"
      ) |>
        dplyr::mutate(drug_approved_indication =
                        dplyr::if_else(
                          is.na(drug_approved_indication),
                          as.logical(FALSE),
                          as.logical(drug_approved_indication)
                        ))
    }else{
      OT_indication_entries$drug_approved_indication <- FALSE
    }
    
    OT_indication <- OT_indication |>
      dplyr::bind_rows(OT_indication_entries)
    
    #cat(i,'\n')
    i <- i + 1
  }
  close.connection(con)
  m <- m + 1
}

OT_drugs_indication <- as.data.frame(
  OT_drugs_moa_final |> 
    dplyr::full_join(OT_indication) |>
    dplyr::select(drug_name, 
                  drug_synonyms,
                  drug_tradenames, 
                  drug_type,
                  drug_moa, 
                  drug_is_approved,
                  drug_description,
                  molecule_chembl_id,
                  child_molecule_chembl_id,
                  parent_molecule_chembl_id,
                  target_genename,
                  target_symbol,
                  target_entrezgene,
                  target_type,
                  target_ensembl_gene_id,
                  disease_efo_id,
                  disease_efo_label,
                  drug_max_phase_indication,
                  drug_max_ct_phase,
                  drug_approved_indication,
                  dplyr::everything())
)

drugs_no_indication <- as.data.frame(
  OT_drugs_indication |> 
    dplyr::filter(is.na(disease_efo_id) & is.na(drug_clinical_source)) |>
    tidyr::separate_rows(child_molecule_chembl_id, sep="&")
)

drugs_indication <- OT_drugs_indication |> 
  dplyr::filter(!is.na(disease_efo_id) & !is.na(molecule_chembl_id))

drugs_no_indication <- as.data.frame(
  drugs_no_indication |> 
    dplyr::anti_join(drugs_indication, 
                     by = c("child_molecule_chembl_id" = "molecule_chembl_id")) |>
    dplyr::group_by_at(dplyr::vars(-c(child_molecule_chembl_id))) |>
    dplyr::summarise(child_molecule_chembl_id =
                       paste(child_molecule_chembl_id,
                             collapse = "&"),
                     .groups = "drop") |>
    dplyr::mutate(child_molecule_chembl_id = dplyr::if_else(
      child_molecule_chembl_id == "NA",
      as.character(NA),
      as.character(child_molecule_chembl_id)
    )) |>
    dplyr::select(drug_name, 
                  drug_synonyms,
                  drug_tradenames, 
                  drug_type,
                  drug_moa, 
                  drug_is_approved,
                  drug_description,
                  molecule_chembl_id,
                  child_molecule_chembl_id,
                  parent_molecule_chembl_id,
                  target_genename,
                  target_symbol,
                  target_entrezgene,
                  target_type,
                  target_ensembl_gene_id,
                  disease_efo_id,
                  disease_efo_label,
                  drug_max_phase_indication,
                  drug_max_ct_phase,
                  drug_approved_indication,
                  dplyr::everything())
)

OT_drugs_indication_final <- 
  drugs_indication |>
  dplyr::bind_rows(drugs_no_indication) |>
  dplyr::arrange(target_symbol)

saveRDS(OT_drugs_indication_final, 
        file = file.path(
          "output",
          paste0(
            "opentargets_drugs_",
            release,".rds"))
)
        