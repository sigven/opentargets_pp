#!/usr/local/bin/Rscript

release <- '2026.06'

OT_datasets <- list()
OT_datasets[['molecule']] <- data.frame()
OT_datasets[['drug_moa']] <- data.frame()
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


####-- MOLECULES --####
options(arrow.skip_nul = TRUE)
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
      'molecule_chembl_id' = molecule_item$id
    )

    # Extract basic fields
    for(e in c('name', 'drugType', 'canonicalSmiles', 'inchiKey',
               'maximumClinicalStage', 'description')){
      df[[e]] <- molecule_item[[e]] %||% NA
    }

    df <- df |>
      dplyr::rename(
        drug_name = name,
        drug_inchi_key = inchiKey,
        drug_canonical_smiles = canonicalSmiles,
        drug_clinical_stage_max = maximumClinicalStage,
        drug_type = drugType,
        drug_description = description)

    # Extract synonyms and trade names
    for(e in c('synonyms','tradeNames')){
      if(!is.null(molecule_item[[e]]) &&
         tibble::is_tibble(molecule_item[[e]][[1]]) &&
         NROW(molecule_item[[e]][[1]]) > 0){
        df[[paste0('drug_', tolower(e))]] <-
          paste(unique(sort(molecule_item[[e]][[1]]$label)), collapse = "|")
      } else {
        df[[paste0('drug_', tolower(e))]] <- NA
      }
    }

    # Parent ID
    df$parent_molecule_chembl_id <- molecule_item[['parentId']] %||% NA
    df$drug_drugbank_xref <- NA
    df$drug_ema_xref <- NA
    df$drug_dailymed_xref <- NA

    # Cross-references (only create columns that have data)
    if(!is.null(molecule_item[['crossReferences']]) &&
       is.data.frame(molecule_item[['crossReferences']][[1]]) &&
       nrow(molecule_item[['crossReferences']][[1]]) > 0) {

      drug_df <- molecule_item[['crossReferences']][[1]]

      xref_data <- drug_df |>
        dplyr::filter(source %in% c('drugbank', 'DailyMed', 'EMA')) |>
        dplyr::mutate(
          id_string = sapply(ids, function(x) paste(x[[1]], collapse = "|"))
        ) |>
        dplyr::distinct(source, .keep_all = TRUE) |>
        dplyr::mutate(source = paste0("drug_", tolower(source), "_xref")) |>
        dplyr::select(source, id_string) |>
        tidyr::pivot_wider(
          names_from = source,
          values_from = id_string
        )

      for(c in c('drug_drugbank_xref',
                   'drug_ema_xref',
                   'drug_dailymed_xref')) {
        if(c %in% colnames(xref_data)) {
          df[[c]] <- xref_data[[c]]
        }
      }

    }

    OT_datasets[['molecule']] <- OT_datasets[['molecule']] |>
      dplyr::bind_rows(df) |>
      dplyr::distinct()

    i <- i + 1
  }

  m <- m + 1
}


phenOncoX_maps <-
  phenOncoX::get_aux_maps(
    cache_dir = "data"
  )

clinical_indication <-
  arrow::read_parquet(
    file.path(
      here::here(), "data",
      release,
      "clinical_indication",
      "clinical_indication.parquet"))

clinical_target <-
  tibble::tibble(
    arrow::read_parquet(
      file.path(
        here::here(), "data",
        release,
        "clinical_target",
        "clinical_target.parquet")) |>
      dplyr::select(-c("diseases")) |>
      dplyr::rename(clinicalReportId = clinicalReportIds) |>
      dplyr::select(-c("id")) |>
      tidyr::unnest(clinicalReportId)
  )

options(arrow.skip_nul = TRUE)

# Keep drugId here so target-free drugs can still be recovered via
# clinical_report as the join anchor below.
clinical_report <-
  tibble::tibble(
    arrow::read_parquet(
      file.path(
        here::here(), "data",
        release,
        "clinical_report",
        "clinical_report.parquet")) |>
      dplyr::select(
        id,
        diseases,
        source,
        type,
        url,
        drugs,
        hasExpertReview,
        clinicalStage) |>
      dplyr::rename(
        clinical_report_type = type,
        clinical_report_source = source,
        clinical_report_expert_review = hasExpertReview,
        clinical_report_url = url,
        drug_clinical_stage_indication = clinicalStage) |>
      dplyr::mutate(
        across(where(is.character), 
               ~stringr::str_remove_all(., "\\x00"))) |>

      tidyr::unnest(diseases) |>
      tidyr::unnest(drugs) |>
      dplyr::filter(!(stringr::str_detect(.data$drugFromSource," \\+ "))) |>
      dplyr::select(-c("drugFromSource")) |>
      dplyr::rename(
        clinicalReportId = id,
        molecule_chembl_id = drugId) |>
      dplyr::distinct()
  )

# Anchor on clinical_report so drugs without a target entry in
# clinical_target are still present (with NA target columns).
result <-
  as.data.frame(
    clinical_report |>
      dplyr::left_join(
        clinical_target |>
          dplyr::rename(
            molecule_chembl_id = drugId,
            target_ensembl_gene_id = targetId
          ) |>
          dplyr::select(-maxClinicalStage),
        by = c("clinicalReportId", "molecule_chembl_id"),
        relationship = "many-to-many"
      )
  ) |>
  dplyr::rename(
    disease_efo_id = diseaseId,
    drug_clinical_id = clinicalReportId
  ) |>
  dplyr::left_join(
    OT_datasets$drug_moa |>
      dplyr::select(
        molecule_chembl_id,
        target_ensembl_gene_id,
        target_type,
        drug_action_type,
        drug_moa,
        drug_moa_references),
    by = c("molecule_chembl_id",
           "target_ensembl_gene_id"),
    relationship = "many-to-many"
  ) |>
  # Join molecule on molecule_chembl_id only so that target-free drugs
  # (where maxClinicalStage from clinical_target is NA) still match.
  # drug_clinical_stage_max is sourced exclusively from the molecule table.
  dplyr::left_join(
    OT_datasets$molecule,
    by = "molecule_chembl_id",
    relationship = "many-to-many"
  ) |>
  dplyr::left_join(
    dplyr::select(
      OT_datasets$drug_target,
      c("target_ensembl_gene_id",
        "target_symbol",
        #"target_biotype",
        "target_name",
        "target_entrezgene")),
    by = "target_ensembl_gene_id",
    relationship = "many-to-many"
  ) |>
  dplyr::rename(
    target_genename = target_name
  ) |>
  dplyr::filter(
    !is.na(disease_efo_id) &
      !is.na(molecule_chembl_id)
  ) |>
  dplyr::mutate(
    disease_efo_id = stringr::str_replace(
      .data$disease_efo_id, "_",":")
  ) |>
  dplyr::left_join(
    phenOncoX_maps$records$efo$efo2name,
    by = c("disease_efo_id" = "efo_id")
  ) |>
  dplyr::rename(
    disease_efo_label = efo_name
  ) |>
  dplyr::filter(!is.na(disease_efo_label))

# For drugs that have at least one target-annotated record, drop the
# NA-target rows produced by clinical_report entries that lacked a
# corresponding clinical_target entry. This prevents mixed target /
# no-target rows for the same drug (e.g. sotorasib) while preserving
# genuinely target-free drugs (e.g. cisplatin).
drugs_with_targets <- result |>
  dplyr::filter(!is.na(target_ensembl_gene_id)) |>
  dplyr::distinct(molecule_chembl_id)

result <- result |>
  dplyr::filter(
    !is.na(target_ensembl_gene_id) |
      !(molecule_chembl_id %in% drugs_with_targets$molecule_chembl_id)
  )

## consider only EMA/FDA/DailyMed approved indications with
## expert review for now
approved_drug_indications <- as.data.frame(
  result |>
    dplyr::filter(
      drug_clinical_stage_max == "APPROVAL",
      drug_clinical_stage_indication == "APPROVAL") |>
    dplyr::group_by(
      drug_name,
      drug_clinical_stage_indication,
      molecule_chembl_id,
      parent_molecule_chembl_id,
      disease_efo_id,
      disease_efo_label) |>
    dplyr::summarise(
      clinical_report_source = paste(
        unique(sort(clinical_report_source)), collapse = "|"),
      .groups = "drop_last"
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      stringr::str_detect(
        .data$clinical_report_source, "EMA|FDA|DailyMed")
    ) |>
    dplyr::select(
      c("molecule_chembl_id",
        "parent_molecule_chembl_id",
        "disease_efo_id")
    ) |>
    dplyr::mutate(
      drug_approved_indication = TRUE) |>

    ## propagate the approved indication status also
    ## to the parent molecule, if applicable
    dplyr::mutate(molecule_chembl_id = paste(
      molecule_chembl_id, parent_molecule_chembl_id, sep="|"
    )) |>
    dplyr::select(-c("parent_molecule_chembl_id")) |>
    tidyr::separate_rows(molecule_chembl_id, sep="\\|") |>
    dplyr::distinct()
)


## consider only EMA/FDA/DailyMed approved indications with
## expert review for now
approved_drugs <- as.data.frame(
  result |>
    dplyr::filter(drug_clinical_stage_max == "APPROVAL") |>
    dplyr::group_by(molecule_chembl_id,
                    parent_molecule_chembl_id) |>
    dplyr::summarise(
      clinical_report_source = paste(
        unique(sort(clinical_report_source)), collapse = "|"),
      .groups = "drop_last"
    ) |>
    dplyr::filter(
      stringr::str_detect(
        .data$clinical_report_source, "EMA|FDA|DailyMed")
    ) |>
    dplyr::select(
      "molecule_chembl_id",
      "parent_molecule_chembl_id"
    ) |>
    dplyr::mutate(drug_is_approved = TRUE) |>

    ## propagate the approved indication status also
    ## to the parent molecule, if applicable
    dplyr::mutate(molecule_chembl_id = paste(
      molecule_chembl_id, parent_molecule_chembl_id, sep="|"
    )) |>
    dplyr::select(-c("parent_molecule_chembl_id")) |>
    tidyr::separate_rows(molecule_chembl_id, sep="\\|") |>
    dplyr::distinct()
)


OT_drugs_complete <- result |>
  dplyr::left_join(
    approved_drug_indications,
    by = c("molecule_chembl_id", "disease_efo_id"),
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(drug_approved_indication = dplyr::if_else(
    is.na(drug_approved_indication),
    as.logical(FALSE),
    as.logical(drug_approved_indication)
  )) |>
  dplyr::left_join(
    approved_drugs,
    by = "molecule_chembl_id",
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(drug_is_approved = dplyr::if_else(
    is.na(drug_is_approved),
    as.logical(FALSE),
    as.logical(drug_is_approved)
  )) |>
  dplyr::select(
    molecule_chembl_id,
    parent_molecule_chembl_id,
    clinical_report_expert_review,
    clinical_report_source,
    clinical_report_type,
    drug_name,
    drug_is_approved,
    drug_approved_indication,
    drug_clinical_id,
    drug_clinical_stage_max,
    drug_clinical_stage_indication,
    drug_description,
    drug_action_type,
    drug_moa,
    drug_moa_references,
    drug_type,
    drug_canonical_smiles,
    drug_inchi_key,
    drug_tradenames,
    drug_synonyms,
    drug_drugbank_xref,
    drug_ema_xref,
    drug_dailymed_xref,
    target_ensembl_gene_id,
    target_symbol,
    target_type,
    target_genename,
    target_entrezgene,
    disease_efo_id,
    disease_efo_label
  )

saveRDS(OT_drugs_complete,
        file = file.path(
          "output",
          paste0(
            "opentargets_drugs_",
            release,".rds"))
)
