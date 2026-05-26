release <- '2026.03'

####---- TARGETS ----####

# wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.04/output/etl/json/targets

basepath <- file.path(here::here(), "data", 
                      release, "pharmacogenomics")
parquet_files <- sort(
  list.files(basepath, pattern = ".parquet", 
             all.files = T, full.names = T))

otp_pgx <- data.frame()

m <- 1
for(parquet_fname in parquet_files){
  cat(paste0('Chunk - ',m))
  cat('\n')
  pgx_data <- 
    arrow::read_parquet(parquet_fname)
  i <- 1
  
  while(i <= NROW(pgx_data)){
    pgx_item <- pgx_data[i,]
    
    df <- data.frame(
      'datasourceId' =  pgx_item$datasourceId,
      stringsAsFactors = F)
    for(e in c('datasourceVersion',
               'datatypeId',
               'directionality',
               'evidenceLevel',
               'genotype',
               'genotypeAnnotationText',
               'haplotypeFromSourceId',
               'haplotypeId',
               'variantRsId',
               'targetFromSourceId',
               'variantFunctionalConsequenceId',
               'genotypeId',
               'phenotypeText',
               'pgxCategory',
               'studyId',
               'isDirectTarget')){
      if(!is.null(pgx_item[[e]])){
        df[,e] <- pgx_item[[e]]
      }else{
        df[,e] <- NA
      }
    }
    
    df$literature <- NA
    if(!is.null(pgx_item$literature)){
      if(is.list(pgx_item$literature)){
        df$literature <- paste(
          unique(sort(pgx_item$literature[[1]])), collapse = "|")
      }
    }
    
    molecule_chemb_id <- c()
    drug_name <- c()
    if(is.list(pgx_item$drugs)){
      if(nrow(pgx_item$drugs[[1]]) > 0){
        for(k in 1:nrow(pgx_item$drugs[[1]])){
          if("drugId" %in% names(pgx_item$drugs[[1]])){
            molecule_chemb_id <- c(
              molecule_chemb_id, pgx_item$drugs[[1]][k,]$drugId)
          }
          if("drugFromSource" %in% names(pgx_item$drugs[[1]])){
            drug_name <- c(
              drug_name, pgx_item$drugs[[1]][k,]$drugFromSource)
          }
        }
      }
    }
    df$molecule_chembl_id <- paste(molecule_chemb_id, collapse = "|")
    df$drug_name <- paste(drug_name, collapse = "|")
    
    df <- janitor::clean_names(df)
    
    otp_pgx <- otp_pgx |>
      dplyr::bind_rows(df)
    
    i <- i + 1
  }
  m <- m + 1
    
}

saveRDS(otp_pgx, file = file.path(
  here::here(), "output", 
  paste0("opentargets_pgx_", release,".rds")))
