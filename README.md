## Scripts to download and preprocess data from Open Targets Platform

This repository contains scripts (in R) to establish global datasets from the Open Targets Platform, including:

-   gene-disease associations
-   drug-target associations
-   drug-disease indications
-   other functional target information (tractability data, function description, cancer hallmarks information etc.)

### Overview

-   Code to download/rsync raw [Parquet files](https://platform.opentargets.org/downloads) from the Open Targets Platform and 
preprocess them (pull out key information, and making them available as flat text records) is available in [code](code)

#### Overall procedure

1.  Download raw Parquet files (using [code/rsync_parquet_data.sh](rsync_parquet_data.sh) towards a `data` folder
2.  Parse the raw Parquet files using R scripts found in [code](code) to produce `.rds` files in an `output` folder
    -   `parse_opentargets_targets_json.R`
    -   `parse_opentargets_drugs_json.R`
    -   `parse_opentargets_associations_json.R`
