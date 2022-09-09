## Scripts to pre-process JSON files from Open Targets Platform

This repository contains scripts (R) to establish global datasets from the Open Targets Platform, including:

-   gene-disease associations
-   drug-target associations
-   drug-disease indications
-   other functional target information (tractability data, function description, cancer hallmarks information etc.)

### Overview

-   Code to download/rsync raw JSON files from the Open Targets Platform and preprocess them is available in [code](code)

#### Procedure

1.  Download raw JSON files (using [code/rsync_json_data.sh](rsync_json_data.sh) towards a `data` folder
2.  Parse the raw JSON files using R scripts found in [code](code) to produce `.rds` files in an `output` folder
    -   `parse_opentargets_targets_json.R`
    -   `parse_opentargets_drugs_json.R`
    -   `parse_opentargets_associations_json.R`
