# Code

1. Get raw JSON data with `rsync_parquet_data.sh`

2. Run individual scripts to pull out drug information, disease associations etc.
   * `parse_opentargets_targets_parquet.R`
   * `parse_opentargets_drugs_parquet.R`
   * `parse_opentargets_associations_parquet.R`