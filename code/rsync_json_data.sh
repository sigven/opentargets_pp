#!/bin/sh

release="22.06"

rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/targets .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/molecule .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/associationByOverallDirect .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/associationByOverallIndirect .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/associationByDatasourceDirect .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/associationByDatatypeDirect .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/associationByDatatypeIndirect .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/mechanismOfAction .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/indication .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/etl/json/diseases .