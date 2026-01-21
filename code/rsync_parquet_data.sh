#!/bin/sh

release="25.12"

rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/reactome .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/go .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/association_overall_indirect .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/association_overall_direct .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/association_by_datasource_direct .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/association_by_datatype_direct .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/association_by_datatype_indirect .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/drug_mechanism_of_action .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/drug_indication .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/drug_molecule .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/known_drug .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/disease .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/disease_phenotype .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/expression .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/target .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/target_essentiality .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/target_prioritisation .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/pharmacogenomics .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/$release/output/literature .
