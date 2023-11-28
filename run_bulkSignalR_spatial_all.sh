#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh

conda activate r43_spatial


Rscript run_bulkSignalR_spatial.R YS15-33091
Rscript run_bulkSignalR_spatial.R YS17-23473
Rscript run_bulkSignalR_spatial.R YS17-26136
Rscript run_bulkSignalR_spatial.R YS19-34318
Rscript run_bulkSignalR_spatial.R YS20-17093
Rscript run_bulkSignalR_spatial.R YS22-14579
Rscript run_bulkSignalR_spatial.R YS16-23837
Rscript run_bulkSignalR_spatial.R YS21-10131


conda deactivate

