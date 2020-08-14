#!/bin/bash

awsdir="s3://ndmg-data/fm2g-edgelists/"
outdir="/data/corr/m2g/fmri/"
mkdir -p $outdir
conndir="_mask_desikan_space-MNI152NLin6_res-2x2x2_mask_file_..m2g_atlases..atlases..label..Human..desikan_space-MNI152NLin6_res-2x2x2.nii.gz"

for dataset in $(aws s3 ls $awsdir | awk '{print $2}' | sed '$ d'); do
    mkdir -p $outdir$dataset
    dsdir=$outdir$dataset
    for file in $(aws s3 ls $awsdir$dataset$conndir/NEW/ | awk '{if ($3 >0) print $4}'); do
        aws s3 cp $awsdir$dataset$conndir/NEW/$file $dsdir
    done
done