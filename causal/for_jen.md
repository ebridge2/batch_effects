# Full Analysis

[Docker container with dependencies](https://hub.docker.com/r/neurodata/batch_effects)

You will want version 0.0.2.

```
docker pull neurodata/batch_effects:0.0.2
```

# Requirements

1. Connectomes stored as edge lists, with the following file organization:

```
<path>/<to>/<base>/AAL/<batch name>/sub-<sub-id>_ses-<ses-id>_<anything>.csv
```

It is assumed that all connectomes will be the AAL parcellation, appropriately numbered from 1 to 116.

2. A phenotypic file stored at:

```
<inputs>/<dir>/phenotypic/<pheno-name>.csv
```

Which is organized as follows according to the CoRR data legend:

```
SUBID,SITE,SESSION,SEX,AGE_AT_SCAN_1,HANDEDNESS,RETEST_DESIGN,RETEST_DURATION,RETEST_UNITS,PRECEDING_CONDITION,VISUAL_STIMULATION_CONDITION,RESTING_STATE_INSTRUCTION
3001,BMB_1,Baseline,2,25.13,R,1,#,#,0,2,Relax and keep your eyes open/closed.
3001,BMB_1,Retest_1,2,25.13,R,1,10,m,0,2,Relax and keep your eyes open/closed.
3002,BMB_1,Baseline,1,23.96,R,1,#,#,0,2,Relax and keep your eyes open/closed.
3002,BMB_1,Retest_1,1,23.96,R,1,10,m,0,2,Relax and keep your eyes open/closed.
3004,BMB_1,Baseline,2,31.15,R,1,#,#,0,2,Relax and keep your eyes open/closed.
3004,BMB_1,Retest_1,2,31.15,R,1,10,m,0,2,Relax and keep your eyes open/closed.
3006,BMB_1,Baseline,1,23,R,1,#,#,0,2,Relax and keep your eyes open/closed.
3006,BMB_1,Retest_1,1,23,R,1,10,m,0,2,Relax and keep your eyes open/closed.
3007,BMB_1,Baseline,2,43.53,R,1,#,#,0,2,Relax and keep your eyes open/closed.
...
```

The key to be aware of is Male = 2 and Female = 1 for the SEX column. 

# Usage

1. Clone this repo locally, to `<repo>/<dir>`.

2. Change the following lines of the file `causal_driver_pdcorr.R`:

```
# Line 23
pheno.name <- "<pheno-name>"
# Line 28
cohort <- "ABCD"
#Line 31
datasets <- c("<batch-1>", "<batch-2>", ...) # whatever batch names you want to include
```

and possibly more depending on whether you did the graph reading/phenotypic data reading step using my code.

Change the following lines in `edgewise_ana_ss.R`:

```
# Line 8
cohort <- "ABCD"
```

and `edgewise_ana_all.R`:

```
# Line 8
cohort <- "ABCD"
```

3. Run the following command to allow the docker container to save outputs:

```
chmod -R 707 /<repo>/<dir>/batch_effects/data/dcorr/
```

4. I usually enter the docker container in a tmux session, and then use the following command:

```
docker run -ti --entrypoint /bin/bash -v /<inputs>/<dir>/:/data -v /<repo>/<dir>/batch_effects/:/base neurodata/batch_effects:0.0.2
```

5. Be aware that the parsed covariate data should be in 1:1 order with the associated graphs. You can verify this by manually copying/pasting the first 122 lines of code from `causal_driver_pdcorr.R` into the terminal session, and ensuring manually that the object `cov.dat` has the subject/session/datasets organized properly. `cov.dat[1,]` should correspond to `gr.dat[1,]` and vice-versa for both objects. If this is not the case, the result will be nonsensical. I would check this manually (literally checking the split name of the output graph in `gr.dat[1,]` which is found in `spl.names[[1]]` and making sure the covariates are correct in `cov.dat[1,]`). I say this because I've never used this for anything but CoRR and this would be the most catastrophic failure that could potentially occur and go unnoticed if we aren't careful.

5. cd and run scripts sequentially. You can do this by either executing the following commands sequentially:

```
cd /base/causal/
Rscript causal_driver_pdcorr.R
Rscript edgewise_ana_subset.R  # this takes forever
Rscript edgewise_ana_all.R  # this does too
```

Or by copying and pasting the lines of each file into the terminal directly. I tend to do the latter, because if there are any issues before the "main" driver, I can catch them without losing the entire session.

The second to last script takes about 3 hours on 40 cores with about 300 connectomes in the overlapping region, and the latter took somewhere around 7 days on 40 cores, with 3500 connectomes. You will need a distance matrix for each of the 116*(115)/2 features of the AAL connectomes, for n subjects, and the distance matrix has rooughly quadratic complexity, so if there are more than 3500 connectomes, I would expect fairly crazy runtimes. If this is the case we should talk about how to speed it up because I was not shooting for efficiency, but if there are like anywhere near 10k subjects in ABCD we might want to explore that.

7. Send me the files `pdcorr_outputs_fMRI_AAL.rds` and `outputs_edgewise_ss.rds` and `outputs_edgewise_all.rds` which will be in the directory `<repo-dir>/batch_effects/data/dcorr/`.
