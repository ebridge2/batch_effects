# Full Analysis

[Docker container with dependencies](https://hub.docker.com/r/neurodata/batch_effects)

You will want version 0.0.2.

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
2. Change the following lines:

```

datasets <- c("<batch-1>", "<batch-2>", ...) # whatever batch names you want to include
pheno.name <- "<pheno-name>"
```

and possibly more depending on whether you did the graph reading/phenotypic data reading step using my code.

3. I usually enter the docker container in a tmux session, and then use the following command:

```
docker run -ti --entrypoint /bin/bash -v /<inputs>/<dir>/:/data -v /<repo>/<dir>/batch_effects/:/base neurodata/batch_effects:0.0.2
```

4. cd and run scripts sequentially:

```
cd /base/causal/
Rscript causal_driver_pdcorr.R
Rscript edgewise_ana_subset.R  # this takes forever
Rscript edgewise_ana_all.R  # this does too
```

The second to last script takes about 3 hours on 40 cores with about 300 connectomes in the overlapping region, and the latter took somewhere around 7 days on 40 cores, with 3500 connectomes. You will need a distance matrix for each of the 116*(115)/2 features of the AAL connectomes, for n subjects, and the distance matrix has rooughly quadratic complexity, so if there are more than 3500 connectomes, I would expect fairly crazy runtimes. If this is the case we should talk about how to speed it up because I was not shooting for efficiency, but if there are like anywhere near 10k subjects in ABCD we might want to explore that.

5. Send me the files `pdcorr_outputs_fMRI_AAL.rds` and `outputs_edgewise_ss.rds` and `outputs_edgewise_all.rds` which will be in the directory `<repo-dir>/batch_effects/data/dcorr/`.
