#!/bin/sh

module load R/3.6.0

# create batchtools.lsf.tmpl in the working directory
cat <<'TEMPLATE' > batchtools.lsf.tmpl
## Default resources can be set in your .batchtools.conf.R by defining the variable
## 'default.resources' as a named list.

#BSUB-J <%= job.name %>                             # Name of the job
#BSUB-o <%= log.file %>                             # Output is sent to logfile, stdout + stderr by default
#BSUB-q <%= resources$queue %>                      # Job queue
#BSUB-W <%= round(resources$walltime / 60, 1) %>    # Walltime (LSF requires minutes, batchtools uses seconds)
#BSUB-M <%= resources$memory %>                     # Memory requirements in KBytes; depends on setting LSF_UNIT_FOR_LIMITS in lsf.conf
#BSUB-P <%= resources$allocation %>

## Export value of DEBUGME environemnt var to slave
export DEBUGME=<%= Sys.getenv("DEBUGME") %>

Rscript -e 'batchtools::doJobCollection("<%= uri %>")'
TEMPLATE


# create master script
cat <<'LMROBSCRIPT' > mvb_master.R
.libPaths(c("/hpc/users/heissj01/.Rlib36",.libPaths()))
library(future.batchtools)
library(listenv)
library(data.table)

set.seed(1220399177L)
SEED = .Random.seed

load("transfer_phth.Rdata")

h = plan(batchtools_lsf,template="batchtools.lsf.tmpl",resources=
list(queue="premium",walltime=3600*5,memory=500,allocation="acc_envepi"))

n_cores = 200

# break the matrix into chunks
beta = t(beta)

chunks = ncol(beta)
chunks = split(1:chunks, rep(1:n_cores,each=ceiling(chunks/n_cores))[1:chunks])


pvals = listenv()

for(PHTH in exposures){

 	FRML = update(frml,paste0("~.+",PHTH))
  
 	pp = listenv()

	## futures
	for(core in 1:n_cores){

	  pp[[ core ]] %<-% {

			library(robustbase)
			library(purrr)

			ff = function(cpg){
				if(sum(!is.na(cpg)) < 100) return(NA_real_)
				PRGS$cpg = cpg[PRGS$j]
				m = lmrob(FRML,PRGS,control=lmrob.control(refine.tol=1e-4,setting="KS2014",seed=SEED))
				coef(summary(m))[PHTH,4]
			}

			ff = possibly(ff,otherwise=NA_real_)

			apply(CHNK,2,ff)


		} %globals% list(
			 PHTH = PHTH
			,PRGS = prgs
			,FRML = FRML
			,CHNK = beta[,chunks[[core]]]
			,SEED = SEED
		) %label% paste0(PHTH,"+",core)
	}

	pvals[[ PHTH ]] = pp
}

pvals = as.list(pvals)
pvals = lapply(pvals,unlist)

saveRDS(pvals,file="pvals_phth.rds")

LMROBSCRIPT

bsub -J mstr -P acc_envepi -q premium -W 1000 -M 20000 -e mstr.err -o mstr.out "Rscript --vanilla mvb_master.R"