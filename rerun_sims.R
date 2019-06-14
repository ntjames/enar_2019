# determine simulations that timed out on original run and rerun on ACCRE
library(tidyverse)
wdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/ENAR2019/code")

# continuous-binary sims
cb_dir <- file.path(wdir,"sims","cb")
cb_files<-dir(cb_dir)
started_cb <- str_subset(cb_files,".out") %>% str_replace(".out","") %>% str_replace("cbmod_sim_","")
comp_cb <- str_subset(cb_files,".RData") %>% str_replace(".RData","") %>% str_replace("cb_sim_","")
timedout_cb<-started_cb[!started_cb %in% comp_cb]

prev_timedout<-sort(as.numeric(timedout_cb))

# binary-binary sims
bb_dir <- file.path(wdir,"sims","bb")
bb_files <- dir(bb_dir)

bb_runs <-1:10800 # total number of sims
started_bb <- str_subset(bb_files,".out") %>% str_replace(".out","") %>% str_replace("bbmod_sim_","")
comp_bb <- str_subset(bb_files,".RData") %>% str_replace(".RData","") %>% str_replace("bb_sim_","")
timedout_bb<-started_bb[!started_bb %in% comp_bb] # started but no .RData output
never_started<-bb_runs[!bb_runs %in% started_bb] # never started (no .out file)

bb_sim_rerun<-sort(c(as.numeric(timedout_bb),as.numeric(never_started)))

# copy and paste this into array list
write.table(matrix(bb_sim_rerun,nrow=1), sep=",",
            row.names=FALSE, col.names=FALSE)

