# Read in all chim1 ASV tables as RDS files, and
# combine all plates run through the dada2 pipeline, and then filter them to remove 
# low abundance reads, using the threshold of keeping only ASVs with >0.0001% per sample. 
# Then run the collapse no mismatch step on the combined ASVs tables.
# The collapse no mismatch step will compare reads and collapse any reads that are 100% identical
# into one read, using the more dominant of the two collapsed ASVs as the final ASV and adding their abundances.

#further testing required

library(dada2)
library(tidyverse)

rds.list <- list() 

rds.list <- list.files(pattern='RDS')

merger <- list() # create empty list

for (i in 1:65){
  merger[[i]] <-readRDS(paste0( rds.list[i]))
  }


seqtabnames<- as.data.frame(rds.list) %>% separate(rds.list, c('flow_id', NA), sep='seqtab', remove=F)

names(merger) <- seqtabnames$flow_id

seqCounts<-as.data.frame(unlist(lapply(merger, function(x) rowSums(x))))

seqCounts<- rownames_to_column(seqCounts, var='flowIDsample' )

seqCounts<- separate(seqCounts, flowIDsample, c('flow_id', 'code'), sep='\\.', remove=F)

colnames(seqCounts)[4]<-'nseqs'


#merger_nodupe<-lapply(merger, function(x) x[!rownames(x) %in% (to.drop %>% filter(flow_id %in% names(x)) %>% select(code)),])

st.all <- mergeSequenceTables(merger[[1]], merger[[2]],merger[[3]],merger[[4]],merger[[5]],merger[[6]],merger[[7]],merger[[8]],merger[[9]],
merger[[10]],merger[[11]],merger[[12]],merger[[13]],merger[[14]],merger[[15]],merger[[16]],merger[[17]],merger[[18]],merger[[19]],
merger[[20]],merger[[21]],merger[[22]],merger[[23]], merger[[24]], merger[[25]],merger[[26]],merger[[27]],merger[[28]],merger[[29]],
merger[[30]],merger[[31]],merger[[34]], merger[[35]], merger[[36]],merger[[37]],merger[[38]],merger[[39]],
merger[[40]],merger[[41]],merger[[42]],merger[[43]], merger[[44]], merger[[45]],merger[[46]],merger[[47]],merger[[48]],merger[[49]],
merger[[50]],merger[[51]],merger[[52]],merger[[53]], merger[[54]], merger[[55]],merger[[56]],merger[[57]],merger[[58]],merger[[59]], 
merger[[60]],merger[[61]],merger[[62]],merger[[63]], merger[[64]], merger[[65]])


st.all.4.ns<- st.all.4[,colSums(st.all.4) >1]
saveRDS(st.all.4.ns, "st.all.ns.chim4.RDS")
#  3289 930855
st.all.4.10<- st.all.4[,colSums(st.all.4) > 9]
saveRDS(st.all.4.10, "st.all.10.chim4.RDS")
#   3289 958133
st.all.4.100<- st.all.4[,colSums(st.all.4) > 99]
saveRDS(st.all.4.100, "st.all.100.chim4.RDS")


#### edited collapseNoMismatch script to assess the importance of collapsing identical sequences of slightly differing lengths

myCNMM<-function (seqtab, minOverlap = 200, orderBy = "abundance", identicalOnly = FALSE, 
    vec = TRUE, band = -1, verbose = FALSE) 
{
    dupes <- duplicated(colnames(st.all.nd.4))
    if (any(dupes)) {
        st <- st.all.nd.4[, !dupes, drop = FALSE]
        for (i in which(dupes)) {
            sq <- colnames(st.all.nd.4)[[i]]
            st[, sq] <- st[, sq] + st.all[, i]
        }
        st.all <- st
    }
    if (identicalOnly) {
        return(st.all)
    }
    unqs.srt <- sort(getUniques(st.all), decreasing = TRUE)
    seqs <- names(unqs.srt)[1:10000] # This addition limits the loop to the top x sequences by abundance
    seqs.out <- character(0)
    collapsed <- matrix(0L, nrow = nrow(seqtab), ncol = ncol(seqtab))
    colnames(collapsed) <- colnames(seqtab)
    rownames(collapsed) <- rownames(seqtab)
    for (query in seqs) {
        added = FALSE
        prefix <- substr(query, 1, minOverlap)
        for (ref in seqs.out) {
            prefix.ref <- substr(ref, 1, minOverlap)
            if (grepl(prefix, ref, fixed = TRUE) || grepl(prefix.ref, 
                query, fixed = TRUE)) {
                if (nwhamming(query, ref, vec = TRUE, band = -) == 
                  0) {
                  collapsed[, ref] <- collapsed[, ref] + seqtab[, 
                    query]
                  added = TRUE
                  break
                }
            }
        }
        if (!added) {
            collapsed[, query] <- seqtab[, query]
            seqs.out <- c(seqs.out, query)
        }
    }
    collapsed <- collapsed[, colnames(collapsed) %in% seqs.out, 
        drop = FALSE]
    if (!is.null(orderBy)) {
        if (orderBy == "abundance") {
            collapsed <- collapsed[, order(colSums(collapsed), 
                decreasing = TRUE), drop = FALSE]
        }
        else if (orderBy == "nsamples") {
            collapsed <- collapsed[, order(colSums(collapsed > 
                0), decreasing = TRUE), drop = FALSE]
        }
    }
    collapsed <- collapsed[, order(colSums(collapsed), decreasing = TRUE), 
        drop = FALSE]
    if (verbose) 
        message("Output ", ncol(collapsed), " collapsed sequences out of ", 
            ncol(seqtab), " input sequences.")
    collapsed
}

col.st <- collapseNoMismatch(st.all, minOverlap = 20,orderBy="abundance", verbose = T)
saveRDS(col.st, "/shared/c3/bio_db/BPA/amplicons/16s/seqtab/collapsed_seqtab.rds")

***

alternate configuration for Euk files



# read in all sequence tables from plates (there are 24)
merger <- list() # create empty list
# read in data with loop
for(i in seq_along(rds.list)){
merger[[i]] <- readRDS(paste0(rds.list[[i]]))
}
# merge all tables
st.all <- mergeSequenceTables(merger[[1]],merger[[2]],merger[[3]],merger[[4]],merger[[5]],merger[[6]],merger[[7]],merger[[8]],merger[[9]],
merger[[10]],merger[[11]],merger[[12]],merger[[13]],merger[[14]],merger[[15]],merger[[16]],merger[[17]],merger[[18]],merger[[19]],
merger[[20]],merger[[21]],merger[[22]],merger[[24]], merger[[25]],merger[[26]],merger[[27]],merger[[28]],merger[[29]],
merger[[30]],merger[[31]],merger[[32]])



saveRDS(st.all, "./Objects/all_seqtab.rds") # save
# collapse ASVs with identical sequence
col.st <- collapseNoMismatch(st.all, minOverlap = 20, verbose = T)
saveRDS(col.st, "./Objects/collapsed_seqtab.rds")


st.all.ns<- st.all[,colSums(st.all) !=1]
saveRDS(st.all.ns, "st.all.ns.RDS")
#  3289 930855
st.all.10<- st.all[,colSums(st.all) !=10]
saveRDS(st.all.10, "st.all.10.RDS")
#   3289 958133
st.all.100<- st.all[,colSums(st.all) !=100]
saveRDS(st.all.100, "st.all.100.RDS")
saveRDS(st.all, "seqtab_v1.rds")










myCNMM<-function (seqtab, minOverlap = 200, orderBy = "abundance", identicalOnly = FALSE, 
    vec = TRUE, band = -1, verbose = FALSE) 
{
    dupes <- duplicated(colnames(seqtab))
    if (any(dupes)) {
        st <- seqtab[, !dupes, drop = FALSE]
        for (i in which(dupes)) {
            sq <- colnames(seqtab)[[i]]
            st[, sq] <- st[, sq] + seqtab[, i]
        }
        seqtab <- st
    }
    if (identicalOnly) {
        return(seqtab)
    }
    unqs.srt <- sort(getUniques(seqtab), decreasing = TRUE)
    seqs <- names(unqs.srt)
    seqs.out <- character(0)
    collapsed <- matrix(0L, nrow = nrow(seqtab), ncol = ncol(seqtab))
    colnames(collapsed) <- colnames(seqtab)
    rownames(collapsed) <- rownames(seqtab)
    for (query in seqs[1:100000]) {
        added = FALSE
        prefix <- substr(query, 1, minOverlap)
        for (ref in seqs.out) {
            prefix.ref <- substr(ref, 1, minOverlap)
            if (grepl(prefix, ref, fixed = TRUE) || grepl(prefix.ref, 
                query, fixed = TRUE)) {
                if (nwhamming(query, ref, vec = vec, band = band) == 
                  0) {
                  collapsed[, ref] <- collapsed[, ref] + seqtab[, 
                    query]
                  added = TRUE
                  break
                }
            }
        }
        if (!added) {
            collapsed[, query] <- seqtab[, query]
            seqs.out <- c(seqs.out, query)
            	 print(length(seqs.out), progress.bar=T)
        }
    }
    collapsed <- collapsed[, colnames(collapsed) %in% seqs.out, 
        drop = FALSE]
    if (!is.null(orderBy)) {
        if (orderBy == "abundance") {
            collapsed <- collapsed[, order(colSums(collapsed), 
                decreasing = TRUE), drop = FALSE]
        }
        else if (orderBy == "nsamples") {
            collapsed <- collapsed[, order(colSums(collapsed > 
                0), decreasing = TRUE), drop = FALSE]
        }
    }
    collapsed <- collapsed[, order(colSums(collapsed), decreasing = TRUE), 
        drop = FALSE]
    if (verbose) 
        message("Output ", ncol(collapsed), " collapsed sequences out of ", 
            ncol(seqtab), " input sequences.")
    collapsed
}

