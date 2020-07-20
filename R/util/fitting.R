## Signature fitting function for mmsig

#' @title convert_signature_definitions_mmsig
#'
#' @description
#' Re-format signature matrix in the mmsig format.
#'
#' @param sig_df String to split [data frame]
#'
#' @return mmsig formatted signatures data frame [data frame]
#'
#' @export
convert_signature_definitions_mmsig <- function(sig_df, reverse = FALSE){
    if (any(grepl(pattern = ">", x = colnames(sig_df)))) {
        sig_df <- data.frame(t(sig_df))
    }
    if (!reverse) {
        sig_df$Substitution.Type <- substr(rownames(sig_df), 3,5)
        sig_df$Trinucleotide <- gsub(pattern = "[^A-Z]", replacement = "", x = gsub(pattern = ">[A-Z]", replacement = "", x = rownames(sig_df)))
        sig_df_out <- sig_df %>% dplyr::select(Substitution.Type, Trinucleotide, everything())
    }
    
    if (reverse) {
        sig_df_out <- sig_df
        rownames(sig_df) <- paste0(substr(rownames(sig_df), 1, 1),
                                   "[",
                                   substr(rownames(sig_df), 2, 2),
                                   ">",
                                   substr(rownames(sig_df), 6, 6),
                                   "]",
                                   substr(rownames(sig_df), 3, 3))   
        sig_df_out <- sig_df
    }

    return(sig_df_out)
}

#' @title calculate_mutation_probabilities
#'
#' @description Calculates probability of a signature resulting in a specific mutation type
#'
#' @details Given a signatures data frame and a weights data frame, calculates the probability of each
#' signature generating each trinucleotide context in the sample.
#'
#' @param signatures_df trincucleotide contexts as rows and signatures as columns [data frame]
#' @param exposures_df samples as rows and signatures as columns [data frame]
#' @return mmsig output [data frame]
#' @export
calculate_mutation_probabilities <- function(signatures_df, exposures_df){
    
    if (any(grepl(pattern = ">", x = colnames(signatures_df)))) {
        signatures_df <- data.frame(t(signatures_df))
    }
    
    if (any(grepl(pattern = "SBS|Signature", x = colnames(exposures_df)))) {
        exposures_df <- t(exposures_df)
    }
    
    signatures_used <- as.matrix(signatures_df[,intersect(rownames(exposures_df), 
                                                          colnames(signatures_df))])
    exposures_df <- exposures_df[intersect(rownames(exposures_df), 
                                           colnames(signatures_df)),]
    
    # this will approximate the mutation matrix I input for each sample
    genomes <- signatures_used %*% exposures_df
    
    boop <- lapply(1:ncol(exposures_df), FUN = function(i){
        M <- genomes[,i,drop = TRUE]
        tmp_probs <- sweep(signatures_used, 2, exposures_df[,i], "*")
        probs <- data.frame(tmp_probs/M)
        probs$Sample.Names <- colnames(exposures_df)[i]
        probs$MutationTypes <- rownames(probs)
        return(probs)
    })
    
    result <- do.call(rbind, boop)
    result <- result %>% dplyr::select(Sample.Names, MutationTypes, everything())
    return(result)
    
}

fit_signatures = function(samples.muts,
                          consigts.defn,
                          sigt.profs, 
                          cos_sim_threshold,
                          dbg=FALSE) {
    "
    samples.muts = 96 classes mutational profile for samples
    consigts.defn = 96 classes mutational profile for mutational signature reference
    sigt.profs = which signatures to fit for each sample
    "
    max.em.iter=2000
    consigt.names <- colnames(consigts.defn)
    samples <- colnames(samples.muts)
    
    # Mutational signature fitting procedure for each infividual sample
    sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
    rownames(sigt.fraction) <- consigt.names
    colnames(sigt.fraction) <- samples
    
    for (j in 1:length(samples)) {
        spit(dbg, "sample number: %d", j)
        sample.mut.freqs = as.numeric(samples.muts[,j])
        sample.mut.freqs[is.na(sample.mut.freqs)] = 0
        sample.sigts <- unique(sigt.profs[[ samples[j] ]])
        sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)] 
        sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
        spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
        alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        spat(dbg, "alpha", alpha)
        sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
        sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]  
        
        if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
        
        spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
        reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs) 
        sample.cos.sim.meas <- MutationalPatterns::cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
        spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
        
        rem.alpha <- sampleAlpha                     # holds the final result
        rem.sample.consigts.defn <- sample.consigts.defn
        spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
        reducing = TRUE
        
        # Signature profile shrinkage by cosine similarity (removing signatures that are not necessary to explain profile)
        while (reducing) { 
            spat(dbg, "in the while, rem.alpha: ", rem.alpha)
            cosReduction <- NULL
            if(any(grepl(pattern = "SBS", names(rem.alpha)))){
                rem.names <- setdiff(names(rem.alpha),c("SBS1","SBS5"))
                if(length(rem.names) == 0){ ## Avoiding script crash when only SBS1 and 5 are present.
                    spit(dbg, "removed all signatures except SBS1 and SBS5: exiting while...")
                    break
                }
            }
            if(any(grepl(pattern = "Signature.", names(rem.alpha)))){
                rem.names <- setdiff(names(rem.alpha),c("Signature.1","Signature.5"))
                if(length(rem.names) == 0){ ## Avoiding script crash when only SBS1 and 5 are present.
                    spit(dbg, "removed all signatures except Signature.1 and Signature.5: exiting while...")
                    break
                }
            }
            
            for(c in rem.names){
                spit(dbg, "doing c: %s", c)
                red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c,drop=FALSE]
                red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
                red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
                red.cos.sim.meas <- MutationalPatterns::cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
                cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
                }
            names(cosReduction) <- rem.names
            if (min(cosReduction) < cos_sim_threshold) {
                spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
                rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction))),drop=FALSE]
                rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
                reducing = TRUE
                } 
            else {
                spit(dbg, "exiting while...")
                reducing = FALSE
                }
        }
        
        spit(dbg,"... while exited")
        rem.alpha.names <- names(rem.alpha)
        
        for (n in 1:length(consigt.names)) {
            if (consigt.names[n] %in% rem.alpha.names) {
                sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
                }
            else {
                sigt.fraction[n,j] <- 0
            }
        }
    }
    spat(dbg, "sigt.fraction", sigt.fraction)
    tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
    colsums.samples.muts <- colSums(samples.muts)
    sig <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)

    
    consigts.defn.tweaked <- convert_signature_definitions_mmsig(consigts.defn, reverse = TRUE)
    mut_prob <- calculate_mutation_probabilities(signatures_df = consigts.defn.tweaked, exposures_df = tdf.sigt.fraction)
    #mut_prob <- calculate_mutation_probabilities(signatures_df = consigts.defn, exposures_df = tdf.sigt.fraction)
    
    out <- list(sig, mut_prob)
    names(out) <- c("weights", "mutation_probability")
    
    return(out)
}



