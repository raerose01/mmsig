## Bootstrapping function for mmsig

bootstrap_fit_signatures <- function(samples.muts, 
                                     consigts.defn,
                                     cos_sim_threshold,
                                     sigt.profs=sigt.profs, 
                                     iterations = 1000){
    
    "
    Bootstrapping function for mm_fit_signatures
    Draw b = iterations mutational profiles for each sample from its multinomial distribution
    Performs signature fitting independently for each mutational profile
    For each mutational signature, returns its relative contribution as point estimate and bootstrapping mean with 95 % CI
    Can take a sample.sigt.profs argument that is passed directly to mm_fit_signatures. 
    "
    
    # Setup
    samples <- names(samples.muts)
    classes <- row.names(samples.muts)
    
    # List to populate with signatures
    mutSigs <- list()
    mutProbs <- list()
    
    # Setup progress bar
    pbar <- plyr::create_progress_bar('text')
    pbar$init(length(samples))
    
    for(i in 1:length(samples)){
        # Loop through samples, generating a data frame of signature contributions for each
        sub <- as.integer(samples.muts[classes,i])
        total <- sum(sub)
        # sample new 96-classes profiles from the multinomial distribution
        bootMat <- data.frame(rmultinom(n = iterations, size = total, prob = sub/total))
        row.names(bootMat) <- classes
        
        # prepare the signatures to fit for each sample
        sig.prof <- list()
        
        for(s in 1:ncol(bootMat)){
            sig.prof[[names(bootMat)[s]]] <- sigt.profs[[samples[i]]]
        }
        
        ### Run mmSig
        sig_out <- fit_signatures(samples.muts=bootMat, 
                                  consigts.defn=consigts.defn,
                                  sigt.profs=sig.prof, 
                                  cos_sim_threshold=cos_sim_threshold,
                                  dbg=FALSE)
        
        mutSigs[[i]] <- sig_out[[1]]
        mutProbs[[i]] <- sig_out[[2]]
        pbar$step()
    }
    
    names(mutSigs) <-  names(mutProbs) <- samples

    # Generate final summary data frame
    
    ## Summary statistics
    my_summary <- function(x){
        c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
    }
    
    mutSigsSummary <- list()
    for(i in 1:length(mutSigs)){
        s <- names(mutSigs)[i]
        temp <- mutSigs[[i]]
        out <- data.frame(t(sapply(temp[names(temp) != "mutations"], my_summary)))
        names(out) <- c('mean', 'CI025', 'CI975')
        out$signature <- row.names(out)
        out$sample <- s
        out <- out[c('sample', 'signature', 'mean', 'CI025', 'CI975')]
        mutSigsSummary[[i]] <- out
    }
    
    mutSigsSummary <- bind_rows(mutSigsSummary)
    
    
    mutProbsSummary <- list()
    for(i in 1:length(mutProbs)){
        s <- names(mutProbs)[i]
        temp <- mutProbs[[i]]
        boop <- lapply(unique(temp$MutationTypes), FUN = function(x){
            tmp <- temp[which(temp$MutationTypes == x),]
            out <- data.frame(t(sapply(tmp[!names(tmp) %in% c("Sample.Names", "MutationTypes")], my_summary)))
            names(out) <- c('mean', 'CI025', 'CI975')
            out$sample <- s
            out$MutationTypes <- x
            out$signature <- rownames(out)
            return(out)
        })
        #out <- out[c('sample', 'signature', 'mean', 'CI025', 'CI975')]
        out <- do.call(rbind, boop)
        rownames(out) <- NULL
        mutProbsSummary[[i]] <- out
    }
    
    names(mutProbsSummary) <- names(mutProbs)
    
    final_out <- list(mutSigsSummary, mutProbsSummary)
    names(final_out) <- c("mutSigsSummary", "mutProbsSummary")
    return(final_out)
}
