
goana.default <- function(de, universe = NULL, species = "Hs", prior.prob = NULL, covariate=NULL, plot=FALSE, ...)
  #	Gene ontology analysis of DE genes
  #	Gordon Smyth and Yifang Hu
  #	Created 20 June 2014.  Last modified 23 June 2016.
{
  #	Ensure de is a list
  if(!is.list(de)) de <- list(DE = de)
  
  #	Stop if components of de are not vectors
  if(!all(vapply(de,is.vector,TRUE))) stop("components of de should be vectors")
  
  #	Ensure gene IDs are of character mode
  de <- lapply(de, as.character)
  if(!is.null(universe)) universe <- as.character(universe)
  
  #	Ensure all gene sets have unique names
  nsets <- length(de)
  names(de) <- trimWhiteSpace(names(de))
  NAME <- names(de)
  i <- which(NAME == "" | is.na(NAME))
  NAME[i] <- paste0("DE",i)
  names(de) <- makeUnique(NAME)
  
  #	Fit trend in DE with respect to the covariate, combining all de lists
  if(!is.null(covariate)) {
    covariate <- as.numeric(covariate)
    if(length(covariate) != length(covariate)) stop("universe and covariate must have same length")
    isDE <- as.numeric(universe %in% unlist(de))
    o <- order(covariate)
    prior.prob <- covariate
    span <- approx(x=c(20,200),y=c(1,0.5),xout=sum(isDE),rule=2)$y
    prior.prob[o] <- tricubeMovingAverage(isDE[o],span=span)
    if(plot) barcodeplot(covariate, index=(isDE==1), worm=TRUE, span.worm=span)
  }
  
  #	Get access to package of GO terms
  suppressPackageStartupMessages(OK <- requireNamespace("GO.db",quietly=TRUE))
  if(!OK) stop("GO.db package required but not installed (or can't be loaded)")
  
  #	Get access to required annotation functions
  suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly=TRUE))
  if(!OK) stop("AnnotationDbi package required but not installed (or can't be loaded)")
  
  #	Load appropriate organism package
  orgPkg <- paste0("org.",species,".eg.db")
  suppressPackageStartupMessages(OK <- requireNamespace(orgPkg,quietly=TRUE))
  if(!OK) stop(orgPkg," package required but not not installed (or can't be loaded)")
  
  #	Get GO to Entrez Gene mappings
  obj <- paste0("org.",species,".egGO2ALLEGS")
  egGO2ALLEGS <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
  if(is.logical(egGO2ALLEGS)) stop("Can't find gene ontology mappings in package ",orgPkg)
  
  #	Convert gene-GOterm mappings to data.frame and remove duplicate entries
  if(is.null(universe)) {
    EG.GO <- AnnotationDbi::toTable(egGO2ALLEGS)
    d <- duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])
    EG.GO <- EG.GO[!d, ]
    universe <- unique(EG.GO$gene_id)
    universe <- as.character(universe)
  } else {
    
    universe <- as.character(universe)
    
    dup <- duplicated(universe)
    if(!is.null(prior.prob)) {
      if(length(prior.prob)!=length(universe)) stop("length(prior.prob) must equal length(universe)")
      prior.prob <- rowsum(prior.prob,group=universe,reorder=FALSE)
      n <- rowsum(rep_len(1L,length(universe)),group=universe,reorder=FALSE)
      prior.prob <- prior.prob/n
    }
    universe <- universe[!dup]
    
    m <- match(AnnotationDbi::Lkeys(egGO2ALLEGS),universe,0L)
    universe <- universe[m]
    if(!is.null(prior.prob)) prior.prob <- prior.prob[m]
    
    AnnotationDbi::Lkeys(egGO2ALLEGS) <- universe
    EG.GO <- AnnotationDbi::toTable(egGO2ALLEGS)
    d <- duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])
    EG.GO <- EG.GO[!d, ]
  }
  
  Total <- length(unique(EG.GO$gene_id))
  if(Total<1L) stop("No genes found in universe")
  
  #	Check prior probabilities
  if(!is.null(prior.prob)) {
    if(length(prior.prob)!=length(universe)) stop("length(prior.prob) must equal length(universe)")
  }
  
  #	Overlap with DE genes
  isDE <- lapply(de, function(x) EG.GO$gene_id %in% x)
  TotalDE <- lapply(isDE, function(x) length(unique(EG.GO$gene_id[x])))
  nDE <- length(isDE)
  
  if(length(prior.prob)) {
    #	Probability weight for each gene
    m <- match(EG.GO$gene_id, universe)
    PW2 <- list(prior.prob[m])
    X <- do.call(cbind, c(N=1, isDE, PW=PW2))
  } else
    X <- do.call(cbind, c(N=1, isDE))
  
  group <- paste(EG.GO$go_id, EG.GO$Ontology, sep=".")
  S <- rowsum(X, group=group, reorder=FALSE)
  
  P <- matrix(0, nrow = nrow(S), ncol = nsets)
  
  if(length(prior.prob)) {
    
    #		Calculate average prior prob for each set
    PW.ALL <- sum(prior.prob[universe %in% EG.GO$gene_id])
    AVE.PW <- S[,"PW"]/S[,"N"]
    W <- AVE.PW*(Total-S[,"N"])/(PW.ALL-S[,"N"]*AVE.PW)
    
    #		Wallenius' noncentral hypergeometric test
    if(!requireNamespace("BiasedUrn",quietly=TRUE)) stop("BiasedUrn package required but is not installed (or can't be loaded)")
    for(j in 1:nsets) for(i in 1:nrow(S)) 
      P[i,j] <- BiasedUrn::pWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i],lower.tail=FALSE) + BiasedUrn::dWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i])
    S <- S[,-ncol(S)]
    
  } else {
    
    #	Fisher's exact test
    for(j in 1:nsets)
      P[,j] <- phyper(q=S[,1+j]-0.5,m=TotalDE[[j]],n=Total-TotalDE[[j]], k=S[,"N"],lower.tail=FALSE)
    
  }
  
  #	Assemble output
  g <- strsplit2(rownames(S),split="\\.")
  TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,keys=g[,1],columns="TERM"))
  Results <- data.frame(GOID = TERM[[1]], Term = TERM[[2]], Ont = g[,2], S, P, stringsAsFactors=FALSE)
  rownames(Results) <- g[,1]
  
  #	Name p-value columns
  colnames(Results)[4+nsets+(1L:nsets)] <- paste0("P.", names(de))
  
  return(list(tbl = as_tibble(Results), gns = inner_join(tibble(gene_id = de$DE), tibble(EG.GO))))
}

