# Infer Logics between triplets from gene_expression eset datasets
# Colbert Sesanker Jan 2013

# Load Packages
library(MASS)
source("partial.R")
tryCatch(
    library(Biobase), 
    error=function(e) {
    message(e)
    message("Now installing Biobase")
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biobase")
    }
) 
geneFile="GeneData_ANRbWBT__Oct_20_2009.Robj"
timePts = 1:7
times = length(timePts)
load(geneFile)
channel=channelNames(GeneData)

# This corrects the discrepancy between featureNames in featureData and assayData
featureNames(featureData(GeneData)) <- featureNames(assayData(GeneData)) 
# Functions for extracting p-values and r2 values
gene_pair_correlation <- function(target, source){
	modelResults <- NA
	modelResults <- try(lm(target  ~ source), silent=T)
	return(modelResults)
}

extract_r_squared <- function(target, source){
	output <- NA	
    cor    <- gene_pair_correlation(target, source)
    output <- if (class(cor) == "lm") summary(cor)$r.squared	    
	return(output)
}

p_value <- function(model){
    coef  <- summary(model)$fstatistic
    p_val <- pf(coef[1],coef[2], coef[3], lower.tail=F)
    return(p_val)
}

extract_p_val <- function(target, source){
	p_val <- NA	
    cor   <- gene_pair_correlation(target, source)
    if (class(cor) == "lm") {
      p_val <- p_value(cor)
    }   
	return(p_val)
}

# Symmetric Logic Functions
AND <- function(source_1, source_2){
    return(unlist(Map("*", source_1, source_2)))
}

OR <- function(source_1, source_2){
    return(unlist(Map("+", source_1, source_2)))
}

MIN <- function(source_1, source_2){
    return(unlist(Map("min", source_1, source_2)))
}

MAX <- function(source_1, source_2){
    return(unlist(Map("max", source_1, source_2)))
}

logics <- c("AND"= AND, "OR" = OR, "MIN" = MIN, "MAX" = MAX)
logic_indices <- list("AND"= 1, "OR" = 2, "MIN" = 3, "MAX" = 4)
# Various measures of correlation significance
significance_measures <- list("p-val" = extract_p_val, "R2"= extract_r_squared)

test_logic <- function(A, B, C, significance_measure="R2", sample=3) {
  sm <- significance_measures[[significance_measure]]
  GeneData <- GeneData[,pData(phenoData(GeneData))$Experiment=="TBC" & pData(phenoData(GeneData))$DNA_type=="cDNA"] 
  logicCorrelations <- matrix(vector(mode="list", length=times), nrow=times, ncol=length(logics))  
	  for (time in timePts){	
        print(paste("time =", time))
        GeneDataTime <- GeneData[,pData(phenoData(GeneData))$TimePt %in% time] #  all genes, all families for a given time point
        source_1 <- unlist(assayData(GeneDataTime)[[channel]][A,]); pts <- length(source_1); 
        sp <- function(x) spline(1:pts, x, n= sample*pts)$y # cubic spline: point density increased by factor of 'sample'
        source_1 <- sp(source_1)
        source_2 <- unlist(assayData(GeneDataTime)[[channel]][B,])
        source_2 <- sp(source_2)
	    target   <- unlist(assayData(GeneDataTime)[[channel]][C,])  
        target   <- sp(target)
        fn       <- featureNames(GeneData)       
        data     <- t(as.data.frame((assayData(GeneDataTime)[[channel]][fn != A & fn != B & fn != C,]))) # all genes except A, B and C       
        data     <- split(data, c(col(data))) # split data frame into a list of rows, each row is now a list and corresponds to the gene families
        data     <- Map(sp, data)
        for (l in 1: length(logics)) {# Test each logic for each time point
          logic <- logics[[l]]; logic_name <- names(logics[l])
          print(paste("logic =", logic_name,",", "time=",time))
          map_logic <- function(source) Map(partial(logic, source),  data)  # maps logic(source, x) to every x in list data       
          # The null model for A
          null_A <- map_logic(source_1)       
          null_A <- Map("unlist", null_A)
		  null_A_model <- Map(partial(sm, target), null_A)          
          # The null model for B
          null_B <- map_logic(source_2)                  
          null_B <- Map("unlist", null_B)   
          null_B_model <- Map(partial(sm, target), null_B)
         
          null_model <- unlist(Filter(complete.cases, c(null_A_model, null_B_model))) # filter out NAs from null model       
          null_mean  <- mean(null_model) 
          net_source <- logic(source_1, source_2)
          model_r_squared <- sm(target, net_source ) # compute actual r squared for proposed logic
          logic_measure <- model_r_squared - null_mean		
		  	 
		  results <- list("null_model"=null_model, "null_mean"=null_mean, "significance"=significance_measure, 
                          "model_r_squared" = model_r_squared, "logic"=logic_name, "logic_measure" = logic_measure, 
                          "time"=time, "source_1"=source_1, "source_2"=source_2, "net_source" = net_source, 
                          "target"=target, "families"=pts, "samples"=samples, "A"=A, "B"=B, "C"=C)
		  logicCorrelations[[time, l]] <- results      # matrix of lists, rows are time, cols are logics          
      }
	}
  return(logicCorrelations) 
} 

# Histogram with mean in red and trio in green, breakpoints = 100
logic_hist <- function(logicCorrelations, time, logic) {    
    logic_index <- logic_indices[[logic]]
    data  <- logicCorrelations[[time, logic_index]]
    logic_trio  <- paste(data$A, logic, data$B," -> ", data$C)
    hist(data$null_model,
         xlab = paste(data$significance,"values for background logic Corellations"),
         main = paste(logic_trio,":","Distribution of", logic, "logic at time", time), 
         breaks = 100)
         abline(v=data$null_mean, col ="red"); abline(v=data$model_r_squared, col = "green")
         legend("topright", lty=c(1,1), col = c("red", "green"), 
               legend = c("Background Mean Logic Correlation", "Actual Logic Correlation"))         
}    

visualize_logic <- function(logicCorrelations, time, logic, dot="l", normalize=F) {
    logic_index <- logic_indices[[logic]]
    data  <- logicCorrelations[[time, logic_index]]
    lim   <- c(data$source_1, data$source_2, data$net_source, data$target)
    x     <- data$families
    logic_trio  <- paste(data$A, logic, data$B," -> ", data$C)
    plot(data$source_1,type=dot, col ="blue", ylim=range(lim),xlab="",ylab=""); par(new=TRUE);
    plot(data$source_2,type=dot, col ="red", ylim=range(lim),xlab="",ylab=""); par(new=TRUE);
    plot(data$net_source,type=dot, col ="green", ylim=range(lim),xlab="",ylab=""); par(new=TRUE); 
    plot(data$target,type=dot, col="black", ylim=range(lim), xlab= "Families", ylab="Expression Level",
        main=paste(logic_trio,":",logic, "logic at time", time))
    legend("topright", lty=c(1,1), col = c("blue", "red", "green", "black"), 
          legend = c(data$A, data$B,"logic of sources", data$C))
}

scatter_logic <- function(logicCorrelations, time, logic) {
  logic_index <- logic_indices[[logic]]
  data        <- logicCorrelations[[time, logic_index]]  
  logic_trio  <- paste(data$A, logic, data$B," -> ", data$C)
  normalized_source <- data$net_source/max(data$net_source); 
  normalized_target <- data$target/max(data$target)   
  lm                <- gene_pair_correlation(normalized_target, normalized_source)   
  r_squared         <- signif(summary(lm)$r.squared, digits = 4)
  p_val             <- signif(p_value(lm), digits = 4)
    plot(normalized_source, normalized_target, type="p", 
         main = paste(logic_trio, "prediction vs actual", data$C, "expression at time", time), 
         xlab=paste(data$A, logic, data$B), ylab=data$C); par(new=TRUE); abline(lm, col="green"); abline(0,1, col="black"); 
    legend("bottomright",lty=c(1,1), col = c("green", "black", NA), title = paste("p-value:", p_val),
           legend = c("best fit", "y = x",paste("R2:", r_squared)))
}

