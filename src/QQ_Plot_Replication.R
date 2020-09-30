#Function to generate qqplots

qqunif = function(p,BH=T,CI=T,BF=T,...)
{
  par(mar=c(5,6,4,1)+.1)
  nn = length(p)
  xx =  -log10((1:nn)/(nn+1))
  plot( xx,  -sort(log10(p)),
        #main='Example Title',
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        pch=16,
        cex.lab=1.4, cex.main=2.0,  mgp=c(3,1,0),
        ... )
  abline(0,1,col='gray')
  if(BH)
  {
    abline(-log10(0.05),1, col='red',lty=1)
    abline(-log10(0.10),1, col='orange',lty=1)
    abline(-log10(0.25),1, col='green',lty=1)
    #legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
    #legend('topleft', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
    #col=c('red','orange','green'),lty=c(1,1,1))
  }
  if (BF)
  {
    abline(h=-log10(0.05/nn)) ## bonferroni
  }
  
  #if(CI)
  {
    ## create the confidence intervals
    #c95 <- rep(0,nn)
    #c05 <- rep(0,nn)
    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)
    ## this portion was posted by anonymous on
    ## http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
    
    #for(i in 1:nn)
    {
      #c95[i] <- qbeta(0.95,i,nn-i+1)
      #c05[i] <- qbeta(0.05,i,nn-i+1)
    }
    
    #lines(xx,-log10(c95),col='gray')
    #lines(xx,-log10(c05),col='gray')
  }
}

#For each results file
trait <- read.csv(file= '#Absolute path to replication results')
trait_values <- trait$pvalue
trait_plot <- qqunif(trait_values)