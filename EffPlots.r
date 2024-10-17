# function eff.plot
#
# Typical usage: eff.plot(eff), where eff <- effect('A:B:C',mod1). That is, eff
# needs to be created using the effect function from package effects.
# It works just as well to write: eff.plot(effect('A:B:C',mod1))
#
# Run library(effects) and source('EffPlots.r') before running eff.plot()
#
# The below arguments are for eff.plot, and not for effect, i.e. if you want to 
# add a label for the y-axis, write eff.plot(effect('A:B:C',mod1),ylab='My label')
#
# However, also check the help for effect. For example, the xlevels attribute of 
# that function is very useful, e.g. eff.plot(effect('my.covariate:my.Factor',my.model)) 
# may give an ugly, unsmooth curve, whereas eff.plot(effect('my.covariate:my.Factor',my.model,xlevels=20)) 
# makes it smoother (50 even smoother).
#
# - Arguments
# eff       an object of class eff created by the effect function
#
# ylab      a text to use for y-label. The default is to use the response name +
#           a reminder that it is plotted on linear scale (even if e.g. link=log is used)
#
# xlab      a text used for x-label. Will be left blank by default.
#
# ylim      limits for the y-axis. Will by default find appropriate values from
#           the eff object.
#
# xlim      limits for x-axis. will by default find values from eff object.
#
# main      title for the graph. Will be left blank by default.
#
# lwd       linewidth to use for error bars and regression line. Default is 1.5 for
#           error bars and 2 for regression lines.
#
# pch       symbol for predicted value. Default is 19.
#
# cex       size of symbol for predicted value. Default is 1.5.
#
# pch.data  symbol for data points. Default is 20.
#
# cex.data  size of symbol for data point. Default is 1.5.
#
# col.data  color of data points. This can be any valid R-color, including colors
#           created by rgb(). For ancova-type with overlay=T different colors 
#           can be used for the different group-levels if the color specification 
#           passed is a vector with more than one color. These colors will be repeated
#           as needed. For all other cases only the first color is used. Default is 
#           c(rgb(.25,.15,.15,.25),rgb(.95,.15,.15,.25),rgb(.15,.15,.95,.25)), 
#           i.e. the first (and only one used in most cases) color is a warm grey, the 
#           secod is light red, and the third is light blue. An alpha of .25 is used to 
#           make the poinst semi-transparent. 
#           If you only want to use a single colour (which is what you would often want
#           unless you have an ANCOVA-type model) then you can use one of the default colours
#           (grey, red, or blue, as above) use e.g. col.data=2, or pass along some valid 
#           colour specification.
#           See grp below.
#
# col.err   as col.data, but for the error bands in ANCOVA
#
# length    length of error bar ends. Default is .05.
#
# cex.lab   size of label text (?). Default is 1.
#
# cex.axis  size of text for x-values. Default is 1.
#
# las.x1    orientation of x-labels, for first first factor. Currently, only implemented for
#           ANOVA-type main effects. See las of par for more info.
#
# ltyp      the order of line types to be cycled through when drawing regression
#           lines. The default order is c(1,5,3,4,2,6). These lines will be repeated if more 
#           than 6 factor levels.
#
# plotdata  should the raw data be plotted? True or False. Default is F
#
# grp       Group data points according to a (random) factor. This is useful for e.g. mixed
#           models, if you want the data points to be grouped along the x-axis according to 
#           their (random) factor levels. This can be combined with col.data. You can use a 
#           factor from your data set, or one from your model object, e.g. model$frame$grouping.factor
#
# ylim.data should ylims be set by the data range, rather than by the ylim argument?
#           True or False. If ylim.data==T the y-axis range will be made to fit data, even if 
#           ylim is used (i.e. ylim.data overrides ylim). Default is F
#
# overlay   for ancova-type interactions: should the separate regression lines
#           be plotted on top of eachother (overlay=T) or on separate panels
#           (overlay=F)? The latter gives grey error bands, rather than error lines.
#           overlay=F gives grey bands also in the simple regreesion case.
#
# factor    used for jittering raw-data when plotted. See documentation for
#           jitter() for information on how to use.
#
# trans     should the results be back transformed? Default is 'lin', which
#           corresponds to no tranfromation. Allowed alternatives are 'sqrt', asr,
#           'log', 'logx', 'logy', 'logxy', depedning on which of the variables
#           have been transformed. If data has been log-transformed before analysis (that is
#           not if a log-link has been used in a glm) trans='log' can be used to
#           get the scale of the plot and the results back on a linear scale. 
#
# add       was a value added to data before transformation? Default is 0. For
#           example is the model was run on log(y+1), the use trans='logy', add=1.
#
# rescale.x if x has been rescaled, then pass center and scale of the rescaled variable to
#           back-scale it, e.g. rescale.x=c(attr(X$var,'scaled:center'),attr(X$var,'scaled:scale'))
#
# xd        used to nudge the x-values by xd (e.g. xd=-.1 or xd=.1) in main-effects plots
#           when two (or more) different values of a non-focal predictor are shown. 
#           Example:
#             eff.plot(effect('ald',m2a,fixed.predictors=list(given.values=c(vild=1))),ylim=c(0,.4),xd=-.1)
#             par(new=T)
#             eff.plot(effect('ald',m2a,fixed.predictors=list(given.values=c(vild=0))),ylim=c(0,.4),xd=.1)
#             par(new=F)
#
#
# Function to get nicer plots for effects from (almost?) any linear model.
# Works for 2 and 3-way interactions and main effects of factors (ANOVA-type models), 
# and also for main effects of continuous predictors (simple regression type) and 
# ANCOVA-type interactions (one continus x one factor).
#
# Should also work for interactions between two or more continuous predictors
# (multiple regression type). However, this has not been thoroughly tried and tested.
#
# This function uses the eff-object created by the function effect in the effects
# package. The effects package also has a plot-function, which differs from this
# one only in layout and design. Please, see the documentation for effect.
#
# For plotting interaction effects note that the number of levels will be used
# to arrange the plot and factor level labels used to label x-axis. Note that
# the order of factors in the lm-command (glm/lmer...) will detrmine the
# appearance of the plot. The inner grouping should be first and the outer last.
#
# For example:
# Inner factor:   123  123    123  123
# Mid factor:      A    B      A    B
# Outer factor:      X           Y
#
# For an ANCOVA-type interaction the order of predictors must be 'Factor:Covariate',
# both in the model (lm, glmer or whatever) specification and in the effect specification.
#
# For an interaction between two continuous variables something like the following
# code will work and give three curves for the levels 0, 100, and 1000:
# eff <- effect('Area_m2:Connectivity',model,xlevels=list(Area_m2=c(0,100,1000), Connectivity=50))

# *****************************************************************************
# *** Conditions for using this code ***
# You are free to use this in your work and produce publishable graphs with it
# freely. You might acknowledge me, but you MUST properly cite Fox' effects
# package.
#
# This script is still under development. Therefore, whenever you use it make
# sure you have the latest version, by asking me.
#
# * DO NOT * redistribute it, but ask me for a new copy if someone else wants to
# use it. This is only so that old bugs will not be spread around, and so that I
# can keep some kind of track of where it goes.
#
# If you have suggestions for improvement just let me know, and we might be able
# to include more functionalities.
#
# *****************************************************************************

# _OO_ 2012-12-01 - 2019-05-23



# *** Note for developing the code:
# plotdata fungerar f?r samtliga fall, utom ifall responsen i en binomial modell
# ?r en n-by-2 skapad med cbind(). Detta fungerar f?r ancova och regressions-
# fallen, men inte f?r anova-fallen. Hur man ska g?ra framg?r nog av de givna
# fallen, men man m?ste ocks? justera vilka kolumner som ska anv?ndas f?r att
# skapa xdata.

# function head ####  

eff.plot <- function(eff, ylab=paste(eff$response, ''),  
        ylim=c(min(el),max(eu)), xlim=c(min(x), max(x)), xlab=' ', main=NULL, 
        trans='lin', add=0, log='', factor=0, rescale.x=NULL,
        lwd=2, pch=19, cex=2, length=.05, cex.lab=1, cex.axis=1,las.x1=0,
        plotdata=F, overlay=T, ltyp=c(1,5,3,4,2,6),
        ylim.data=F,col.data=c(rgb(.25,.15,.15,.25),rgb(.95,.15,.15,.25),rgb(.15,.15,.95,.25)),
        col.err=grey(.8),
        pch.data=20,cex.data=1.5,grp=NULL,
        xaxt='s',yaxt='s',bty='l',xd=0,col=1,line=2.5,tick=T,padj=NA,hadj=NA,lwd.ax=1,
        abl=NULL
        ){
  
  
  #'(linear scale)'),

if (trans=='lin'){
    ee <- matrix(summary(eff)$effect)
    el <- matrix(summary(eff)$lower)
    eu <- matrix(summary(eff)$upper)
    ydata <- eff$data[,1]
} else if (trans=='sqrt'){
    ee <- matrix(summary(eff)$effect^2)
    el <- matrix(summary(eff)$lower^2)
    eu <- matrix(summary(eff)$upper^2)
    ydata <- eff$data[,1]^2
} else if (trans=='log' | trans=='logy' | trans=='logxy'){
    ee <- matrix(exp(summary(eff)$effect)-add)
    el <- matrix(exp(summary(eff)$lower)-add)
    eu <- matrix(exp(summary(eff)$upper)-add)
    ydata <- exp(eff$data[,1])-add
} else if (trans=='logit'){
  ee <- matrix(plogis(summary(eff)$effect))
  el <- matrix(plogis(summary(eff)$lower))
  eu <- matrix(plogis(summary(eff)$upper))
  ydata <- plogis(eff$data[,1])
} else if (trans=='asr'){
  ee <- matrix(sin(summary(eff)$effect)^2)
  el <- matrix(sin(summary(eff)$lower)^2)
  eu <- matrix(sin(summary(eff)$upper)^2)
  ydata <- sin(eff$data[,1])^2
}

col.line <- col.data
for (i in 1:length(col.line)){
  col.line[i] <- paste0(substr(col.line[i],1,7),'FF')
}

  
cd <- c(rgb(.25,.15,.15,.25),rgb(.95,.15,.15,.25),rgb(.15,.15,.95,.25))
if (is.numeric(col.data)){
  if (col.data<4){
    col.data <- cd[col.data]
  }
}
                   
# ANOVAs ####
if (all(sapply(eff$variables, function(x) x$is.factor))){  # if all variables in effect are factors

# 3-way interaction ####
if (length(dim(summary(eff)$effect)) ==3){
    x <- NA
    xdata <- eff$data[,1]*0
    l <- m <- 0
    for (i in 1:length(levels(eff$x[,3]))){
        for (j in 1:length(levels(eff$x[,2]))){
            for (k in 1:length(levels(eff$x[,1]))){
                l <- l+1
                m <- m+1
                x[l] <- m
                # Detta kanske inte funkar i alla l?gen! Kolla 2-way interaction f?r att fatta hur det ska ?ndras. 
                # Det ?r ?ndrat nu, s? b?r funka. Originalkoden st?r kvar p? 2-way int.
                xdata[eff$data[names(eff$x)[1]]==levels(eff$x[,1])[k] & 
                   eff$data[names(eff$x)[2]]==levels(eff$x[,2])[j] & 
                   eff$data[names(eff$x)[3]]==levels(eff$x[,3])[i]] <- m
                
            }
            m <- m+1
        }
        m <- m+2
    }

    if (plotdata==T){
       if (ylim.data==T){
         ylim=c(min(ydata[!is.nan(ydata)]), max(ydata[!is.nan(ydata)]))
       }
      # Not tested
      if (!is.null(grp)){
        grp <- as.numeric(grp)
        nd <- unique(diff(grp))
        nd <- min(nd[nd>0])
        xdata <- xdata + (grp-mean(unique(grp)))*nd/7
      }
      
       plot(jitter(ydata)~jitter(xdata),pch=pch.data,cex=cex.data,col=col.data,
           ylim=ylim,bty=bty, xlim=c(min(x)-1, max(x)+1), xaxt='n',xlab='',
           ylab=ylab,main=main, log=log, cex.lab=cex.lab)
       points(ee~x, pch=pch,cex=cex,col=col)

    } else {
        plot(ee~x,ylim=ylim,bty=bty,pch=pch,cex=cex,xaxt='n',xlab='',
            xlim=c(min(x)-1, max(x)+1), ylab=ylab,main=main, log=log, cex.lab=cex.lab,
            col=col)
    }
    arrows(x,el,x,eu,length=length,lwd=lwd,angle=90,code=3,col=col)
    axis(side=1, at=x, labels=eff$x[,1])
    axis(side=1, at=tapply(x,eff$x[,2]:eff$x[,3],mean),
        labels=rep(levels(eff$x[,2]),each=length(levels(eff$x[,3]))),padj=2,tick=F)
    axis(side=1, at=tapply(x,eff$x[,3],mean),labels=levels(eff$x[,3]),padj=3,tick=F)

} else if (all(sapply(eff$variables, function(x) x$is.factor)) & length(dim(summary(eff)$effect)) ==2){

# 2-way interactions ####
    x <- NA
    l <- m <- 0
    xdata <- eff$data[,1]*0

  #  print(levels(eff$x[,1]))
  #  print(levels(eff$x[,2]))
        for (j in 1:length(levels(eff$x[,2]))){
            for (k in 1:length(levels(eff$x[,1]))){
                l <- l+1
                m <- m+1
                x[l] <- m
                #xdata[eff$data[2]==levels(eff$x[,1])[k] & 
                #   eff$data[3]==levels(eff$x[,2])[j]] <- m
                xdata[eff$data[names(eff$x)[1]]==levels(eff$x[,1])[k] & 
                        eff$data[names(eff$x)[2]]==levels(eff$x[,2])[j]] <- m
            }
            m <- m+2
        }
  #  print(xdata)
    
    if (plotdata==T){
       if (ylim.data==T){
         ylim=c(min(ydata[!is.nan(ydata)]), max(ydata[!is.nan(ydata)]))
       }
      # Not tested
      if (!is.null(grp)){
        grp <- as.numeric(grp)
        nd <- unique(diff(grp))
        nd <- min(nd[nd>0])
        xdata <- xdata + (grp-mean(unique(grp)))*nd/7
      }
      
       plot(jitter(ydata)~jitter(xdata),pch=pch.data,cex=cex.data,col=col.data,
           ylim=ylim,bty=bty, xlim=c(min(x)-.5, max(x)+.5), xaxt='n',xlab='',
           ylab=ylab,main=main, log=log, cex.lab=cex.lab)
       points(ee~x, pch=pch,cex=cex,col=col)
       
    } else {
    plot(ee~x,ylim=ylim,bty=bty,pch=pch,cex=cex,
        xlim=c(min(x)-.5, max(x)+.5), xaxt='n',xlab='',ylab=ylab,main=main, log=log, 
	      cex.lab=cex.lab,col=col)
    }
    arrows(x,el,x,eu,length=length,lwd=lwd,angle=90,code=3,col=col)
    axis(side=1, at=x, labels=eff$x[,1],cex.axis=cex.axis)
    axis(side=1, at=tapply(x,eff$x[,2],mean),labels=levels(eff$x[,2]),padj=2,tick=F)


} else if (length(dim(summary(eff)$effect)) ==1){
  # Main effects ####
    x <- 1:length(ee)
    x0 <- x
    xlm <- c(min(x)-.5, max(x)+.5)
    x <- x + xd
    
    if (plotdata==T){
       xdata <- as.numeric(eff$data[,eff$term])
       if (ylim.data==T){
         ylim=c(min(ydata[!is.nan(ydata)]), max(ydata[!is.nan(ydata)]))
       }
       
       if (!is.null(grp)){
         grp <- as.numeric(grp)
         nd <- unique(diff(grp))
         nd <- min(nd[nd>0])
         xdata <- xdata + (grp-mean(unique(grp)))*nd/7
       }

#       plot(jitter(ydata)~jitter(xdata),pch=pch.data,cex=cex.data, col=col.data,
#           ylim=ylim,bty=bty, xlim=xlm, xaxt='n',xlab='',
#           ylab=ylab,main=main, log=log, cex.lab=cex.lab,cex.axis=cex.axis)
       plot(jitter(ydata)~jitter(xdata),pch=pch.data,cex=cex.data, col=col.data,
            ylim=ylim,bty=bty, xlim=xlm, xaxt='n',yaxt='n',
            xlab='',ylab='',main=main, log=log)
       
       points(ee~x, pch=pch,cex=cex,col=col)
    } else {
    #    plot(ee~x,ylim=ylim,bty=bty,pch=pch,cex=cex,xaxt='n',xlab=xlab,ylab=ylab,main=main,
    #        xlim=xlm, log=log, cex.lab=cex.lab,col=col,cex.axis=cex.axis,yaxt=yaxt)
      plot(ee~x,ylim=ylim,bty=bty,pch=pch,cex=cex,xaxt='n',yaxt='n',main=main,
           xlim=xlm, log=log, col=col,xlab='',ylab='')
    }
    arrows(x,el,x,eu,length=length,lwd=lwd,angle=90,code=3,col=col) 
    axis(side=1, at=x0, labels=eff$x[,1], cex.axis=cex.axis,line=NULL,xaxt=xaxt,
         tick=tick,padj=padj,lwd=lwd.ax,lwd.ticks=lwd.ax,las=las.x1)
    axis(side=2, cex.axis=cex.axis,line=NULL,yaxt=yaxt,tick=tick,padj=-padj,lwd=lwd.ax,lwd.ticks=lwd.ax)
    mtext(xlab,side=1,line=line,cex=cex.lab)
    mtext(ylab,side=2,line=line,cex=cex.lab)

}

# Simple regression effect ####
} else if (!any(sapply(eff$variables, function(x) x$is.factor)) & length(eff$variables)==1){
    x <- sapply(eff$variables, function(x) x$levels)
    if (trans=='logx' | trans=='logxy'){
       x <- exp(sapply(eff$variables, function(x) x$levels)-add)
    } 
# Trying to fix "back-scaling", but seems difficult as scale attributes are not passed
    if (!is.null(rescale.x)){
      x <- x*rescale.x[2]+rescale.x[1]
    }
    
    
    if (plotdata==T){
        xdata <- eff$data[,eff$variables[[1]]$name]  # Den h?r ?r inte fixad f?r trans='logx' --- Jo, nu ska den vara det
        if (trans=='logx' | trans=='logxy'){
          xdata <- exp(xdata)
        }
        if (!is.null(rescale.x)){
          xdata <- xdata*rescale.x[2]+rescale.x[1]
        }
        if (is.matrix(ydata)==T){                   # if dependent a cbind from a binomial fit
            ydata <- ydata[,1]/(ydata[,1]+ydata[,2])
        }
        if (ylim.data==T){
          ylim=c(min(ydata[!is.nan(ydata)]), max(ydata[!is.nan(ydata)]))
        }
    }
    if (xlab==' '){ xlab <- eff$variables[[1]]$name}
    plot(ee~x, type='n', ylim=ylim, xlim=xlim, bty=bty, xlab='', ylab='', main=main, log=log, 
         xaxt='n',yaxt='n')
    if (overlay==T){
       lines(el~x)
       lines(eu~x)
       if (lwd==1.5){lwd=2}
    }else{
        if (col==1){col=gray(0.8)}
        polygon(c(x, rev(x)), c(eu, rev(el)), col=col.err, density=-1, border=NA)
    }
    lines(ee~x, lwd=lwd)

    if (plotdata==T){
            points(jitter(ydata, factor=factor)~jitter(I(xdata+add)),pch=pch.data,cex=cex.data, col=col.data)
    }
    axis(side=1, cex.axis=cex.axis,line=NULL,xaxt=xaxt,tick=tick,padj=padj,lwd=lwd.ax,lwd.ticks=lwd.ax)
    axis(side=2, cex.axis=cex.axis,line=NULL,yaxt=yaxt,tick=tick,padj=-padj,lwd=lwd.ax,lwd.ticks=lwd.ax)
    mtext(xlab,side=1,line=line,cex=cex.lab)
    mtext(ylab,side=2,line=line,cex=cex.lab)
    

# "ANCOVA-effect" ####
} else if (any(sapply(eff$variables, function(x) x$is.factor)) & length(eff$variables)==2){
                            
    fact <- which(sapply(eff$variables, function(x) x$is.factor))
    covar <- which(sapply(eff$variables, function(x) x$is.factor)==F)
    factdat <- eff$data[, eff$variables[[fact]]$name]
    if (length(ltyp)<length(levels(factdat))){ltyp <- rep(ltyp,length.out=length(levels(factdat)))}
    
    if (length(col.data)<length(levels(factdat))){col.data <- rep(col.data,length.out=length(levels(factdat)))}
    col.line=col.data
    
    if (length(col.err)<length(levels(factdat))){col.err <- rep(col.err,length.out=length(levels(factdat)))}
    
    if (trans=='logx' | trans=='logxy'){
        x <- exp(sapply(eff$variables, function(x) x$levels)[[covar]]-add)
    }else{
        x <- sapply(eff$variables, function(x) x$levels)[[covar]]
        if (!is.null(rescale.x)){
          x <- x*rescale.x[2]+rescale.x[1]
        }
    }

    i <- 0
    mainauto <- F
    for (L in levels(factdat)){
        i <- i+1
        
        #ydata <- plogis(eff$data[,1])
        if (plotdata==T){
            ydata <- eff$data[factdat==L,1]
            xdata <- eff$data[factdat==L, eff$variables[[2]]$name]
            cdata <- col.data[factdat==L]
            if (trans=='logx' | trans=='logxy' | trans=='log'){
              xdata <- exp(eff$data[factdat==L, eff$variables[[2]]$name])
            }
            if (!is.null(rescale.x)){
              xdata <- xdata*rescale.x[2]+rescale.x[1]
            }
            
            if (trans=='log' | trans=='logy' | trans=='logxy'){
              ydata <- exp(ydata)-add
            } else if(trans=='logit'){
              ydata <- plogis(ydata)
            } else if(trans=='sqrt'){
              ydata <- ydata^2
            }
            if (is.matrix(ydata)==T){                   # if dependent a cbind from a binomial fit
                ydata <- ydata[,1]/(ydata[,1]+ydata[,2])
            }
            if (ylim.data==T){
              ylim=c(min(ydata[!is.nan(ydata)]), max(ydata[!is.nan(ydata)]))
            }
        }

        
        if (trans=='log' | trans=='logy' | trans=='logxy'){
          ee <- matrix(exp(summary(eff)$effect[i,])-add)
          el <- matrix(exp(summary(eff)$lower[i,])-add)
          eu <- matrix(exp(summary(eff)$upper[i,])-add)
        }else if(trans=='logit'){
          ee <- matrix(plogis(summary(eff)$effect[i,]))
          el <- matrix(plogis(summary(eff)$lower[i,]))
          eu <- matrix(plogis(summary(eff)$upper[i,]))
        } else if(trans=='sqrt'){
          ee <- matrix(summary(eff)$effect[i,]^2)
          el <- matrix(summary(eff)$lower[i,]^2)
          eu <- matrix(summary(eff)$upper[i,]^2)
        } else {
          ee <- summary(eff)$effect[i,]
          eu <- summary(eff)$upper[i,]
          el <- summary(eff)$lower[i,]
        }
        
        
        if (xlab==' '){ xlab <- eff$variables[[2]]$name}
        if ((is.null(main) | mainauto==T) & overlay==F){ 
           main[i] <- L
           mainauto <- T
        }
        if (overlay==T){
          #main <- L
           if (i==1){
              if (lwd==1.5){lwd=2}
              plot(ee~x, type='l', ylim=ylim, xlim=xlim, bty=bty, xlab=xlab, ylab=ylab, 
                main=main[i], log=log, cex.lab=cex.lab, lty=ltyp[i], lwd=lwd,col=col.data[i])
           }else{
              lines(ee~x, lwd=lwd, lty=ltyp[i],col=col.line[i])
           }
          polygon(c(x, rev(x)), c(eu, rev(el)), col=col.err[i], density=-1, border=NA)
          if (plotdata==T){
            points(jitter(ydata, factor=factor)~jitter(I(xdata+add)),pch=pch.data,cex=cex.data,col=col.data[i])
          }
           lines(ee~x, lwd=lwd, lty=ltyp[i],col=col.line[i])
        } else {   # not overlay
          plot(ee~x, type='n', ylim=ylim, xlim=xlim, bty=bty, xlab=xlab, ylab=ylab, 
          main=main[i], log=log, cex.lab=cex.lab,cex.axis=cex.axis)
          polygon(c(x, rev(x)), c(eu, rev(el)), col=col.err[i], border=0) #gray(0.8)
          if (plotdata==T){
            points(jitter(ydata, factor=factor)~jitter(I(xdata+add)),pch=pch.data,cex=cex.data,col=cdata)
          }
          lines(ee~x, lwd=lwd,col=col.line[i])
          
        }
        
#        if (plotdata==T){  
#          if (overlay==T){
#              points(jitter(ydata, factor=factor)~jitter(I(xdata+add)),pch=pch.data,cex=cex.data,col=col.data[i])
#          } else {
#            points(jitter(ydata, factor=factor)~jitter(I(xdata+add)),pch=pch.data,cex=cex.data,col=col.data[1])
#          }
#        }
        
        if (!is.null(abl)){
          abline(abl,lty=2)
        }
        
    }
} else if(!any(sapply(eff$variables, function(x) x$is.factor)) & length(dim(summary(eff)$effect)) ==2){
  # Regression, 2 variables with interaction ----
      #fact <- which(sapply(eff$variables, function(x) x$is.factor))
      #covar <- which(sapply(eff$variables, function(x) x$is.factor)==F)
      fact <- which.min(c(length(eff$variables[[1]]$levels),length(eff$variables[[2]]$levels)))
      covar <- which.max(c(length(eff$variables[[1]]$levels),length(eff$variables[[2]]$levels)))
      
      factdat <- eff$data[, eff$variables[[fact]]$name]
      
      factLevs <- length(eff$variables[[fact]]$levels)
      
      if (length(ltyp)<factLevs){ltyp <- rep(ltyp,length.out=factLevs)}
      #if (length(ltyp)<length(levels(factdat))){ltyp <- rep(ltyp,length.out=length(levels(factdat)))}
      #if (length(col.data)<length(levels(factdat))){col.data <- rep(col.data,length.out=length(levels(factdat)))}
      
      #if (length(col.err)<factLevs){col.err <- rep(col.err,length.out=factLevs)}
      col.err=rgb(.75,.75,.75,.5)
      
      
      if (trans=='logx' | trans=='logxy'){
        x <- exp(sapply(eff$variables, function(x) x$levels)[[covar]]-add)
      }else{
        x <- sapply(eff$variables, function(x) x$levels)[[covar]]
        if (!is.null(rescale.x)){
          x <- x*rescale.x[2]+rescale.x[1]
        }
      }
      
      i <- 0
      mainauto <- F
      for (L in 1:factLevs){
        i <- i+1
        #ydata <- plogis(eff$data[,1])
        if (plotdata==T){
          #ydata <- eff$data[factdat==L,1]
          ydata <- eff$data[,1]
          xdata <- eff$data[, eff$variables[[covar]]$name]
          #cdata <- col.data[factdat==L]
          cdata <- rgb(.75,.75,.75,.5)
          if (trans=='logx' | trans=='logxy' | trans=='log'){
            xdata <- exp(eff$data[factdat==L, eff$variables[[2]]$name])
          }
          if (!is.null(rescale.x)){
            xdata <- xdata*rescale.x[2]+rescale.x[1]
          }
          
          if (trans=='log' | trans=='logy' | trans=='logxy'){
            ydata <- exp(ydata)-add
          } else if(trans=='logit'){
            ydata <- plogis(ydata)
          } else if(trans=='sqrt'){
            ydata <- ydata^2
          }
          if (is.matrix(ydata)==T){                   # if dependent a cbind from a binomial fit
            ydata <- ydata[,1]/(ydata[,1]+ydata[,2])
          }
          if (ylim.data==T){
            ylim=c(min(ydata[!is.nan(ydata)]), max(ydata[!is.nan(ydata)]))
          }
        }


        if (fact==1){
          if (trans=='log' | trans=='logy' | trans=='logxy'){
            ee <- matrix(exp(summary(eff)$effect[i,])-add)
            el <- matrix(exp(summary(eff)$lower[i,])-add)
            eu <- matrix(exp(summary(eff)$upper[i,])-add)
          }else if(trans=='logit'){
            ee <- matrix(plogis(summary(eff)$effect[i,]))
            el <- matrix(plogis(summary(eff)$lower[i,]))
            eu <- matrix(plogis(summary(eff)$upper[i,]))
          } else if(trans=='sqrt'){
            ee <- matrix(summary(eff)$effect[i,]^2)
            el <- matrix(summary(eff)$lower[i,]^2)
            eu <- matrix(summary(eff)$upper[i,]^2)
          } else {
            ee <- summary(eff)$effect[i,]
            eu <- summary(eff)$upper[i,]
            el <- summary(eff)$lower[i,]
          }
        }else{
          if (trans=='log' | trans=='logy' | trans=='logxy'){
            ee <- matrix(exp(summary(eff)$effect[,i])-add)
            el <- matrix(exp(summary(eff)$lower[,i])-add)
            eu <- matrix(exp(summary(eff)$upper[,i])-add)
          }else if(trans=='logit'){
            ee <- matrix(plogis(summary(eff)$effect[,i]))
            el <- matrix(plogis(summary(eff)$lower[,i]))
            eu <- matrix(plogis(summary(eff)$upper[,i]))
          } else if(trans=='sqrt'){
            ee <- matrix(summary(eff)$effect[,i]^2)
            el <- matrix(summary(eff)$lower[,i]^2)
            eu <- matrix(summary(eff)$upper[,i]^2)
          } else {
            ee <- summary(eff)$effect[,i]
            eu <- summary(eff)$upper[,i]
            el <- summary(eff)$lower[,i]
          }
        }
        
        if (xlab==' '){ xlab <- eff$variables[[2]]$name}
        if ((is.null(main) | mainauto==T) & overlay==F){ 
          main[i] <- L
          mainauto <- T
        }
        if (overlay==T){
          #main <- L
          if (i==1){
            if (lwd==1.5){lwd=2}
            plot(ee~x, type='l', ylim=ylim, xlim=xlim, bty=bty, xlab=xlab, ylab=ylab, 
                 main=main[i], log=log, cex.lab=cex.lab, lty=ltyp[i], lwd=lwd,col=col.data[i])
          }else{
            lines(ee~x, lwd=lwd, lty=ltyp[i],col=col.line[i])
          }
          #polygon(c(x, rev(x)), c(eu, rev(el)), col=col.err[i], density=-1, border=NA)
          polygon(c(x, rev(x)), c(eu, rev(el)), col=col.err, density=-1, border=NA)
          if (plotdata==T){
            points(jitter(ydata, factor=factor)~jitter(I(xdata+add)),pch=pch.data,cex=cex.data,col=col.data[i])
          }
          lines(ee~x, lwd=lwd, lty=ltyp[i],col=col.line[i])
        } else {   # not overlay
          plot(ee~x, type='n', ylim=ylim, xlim=xlim, bty=bty, xlab=xlab, ylab=ylab, 
               main=main[i], log=log, cex.lab=cex.lab,cex.axis=cex.axis)
          polygon(c(x, rev(x)), c(eu, rev(el)), col=col.err[i], border=0) #gray(0.8)
          if (plotdata==T){
            points(jitter(ydata, factor=factor)~jitter(I(xdata+add)),pch=pch.data,cex=cex.data,col=cdata)
          }
          lines(ee~x, lwd=lwd,col=col.line[i])
          
        }
        
        if (!is.null(abl)){
          abline(abl,lty=2)
        }
        
    }

# Multiple regression with two x-varialbles ####
# assuming that the 1st variable has 3 levels, and the 2nd is the one to be used 
# on the x-axis
} else if (!any(sapply(eff$variables, function(x) x$is.factor)) & length(eff$variables)==2){
    #fact <- which(sapply(eff$variables, function(x) x$is.factor))
    #covar <- which(sapply(eff$variables, function(x) x$is.factor)==F)
    #factdat <- eff$data[, eff$variables[[fact]]$name]
    xlevdat <- sapply(eff$variables, function(x) x$levels)[[1]]
    
    if (trans=='logx' | trans=='logxy'){
        x <- exp(sapply(eff$variables, function(x) x$levels)[[2]]-add)
    }else{
        x <- sapply(eff$variables, function(x) x$levels)[[2]]
    }

    i <- 0
    mainauto <- F
    for (L in xlevdat){
        i <- i+1

        if (plotdata==T){
            ydata <- eff$data[xlevdat==L,1]
            xdata <- eff$data[xlevdat==L, eff$variables[[2]]$name]
            if (is.matrix(ydata)==T){                   # if dependent a cbind from a binomial fit
                ydata <- ydata[,1]/(ydata[,1]+ydata[,2])
            }
            if (ylim.data==T){
              ylim=c(min(ydata[!is.nan(ydata)]), max(ydata[!is.nan(ydata)]))
            }
            
        }

        #print(summary(eff)$effect)
        if (trans=='log' | trans=='logy' | trans=='logxy'){
            ee <- matrix(exp(summary(eff)$effect[i,])-add)
            el <- matrix(exp(summary(eff)$lower[i,])-add)
            eu <- matrix(exp(summary(eff)$upper[i,])-add)
        } else {
            ee <- summary(eff)$effect[i,]
            eu <- summary(eff)$upper[i,]
            el <- summary(eff)$lower[i,]
        }
        
        if (xlab==' '){ xlab <- eff$variables[[2]]$name}
        #if ((main==' ' | mainauto==T) & overlay==F){ 
           main <- paste(eff$variables[[1]]$name, ': ', L, sep='')
           #mainauto <- T
        #}
        if (overlay==T){
           if (i==1){
              if (lwd==1.5){lwd=2}
              plot(ee~x, type='l', ylim=ylim, xlim=xlim, bty=bty, xlab=xlab, ylab=ylab, 
                main=main, log=log, cex.lab=cex.lab, lty=ltyp[i], lwd=lwd)
           }else{
              lines(ee~x, lwd=lwd, lty=ltyp[i])
           }
           lines(el~x, lty=ltyp[i])
           lines(eu~x, lty=ltyp[i])
        } else {
          if (i==1){
            if (lwd==1.5){lwd=2}
          plot(ee~x, type='n', ylim=ylim, xlim=xlim, bty=bty, xlab=xlab, ylab=ylab, 
              main=main, log=log, cex.lab=cex.lab, lty=ltyp[i],lwd=lwd)                 # row.names(summary(eff)$effect)[i]
          }else{
            lines(ee~x, lwd=lwd, lty=ltyp[i])
          }
          polygon(c(x, rev(x)), c(eu, rev(el)), col=rgb(.5,.5,.5,.25), density=-1, border=NA)
          lines(ee~x, lwd=lwd, lty=ltyp[i])
        }
        
        if (plotdata==T){  
            points(jitter(ydata, factor=factor)~jitter(I(xdata+add)),pch=pch.data,cex=cex.data,col=col.data[1])
        }
    }
}


}  # end function
