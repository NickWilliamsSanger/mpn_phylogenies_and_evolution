require("RColorBrewer")
require("scales")
require("tidyverse")


##  Functions that extend the tree plotting available in the rsimpop package


add_annotation=function(pdx,tree,annot_function,control=NULL){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(pdx,tree,tree$edge[i,2],control=control))
}



add_defaults=function(control,defaults,mandatory_fields=c()){
  if(is.null(control)){
    control=list()
  }
  for(field in mandatory_fields){
    if(is.null(control[[field]])){
      stop(sprintf("Required parameter %s not supplied in control list",field))
    }
  }
  
  for(field in names(defaults)){
    if(is.null(control[[field]])){
      control[[field]]=defaults[[field]]
    }
  }
  control
}

add_binary_proportion=function(pdx,##<< PDX object or list including details matrix
                               tree,##<< enhanced phylo returned from plot_tree
                               node,##<< Tree node - maps to pdx$details$node
                               control,##<< Control parameters.. Requires bfield
                               ...
                               ){
  control=add_defaults(control,defaults=list( b.add.line=TRUE,
                                              b.add.text=FALSE),
                       mandatory_fields="bfield")

  bfield=control$bfield
  #browser()
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(pdx,tree,node)
  bdat=pdx$details[[bfield]][info$idx]
  if(is.null(bdat) || class(bdat)!="logical"){
    stop("Error in provided bfield (does it exist and is it boolean?)")
  }
  pass=sum(bdat,na.rm=TRUE)
  fail=sum(!bdat,na.rm=TRUE)
  tot=pass+fail
  ycross=info$yb+(fail/tot)*(info$yt-info$yb)
  ##Could add in a third category NA
  #missing=sum(is.na(bdat))
  if(control$b.add.line){
    arrows(y0=info$yb,y1=ycross,x0=info$x,length = 0,col="red",lend=1,...)
    arrows(y0=ycross,y1=info$yt,x0=info$x,length = 0,col="blue",lend=1,...)
  }
  if(control$b.add.text){
    text(y=ycross,x=info$x,label=pass,pos = 4,offset = 0,...)
  }
}


##Chromosome plot
test_chr=function(col_ctr,ptarget){
  tot=sum(col_ctr)
  if(tot>20){
    chires=chisq.test(col_ctr,p=ptarget,simulate.p.value = TRUE,B=500)
    pval=chires$p.value
    if(is.na(pval)){
      browser()
    }
    idx=which.max(abs(chires$stdres))
    stdres=chires$stdres[idx]
  }else{
    pval=1
    idx=-1
    stdres=0
  }
  list(idx.max.deviation=idx,pval=pval,stdres=stdres)
}

get_chr_cols=function(){
  cc=hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0,direction = 1)(24)
  cc[as.vector(rbind(1:12,13:24))]##Alternate
}

add_chromosome_decomp=function(pdx,##<< PDX object or list including details matrix
                               tree,##<< enhanced phylo returned from plot_tree
                               node,##<< Tree node - maps to pdx$details$node
                               control,##<< Control parameters.. Requires bfield
                               ...)
{
                               ##tree,node,res,ptarget){
  info=get_edge_info(pdx,tree,node)
  control=add_defaults(control,defaults=list(pcut=0.01),
                       mandatory_fields=c("targetprop"))
  ptarget=control$targetprop

  legend=FALSE
  if(is.null(tree$direction)){
    stop("Need to add to existing plotted tree")
  }
  chrs=sprintf("%s",c(1:22,"X","Y"))
  chr=match(pdx$details$Chrom[info$idx.in.detail],chrs)
  col_ctr=tabulate(chr,nbins=24)
  cols=get_chr_cols()
  tot=sum(col_ctr)
  y0=info$yb
  y1=info$yt
  x=info$x
  mh=y1-y0
  start=y1-0.05*mh
  #browser()
  unit=0.9*mh/tot
  at=c()

  chrtest=test_chr(col_ctr,ptarget)

  for(i in 1:length(col_ctr)){
    if(chrtest$idx.max.deviation==i){
      y0=start
      y1=start-unit*col_ctr[i]
      mh=y1-y0
    }
    if(col_ctr[i]>0){
      rect(xleft=x-0.25,xright=x+0.25,ytop = start,ybottom=start-unit*col_ctr[i],col=cols[i],border="black")
      start=start-unit*col_ctr[i]

    }

  }

  if(chrtest$pval<control$pcut){
    ##idx=which.max(abs(chires$stdres))
    ccol="blue"
    voff=(runif(1)-0.5)*0.5
    yy=y0+mh*0.5+voff*mh
    txtbox(x=x-0.25,y=yy,txt = sprintf("%schr%s\nP=%3.2g",ifelse(chrtest$stdres>0,"+","-"),chrs[chrtest$idx],chrtest$pval),ccex = 1,pos = NULL,offset = 0,col=ccol)
    points(x,yy,cex=1,pch=4,col=ccol)
  }
  col_ctr
}


get_consequence_scheme=function(allowed.value=c("missense","nonsense","ess_splice","frameshift","inframe","cna","loh","start_lost","nc_ess_splice")){
  mut.scheme=read.table("consequence_col_scheme.txt",head=T,stringsAsFactors = FALSE)
  colnames(mut.scheme)=c("value","col","pch")
  if(length(allowed.value)==0){
    mut.scheme
  }else{
    mut.scheme %>% filter(value %in% allowed.value)
  }
}


add_simple_labels=function(
                    pdx,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                    tree,##<< enhanced phylo returned from plot_tree
                    node,##<< Node (see details)
                    control,
                    ... ##<< paremeters for points (not color)
){
  control=add_defaults(control,defaults=list( query.field="VC",
                                              label.field="GENE",
                                              query.allowed.df=get_consequence_scheme(),
                                              cex.label=1,
                                              b.add.label=TRUE,
                                              b.add.marker=TRUE)
                       )
  info=get_edge_info(pdx,tree,node)
  ###browser()
  details=pdx$details
  query.field=control$query.field
  query.allowed.df=control$query.allowed.df
  label.field=control$label.field
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>50){
    stop("too many variants to annotate")
  }
  query.value=details[[query.field]][idx]
  idx.match=match(query.value,query.allowed.df$value)
  cols=query.allowed.df$col[idx.match]
  pch=query.allowed.df$pch[idx.match]

  vlabels=details[[label.field]][idx]
  ## spread out
  N=length(idx)
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  if(N>0){
    yd=info$yt-info$yb
    if(N==1){
      y=0.5*(info$yb+info$yt)
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    if(control$b.add.marker){
      points(rep(info$x,N),y,col=cols,pch=pch,...)
    }
    if(control$b.add.label){
      text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=control$cex.label)
    }
  }
  list(node=node,value=query.value)
}


add_specified_labels=function(
  pdx,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
  tree,##<< enhanced phylo returned from plot_tree
  node,##<< Node (see details)
  control,
  ... ##<< paremeters for points (not color)
){
  info=get_edge_info(pdx,tree,node)
  thisnode=node
  if(is.null(control$vars)){
    return(NULL)
  }
  vars=control$vars %>% filter(node==thisnode)
  ###browser()
  if(dim(vars)[1]==0){
    return(NULL)
  }
  if(control$b.annot.just.shared && thisnode<=length(tree$tip.label)){
    return(NULL)
  }
  cols=vars$col
  pch=vars$pch
  vlabels=vars$label
  #vlabels=details[[label.field]][idx]
  ## spread out
  N=length(vlabels)
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  y=c()
  if(N>0){
    yd=info$yt-info$yb
    ##browser()
    if(N==1){
      y=info$yb+rbeta(1,shape1 = 10,shape2 = 10)*yd
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    if(control$b.add.marker){
      points(rep(info$x,N),y,col=cols,pch=pch,bg=cols,...)
    }
    if(control$b.add.label){
      text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=control$cex.label,col=control$col.label)
    }
  }
  list(node=node,value=N,y=y,x=info$x)
}

add_length=function(
  pdx,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
  tree,##<< enhanced phylo returned from plot_tree
  node,##<< Node (see details)
  control,
  ... ##<< paremeters for points (not color)
){
  info=get_edge_info(pdx,tree,node)
  thisnode=node
  #browser()
  ll=tree[[control$el_param]]
  if(is.null(ll)){
    stop(sprintf("tree$%s does not exist",control$el_param))
  }
  len=ll[which(tree$edge[,2]==node)]
  y=info$yb+(0.2+0.6*runif(1))*(info$yt-info$yb)
  if(len<control$cutoff && length(info$samples)>=control$min.shared){
  text(info$x,y,labels = sprintf("%3.1f",len),
       pos = 2,offset = 0,
       cex=control$cex.annot.label)
  }
  list(node=node,value=len)
}





add_vaf=function(pdx,##<< PDX object or list including details matrix
                 tree,
                 node,
                 control,
                 ...
){
  
  control=add_defaults(control,defaults=list( samples=c(),
                                              b.plot.bars=TRUE,
                                              lwd.rect=1,
                                              min.depth=1,
                                              filter.on=NULL,
                                              cex.label=1,
                                              min.mean.vaf=0.35,
                                              mode="recap"))
  
  #browser()
  if(control$mode=="recap"){
    barcol="red"
    txtcol="blue"
  }else{
    barcol="darkgrey"
    txtcol="black"
  }
  
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(pdx,tree,node)
  #browser()
  
  if(length(info$idx)==0){
    return(NULL)
  }
  
  if(!is.null(control$filter.on)){
    # browser()
    info$idx=info$idx[which(pdx$details[[control$filter.on]][info$idx])]
    if(length(info$idx)==0){
      return(NULL)
    }
  }
  
  
  if(control$b.plot.bars){
    plotF=plotBars
  }else{
    plotF=plotPie
  }
  samples=control$samples
  if(is.null(samples) || length(samples)==0){
    samples=info$samples
  }
  
  ##browser()
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
    adep=rowSums(pdx$dep[,samples])
  }else{
    df=data.frame(mtr=pdx$mtr[info$idx,samples],
                  dep=pdx$dep[info$idx,samples],stringsAsFactors = FALSE)
    adep=pdx$dep[,samples]
  }
  
  
  mdep=mean(df$dep)
  depth_zscore=sqrt(length(df$dep))*(mdep-mean(adep))/sd(adep)
  
  df=cbind(df,pdx$details[info$idx,])
  df=df[which(df$dep>=control$min.depth),]
  
  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  N=dim(df)[1]
  #cat(N,"\n")
  if(N==0){
    return(df)
  }
  
  
  df=df[order(df$vaf),]
  yd=info$yt-info$yb
  
  
  if(N==1){
    y=0.5*(info$yb+info$yt)
    width=yd
  }else{
    y1=seq(info$yb,info$yt,length.out = N+2)
    #Don't use the ends..
    y=y1[2:(N+1)]
    width=y[2]-y[1]
  }
  
  if(!control$b.plot.bars){
    r=0.8  ##r>0.5 will cause likely overlap problems
  }else{
    r=0.4
  }
  for(i in 1:N){
    vaf=min(df$vaf[i],0.999)
    #if(!is.na(df$label[i]) && nchar(df$label[i])>0){
    #  text(x=info$x-r,y=y[i],labels = df$label[i],cex=control$cex.label,pos=2,offset=0)
    #}
    if(is.na(vaf)){
      plotF(x=info$x,y = y[i],radius=r,col=c("lightgray","lightgray"),prop=c(0,1),border="lightgray",width=width)
    }else{
      plotF(x=info$x,y = y[i],radius = r,col=c(barcol,"white"),prop = c(vaf,1-vaf),width = width)
    }
  }
  for(i in 1:N){
    if(!is.null(df$label) && !is.na(df$label[i]) && nchar(df$label[i])>0){
      text(x=info$x-r,y=y[i],labels = df$label[i],cex=control$cex.label,pos=2,offset=0)
      points(x=info$x,y=y[i],pch=df$pch[i],col=df$col[i],bg=df$col[i])
      #cat(df$label[i],":",df$vaf[i],": rank=",length(which(df$vaf[i]>df$vaf))/length(df$vaf),":",df$vaf,"\n")
    }
  }
  
  
  if( !control$b.plot.bars){
    return(df)
  }
  ##Now check to see if we need to highlight
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  if(control$mode=="recap"){
    rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,border="darkgrey",lwd=control$lwd.rect)
    if(MTR/DEP>0.01){
      txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
      text(txt,x=info$x,y=info$yb+0.3*(info$yt-info$yb),col=txtcol,cex=0.6)
    }
    
  }else{
    min.mean.vaf=control$min.mean.vaf
    z=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
    z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
    z$p.value=max(z$p.value,z2$p.value)
    txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
    if(z$p.value<0.05){
      if(z$p.value<0.05/dim(tree$edge)[1]){
        border.color="red"
      }else{
        border.color="blue"
      }
    }else{
      border.color="darkgrey"
    }
    
    rect(xleft=info$x-r,xright=info$x+r,ybottom=y[1]-width/2,ytop=y[N]+width/2,border=border.color,lwd=control$lwd.rect)
    if(border.color!="darkgrey"){
      text(txt,x=info$x,y=y[1]+0.3*(y[N]-y[1]),col="black",cex=0.6)
    }
    
    if(!is.na(depth_zscore)){
      p=pnorm(depth_zscore,lower.tail=T)
      if(p<0.01){
        text(sprintf("%3.1f",mdep),x=info$x,y=y[1]+0.6*(y[N]-y[1]),col="red",cex=0.6)
        #text(sprintf("%2.1f",z),x=coords$b1,y=yy[1]+0.8*(yy[A]-yy[1]),col="red",cex=0.6)
      }
    }
  }
  arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)
  df
}

add_colored_vaf=function(pdx,##<< PDX object or list including details matrix
                         tree,
                         node,
                         control,
                         ...
){
  control=add_defaults(control,defaults=list( samples=c(),
                                              lwd.rect=1,
                                              min.depth=1,
                                              filter.on=NULL,
                                              mean.type="agg",
                                              cex.label=1,
                                              mode="recap",
                                              b.add.vaf.label=TRUE))
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(pdx,tree,node)
  ##browser()
  
  if(length(info$idx)==0){
    return(NULL)
  }
  
  if(!is.null(control$filter.on)){
    info$idx=info$idx[which(pdx$details[[control$filter.on]][info$idx])]
    if(length(info$idx)==0){
      return(NULL)
    }
  }
  
  samples=control$samples
  if(is.null(samples) || length(samples)==0){
    samples=info$samples
  }
  
  ##browser()
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(pdx$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(pdx$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
  }else{
    df=data.frame(mtr=pdx$mtr[info$idx,samples],
                  dep=pdx$dep[info$idx,samples],stringsAsFactors = FALSE)
  }
  df=cbind(df,pdx$details[info$idx,])
  if(control$mean.type!="agg"){
    df=df[which(df$dep>=control$min.depth),]
  }
  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  vaf=mean(df$vaf)
  
  
  N=dim(df)[1]
  if(N==0){
    return(df)
  }
  
  # cat("Fix me in add_color_vaf\n")
  r=0.3 ##0.25*length(tree$tip.label)/20 #0.3
  #rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,
  #     col=thiscol)
  
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  
  if(control$mean.type=="agg"){
    vaf=MTR/DEP
  }
  thiscol=getVafColor(vaf)$col
  if(control$mean.type=="agg" & DEP<control$min.agg.depth){
    thiscol="grey"
  }
  #thiscol2=getVafColor(MTR/DEP)$col
  
  rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,
       col=thiscol)
  if(control$mode=="recap"){
    if(MTR/DEP>0.01){
      txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
      if(control$b.add.vaf.label){
        text(txt,x=info$x,y=info$yb+0.3*(info$yt-info$yb),col="black",cex=0.6)
      }
    }
    rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,border="darkgrey",lwd=control$lwd.rect)
  }else{
    min.mean.vaf=0.45
    z=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
    z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
    z$p.value=max(z$p.value,z2$p.value)
    txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
    if(z$p.value<0.05){
      if(z$p.value<0.05/dim(tree$edge)[1]){
        border.color="red"
      }else{
        border.color="blue"
      }
    }else{
      border.color="darkgrey"
    }
    
    rect(xleft=info$x-r,xright=info$x+r,ybottom = info$yb, ytop=info$yt,border=border.color,lwd=control$lwd.rect)
    if(border.color!="darkgrey"){
      text(txt,x=info$x,y=info$yb+0.3*(info$yt-info$yb),col="black",cex=0.6)
    }
  }
  
  
  
  
  #arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)
  df
}











plotBars=function(x,y,radius,col,prop,border="black",width=1){
  #cat(prop,"\n")
  if(width<2){
    arrows(x0 = x-radius,y0=y,x1=x-radius+2*radius*prop[1],col="darkgrey",lend=2,length=0)
    arrows(x0 = x-radius+2*radius*prop[1],y0=y,x1=x-radius+2*radius,col=rgb(0.98,0.98,0.98),lend=2,length=0)
  }else{
    rect(xleft = x-radius,xright =x-radius+2*radius*prop[1],ybottom = y-width/2,ytop=y+width/2,border = NA,col="darkgrey")
    rect(xleft =  x-radius+2*radius*prop[1],xright =x+radius,ybottom = y-width/2,ytop=y+width/2,border = NA,col=rgb(0.98,0.98,0.98))
  }
  1
}

plotBars2=function(x,y,radius,col,prop,border="black",width=1){

    rect(xleft = x-radius,xright =x-radius+2*radius*1,ybottom = y-width/2,ytop=y+width/2,border = NA,col=rgb(red = prop,green=(1-prop),blue = 0.01))

  1
}


plotPie=function(x,y,radius,col,prop,llwd=0.5,border="black",width=NA){
  lims=par("usr")
  as=dev.size()
  asr=as[1]/as[2]
  yscale=asr*(lims[4]-lims[3])/(lims[2]-lims[1])
  prop=prop/sum(prop)
  cutpoint=c(0,cumsum(prop)*2*pi)

  N=2*pi/0.05
  n=ceiling(N*diff(cutpoint)/(2*pi))
  d=diff(cutpoint)/n
  if(length(prop)>1){
    for(i in 2:length(cutpoint)){
      polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              border=border,col=col[i-1],lwd=llwd)
    }
  }else{
    i=2
    polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            border=border,col=col[1],lwd=llwd)
  }
  yscale
}


##Gets unique colour pch combos and returns in dataframe with columns "col" and "pch"
get_color_pch_df=function(n){
  pch.list=c(18,17,16,15,0:6)
  if(n>length(pch.list)*8){
    stop("Too many colours requested")
  }
  cols=rep(RColorBrewer::brewer.pal(8,"Set1"),times=length(pch.list))
  pch=rep(pch.list,each=8)
  data.frame(col=cols,pch=pch,stringsAsFactors = FALSE)[1:n,]

}

get_qdf=function(values){
  if(length(values)>length(unique(values))){
    stop("get_qdf: please provide values without duplication")
  }
  cbind(data.frame(value=values,stringsAsFactors = FALSE),
        get_color_pch_df(length(values)))
}

plot_tree_labels_consequence=function(tree,details,consequences,
                                      query.allowed.df=get_qdf(consequences),
                          query.field="VC",
                          label.field="GENE",
                          cex.label=1){
  plot_tree_labels(tree,details,
                   query.allowed.df = query.allowed.df,
                   query.field=query.field,
                   label.field=label.field,
                   cex.label=cex.label)
}


plot_tree_labels=function(tree,pdx,
                          query.field="VC",
                          query.allowed.df=get_consequence_scheme(),
                          label.field="GENE",
                          cex.label=1){

  ##annot_function(pdx,tree,tree$edge[i,2],control=control)
  control=add_defaults(list(),defaults=list( query.field=query.field,
                                              label.field=label.field,
                                              query.allowed.df=get_consequence_scheme(),
                                              cex.label=cex.label,
                                              b.add.label=TRUE,
                                              b.add.marker=TRUE)
  )
  res=add_annotation(pdx,
                     tree,
                     add_simple_labels,control=control)

  with(control$query.allowed.df,legend("topleft",legend=value,col=col,pch=pch))
}

add_vertical_labels=function(
  pdx,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
  tree,##<< enhanced phylo returned from plot_tree
  node,##<< Node (see details)
  control,
  ... ##<< paremeters for points (not color)
){
  info=get_edge_info(pdx,tree,node)
  thisnode=node
  if(is.null(control$vars)){
    return(NULL)
  }
  vars=control$vars %>% filter(node==thisnode)
  #browser()
  if(dim(vars)[1]==0){
    return(NULL)
  }
  scheme=get_driver_scheme2()
  vars$group=gsub(":.+","",vars$label2)
  #col=scheme$colour[match(vars$label2,scheme$driver)]
  col=scheme$colour[match(vars$group,scheme$group)]
  pch=vars$pch
  
  vlabels=vars$label
  if(control$b.show.just.group){
    vlabels=unique(vars$group)
  }
  N=length(vlabels)
  
  col1=col
  col=unique(col)
  N1=length(col)
  if(N1>4){
    stop("too many distinct driver types <= 4. Please fix code..")
  }
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  ww=min(0.25,0.25*length(tree$tip.label)/20)   ##half Rect width   ## Or less if very small tree
  offsetodd=c(0,-0.5,0.5)*ww/0.25
  offseteven=c(0.25,-0.25,0.75,-0.75)*ww/0.25
  #idx2=info$idx.in.details
  
  if(N1 %% 2){
    offset=offsetodd[1:N1]
    
  }else{
    offset=offseteven[1:N1]
    ##oidx=(offset)
  }
  offset=sort(offset)
  #browser()
  
  if(N>0){
    yd=info$yt-info$yb
    if(N==1){
      #y=0.5*(info$yb+info$yt)
      y=info$yb+rbeta(1,shape1 = 5,shape2 = 5)*yd
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    
    if(control$b.add.label){
      text(rep(info$x+min(offset)-ww,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=control$cex.label,col=col1)
    }
    for(i in 1:N1){
      rect(xleft=info$x+offset[i]-ww,xright=info$x+offset[i]+ww,
           ytop =info$yt,ybottom = info$yb,col=col[i],border=NA)
    }
  }
  list(node=node,value=N)
}

plot_labelled_tree=function(tree,
                            pdx,
                            cv=c("missense","nonsense","ess_splice","frameshift","inframe","loh","cna"),
                            cex.annot.label=0.5,legpos="topleft",style="classic",genes=NULL,seed=-1,b.annot.just.shared=FALSE,col.label="black",b.show.just.group=FALSE){
  if(seed>0){
    set.seed(seed)
  }else{
    seed=sample(100000,1)+100000
    set.seed(seed)
  }
  muts=get_all_tree_drivers(pdx,cv = cv,genes=genes)
  print(muts)
  
  scheme=get_consequence_scheme(allowed.value=cv)
  colnames(scheme)[1]="cv"
  vars=muts %>% left_join(scheme,by="cv")
  ##browser()
  control=add_defaults(list(),defaults=list(vars=vars,
                                            cex.label=cex.annot.label,
                                            b.add.label=cex.annot.label>0,
                                            b.add.marker=TRUE,
                                            b.annot.just.shared=b.annot.just.shared,
                                            col.label=col.label,
                                            b.show.just.group=b.show.just.group))
  if(style=="classic"){
    annotfun=add_specified_labels
  }else{
    annotfun=add_vertical_labels
  }
  res=add_annotation(pdx,
                     tree,
                     annotfun,control=control)
  
  if(!is.null(legpos) && !is.null(vars) && dim(vars)[1]>0){
    ##leg=legend(legpos,mut.scheme$cv,col = mut.scheme$col,pch=mut.scheme$pch,pt.bg=mut.scheme$col,pt.cex = cex.marker)$rect
    ##browser()
    leg=with(unique(control$vars[,c("cv","col","pch")]),legend(legpos,legend=cv,col=col,pch=pch,pt.bg=col)$rect)
  }
  tree$seed=seed
  tree
}

plot_chromosome_annotated_tree=function(tree,pdx,pthreshold=0.01){

  chrs=sprintf("%s",c(1:22,"X","Y"))
  chr=match(pdx$details$Chrom,chrs)
  col_ctr=tabulate(chr,nbins=24)

  ptarget=col_ctr/sum(col_ctr)
  ptarget=(ptarget+0.0001)/sum(ptarget+0.0001)
  control=list(targetprop=ptarget,pcut=pthreshold)
  res=add_annotation(pdx,
                     tree,
                     add_chromosome_decomp,control=control)

  leg=legend("topleft",sprintf("chr%s",c(1:22,"X","Y")),col = get_chr_cols(),pch=15,
             pt.cex=2)$rect
}


#add_vaf_bar
#add_vof_pie
#add_label

##Need to hack at get_all_tree_drivers to make a satisfactory interface for specifying
## what drivers to annotate the tree with!
get_all_tree_drivers=function(pdx,drivers=get_drivers("drivers_red.bed"),
                              cv=c("missense","nonsense","ess_splice","frameshift","inframe","loh","cna"),genes=NULL){
  drivers=get_tree_drivers(pdx,drivers=drivers,cv = cv,genes=genes)

  drivers$label=with(drivers,sprintf("%s:%s",GENE,HGVS_PROTEIN))
  ##map copy number
  cna=rbind(get_cna_profile(pdx),get_loh_profile(pdx))[,c("label","profile","cv")]
  drivers=rbind(drivers[,c("label","profile","cv")],cna)
  #drivers$node=pdx$tree_ml$edge[pdx$summary$edge_ml[match(drivers$profile,pdx$summary$profile)],2]
  drivers$node=get_node_for_profile(pdx,drivers$profile)
  
  drivers$label2=gsub("_[A-Z]$","",drivers$label)
  drivers$label2=gsub("_[a-z]$","",drivers$label2)
  drivers
}

get_node_for_profile=function(pdx,profile){
  ##df and the tree are aligned.
  df=pdx$df$df
  if(any(pdx$df$samples!=pdx$tree_ml$tip.label)){
    stop("profile df is inconsistent with tree!")
  }
  df$node=pdx$tree_ml$edge[,2]
  df$node[match(profile,df$profile)]
}


get_tree_drivers=function(pdx,b_just_drivers=FALSE,genes=c("EIF6","SBDS","TP53","RPL5"),
                          cv=c("missense","nonsense","ess_splice","frameshift",
                               "inframe","start_lost","nc_ess_splice"),
                          drivers=drivers){
  res=pdx
  res$dat$details$ID=with(res$dat$details,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))

  defdrivers=drivers
  if(is.null(res$dat$details$color)){
    res$dat$details$color="black"
  }
  ##Restrict to sites matching drivers
  res$dat$details$profile=res$df$df$profile[res$summary$edge_ml]
  allprofiles=unique(res$dat$details$profile)
  ##browser()
  if(b_just_drivers){
    idx.keep=which(res$dat$details$ID %in% defdrivers$ID)
  }else{
    if(length(genes)>0){
      idx.keep=which(res$dat$details$GENE %in% genes  & res$dat$details$VC %in% cv)
    }else{
      idx.keep=which(res$dat$details$VC %in% cv)
    }
  }
  if(length(idx.keep)==0){
    cat("No drivers to add to tree!\n")
    return(data.frame(CHROM=character(),
                      POS=integer(),
                      REF=character(),
                      ALT=character(),
                      cv=character(),
                      GENE=character(),
                      HGVS_PROTEIN=character(),
                      VC=character(),
                      profile=character()))
  }
  #browser()
  drivers=res$dat$details[idx.keep, c("Chrom","Pos","Ref","Alt","color","GENE","VC","HGVS_PROTEIN","profile")]#
  colnames(drivers)[c(1:4,7)]=c("CHROM","POS","REF","ALT","cv")#
  drivers$label=with(drivers,sprintf("%s:%s",GENE,HGVS_PROTEIN))
  #drivers$node=res$tree_ml$edge[res$summary$edge_ml[match(drivers$profile,res$summary$profile)],2]
  drivers$node=get_node_for_profile(res,drivers$profile)
  drivers$label2=gsub("_[A-Z]$","",drivers$label)
  drivers$label2=gsub("_[a-z]$","",drivers$label2)

  drivers
}

get_cna_profile=function(pdx){
  if(length(pdx$meta$CNA)==0){
    return(NULL)
  }
  loh=sapply(pdx$meta$CNA,function(x){zeros=rep(0,length(pdx$df$samples))
  zeros[match(x$samples,pdx$df$samples)]=1;
  c(label=x$LABEL,profile=paste(zeros,collapse=""))
  })
  loh=as.data.frame(t(loh),stringsAsFactors=FALSE)
  loh$cv="cna"
  loh
}

get_loh_profile=function(pdx){
  if(length(pdx$meta$LOH)==0){
    return(NULL)
  }
  loh=sapply(pdx$meta$LOH,function(x){zeros=rep(0,length(pdx$df$samples))
  zeros[match(x$samples,pdx$df$samples)]=1;
  c(label=x$LABEL,profile=paste(zeros,collapse=""))
  })
  loh=as.data.frame(t(loh),stringsAsFactors=FALSE)
  ###browser()
  loh=loh[grep("1",loh$profile),]
  loh$label=gsub("^LOH_","",loh$label)
  loh$label=gsub("chr9UPD","9pUPD",loh$label)
  loh$cv="loh"
  loh
}

get_drivers=function(driver_file){
  defdrivers=read.table(driver_file,head=FALSE,stringsAsFactors = FALSE)
  colnames(defdrivers)=c("CHROM","POS","REF","ALT")
  defdrivers$ID=sprintf("%s:%s:%s:%s",defdrivers$CHROM,defdrivers$POS,defdrivers$REF,defdrivers$ALT)
  ##defdrivers=defdrivers[which(defdrivers$CHROM=="ZZZ"),] ## HACK to ignore driver file for now

  defdrivers
}

##Invokes high level pdx....
plot_basic_tree=function(pdx,label,style="classic",cex.annot.label=1,legpos="topleft",genes=NULL,cv=c(),seed=-1,b.annot.just.shared=FALSE,col.label="black",b.show.just.group=FALSE,...){
  tree=plot_tree(pdx$tree_ml,...);title(label)
  tree=plot_labelled_tree(tree,pdx,style=style,cex.annot.label=cex.annot.label,legpos=legpos,genes=genes,cv=cv,seed=seed,b.annot.just.shared=b.annot.just.shared,col.label=col.label,
                          b.show.just.group = b.show.just.group)
  tree
}


txtbox=function(xt, yt, txt,ccex=1,bcol="white",...){
  sw   <- strwidth(txt,cex=ccex)
  sh   <- strheight(txt,cex=ccex)
  frsz <- 0.3
  rect(
    xt - 2*sw*(0.5+ frsz),
    yt - sh*(0.5+ frsz),
    xt ,
    yt + sh*(0.5+ frsz),col=bcol
  )
  text(x=xt -sw*(0.5+ frsz),y=yt,label=txt,cex=ccex,...)
}

plot_chromosome_annotated_tree_pdx=function(pdx,label){
tree=plot_tree(pdx$tree_ml);title(label)
plot_chromosome_annotated_tree(tree,pdx$dat)
}


plot_vaf_tree=function(pdx,label,samples=NULL,cex.label=1,b.add.drivers=TRUE,filter.on=NULL,min.depth=1,mode="qc",b.plot.bars=TRUE,lwd=1,min.mean.vaf=0.35){
  tree=plot_tree(pdx$tree_ml,cex.label = cex.label,mar=c(1,1,2,3)+0.1)
  title(label)
  ##mtext(label,side=3,xpd=TRUE)
  idx=c()
  if(b.add.drivers){
    
    drivers=get_tree_drivers(pdx,drivers=get_drivers("drivers_red.bed"))
    colnames(drivers)[1:4]=c("Chrom", "Pos", "Ref", "Alt")
    scheme=get_consequence_scheme()
    colnames(scheme)[1]="cv"
    drivers=drivers %>% left_join(scheme,by="cv")
    drivers$ID=with(drivers,sprintf("%s:%s:%s:%s",Chrom,Pos,Ref,Alt))
    #browser()
    if(dim(drivers)[1]>0){
      pdx$dat$details=pdx$dat$details %>% left_join(drivers[,c("Chrom", "Pos", "Ref", "Alt","node","label", "label2","col","pch")],
                                                    by=c("Chrom", "Pos", "Ref", "Alt","node"))
      
      
    }
    #control=list(cex.label=0.5,samples=samples,filter.on=filter.on,min.depth=min.depth)
    #}else{
    # control=list(filter.on=filter.on,min.depth=min.depth)
  }else{
    idx=which(pdx$dat$details$keep_embryonic==1)
    pdx$dat$details$label=NA
    if(length(idx)>0){
      #pdx$dat$details$label[idx]=sprintf("emb%s",1:length(idx))
    }
  }
  control=list(cex.label=cex.label,samples=samples,filter.on=filter.on,min.depth=min.depth,mode=mode,b.plot.bars=b.plot.bars,lwd.rect=lwd,min.mean.vaf=min.mean.vaf)
  #browser()
  res=add_annotation(pdx$dat,tree,add_vaf,control=control)
  ##tree=plot_labelled_tree(tree,pdx,style="classic",cex.annot.label=1,legpos="topleft")
  tree$embryonic_idx=idx
  return(tree)
}

plot_color_vaf_tree=function(pdx,label,samples=NULL,cex.label=1,b.add.drivers=TRUE,
                             filter.on=NULL,min.depth=1,min.agg.depth=50,mean.type="agg",mode="recap",b.add.vaf.label=TRUE,legpos.vaf="left",...){
  tree=plot_tree(pdx$tree_ml,cex.label = cex.label,mar=c(1,1,2,3)+0.1)
  title(label)
  control=list(cex.label=0.5,
               samples=samples,
               filter.on=filter.on,
               min.depth=min.depth,
               mean.type=mean.type,
               min.agg.depth=min.agg.depth,
               mode=mode,
               b.add.vaf.label=b.add.vaf.label
  )
  #browser()
  res=add_annotation(pdx$dat,tree,add_colored_vaf,control=control)
  if(b.add.drivers){
    tree=plot_labelled_tree(tree,pdx,style="classic",cex.annot.label=0.6,...)
  }
  add_vaf_legend(legpos.vaf)
}

add_vaf_legend=function(legpos){
  if(!is.null(legpos)){
    vafc=getVafColor()
    legend(legpos,vafc$strVaf,col=vafc$colR,pch=15,pt.cex =3,title="Mean VAF")
  }
}

plot_len=function(pdx,label,edge.len.param="edge.length",cex.label=1,cex.annot.label=0.5){
  tree=plot_tree(pdx$tree_ml,cex.label = cex.label,mar=c(1,1,2,3)+0.1)
  title(label)
  control=list(el_param=edge.len.param,cex.annot.label=cex.annot.label)
  #browser()
  res=add_annotation(pdx$dat,tree,add_length,control=control)
  tree
}



getVafColor=function(vaf=-1,brewer.palette="YlOrRd"){
  colRange=brewer.pal(8,brewer.palette)
  #colorRampPalette(c("yellow", "red"),bias=10)(5)
  vafRange=c(0,0.02,0.05,seq(0.1,0.5,length.out = 5))
  strVaf=sprintf("%3.2f-%3.2f",vafRange[-length(vafRange)],vafRange[-1])
  if(!is.na(vaf) && vaf>=0){
    thiscol=colRange[findInterval(vaf,vafRange,rightmost.closed = TRUE,all.inside = TRUE)]
  }else{
    thiscol="grey"
  }
  list(colR=colRange,vafR=vafRange,strVaf=strVaf,col=thiscol)
}


plot_pooled_vaf_tree=function(pdx,label,cex.label=1,b.add.drivers=TRUE,lwd=1){
  plot_vaf_tree(pdx,label,samples=NULL,cex.label=cex.label,b.add.drivers=b.add.drivers,lwd=lwd)
}

plot_all_vaf=function(pdx,label,vtype="SNV",samples=colnames(pdx$dat$mtr),expected_vaf=0.35){
  pdx$dat$details$is_vtype=pdx$dat$details$TYPE %in% vtype
  ##Reorder so that samples are in tip order...
  samples=intersect(c(setdiff(samples,pdx$tree_ml$tip.label),setdiff(pdx$tree_ml$tip.label,"zeros")),samples)

  for(x in samples){
    plot_vaf_tree(pdx,
                  sprintf("%s: Annotated with VAF from %s\n Mean Depth=%3.2f",label,x,mean(pdx$dat$dep[,x],na.rm=T)),
                  samples=x,filter.on="is_vtype",min.depth=8,min.mean.vaf = expected_vaf)
  }
}

add_heatmap=function(tree,heatmap,heatvals=NULL,border="black",cex.label=2,pos=4,force.count=-1,txtcols=rep("black",dim(heatmap)[1])){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.05*ymax
  if(force.count<0){
    force.count=dim(heatmap)[1]
  }
  gap=ymax*tree$vspace.reserve/force.count
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-gap
    if(!is.na(border) && border=="none"){
      ##collapse 
      cc=unique(heatmap[i,])
      for(ccc in cc){
        idxx=which(heatmap[i,]==ccc)
        start=idxx[1]
        end=start
        for(j in idxx){
          if(j>(end+1)){
            rect(xleft=start-0.5,xright=end+0.5,ybottom = bot,ytop=top,col = heatmap[i,end],border=NA)
            start=j
            end=j
          }else{
            end=j
          }
        }
        rect(xleft=start-0.5,xright=end+0.5,ybottom = bot,ytop=top,col = heatmap[i,end],border=NA)
      }
    }else{
      rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd=0,ljoin=2)
    }
    
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      text(labels[i],x=par("usr")[1]+1,y=0.5*(top+bot),pos = pos,cex = cex.label,col=txtcols[i])
    }
    top=bot
  }
  tree
}

