library("rtreefit")
library("metafor")
source("plot_tree_annots_extra.R")
MEAN_AGE_AT_DELIVERY_DAYS=40*7 - 14 #Median gestational age=40 and subrtract 14 for LMP->fertilisation

## Call the rtreefit on the "PD" data structure and place results in the data structure in PD$fit$<method>$<nullmodel|altmodel>
wraptreefit=function(PD,niter=20000,b.fit.null=FALSE,method=c("poisson_tree"),cores=4,b_use_cache=TRUE,stan_control=list(adapt_delta=0.95)){
  res=refit(PD,niter=niter,b.fit.null=b.fit.null,method=method,cores=cores,b_use_cache=b_use_cache,stan_control=stan_control)
  if(is.null(PD$fit)){
    PD$fit=list()
    PD$fit[[method]]=list()
  }
  if(b.fit.null){
    PD$fit[[method]]$nullmodel=res
  }else{
    PD$fit[[method]]$altmodel=res
  }
  PD
}

## Calls rtreefit on the "PD" data structure and returns summary information on the timing of branches specified in the PD$nodes data structure.
## For b.fit.null=FALSE the alternative model is fitted for nodes with PD$nodes$status>=0, but timing is reported for all nodes specified in PD$nodes$node. 
refit=function(PD,niter=20000,b.fit.null=FALSE,method=c("poisson_tree"),cores=4,b_use_cache=TRUE,stan_control=list(adapt_delta=0.95)){
  tree=get_tree_for_fit(PD)
  if(b.fit.null){
    switch_nodes=c()
    if(dim(PD$nodes)[1]==0){
      split_nodes=c()
    }else{
      split_nodes=PD$nodes$node
    }
    
  }else{
    switch_nodes=PD$nodes$node[which(PD$nodes$status>=0)]
    split_nodes=PD$nodes$node
  }
  ptree=fit_tree(tree=tree,switch_nodes = switch_nodes,xcross = c(),niter = niter,cores=cores,model  = method,split_nodes = split_nodes,stan_control = stan_control)
  fields=c("ultratree","lambda","split_nodes","upper_node_lims","lower_node_lims","fullres")
  res=ptree[fields]
  res$summary=with(res,
                   cbind(PD$nodes,chknode=split_nodes,
                         as.data.frame(lower_node_lims) %>% (function(x){colnames(x)=sprintf("lower_%s",colnames(x));x}),
                         as.data.frame(upper_node_lims) %>% (function(x){colnames(x)=sprintf("upper_%s",colnames(x));x})
                   ))
  res$param=list(niter=niter,switch_nodes=switch_nodes,split_nodes=split_nodes,tree=tree)
  res
}

## Retrieves the mutation counte per branch filtered for e.g. excluding CN/LOH regions and also returns the branch specific estimated sensitivity and 
## copies a dataframe (agedf) that has the age at sampling for each tip.  The returned extended Phylo object will be the input tree for rtreefit::fit_tree.
get_tree_for_fit=function(PD,type="local",suffix="",atype="hybrid.multi"){
  thistree=PD$pdx$tree_ml
  thistree$agedf=PD$pdx$agedf
  if(type %in% c("global","local")){
    thistree$edge.length=thistree[[sprintf("el.snv.%s.filtered%s",type,suffix)]]
    thistree$sensitivity=thistree[[sprintf("per.branch.sensitivity.%s",atype)]]
  }else{
    stop("unsupported")
    thistree$edge.length=thistree[[sprintf("el.snv%s",suffix)]]
    thistree$sensitivity=thistree[[sprintf("sensitivity.snv%s",suffix)]]
  }
  thistree$sensitivity=ifelse(thistree$sensitivity>0.99,0.99,thistree$sensitivity)
  thistree$agedf$age=thistree$agedf$age_at_sample_pcy
  thistree
}


get_stan_result_df=function(PD,model="poisson_tree",bIsNull=FALSE){
  if(!bIsNull){
    fitres=PD$fit[[model]]$altmodel
    tree=get_colored_markup_tree2(PD$fit[[model]]$altmodel$ultratree,PD$nodes[PD$nodes$status>=0,])
    nodes=rbind(data.frame(node=-1,driver="WT",status=1,child_count=-1,driver2="WT",driver3="WT",stringsAsFactors = FALSE),
                PD$nodes %>% dplyr::select( node, driver ,status, child_count, driver2 ,driver3))
    nodes=nodes[which(nodes$status>=0),]
    tipstatus=tree$rate[which(tree$edge[,2]<=length(tree$tip.label))]
    nodes$patient=PD$patient
    nodes$type="local"
    M=dim(nodes)[1]
    nodes$colony_count=sapply(1:M,function(i){n=length(which(tipstatus==i));if(i==1){n-1}else{n}})
  }else{
    fitres=PD$fit[[model]]$nullmodel
    nodes=data.frame(node=-1,driver="Global",status=1,child_count=-1,driver2="Global",driver3="Global",stringsAsFactors = FALSE)
    nodes$patient=PD$patient
    nodes$type="local"
    nodes=data.frame(node=-1,driver="",status=1,child_count=-1,driver2="",driver3="",stringsAsFactors = FALSE)
  }
  for(x in c("mean","sd","lb","ub","median")){
    nodes[[x]]=fitres$lambda[[x]]
  }
  type="local"
  for(x in c("mean","sd","lb","ub","median")){
    if(type!="snv"){
      nodes[[sprintf("%s_rescaled",x)]]=fitres$lambda[[x]]*PD[[sprintf("%sx.correction2",type)]]
      nodes[["correction"]]=PD[[sprintf("%sx.correction2",type)]]
    }else{
      stop("unsupported path")
    }
  }
  post.dist=rstan::extract(fitres$fullres,par="lambda")
  N=dim(post.dist$lambda)[1]
  nodes[["p_lt_wt"]]=sapply(1:dim(nodes)[1],function(i) mean(post.dist$lambda[,i]<post.dist$lambda[,1]))
  nodes[["p_lt_wt"]]=ifelse(nodes[["p_lt_wt"]]<1/N,1/N,nodes[["p_lt_wt"]])
  nodes[["rate_diff"]]=sapply(1:dim(nodes)[1],function(i) mean(post.dist$lambda[,i]-post.dist$lambda[,1]))
  if(!is.null(fitres$k) && length(fitres$k$mean)>0){
    for(x in c("mean","sd","lb","ub","median")){
      nodes[[sprintf("k_%s",x)]]=fitres$k[[x]]
    }
  }
  if(!is.null(fitres$p) && length(fitres$p$mean)>0){
    for(x in c("mean","sd","lb","ub","median")){
      nodes[[sprintf("p_%s",x)]]=fitres$p[[x]]
    }
  }
  nodes
}

## Extra visualisation code
get_colored_markup_tree=function(tree,nodes){
  cols=RColorBrewer::brewer.pal(10,"Paired")
  tmp=markup_tree(tree,nodes)
  tree$color=cols[tmp$type+1]
  tree
}
get_colored_markup_tree2=function(tree,nodeinfo){
  cols=RColorBrewer::brewer.pal(10,"Paired")
  nodes=nodeinfo$node[which(nodeinfo$status>=0)]
  ##The following returns a list
  tmp=markup_tree(tree,nodes)
  excluded=which(nodeinfo[which(nodeinfo$status>=0),"status"]==0)
  tree$color=ifelse(tmp$type %in% c(excluded),"grey",cols[tmp$type+1])
  tree$rate=tmp$type+1
  tree
}


markup_tree=function(tree,
         switch_nodes##<< child node ids of edges with a switch in lambda.
){
  ###
  parentidx=match(tree$edge[,1],tree$edge[,2])-1
  parentidx=ifelse(is.na(parentidx),-1,parentidx)
  N=length(tree$tip.label)
  #Index of tips
  tipidx=match(1:N,tree$edge[,2])
  edges=tree$edge
  m=tree$edge.length
  ##Parents of every edge including the edge itself [order from root to tip]
  parentlist=lapply(edges[,2],function(node) {
    ##We need to make sure that parent edges are ordered from root to tip
    parents=rev(get_parents(node,edges))##
    if(length(parents)>0){
      #parents
      match(parents,edges[,2])
    }else{
      parents
    }
  }
  )
  ##Index of applicable lambda rate...
  rates=rep(0,length(m))
  if(length(switch_nodes)>0){
    for(i in 1:length(switch_nodes)){
      snode=switch_nodes[i]
      ##We want to include the parent node as well
      nodes=c(snode,get_all_node_children(node = snode,tree = tree))
      ##Check that nodes in correct order.
      if(length(unique(rates[match(nodes,tree$edge[,2])]))>1){
        stop("switch nodes provided in wrong order!")
      }
      ##browser()
      idxanc=match(nodes,tree$edge[,2])
      if(length(idxanc)>0){
        rates[idxanc]=i
      }
    }
  }
  #Rates of parent node.
  ratesp=rates[match(tree$edge[,1],tree$edge[,2])]
  #If node has no parent assume the wild type rate
  ratesp[which(is.na(ratesp))]=0
  #Index of edges where rate switches occur.
  #if(b_pool_rates){
  nlambda=max(rates)+1
  list(type=rates,ntype=nlambda)
}
get_parents=function(node,mut_table,exclude_root=TRUE){
  idx=which(mut_table[,2]==node)
  parents=node##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=mut_table[idx,1]
    parents=c(parents,parent)
    idx=which(mut_table[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)]
  }else{
    parents
  }
}

### Call plot fitted model 
do_lambda_plot=function(fl,PDD,ftype="local",p.low.or.high="high"){
  res=get_lambda_plot(fl,PDD,ftype,p.low.or.high)
  res$fl
}

##  Do forest plot - works for single samples or multiple samples 
get_lambda_plot=function(fl,# Result of get_stan_result_df on PD 
                         PDD,# List of per patient PD objects
                         ftype="local",p.low.or.high="high",b.do.plot=TRUE,b.add.telomeres=TRUE){
  N=dim(fl)[1]
  ##Override driver column
  fl$driver=fl$driver3
  coldf=getcolor_df(data.frame(patient=sapply(PDD,function(x) x$patient) ))
  
  ##Meta-Analysis
  tmp=fl %>% filter(driver=="WT")
  wtrma=summary(rma(tmp$mean_rescaled,tmp$sd_rescaled))
  if(dim(tmp)[1]>1){
    print(wtrma)
    rate_meta_summary=sprintf("%3.2f (%3.2f-%3.2f)",wtrma$beta,wtrma$ci.lb,wtrma$ci.ub)
  }else{
    rate_meta_summary=sprintf("%3.2f (%3.2f-%3.2f)",tmp$mean_rescaled,tmp$lb_rescaled,tmp$ub_rescaled)
  }
  fl=fl %>% left_join(coldf,by="patient")
  fl$patient2=fl$patient
  fl=fl %>% mutate(patient=ifelse(driver=="WT",patient,""))
  MUSTARD="#ffdb58"
  N2=length(which(fl$driver!="WT"))
  fl=fl %>% mutate(col=ifelse(grepl("JAK2",driver),"red",ifelse(grepl("WT",driver),"grey",MUSTARD))) %>%
    mutate(col2=ifelse(grepl("JAK2",driver),"red",ifelse(grepl("WT",driver),"grey","black")))
  if(p.low.or.high=="high"){
    pp=ifelse(fl$patient!="","",ifelse(fl$p_lt_wt<0.025,ifelse(fl$p_lt_wt<0.025/N2,"**","*"),""))
    extra="Posterior Prob. > Wild type rate"
    label2="Prob. > WT"
  }else{
    pp=ifelse(fl$patient!="","",ifelse(fl$p_lt_wt<0.001,"<0.001",sprintf("%4.3f",fl$p_lt_wt)))
    extra="Posterior Prob. <= Wild type rate"
    label2="Prob. < WT"
  }
  fl$p_lt_wt[grepl("WT",fl$driver)]=NA
  res=list(fl=fl,rmeta=wtrma,rate_meta_summary=rate_meta_summary)
  
  
  if(!b.do.plot){
    return(res)
  }
  
  bwidth=0.12
  telomere_panal_start=30
  tbs=telomere_panal_start
  if(!b.add.telomeres){
    tpl=0 ##telomere width
    tbs=35
  }else{
    tpl=12
    ##Pull in all the telomere data..
    ##Get sample per timepoint info.
    df=do.call("rbind",lapply(PDD,function(PD) PD$pdx$agedf %>% filter(tip.label!="zeros")))
    df$driver=df$driver3
    dfa=df %>% group_by(patient,age_at_sample_pcy) %>% summarise(N=n())
    dfa=dfa[order(dfa$patient,dfa$age_at_sample_pcy),]
    dfa$sample_idx=do.call("c",lapply(unique(dfa$patient),function(x) 1:length(which(dfa$patient==x))))
    dfx=df %>% inner_join(dfa[,-3],by=c("patient","age_at_sample_pcy"))
    #browser()
    dfx=dfx %>% left_join(dfx %>% group_by(patient) %>% summarise(offset=median(unique(sample_idx-1))),by=c( "patient"))
    df2=df %>% group_by(patient,driver) %>% summarise(mtelo=mean(telo_mean_length))
  }
  
  plot(NULL,xlim=c(-10,tbs+3*tpl),ylim=c(0,N+2),xlab="",yaxt="n",ylab="",main="",xaxt="n",bty="n")
  
  with(fl,rect(xleft = lb_rescaled,xright=ifelse(ub_rescaled>tbs,tbs,ub_rescaled),ybottom = (N:1)-bwidth,ytop =(N:1)+bwidth, col=col,border=NA))
  with(fl,points(y=N:1,x=median_rescaled,cex=1,pch=15))
  with(fl,text(y=N:1,x=rep(7,N),labels = sprintf("%s",driver),col="black",pos = 2))
  with(fl,text(y=N:1,x=rep(-10,N),labels = patient,col="black",pos = 4,offset=0.5))
  with(fl,text(y=N:1,x=rep(tbs,N),labels = pp,col="black",pos = 2))
  with(fl,text(y=N:1,x=rep(10,N),labels = colony_count,col="black",pos = 2,cex=1))
  abline(v=7)
  axis(side=1,at=c(-10,tbs+3*tpl),labels=rep("",2),lwd.ticks = 0)
  axis(side = 3,at=c(-10,tbs+3*tpl),labels=rep("",2),lwd.ticks = 0)
  abline(v=-10)
  
  axis(side = 1,at=seq(10,tbs-5,5))
  segments(x0=wtrma$b,y0=-10,y1=N+0.5,lty="dotted")
  text(x=wtrma$b,y = N+0.5,labels =rate_meta_summary ,pos = 3,cex=1,offset=0.2)
  segments(x0=tbs+tpl*0,y0=-10,y1=N+0.5)
  pchh=c(5,2,6)
  pchm=c(23:25)
  abline(v=10)
  abline(v=tbs)
  text(y=N+2,x=20,labels="Mutation Rate (SNV/Year)",pos=NULL,cex=1.2)
  text(y=N+2,x=8.5,labels="N",pos=NULL,cex=1.2)
  devnull=sapply(which(nchar(fl$patient)>1),function(i) segments(x0=-10,x1=tbs+3*tpl,y0=N-i+1+0.5,lwd=0.5))
  if(b.add.telomeres){
    abline(v=tbs+3*tpl)
    segments(x0=tbs+tpl*1,y0=-10,y1=N+0.5)
    for(k in 0:2){
      segments(x0=tbs+tpl*k,y0=-10,y1=N+0.5)
      text(y=N+1,x=tbs+tpl*k+1,labels=sprintf("Timepoint %s",k+1),col="black",pos=4)
      for(i in 1:N){
        tmp=dfx %>% filter(patient==fl$patient2[i],driver==fl$driver[i],sample_idx==k+1)
        points(x=tbs+tpl*k+(tpl/8)*tmp$telo_mean_length/1000,y=N-i+1+0.25*(runif(dim(tmp)[1])-0.5),pch=19,cex=0.5,col=fl$col[i])
        points(x=tbs+tpl*k+(tpl/8)*mean(tmp$telo_mean_length/1000),y=N-i+1,cex=1.2,pch=23,col="black",bg=fl$col[i])
      }
      tscale=seq(0,6,2)
      axis(side = 1,at=tbs+tpl*k+tscale*(tpl/8),labels=tscale)
    }
    ps=tbs+3*tpl+1
    ns=ps+tpl
    text(y=N+2,x=tbs+(tpl/2)*3,labels="Telomere Length(kb)",pos=NULL,cex=1.2)
  }
  return(res)
}


collate_driver_info=function(PDD,treemodel="poisson_tree",b.is.null=FALSE){
  ##Does the messy collation of node timings
  type=ifelse(b.is.null,"nullmodel","altmodel")
  #cat("WARNING: Review the labelling of the output\n")
  inf1=do.call("rbind",lapply(PDD,function(x){
    out=x$fit[[treemodel]][[type]]$summary %>% mutate(patient=x$patient)
    out=out %>% left_join(add_cumulative_mutcount(x,x$fit[[treemodel]][[type]]$summary$node))
    idx=match("patient",colnames(out))
    cbind(data.frame(patient=out[,idx],stringsAsFactors = FALSE),out[,-idx])
    ##Add
  }
  ))
  
  ##Make driver correspond to events on the current node
  inf1$driver=gsub(":.*","",inf1$driver3)
  ##Rename 9pUPD_A etc
  inf1$driver=gsub("_[A-Z]$","",inf1$driver)
  inf1$driver=gsub("_[a-z]$","",inf1$driver)
  #Filter to non-private branches that are of interest (non-nuisance)
  inf2=inf1 %>% filter(status==1 | child_count>1)
  
  #The driver scheme maps specific mutations to drivers..
  ds=get_driver_scheme()
  ##The following is non-missing for copy number events
  inf2$group=ds$group[match(inf2$driver,ds$group)]
  ##Specific genes are filled in order or priority
  ## The driver column can have multiple drivers - so we prioritise
  inf2$group[grep("CBL",inf2$driver)]="CBL"
  inf2$group[grep("TET2",inf2$driver)]="TET2"
  inf2$group[grep("PPM1D",inf2$driver)]="PPM1D"
  inf2$group[grep("DNMT3A",inf2$driver)]="DNMT3A"
  inf2$group[grep("JAK2",inf2$driver)]="JAK2"
  
  pti2=PD$pdx$meta$COLONY_TIMEPOINTS
  rti=PD$pdx$meta$RECAP_TIMEPOINTS
  
  pti=rbind(as.data.frame(pti2[,c("patient","age_at_sample_pcy")]),as.data.frame(rti[,c("patient","age_at_sample_pcy")]))
  df=pti %>% group_by(patient) %>% summarise(N=n(),max_age_at_sample=max(age_at_sample_pcy),min_age_at_sample=min(age_at_sample_pcy))

  inf= inf2 %>% left_join(df,by="patient") %>% left_join(pti2 %>% group_by(patient,age_at_diagnosis_pcy) %>% summarise(totcolonies=sum(N)),by="patient")
  
  inf=inf[order(inf$max_age_at_sample,inf$patient,inf$upper_median),]
  
  ##The following reinstate 1q+ and CBL as more specific classes (will potentially break with updated driver_scheme!)
  ##Should update driver_scheme instead..
  ds$group[grepl("1q+",ds$driver)]="1q+"
  ds$colour[grepl("1q+",ds$driver)]="brown"
  ds$group[grepl("CBL",ds$driver)]="CBL"
  inf$group[grepl("1q+",inf$driver)]="1q+"
  inf$col=ds$colour[match(inf$group,ds$group)]
  
  inf=inf %>% filter(!is.na(col))
  inf
}

## Plot acquisition time line for single or multiple patients.
do_acquisition_plot2b=function(inf,## Result of calling collate_driver_info(PDD)
                               pti,## Colony timepoints e.g. PD$pdx$meta$COLONY_TIMEPOINTS
                               rti,## Recapture timepoints e.g. PD$pdx$meta$RECAP_TIMEPOINTS
                               extra="",xmax=100,b.add.text=FALSE){
  inf2=inf
  patients=unique(inf$patient)
  N=length(patients)
  ## Gestational Age at Birth
  ab=MEAN_AGE_AT_DELIVERY_DAYS/365.25
  ab2=(MEAN_AGE_AT_DELIVERY_DAYS+14)/365.25
  #manipulate ranges so that we can extend pre-birth period...
  for(field in c("lower_median","lower_lb95","lower_ub95","upper_median","upper_lb95","upper_ub95","max_age_at_sample","age_at_diagnosis_pcy")){
    inf[[field]]=inf[[field]]-ab
    inf[[field]]=ifelse(inf[[field]]<0,xmax*0.1*inf[[field]]/ab2,inf[[field]])
  }
  pti$age_at_sample_pcy=pti$age_at_sample_pcy-ab
  rti$age_at_sample_pcy=rti$age_at_sample_pcy-ab
  
  par(mar=c(6,6,2,2)+0.1)
  plot(NULL,xlim=c(-xmax*0.1,xmax*1.15),ylim=c(0,3*length(patients)),yaxt="n",xlab="",ylab="",bty="n",xaxt="n")
  if(!is.null(extra)){
    title(sprintf("Age of Acquisition of Drivers%s",extra),cex=2)
  }
  axis(side = 1,at = seq(0,90,10),labels = seq(0,90,10))
  axis(side =1, at = -xmax*0.1+xmax*0.1*c(0,26/40,1),labels=c("0","26",""),cex.axis=0.6,las=1)
  axis(side =1, at = -xmax*0.1+xmax*0.1*c(13/40),labels=c("13"),cex.axis=0.6,las=1)
  mtext(side = 1,line = 4,text = "Age\n(Years)",at = xmax/2)
  mtext(side = 1,line = 4,text = "Gestational Age\n(Weeks)",at = -0.1*xmax)
  k=1
  width=0.2
  xgap=0.3##0.5*(0.5-width)
  kk=0
  pch.sample=2
  pch.diagnosis=18
  yb=0
  
  for(patient in patients){
    idx=which(inf$patient==patient)
    if(length(idx)>1){
      gap=3/(length(idx)+1)
    }else{
      gap=1.5
    }
    yb=3*kk+gap
    rect(xleft=0,xright=inf$max_age_at_sample[idx[1]],ybottom=3*kk+xgap,ytop=3*(kk+1)-xgap,col="lightgray",border=NA,lwd=1)
    rect(xleft=-0.1*xmax,xright=0,ybottom=3*kk+xgap,ytop=3*(kk+1)-xgap,col="lightpink",border=NA,lwd=1)
    
    for(j in idx){
      rect(xleft = inf$lower_median[j],xright=inf$upper_median[j],ybottom = yb-width,ytop=yb+width,border=NA,col=inf$col[j])
      segments(x0=inf$lower_lb95[j],x1=inf$lower_ub95[j],y0=yb,col="black",lwd=1)
      segments(x0=inf$upper_lb95[j],x1=inf$upper_ub95[j],y0=yb,col="black",lwd=1)
      if(b.add.text){
        text(x=inf$max_age_at_sample[idx[1]]+5,y=yb,labels = inf$driver3[j],pos = 4)
      }
      yb=yb+gap
    }
    ymid=3*kk+1.5
    mtext(text = patient,side = 2,line=1,at=ymid,las=2)
    points(y=ymid,x=inf$age_at_diagnosis_pcy[j],pch=pch.diagnosis,col="black",cex=1.5)
    tmp=sort(unique(rti$age_at_sample_pcy[which(rti$patient==patient)]))
    if(!is.null(tmp)){
      segments(x0=tmp,y0=rep(ymid,length(tmp))-0.2,y1=rep(ymid,length(tmp))+0.2,col="blue",lwd=2,lend=2)
    }
    tmp=sort(unique(pti$age_at_sample_pcy[which(pti$patient==patient)]))
    segments(x0=tmp,y0=rep(ymid,length(tmp))-0.2,y1=rep(ymid,length(tmp))+0.2,col="red",lwd=2,lend=2)
    kk=kk+1
  }
  ##Add a manual legend
  if(!b.add.text){
    tmp=inf %>% group_by(group,col) %>% summarise(N=n())
    yloc=3*kk/4 +3
    #xxl=90
    xxl=0.85*xmax
    offset=5
    for(i in 1:dim(tmp)[1]){
      rect(xleft=xxl,xright=xxl+offset,ybottom=yloc+i,ytop=yloc+i+2*width,col=tmp$col[i],border=NA)
      text(x=xxl+offset,y=yloc+i+width,labels = tmp$group[i],pos=4,offset=1)
    }
    yloc=3*kk/4
    points(x=xxl+0.5*offset,y=yloc,pch=pch.diagnosis,cex=1.5)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Diagnosis")
    yloc=yloc-2
    segments(x0=xxl+0.5*offset,y0=yloc-0.2,y1=yloc+0.2,col="red",lwd=2,lend=2)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Colony Sampling")
    yloc=yloc-2
    segments(x0=xxl+0.5*offset,y0=yloc-0.2,y1=yloc+0.2,col="blue",lwd=2,lend=2)
    text(x=xxl+offset,y=yloc,pos=4,offset=1,labels="Recapture Sampling")
  }
  abline(v=0,lty="dotted")
  inf2
}

getcolor_df=function(mut_count){
  col.df=data.frame(patient=unique(as.character(mut_count$patient)),stringsAsFactors = FALSE)
  col.df$col=ggplot_color(dim(col.df)[1])#
  col.df
}

ggplot_color=function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

add_cumulative_mutcount=function(PD,nodes){
  tree=PD$pdx$tree_ml
  tree$edge.length=PD$localx.correction2*tree$el.snv.local.filtered/tree$per.branch.sensitivity.hybrid.multi
  nh=nodeHeights(tree)
  mutcount=matrix(nh[match(nodes,tree$edge[,2]),],ncol=2)
  colnames(mutcount)=c("lower_mutcount_adj","upper_mutcount_adj")
  tree$edge.length=tree$el.snv.local.filtered
  nh=nodeHeights(tree)
  mutcount2=matrix(nh[match(nodes,tree$edge[,2]),],ncol=2)
  colnames(mutcount2)=c("lower_mutcount","upper_mutcount")
  cbind(data.frame(node=nodes),as.data.frame(mutcount),as.data.frame(mutcount2))
}

## This is called by the tree annotation code. Note that this version only shows drivers that are listed in the files driver_scheme.txt and driver_groups.txt. 
get_driver_scheme=function(){
  driver.scheme=read.table("driver_scheme.txt",head=T,stringsAsFactors = FALSE,sep="\t")
  driver.group=read.table("driver_groups.txt",head=T,stringsAsFactors = FALSE,sep="\t",comment.char = "")
  driver.scheme=driver.scheme[,c("driver","number")] %>% inner_join(driver.group,by="number")
  driver.scheme
}

## Expand short branches to better visualise the top of the tree.
expand_short_branches=function(pdx, ## e.g. PD$pdx
                               prop=0.05  ## Minimum branch length as a proportion of the height of the tree.
                               ){
  tree=plot_tree(pdx$tree_ml,b_do_not_plot = TRUE)
  minmuts=round(prop*max(tree$coords$a1))
  pdx$tree_ml$color=ifelse(pdx$tree_ml$edge.length<minmuts,"grey","black")
  pdx$tree_ml$edge.length.unexpanded=pdx$tree_ml$edge.length
  pdx$tree_ml$edge.length=ifelse(pdx$tree_ml$edge.length<minmuts,minmuts,pdx$tree_ml$edge.length)
  pdx$minmuts=minmuts
  pdx
}

###  Copy number timing code is in cna.R


