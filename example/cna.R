library("data.table")
library("rstan")
CACHE="."


## Collate ASCAT identified copy number neutral LOH events and attempt to refine the timing of the events based on 
## the counts of heterozygous and homozygous somatic mutations
get_hethom=function(PD,refcolony,b_do_plot=FALSE){
  loh=PD$pdx$meta$LOH[which(sapply(PD$pdx$meta$LOH,function(x) x$ploidy)==2)]
  if(length(loh)==0){
    return(NULL)
  }
  if(exists("CENSOR_GERMLINE")){
  }else{
  baff=get_cnfile(PD,refcolony,type="baf")
  bafref=read_tsv(baff)
  bafrefhet=bafref[!is.na(bafref$BAF) & bafref$BAF>0.2 & bafref$BAF<0.8,]
  }
  ##Here we are only interested in copy number neutral LOH
  message("Reading local mutation model from ",MUTCOUNTBIN," \n")
  bins=readRDS(MUTCOUNTBIN)
  
  hethomres=lapply(loh,function(x){
    
    label=sprintf("%s (%s:%s-%s)",x$LABEL,x$chr,x$start,x$end)
    ##representative sample
    colony=x$samples[1]
    ascatfile=get_cnfile(PD,colony,type="ascat")
    ## Get the filespecs..  Checks for a local copy then looks on the servers..
    if(exists("CENSOR_GERMLINE")){
      Pos1=x$start
      Pos2=x$end
      if(b_do_plot){
        plot_cn(get_cn(cn_summary_file = ascatfile),label = sprintf("%s: Representative of %d colonies : %s",colony,length(x$samples),label))
        warning("Germline info is censored so not plotting BAFs..")
        b_do_plot=FALSE
      }
    }else{
      
      baffile=get_cnfile(PD,colony,type="baf")
      baf=read_tsv(baffile);
      ##  Here we refine the breakpoints
      baf=baf[which(!is.na(baf$BAF) & baf$Chromosome==x$chr),]
      bafhet=bafrefhet[which(bafrefhet$Chromosome==x$chr),]
      baf=baf[baf$Position %in% bafhet$Position,]
      mbaf=frollmean(2*abs(baf$BAF-0.5),20,align = "right")
      mbaf2=frollmean(2*abs(baf$BAF-0.5),20,align = "left")
      Pos1=baf$Position[min(which(mbaf2>0.8 & baf$Position>=x$start))]
      Pos2=baf$Position[max(which(mbaf>0.8 & baf$Position<=x$end))]
      cat("Pos1=",Pos1,",Pos2=",Pos2,"\n")
      x$start=Pos1
      x$end=Pos2
    }
    
    hethom=get_loh_tree(PD,x$node,x)
    hethom$label=x$LABEL
    hethom$node=x$node
    hethom$het.sensitivity=PD$pdx$tree_ml$sensitivity.snv[match(x$node,PD$pdx$tree_ml$edge[,2])]
    hethom$chr=x$chr
    hethom$start=x$start
    hethom$end=x$end
    
    if(b_do_plot){
      kable(hethom$inf, "html") %>% kable_styling("striped")
      plot_cn(get_cn(cn_summary_file = ascatfile),label = sprintf("%s: Representative of %d colonies : %s",colony,length(x$samples),label))
      
      plot(baf$Position/1e6,baf$BAF,pch=19,cex=0.5,xlab="Position(Mb)",ylab="BAF");title(sprintf("%s: BAF: Representative of %d colonies : %s",colony,length(x$samples),label))
      lines(baf$Position/1e6,mbaf,col="magenta")
      abline(v=hethom$inf$Pos[hethom$inf$Chrom==x$chr & hethom$inf$het]/1e6,col="green")
      abline(v=hethom$inf$Pos[hethom$inf$Chrom==x$chr & hethom$inf$hom]/1e6,col="blue")
      abline(v=c(x$start/1e6,x$end/1e6),col="red")
      abline(v=c(Pos1/1e6,Pos2/1e6),col="red",lty="dashed")
      legend("bottomright",c("HET","HOM","LOH Boundary","Rolling Homozygosity (20 SNPs)"),col=c("green","blue","red","magenta"),lwd=2,bg = "white",cex = 0.8)
    }
    idx=which(bins$chr==x$chr & bins$start>=(1e5*floor(x$start/1e5)) & bins$end<=(1e5*floor(x$end/1e5)))
    hethom$kb=length(idx)*1e5
    hethom$count_in_bin=sum(bins$count[idx])
    hethom$count_se=sqrt(sum(bins$count[idx]))
    hethom$pmut=sum(bins$prob[idx])
    ##Estimate of uncertainty in local density
    hethom$pmut_se=hethom$pmut/sqrt(hethom$count_in_bin)
    ##Use the null model estimate for the ultrametric tree
    ut=PD$fit$poisson_tree$nullmodel$ultratree
    T=ut$edge.length[which(ut$edge[,2]==x$node)]
    ##Use null model lambda
    lambda=PD$fit$poisson_tree$nullmodel$lambda$mean
    dat=list(N=hethom$nhom,M=hethom$nhet,phom=0.5,phet=hethom$het.sensitivity,L=lambda*T*hethom$pmut,Lstd=lambda*T*hethom$pmut_se)
    ##Fit stan model...
    stanr=stan("stanmodels/cntime.stan", data=dat,chains = 3,iter=20000,cores=3)
    sumc=rstan::summary(stanr,probs = c(0.025,  0.50, 0.975))$summary
    for( z in colnames(sumc)){
      hethom[[sprintf("x%s",z)]]=sumc["x",z]
    }
    for( z in colnames(sumc)){
      hethom[[sprintf("l%s",z)]]=sumc["l",z]
    }
    hethom
  })
  othernames=c("inf")
  keepnames=setdiff(names(hethomres[[1]]),othernames)
  list(hethom=do.call("rbind",lapply(hethomres,function(x) as.data.frame(x[keepnames]))),
       other=do.call("rbind",lapply(hethomres,function(x) as.data.frame(x[othernames]))))
}

## Collate ASCAT identified duplications and attempt to refine the timing of the events based on 
## counts of somatic variants with VAF=1/3 and 2/3.
get_cnainf=function(PD,refcolony,b_do_plot=FALSE){
  if(length(PD$pdx$meta$CNA)==0){
    return(NULL)
  }
  bins=readRDS(MUTCOUNTBIN)
  if(exists("CENSOR_GERMLINE")){
  }else{
    baff=get_cnfile(PD,refcolony,type="baf")
    bafref=read_tsv(baff)
    bafrefhet=bafref[!is.na(bafref$BAF) & bafref$BAF>0.2 & bafref$BAF<0.8,]
  }
  
  idxx=grep("\\+",sapply(PD$pdx$meta$CNA,function(x) x$LABEL))
  if(length(idxx)==0){
    return(NULL)
  }
  hethomres=lapply(PD$pdx$meta$CNA[idxx],function(x){
    label=sprintf("%s (%s:%s-%s)",x$LABEL,x$chr,x$start,x$end)
    ##representative sample
    colony=x$samples[1]
    ascatfile=get_cnfile(PD,colony,type="ascat")
    if(exists("CENSOR_GERMLINE")){
      if(b_do_plot){
        plot_cn(get_cn(cn_summary_file = ascatfile),label = sprintf("%s: Representative of %d colonies : %s",colony,length(x$samples),label))
        warning("Germline info is censored so not plotting BAFs..")
        b_do_plot=FALSE
      }
    }else{
      baff=get_cnfile(PD,colony,type="baf")
      baf=read_tsv(baff);
      baf=baf[which(!is.na(baf$BAF) & baf$Chromosome==x$chr),]
      bafhet=bafrefhet[which(bafrefhet$Chromosome==x$chr),]
      baf=baf[baf$Position %in% bafhet$Position,]
      mbaf=frollmean(2*abs(baf$BAF-0.5),20,align = "right")
      mbaf2=frollmean(2*abs(baf$BAF-0.5),20,align = "left")
      idx=which(mbaf>0.8)
    }
    hethom=get_duplication_dat(PD,x$node,x)
    if(is.null(hethom)){
      return(NULL)
    }
    hethom$label=x$LABEL
    hethom$node=x$node
    hethom$het.sensitivity=PD$pdx$tree_ml$sensitivity.snv[match(x$node,PD$pdx$tree_ml$edge[,2])]
    hethom$chr=x$chr
    hethom$start=x$start
    hethom$end=x$end
    if(b_do_plot){
      kable(hethom$inf, "html") %>% kable_styling("striped")
      plot_cn(get_cn(cn_summary_file = ascatfile),label = sprintf("%s: Representative of %d colonies : %s",colony,length(x$samples),label))
      xx=cbind(baf$Position/1e6,baf$BAF)
      smoothScatter(xx,nrpoints = 0,xlab="Position(Mb)",ylab="BAF")
      abline(v=hethom$inf$Pos[hethom$inf$Chrom==x$chr & hethom$inf$is_dup]/1e6,col="green")
      abline(v=hethom$inf$Pos[hethom$inf$Chrom==x$chr & !hethom$inf$is_dup]/1e6,col="blue")
      abline(v=c(x$start/1e6,x$end/1e6),col="red")
      legend("bottomright",c("VAF=2/3","VAF=1/3","CNA Boundary"),col=c("green","blue","red"),lwd=2,bg = "white",cex = 0.8)
    }
    
    idx=which(bins$chr==x$chr & bins$start>=(1e5*floor(x$start/1e5)) & bins$end<=(1e5*floor(x$end/1e5)))
    hethom$kb=length(idx)*1e5
    hethom$count_in_bin=sum(bins$count[idx])
    hethom$count_se=sqrt(sum(bins$count[idx]))
    hethom$pmut=sum(bins$prob[idx])
    ##Estimate of uncertainty in local density
    hethom$pmut_se=hethom$pmut/sqrt(hethom$count_in_bin)
    ##Use a global estimate of 20..
    ##Use the null model estimates
    ut=PD$fit$poisson_tree$nullmodel$ultratree
    T=ut$edge.length[which(ut$edge[,2]==x$node)]
    ##Use null model lambda
    lambda=PD$fit$poisson_tree$nullmodel$lambda$mean
    #lambda=20
    dat=list(N=hethom$dupcount,M=hethom$ndupcount,pdup=hethom$het.sensitivity,pndup=hethom$het.sensitivity,L=lambda*T*hethom$pmut,Lstd=lambda*T*hethom$pmut_se)
    ##Fit stan model...
    stanr=stan("stanmodels/cntimedup.stan", data=dat,chains = 3,iter=20000,cores=3)
    sumc=rstan::summary(stanr,probs=c(0.025,0.5,0.975))$summary
    for( z in colnames(sumc)){
      hethom[[sprintf("x%s",z)]]=sumc["x",z]
    }
    for( z in colnames(sumc)){
      hethom[[sprintf("l%s",z)]]=sumc["l",z]
    }
    hethom
  })
  othernames=c("inf")
  keepnames=setdiff(names(hethomres[[1]]),othernames)
  list(hethom=do.call("rbind",lapply(hethomres,function(x) as.data.frame(x[keepnames]))),
       other=do.call("rbind",lapply(hethomres,function(x) x[othernames])))
}


get_cnfile=function(PD,colony,type="ascat"){
  #browser()
  if(!type %in% c("ascat","baf")){
    stop("get_cnfile:unsupported type")
  }
  idx=which(PD$pdx$cfg$SHORT_LABEL==colony)
  sample=PD$pdx$cfg$LABEL[idx]
  project=PD$pdx$cfg$PROJECT[idx]
  if(type=="ascat"){
    suffix="ascat_ngs.summary.csv"
  }else{
    ##baf
    suffix="ascat_ngs.cn.tsv.gz"
  }
  target=sprintf("%s/%s_%s_%s",CACHE,project,sample,suffix)
  if(!file.exists(target)){
    src=sprintf("%s/cancer_ref01/nst_links/live/%s/%s/%s.%s",NFS,project,sample,sample,suffix)
    status=file.copy(src,target)
    if(!status){
      cat(sprintf("local copy: %s -> %s",src,target))
      stop("unable to copy file..")
    }
  }
  target
}
## Get the counts of heterozygous and homozygous variants for the specified branch where a simple copy number neutral LOH has occurred.
get_loh_tree=function(PD,node,loh,exclude=c(),hom.threshold=0.8){
  inf=PD$pdx$dat
  mtr=inf$mtr
  dep=inf$dep
  cs=PD$pdx$meta$clones_short
  kids=c(node,get_all_node_children(node,PD$pdx$tree_ml))
  tips=kids[kids<=length(PD$pdx$tree_ml$tip.label)]
  mt=PD$pdx$tree_ml$tip.label[tips]
  wt=setdiff(cs,mt)
  ##check
  wt=setdiff(wt,exclude)
  mt=setdiff(mt,exclude)
  if(length(union(mt,intersect(loh$samples,cs)))!=length(mt)){
    cat("mt:",mt,"\n")
    cat("loh$samples",setdiff(loh$samples,cs),"\n")
    browser()
    stop("Inconsistency between LOH samples and associated node tips")
  }
  idx=with(inf$details,which(Chrom==loh$chr & Pos>loh$start & Pos<loh$end & TYPE=="SNV" & BGLOD==0))
  idx2=which(inf$details$node[idx]==node)
  cat("Number of muts overlapping loh shared in tree =",length(idx2),"\n")
  tmp=inf$details[idx[idx2],c("Chrom","Pos")]
  
  if(length(mt)==1){
    MT=inf$mtr[idx[idx2],mt]
    DP=inf$dep[idx[idx2],mt]
  }else{
    MT=rowSums(inf$mtr[idx[idx2],mt])
    DP=rowSums(inf$dep[idx[idx2],mt])
  }
  if(length(wt)==1){
    MTW=inf$mtr[idx[idx2],wt]
    DPW=inf$dep[idx[idx2],wt]
  }else{
    MTW=rowSums(inf$mtr[idx[idx2],wt])
    DPW=rowSums(inf$dep[idx[idx2],wt])
  }
  tmp$MT=MT
  tmp$DP=DP
  tmp$VAF=MT/DP
  tmp$MTW=MTW
  tmp$DPW=DPW
  tmp$VAFW=MTW/DPW
  status=classify_2state(MT,DP,probs = c(0.5,1))
  nhet=sum(status==1,na.rm = TRUE)
  nhom=sum(status==2,na.rm = TRUE)
  list(nhet=nhet,nhom=nhom,inf=tmp)
}

classify_2state=function(mtr,depth,probs,epsilon=0.01,threshold=0.95){
  if(length(probs)!=2){
    stop("Need to supply to probs")
  }
  #Maths for bernoulli trial.
  #p(mut)=p(mut|mutant read)p(mutant read)+p(mut|ref read)p(ref read)
  #p(mut)=(1-epsilon)*mu+(epsilon/3)(1-mu)
  probs=probs*(1-epsilon)+(epsilon/3)*(1-probs)
  l1=dbinom(mtr,depth,prob=probs[1] )
  l2=dbinom(mtr,depth,prob=probs[2])
  l1n=l1/(l1+l2)
  ifelse(l1n>threshold,1,ifelse(1-l1n>threshold,2,NA))
}

plot_cn=function(dat,label="",xlab="",yylim=c(-1,4)){
  if(is.null(dat)){
    warning("Empty input to plot_cn")
    return(NULL)
  }
  inf=get_chr_info2()
  p.start=inf$start_pos[match(dat$chr,inf$chr)]+dat$start
  p.end=inf$start_pos[match(dat$chr,inf$chr)]+dat$end
  ###browser()
  plot(NULL,
       pch=19,cex=0.5,xaxt="n",
       ylim=yylim,ylab="",
       yaxt="n",xlab=xlab,xlim=c(0,p.end[length(p.end)]),main=label)
  abline(h=0:10,col="grey")
  axis(side = 2,at = 0:4,las=1)
  rect(xleft = p.start,xright=p.end,
       ybottom = dat$minor-0.1,ytop=dat$minor-0.01,col="green",border="green")
  
  rect(xleft = p.start,xright=p.end,
       ybottom = dat$major+0.01,ytop=dat$major+0.1,col="red",border = "red")
  abline(v=inf$start_pos)
  mtext(at=inf$start_pos+0.5*inf$end,side = 1,cex=0.8,line = 1,text = inf$chr)
  
}

get_cn=function(cn_summary_file){
  if(!file.exists(cn_summary_file)){
    warning(sprintf("%s: does not exist",cn_summary_file))
    return(NULL)
  }
  cn=read.csv(cn_summary_file,header = FALSE)
  cn$start=cn$V3
  cn$end=cn$V4
  cn$chr=cn$V2
  ###V7=total copy number
  ## V8=minor allele copy number
  cn$major=cn$V7-cn$V8
  cn$minor=cn$V8
  cn[,-grep("^V",colnames(cn))]
}

## Get the counts of somatic variants with VAF=1/3 and 2/3 for the specified branch where a simple duplication has occureed.
get_duplication_dat=function(PD,node,loh,exclude=c(),hom.threshold=0.8){
  inf=PD$pdx$dat
  mtr=inf$mtr
  dep=inf$dep
  cs=PD$pdx$meta$clones_short
  kids=c(loh$node,get_all_node_children(loh$node,PD$pdx$tree_ml))
  tips=kids[kids<=length(PD$pdx$tree_ml$tip.label)]
  mt=PD$pdx$tree_ml$tip.label[tips]
  wt=setdiff(cs,mt)
  ##check
  wt=setdiff(wt,exclude)
  mt=setdiff(mt,exclude)
  if(length(union(mt,intersect(loh$samples,cs)))!=length(mt)){
    cat("mt:",mt,"\n")
    cat("loh$samples",setdiff(loh$samples,cs),"\n")
    browser()
    stop("Inconsistency between LOH samples and associated node tips")
  }
  idx=with(inf$details,which(Chrom==loh$chr & Pos>loh$start & Pos<loh$end & TYPE=="SNV" & BGLOD==0))
  idx2=which(inf$details$node[idx]==node)
  cat("Number of muts overlapping loh shared in tree =",length(idx2),"\n")
  #cat("WT counts:\n")
  #print(inf$mtr[idx[idx2],wt])
  tmp=inf$details[idx[idx2],c("Chrom","Pos")]
  #print(inf$details[idx[idx2],1:10])
  ##browser()
  cat("WT:",wt,"\n")
  cat("MT:",mt,"\n")
  
  if(length(mt)==1){
    MT=inf$mtr[idx[idx2],mt]
    DP=inf$dep[idx[idx2],mt]
    OFS=inf$ofs[idx[idx2],mt]==0
  }else{
    if(length(idx[idx2])<2){
      return(NULL)
    }
    MT=rowSums(inf$mtr[idx[idx2],mt])
    DP=rowSums(inf$dep[idx[idx2],mt])
    OFS=rowSums(inf$ofs[idx[idx2],mt]==0)>0
  }
  if(length(wt)==1){
    MTW=inf$mtr[idx[idx2],wt]
    DPW=inf$dep[idx[idx2],wt]
  }else{
    MTW=rowSums(inf$mtr[idx[idx2],wt])
    DPW=rowSums(inf$dep[idx[idx2],wt])
  }
  tmp$MT=MT
  tmp$DP=DP
  tmp$VAF=MT/DP
  tmp$MTW=MTW
  tmp$DPW=DPW
  tmp$VAFW=MTW/DPW
  tmp$OFS=OFS
  status=classify_2state(tmp$MT,tmp$DP,probs = c(1/3,2/3))
  tmp$is_dup=ifelse(status==2,1,ifelse(status==1,0,NA))
  list(dupcount=sum(status==2,na.rm=TRUE),ndupcount=sum(status==1,na.rm=TRUE),inf=tmp)
}




##Get chromosome lengths. TODO replace hardcoding with genome.fa.fai.
get_chr_info2=function(){
  chr_info=data.frame(chr=c(sprintf("%d",1:22),c("X","Y")) ,end=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
                                                                  146364022,141213431,135534747,134996516,133851895,115169878,107349540,
                                                                  102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,
                                                                  155270560,59373566))
  chr_info$start_pos=c(0,cumsum(chr_info$end)[-length(chr_info$end)])
  chr_info
}