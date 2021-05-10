#' starfish_link
#'
#' This function loads a SV dataframe and generates "seed" CGR regions, "connected" CGR regions, and complex SVs.
#'
#'
#' @param sv_file a SV dataframe with 8 columns: "chrom1","end1", "chrom2","end2","svtype" (DEL,DUP,h2hINV,t2tINV,TRA),"strand1" (+/-) and "strand2" (+/-),"sample". Other svtypes like INV, INS, BND are not accepted
#' @param prefix the prefix for all intermediate files, default is none
#' @return a list of files: $interleave_tra_complex_sv contains complex SVs, $mergecall contains "seed" CGR regions, and $starfish_call contains "connected" CGR regions.
#' @export
starfish_link=function(sv_file,prefix=""){

  chrlist=c(as.character(1:22),"X")

  svtotal=unique(sv_file)
  ########## no "chr" in chromosome name #######
  svtotal$chrom1=gsub("Chr|chr","",svtotal$chrom1)
  svtotal$chrom2=gsub("Chr|chr","",svtotal$chrom2)
  ###### pesudo cnv file #######
  cnvtotal=as.data.frame(matrix(nrow=1,ncol=4))
  colnames(cnvtotal)=c("chromosome","start","end","total_cn")
  cnvtotal[1,]=c("X","0","1","0")
  ##############################
  svtotal <-svtotal[(svtotal$chrom1 %in% chrlist & svtotal$chrom2 %in% chrlist),]
  samplelist=unique(svtotal$sample)

  ss6=data.frame()
  ss3=data.frame()
  sv=data.frame()
  interleave_sv_total=data.frame()

  for (i in c(1:length(samplelist))) {


    ############ 6 interleaved SVs ##############

    SV =svtotal[svtotal$sample==samplelist[i],]
    cnv=cnvtotal
    SVSample <- ShatterSeek::SVs(chrom1=as.character(SV$chrom1), pos1=as.numeric(SV$end1),chrom2=as.character(SV$chrom2), pos2=as.numeric(SV$end2),SVtype=as.character(SV$svtype),strand1=as.character(SV$strand1),strand2=as.character(SV$strand2))

    CN_data=ShatterSeek::CNVsegs(chrom=as.character(cnv$chromosome),start=as.numeric(cnv$start),end=as.numeric(cnv$end),total_cn=as.numeric(cnv$total_cn))
    chromothripsis6 <- ShatterSeek::shatterseek(SV.sample=SVSample,seg.sample=CN_data,min.Size=6)
    chrss6=chromothripsis6@chromSummary
    chrss6$sample=samplelist[i]

    ss6=rbind(ss6,chrss6)

    ########## 3 interleaved SVs #################

    chromothripsis3 <- ShatterSeek::shatterseek(SV.sample=SVSample,seg.sample=CN_data,min.Size=3)
    chrss3=chromothripsis3@chromSummary
    chrss3$sample=samplelist[i]
    ss3=rbind(ss3,chrss3)

    ########## interleaved sv #############
    interleave_size=3
    idx=which(sapply(chromothripsis3@detail$connComp,length)>=interleave_size)
    if(length(idx)>0){
      interleave_sv=chromothripsis3@detail$SV[as.numeric(unlist(chromothripsis3@detail$connComp[idx])) ,]
      chrss_interleaved_intrasv=merge(SV,interleave_sv,by.x=c("chrom1","end1","chrom2","end2","strand1","strand2","svtype"),by.y=c("chrom1","pos1","chrom2","pos2","strand1","strand2","SVtype"))

      interleave_sv_total=rbind(interleave_sv_total,chrss_interleaved_intrasv)

    }
  }

  ###### test for 6 interleaved SVs ########
  shattercall6=ss6
  ft=shattercall6$pval_fragment_joins
  ent=shattercall6$chr_breakpoint_enrichment
  ext=shattercall6$pval_exp_chr
  exct=shattercall6$pval_exp_cluster

  shattercall6$ft_fdr=p.adjust(ft, method = "fdr")
  shattercall6$ent_fdr=p.adjust(ent, method = "fdr")
  shattercall6$ext_fdr=p.adjust(ext, method = "fdr")
  shattercall6$exct_fdr=p.adjust(exct, method = "fdr")

  shattercall6$total_intra=shattercall6$number_DEL+shattercall6$number_DUP+shattercall6$number_h2hINV+shattercall6$number_t2tINV
  if (length(samplelist)>50){
  shattercall6$call6=ifelse(((shattercall6$total_intra>5)&(shattercall6$ft_fdr>0.2)&(shattercall6$ent_fdr<=0.2|shattercall6$exct_fdr<=0.2)),"chrss","no")
  } else {
    shattercall6$call6=ifelse(((shattercall6$total_intra>5)&(shattercall6$ft_fdr>0.05)&(shattercall6$ent_fdr<=0.2|shattercall6$exct_fdr<=0.2)),"chrss","no")
  }
  ###### test for 3 interleaved SVs ########
  shattercall3=ss3
  ft=shattercall3$pval_fragment_joins
  ent=shattercall3$chr_breakpoint_enrichment
  ext=shattercall3$pval_exp_chr
  exct=shattercall3$pval_exp_cluster

  shattercall3$ft_fdr=p.adjust(ft, method = "fdr")
  shattercall3$ent_fdr=p.adjust(ent, method = "fdr")
  shattercall3$ext_fdr=p.adjust(ext, method = "fdr")
  shattercall3$exct_fdr=p.adjust(exct, method = "fdr")

  shattercall3$total_intra=shattercall3$number_DEL+shattercall3$number_DUP+shattercall3$number_h2hINV+shattercall3$number_t2tINV

  if (length(samplelist)>50){
    shattercall3$call3=ifelse(((shattercall3$total_intra>2)&(shattercall3$number_TRA>3)&(shattercall3$ft_fdr>0.2)),"chrss","no")
  } else{
    shattercall3$call3=ifelse(((shattercall3$total_intra>2)&(shattercall3$number_TRA>3)&(shattercall3$ft_fdr>0.05)),"chrss","no")
  }
  ######## merge 3 and 6 interleaved SVs #####

  shattercall6$call3= shattercall3$call3
  shattercall3$call6=shattercall6$call6
  keeps=c("chrom","start","end","sample","call3","call6")
  shattercall6=shattercall6[keeps]
  shattercall3=shattercall3[keeps]

  mergecall=rbind(shattercall6[shattercall6$call6=="chrss"&shattercall6$call3=="no",],shattercall3[shattercall3$call3=="chrss",])
  if(nrow(mergecall)==0){
    print("No CGR region is identified!")
  } else {

  colnames(mergecall)=c("chr","start","end","sample","call3","call6")

  ######## call linking groups #######


  chrss=mergecall
  chrss$format_id=paste0(chrss$sample,"_",chrss$chr)
  chrss$size=chrss$end-chrss$start
  chrss$chrss_status="chrss"
  chrss$chrss_chr=chrss$chr
  svlist=svtotal
  svlist$start1=svlist$end1-1
  svlist$start2=svlist$end2-1
  svlist$size=svlist$end2-svlist$end1
  format_id=chrss$format_id

  ss=data.frame()


  chrsslink=function(chrss,formatid=format_id){


    for (i in 1:length(chrss$chr)){

      ##### identify TRA in 10k ##########
      inter1=max(chrss$start[i]-10000,0)
      inter2=chrss$end[i]+10000
      sv=svlist[svlist$sample==chrss$sample[i],]
      svtra=sv[((sv$chrom1==chrss$chr[i]&sv$end1>=inter1&sv$end1<=inter2)|(sv$chrom2==chrss$chr[i]&sv$end2>=inter1&sv$end2<=inter2))&(sv$svtype=="TRA"),]
      if (length(svtra$chrom1)>0){
        svtra$chr=ifelse(svtra$chrom1==chrss$chr[i],svtra$chrom2,svtra$chrom1)
        svtra$tra_pos1=ifelse(svtra$chrom1==chrss$chr[i],svtra$end2,svtra$end1)
        svtra$chrss_chr=ifelse(svtra$chrom1==chrss$chr[i],svtra$chrom1,svtra$chrom2)
        svtra$sample=chrss$sample[i]
        svtra <- plyr::ddply(svtra, .(chr), transform, n = length(chr))
        svtra=svtra[svtra$n>0,]
        if (length(svtra$chrom1)>0){
          svtramin=aggregate(svtra$tra_pos1, by = list(svtra$chr), min)
          colnames(svtramin)=c("chr","start")
          svtramax=aggregate(svtra$tra_pos1, by = list(svtra$chr), max)
          colnames(svtramax)=c("chr","end")
          svtra2=merge(svtramin,svtramax,by.x=("chr"),by.y=("chr"))
          svtra=svtra[!duplicated(svtra[,c('chr')]),]
          svtra3=merge(x = svtra2, y = svtra[ , c("chr","chrss_chr", "sample","n")], by = "chr", all.x=TRUE)
          ss=rbind(ss,svtra3)
        }
      }
    }

    if(nrow(ss)>0){

      ss$format_id=paste0(ss$sample,"_",ss$chr)
      idlist=unique(ss$format_id)
      # requires at least two chrss linked tra
      linkss=ss[1,]
      linkss2=data.frame()
      for (id in idlist){
        ss_id=ss[ss$format_id==id,]
        if(sum(ss_id$n)>1){
          # judge the chromosome is an original chrss or not
          if(!(id %in% formatid)){
            # caculate the region for the link region
            linkss$chr=unique(ss_id$chr)
            linkss$start=min(ss_id$start,ss_id$end)
            linkss$end=max(ss_id$start,ss_id$end)
            linkss$chrss_chr=paste(c(chrss[chrss$format_id==id,]$chrss_chr,paste(c(ss_id$chrss_chr),collapse="_",sep="_")),collapse="_",sep="_")
            linkss$sample=unique(ss_id$sample)

            linkss$format_id=unique(ss_id$format_id)
            linkss$chrss_status="no"
            linkss2=rbind(linkss2,linkss)

          }
          if((id %in% formatid)){
            # caculate the region for the link region
            linkss$chr=unique(ss_id$chr)
            # assign start and end based on TRA position
            linkss$start=min(ss_id$start,ss_id$end)
            linkss$end=max(ss_id$start,ss_id$end)
            linkss$chrss_chr=paste(c(chrss[chrss$format_id==id,]$chrss_chr,paste(c(ss_id$chrss_chr),collapse="_",sep="_")),collapse="_",sep="_")

            linkss$sample=unique(ss_id$sample)

            linkss$format_id=unique(ss_id$format_id)
            linkss$chrss_status="chrss"
            linkss2=rbind(linkss2,linkss)

          }
        }
      }
      # assign start and end based on chrss cluster

      linkss2=merge(linkss2,chrss[,c("format_id","start","end")],by.x=("format_id"),by.y=("format_id"),all.x=T)
      linkss2$start=ifelse(!is.na(linkss2$start.y),ifelse(linkss2$start.x<linkss2$start.y,linkss2$start.x,linkss2$start.y),linkss2$start.x)
      linkss2$end=ifelse(!is.na(linkss2$end.y),ifelse(linkss2$end.x>linkss2$end.y,linkss2$end.x,linkss2$end.y),linkss2$end.x)
      keeps=c("chr","start","end","chrss_chr","sample","format_id","chrss_status")
      linkss2=linkss2[keeps]
      # merge with original chrss
      chrss2=chrss[keeps]
      linkss2=rbind(linkss2,chrss2)
      idlist=unique(linkss2$format_id)

      linkss=linkss2[1,]
      linkss3=data.frame()
      for (id in idlist){
        ss_id=linkss2[linkss2$format_id==id,]
         # caculate the region for the link region
        linkss$chr=unique(ss_id$chr)
        linkss$start=min(ss_id$start,ss_id$end)
        linkss$end=max(ss_id$start,ss_id$end)
        linkss$chrss_chr=paste(c(ss_id$chr,ss_id$chrss_chr),collapse="_",sep="_")
        linkss$sample=unique(ss_id$sample)
        linkss$format_id=unique(ss_id$format_id)
        linkss$chrss_status=unique(linkss2[linkss2$format_id==id,]$chrss_status)

        linkss3=rbind(linkss3,linkss)

      }
    } else {
      linkss3=chrss[c("chr","start","end","chrss_chr","sample","format_id","chrss_status")]
    }

    return(linkss3)
  }



  chrsslink1=chrsslink(chrss)
  chrsslink2=chrsslink(chrsslink1)
  # iteratively search and update the chromothriptic "seed" list
  round=0
  while (!identical(chrsslink1,chrsslink2)){
    chrsslink1=chrsslink(chrsslink2)
    chrsslink2=chrsslink(chrsslink1)
    chrsslink1$chrss_chr <- sapply(strsplit(chrsslink1$chrss_chr, "_", fixed = TRUE), function(x)
      paste(unique(x), collapse = "_"))
    chrsslink2$chrss_chr <- sapply(strsplit(chrsslink2$chrss_chr, "_", fixed = TRUE), function(x)
      paste(unique(x), collapse = "_"))
    round=round+1
    print(paste0("Iteration ",round))
    print("Starfish is connecting chromothriptic regions...")
  }

  # split and unique value in chrss_chr

  chrsslink3=chrsslink1
  chrsslink1$chr=gsub("X","23",chrsslink1$chr)
  chrsslink1$chrss_chr=gsub("X","23",chrsslink1$chrss_chr)
  chrsslink1$chr=gsub("Y","24",chrsslink1$chr)
  chrsslink1$chrss_chr=gsub("Y","24",chrsslink1$chrss_chr)
  chrsslink1$link_chromosome=0
  chrsslink1$chrss_chr=strsplit(as.character(chrsslink1$chrss_chr), "_", fixed = TRUE)



  for (j in 1:length(chrsslink1$sample)){
    id=chrsslink1$sample[j]
    chr=chrsslink1$chr[j]
    chrssid=chrsslink1[chrsslink1$sample==id,]
    # identify chrss_chr contains chr with grepl
    link_chrss=as.character(sort(as.numeric(unique(unlist(c(chrssid[chrssid$chr==chr,]$chrss_chr))))))
    # find intersect rows with link_chrss
    link_tra=as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_chrss) length(intersect(y, x))>0))]))))))
    round=0
    while (!identical(link_chrss,link_tra)){
    link_chrss=link_tra=as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_tra) length(intersect(y, x))>0))]))))))
    link_tra=link_tra=as.character(sort(as.numeric(unique(unlist(c(chrssid$chrss_chr[unlist(lapply(chrssid$chrss_chr, function(x,y=link_chrss) length(intersect(y, x))>0))]))))))
    round=round+1
    }
    link_chr=paste(link_chrss,collapse="_")
    chrsslink1$link_chromosome[j]=link_chr
  }

  chrsslink1$chrss_chr <- sapply(chrsslink1$chrss_chr, function(x) paste(unique(x), collapse = "_"))
  chrsslink1$chr=gsub("23","X",chrsslink1$chr)
  chrsslink1$chr=gsub("24","Y",chrsslink1$chr)

  chrsslink1$chrss_chr=gsub("23","X",chrsslink1$chrss_chr)
  chrsslink1$link_chromosome=gsub("23","X",chrsslink1$link_chromosome)
  chrsslink1$chrss_chr=gsub("24","Y",chrsslink1$chrss_chr)
  chrsslink1$link_chromosome=gsub("24","Y",chrsslink1$link_chromosome)
  chrsslink1$cluster_id=paste(chrsslink1$sample,chrsslink1$link_chromosome,sep="_")
  chrsslink1$chrss_status=gsub("no","link",chrsslink1$chrss_status)

  ########## integrate TRA back to interleave SV cluster ##############

  samplelist=unique(chrsslink1$sample)
  ssmatrix=data.frame()

  for(i in c(1:length(samplelist))){


    svss=svtotal[svtotal$sample==samplelist[i],]
    svss_inter=svss[svss$svtype=="TRA",]
    svss_inter=svss_inter[c("chrom1","end1","chrom2","end2","svtype")]

    svss_intra=interleave_sv_total[interleave_sv_total$sample==samplelist[i],]
    svss_intra=svss_intra[c("chrom1","end1","chrom2","end2","svtype")]

    svss=rbind(svss_inter,svss_intra)
    if(nrow(svss)>0){
      svlist3=chrsslink1[chrsslink1$sample==samplelist[i],]
      svss$complex=0

      for (j in 1:length(svss$chrom1)){
        t1=0
        t2=0
        for (k in 1:length(svlist3$chr)){
          t1=ifelse(svss$chrom1[j]==svlist3$chr[k] & svss$end1[j]>=svlist3$start[k] & svss$end1[j]<=svlist3$end[k],1,0)+t1
        }
        for (k in 1:length(svlist3$chr)){
          t2=ifelse(svss$chrom2[j]==svlist3$chr[k] & svss$end2[j]>=svlist3$start[k] & svss$end2[j]<=svlist3$end[k],1,0)+t2
        }
        svss$complex[j]=t1+t2
      }


      svss2=svss[svss$complex>1,]
      if(nrow(svss2)>0){
        svss2=svss2[c("chrom1","end1","chrom2","end2","svtype","complex")]
        svss2$sample=samplelist[i]

        ssmatrix=rbind(ssmatrix,svss2)
      }
    }

  }


  smatrix2=unique(ssmatrix)

  starfish_list=list("ss6"=ss6,"ss3"=ss3,"interleave_sv"=interleave_sv_total,"interleave_tra_complex_sv"=smatrix2,"shattercall6"=shattercall6,"shattercall3"=shattercall3,"mergecall"=mergecall,"starfish_call"=chrsslink1)

  filename=ifelse(prefix=="","Connected_CGR_event.csv",paste0(prefix,"_connected_CGR_event.csv"))
  write.csv(chrsslink1,filename,row.names = F)

  sv_filename=ifelse(prefix=="","Connected_CGR_complex_SV.csv",paste0(prefix,"_connected_CGR_complex_SV.csv"))
  write.csv(smatrix2,sv_filename,row.names = F)


  return(starfish_list)

}
}




