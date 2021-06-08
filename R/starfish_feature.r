#' starfish_feature
#'
#' This function loads "connected" CGR regions and complex SVs reported by starfish_link, then combines CNV calls and sample gender to construct a CGR region vs. feature matrix for downstream clustering and classification.
#'
#' @param cgr "connected" CGR regions, which is the output of starfish_link_out$starfish_call
#' @param complex_sv complex SVs, which is the output from starfish_link_out$interleave_tra_complex_sv
#' @param cnv_file a CNV dataframe with 5 columns: "chromosome","start","end","total_cn", and "sample". "total_cn" should contain absolute copy numbers
#' @param gender_file a sample table with 2 columns: "sample" and "gender" ("Female, "female","F","f","Male","male","M","m"). If gender is unknown, any other characters could be put here like "unknown"
#' @param prefix the prefix for all intermediate files, default is none
#' @param genome_v which genome assembly was used to call SV and CNV. It should be "hg19" or "hg38", default is "hg19"
#' @param cnv_factor the CN fluctuation beyond or below baseline to identify loss and gain fragments for samples with decimal CN, default is 0.15
#' @param arm_del_rm the logical value of removing arm level deletion or not, default is TRUE
#' @return a list of files: $cluster_feature is the CGR region vs. feature matrix,$cnv_baseline is the cnv file with baseline annotation
#' @export

starfish_feature=function(cgr,complex_sv,cnv_file,gender_file,prefix="",genome_v="hg19",cnv_factor=0.15,arm_del_rm=TRUE){
  print("Starfish is computing the feature matrix, please be patient...")

  if(genome_v=="hg19"){
    chrlength=hg19_chr_length
  } else if(genome_v=="hg38"){
    chrlength=hg38_chr_length
  } else {
    print("Wrong genome version!")
  }


  chrlist=c(as.character(1:22),"X")
  cgr$range_size=cgr$end-cgr$start+1
  cgr$chr=as.character(cgr$chr)
  chrlength$chrom=as.character(chrlength$chrom)
  cgr=merge(cgr,unique(chrlength[c("chrom","chr_size")]),by.x=("chr"),by.y=("chrom"))
  samplelist=unique(cgr$cluster_id)

  complex_sv$end1=complex_sv$pos1
  complex_sv$end2=complex_sv$pos2
  complex_sv$start1=complex_sv$end1-1
  complex_sv$start2=complex_sv$end2-1
  complex_sv=complex_sv[complex_sv$sample %in% cgr$sample,]
  complex_sv$size=complex_sv$end2-complex_sv$end1

  ###### format CNV file ######
  cnv_total=cnv_file[!is.na(cnv_file$total_cn),]
  cnv_total=cnv_total[c("chromosome","start","end","total_cn","sample")]
  cnv_total$chromosome=as.character(cnv_total$chromosome)
  cnv_total$chromosome=gsub("Chr|chr","",cnv_total$chromosome)
  cnv_total$total_cn=round(cnv_total$total_cn,digits=2)
  cnv_total$total_cn=round(cnv_total$total_cn,digits=1)
  cnv_total=cnv_total[cnv_total$chromosome %in% chrlist,]
  cnv_total=cnv_total[cnv_total$start!=cnv_total$end,]
  cnv_total$size=cnv_total$end-cnv_total$start+1
  cnv_total=cnv_total[with(cnv_total,order(sample, chromosome,start)),]
  ############## matrix of features  ###############
  cluster=cgr[c("cluster_id","sample")]
  gender_file$gender=ifelse(gender_file$gender%in% c("female","f","F"),"Female",ifelse(gender_file$gender %in% c("male","m","M"),"Male",gender_file$gender))


  if(nrow(cluster[!cluster$sample %in% gender_file$sample,])>0){
    cluster_no_gender=cluster[!cluster$sample %in% gender_file$sample,]
    print(paste0("There is no gender file for sample ",c(unique(cluster_no_gender$sample))))
  }

  cluster=merge(cluster,gender_file,by.x="sample",by.y="sample")
  cluster=cluster[!duplicated(cluster[,c("cluster_id")]),]

  cluster$brk_10mb <-  0
  cluster$brk_1mb <-  0
  cluster$sv_avg_size <- 0
  cluster$CN_state <- 0
  cluster$cnv_state_mean <- 0
  cluster$max_CNV <- 0
  cluster$del_size=0
  cluster$Loss_size_percentage=0
  cluster$Loss_chr_percentage=0
  cluster$dup_size=0
  cluster$Gain_size_percentage=0
  cluster$Gain_chr_percentage=0
  cluster$range_size_sum=0
  cluster$range_size_mean=0
  cluster$arm_deletion=0
  cluster$telo_deletion_size=0
  cluster$max_telo_loss_percentage=0
  cluster$Telo_loss_ratio_mean=0
  cluster$Loss_number_10Mb=0
  cluster$Gain_number_10Mb=0
  cluster$Chrss_size_ratio=0

  cluster$telo_vs_chrss_del_mean_ratio=0
  cluster$telo_vs_chrss_del_max_ratio=0
  cluster$chrss_chr_number=0
  cluster$Loss_size_mean_ratio=0
  cluster$Gain_size_mean_ratio=0
  cluster$Loss_number_length=0
  cluster$Gain_number_length=0

  cluster$Brk_dispersion_MAD <- 0
  cluster$Brk_dispersion_MAD_mean_total=0

  ####### only calculate feature for samples with cnv #########

  if(nrow(cluster[!cluster$sample %in% cnv_total$sample,])>0){
    cluster_no=cluster[!cluster$sample %in% cnv_total$sample,]
  print(paste0("There is no CNV file for sample ",c(unique(cluster_no$sample))))
  }

  cluster=cluster[cluster$sample %in% cnv_total$sample,]
  row.names(cluster)=1:nrow(cluster)
  cnvtotal=data.frame()


  for (i in 1:nrow(cluster)) {
    print(paste0(format(i/nrow(cluster)*100,digits=1),"% is done..."))
    print(paste0("Event ",cluster$cluster_id[i]))
    chrss_cluster=cgr[cgr$cluster_id==cluster$cluster_id[i],]

    donor_sex=cluster$gender[i]
    if(!donor_sex %in% c("Female","Male")){
      print (paste0("There is no gender information for sample ",cluster$sample[i],", uses Female as default"))
    }
    donor_sex=ifelse((!donor_sex %in% c("Female","Male")),"Female",donor_sex)

    chrn=length(unique(chrss_cluster$chr))
    cluster$chrss_chr_number[i]=chrn
    chrss_size_mean=mean(chrss_cluster$range_size)
    chrss_size_sum=sum(chrss_cluster$range_size)
    cluster$range_size_mean[i]=chrss_size_mean
    cluster$range_size_sum[i]=chrss_size_sum
    range_ratio=chrss_size_sum/sum(chrss_cluster$chr_size)
    cluster$Chrss_size_ratio[i]=range_ratio


    cnv=cnv_total[cnv_total$sample==cluster$sample[i],]
    cnv=cnv[cnv$chromosome %in% chrss_cluster$chr,]
    cnv=na.omit(cnv)
    cnv$ovl=0
    cnv$telo_loss=0
    if(nrow(cnv)>0){

      if(nrow(cnv)>1){
      ######### merge neighbor fragments with the same cnv ######
      for (t in 1:(nrow(cnv)-1)){
        if (cnv$chromosome[t]==cnv$chromosome[t+1]&cnv$total_cn[t]==cnv$total_cn[t+1]){
          cnv$start[t+1]=cnv$start[t]
          cnv$ovl[t]=1

        }

      }
      }
      cnv=cnv[cnv$ovl==0,]
      cnv$size=cnv$end-cnv$start+1
      ###################################
      ##### determine background cnv and deletion #####
      cnv_background=as.data.frame(cnv %>% dplyr::group_by(chromosome,total_cn) %>% dplyr::summarize(size=sum(size)))
      # cnv_background=cnv_background[cnv_background$total_cn>1.5,]
      cnv_background2=as.data.frame(cnv_background %>% dplyr::group_by(chromosome) %>% dplyr::top_n(n=1, wt = size))
      colnames(cnv_background2)=c("chromosome","background_cnv","background_size")

      cnv=merge(cnv,cnv_background2,by.x="chromosome",by.y="chromosome")



      # ####### use the whole CGR baseline CNV as the baseline #######
      # cnv_background=as.data.frame(cnv %>% dplyr::group_by(total_cn) %>% dplyr::summarize(size=sum(size)))
      # cnv_background2=cnv_background[cnv_background$size==max(cnv_background$size),]
      # cnv$background_cnv=cnv_background2$total_cn
      # cnv$background_size=cnv_background2$size
      ###############################


      ########### determine cnv is decimal or integer #######
      cnv_integer=sum(cnv$total_cn %%1)
      cnv_integer=ifelse(cnv_integer>0,"decimal","integer")

      if(cnv_integer=="integer"){
      if (donor_sex=="Male"){
        cnv$copy_loss=ifelse((cnv$total_cn<cnv$background_cnv|(cnv$chromosome!="X"&cnv$total_cn<2)|(cnv$chromosome=="X"&cnv$total_cn<1)),"deletion",ifelse((cnv$total_cn>cnv$background_cnv)&((cnv$chromosome=="X"&cnv$total_cn>1)|(cnv$chromosome!="X"&cnv$total_cn>2)) ,"gain","neutral"))

      } else {
        cnv$copy_loss=ifelse((cnv$total_cn<(cnv$background_cnv)|cnv$total_cn<2),"deletion",ifelse(cnv$total_cn>(cnv$background_cnv)&cnv$total_cn>2,"gain","neutral"))

      }
      } else if(cnv_integer=="decimal") {
        loss_factor=1-cnv_factor
        gain_factor=1+cnv_factor
        if (donor_sex=="Male"){

          cnv$copy_loss=ifelse((cnv$total_cn<(cnv$background_cnv*loss_factor)|(cnv$chromosome!="X"&cnv$total_cn<2*loss_factor)|(cnv$chromosome=="X"&cnv$total_cn<1*loss_factor)),"deletion",ifelse((cnv$total_cn>cnv$background_cnv*gain_factor)&((cnv$chromosome=="X"&cnv$total_cn>1*gain_factor)|(cnv$chromosome!="X"&cnv$total_cn>2*gain_factor)) ,"gain","neutral"))

        } else {

          cnv$copy_loss=ifelse((cnv$total_cn<(cnv$background_cnv*loss_factor)|cnv$total_cn<2*loss_factor),"deletion",ifelse(cnv$total_cn>(cnv$background_cnv*gain_factor)&cnv$total_cn>2*gain_factor,"gain","neutral"))

        }
      }


      ######## merge same copy gain or loss fragment #########
      cnv_ori=cnv
      if (nrow(cnv)>1){
      for (t in 1:(nrow(cnv)-1)){
        if (cnv$chromosome[t]==cnv$chromosome[t+1]&cnv$copy_loss[t]==cnv$copy_loss[t+1]){
          cnv$start[t+1]=cnv$start[t]
          cnv$ovl[t]=1

        }

      }
      }
      cnv=cnv[cnv$ovl==0,]

      #################################

      cnv=merge(cnv,unique(chrlength[c("chrom","chr_size")]),by.x=("chromosome"),by.y=("chrom"))
      cnvtotal=rbind(cnvtotal,cnv)

      ###### get overlapped cnv #####
      # ###### trim 10bp from the both end of chrss range to allow accurate cnv overlap #####
      # chrss_cluster$start=chrss_cluster$start-10
      # chrss_cluster$end=chrss_cluster$end+10

      #############

      chrss_range=makeGRangesFromDataFrame(chrss_cluster[c("chr","start","end","range_size")],keep.extra.columns = T)
      cnv_range=makeGRangesFromDataFrame(cnv,keep.extra.columns = T)
      cnv_ori_range=makeGRangesFromDataFrame(cnv_ori,keep.extra.columns = T)
      chr_range=makeGRangesFromDataFrame(chrlength,keep.extra.columns = T)
      ######## determine arm-level deletion if a fragment is larger than 95% of the arm #######
      chr_ovl=as.data.frame(findOverlapPairs(cnv_range, chr_range,ignore.strand=TRUE))
      chr_ovl$cnv_arm=chr_ovl$first.X.width/chr_ovl$second.arm_size
      chr_ovl$del_arm=ifelse(chr_ovl$first.copy_loss=="deletion"&chr_ovl$cnv_arm>=0.95,"arm_deletion","no_arm_deletion")
      arm_del_chr=as.character(chr_ovl[chr_ovl$del_arm=="arm_deletion",]$first.X.seqnames)
      chr_name=unique(cnv$chromosome)
      arm_deletion=vector()

      for (c in 1:length(chr_name)){

        chr_ovl_chr=chr_ovl[chr_ovl$first.X.seqnames==chr_name[c],]
        #### determine by deletion proportion ####
        arm_deletion_chr=ifelse(chr_ovl_chr$del_arm[1]=="arm_deletion"|chr_ovl_chr$del_arm[nrow(chr_ovl_chr)]=="arm_deletion",1,0)

        ##### determine by top deleted #####

        arm_deletion=c(arm_deletion,arm_deletion_chr)

      }

      ###### don't consider arm deletion if arm_del is not TRUE ######
      if (arm_del_rm!=TRUE) {
      arm_del_chr=0 }

      ######## define telomere loss size ###############
      telo_deletion_size=vector()
      telo_deletion_ratio=vector()
      telo_loss=data.frame()
      for (c in 1:length(chr_name)){
        cnv_ovl_chr=cnv[cnv$chromosome==chr_name[c],]

        telo_deletion_chr_size1=ifelse(((!cnv_ovl_chr$chromosome[1] %in% arm_del_chr) & cnv_ovl_chr$copy_loss[1]=="deletion"),cnv_ovl_chr$end[1]+1,0)

        telo_deletion_chr_size2=ifelse(((!cnv_ovl_chr$chromosome[nrow(cnv_ovl_chr)] %in% arm_del_chr) &cnv_ovl_chr$copy_loss[nrow(cnv_ovl_chr)]=="deletion"),cnv_ovl_chr$chr_size[nrow(cnv_ovl_chr)]-cnv_ovl_chr$start[nrow(cnv_ovl_chr)]+1,0)

        cnv_ovl_chr$telo_loss[1]=ifelse(((!cnv_ovl_chr$chromosome[1] %in% arm_del_chr) & cnv_ovl_chr$copy_loss[1]=="deletion"),1,0)
        cnv_ovl_chr$telo_loss[nrow(cnv_ovl_chr)]=ifelse(((!cnv_ovl_chr$chromosome[nrow(cnv_ovl_chr)] %in% arm_del_chr) &cnv_ovl_chr$copy_loss[nrow(cnv_ovl_chr)]=="deletion"),1,0)
        telo_loss=rbind(telo_loss,cnv_ovl_chr[cnv_ovl_chr$telo_loss==1,])

        telo_deletion_chr_size=ifelse(nrow(cnv_ovl_chr)==1,telo_deletion_chr_size1,telo_deletion_chr_size1+telo_deletion_chr_size2)

        telo_deletion_chr_ratio=telo_deletion_chr_size/(unique(cnv_ovl_chr$chr_size))

        telo_deletion_size=c(telo_deletion_size,telo_deletion_chr_size)
        telo_deletion_ratio=c(telo_deletion_ratio,telo_deletion_chr_ratio)

      }


      arm_deletion=max(arm_deletion)
      telo_deletion_size_sum=sum(telo_deletion_size)
      telo_deletion_size_mean=0
      telo_deletion_size_max=0
      if (telo_deletion_size_sum>0){
        telo_deletion_size_t=telo_deletion_size[telo_deletion_size>0]
        telo_deletion_size_mean=mean(telo_deletion_size_t)
        telo_deletion_size_mean=ifelse(is.na(telo_deletion_size_mean),0,telo_deletion_size_mean)
        telo_deletion_size_max=max(telo_deletion_size_t)
        telo_deletion_size_max=ifelse(is.na(telo_deletion_size_max),0,telo_deletion_size_max)
      }


      telo_deletion_ratio_max=max(telo_deletion_ratio)
      telo_deletion_ratio_mean=mean(telo_deletion_ratio)

      cluster$arm_deletion[i]=arm_deletion
      cluster$telo_deletion_size[i]=telo_deletion_size_sum
      cluster$max_telo_loss_percentage[i]=telo_deletion_ratio_max
      cluster$Telo_loss_ratio_mean[i]=telo_deletion_ratio_mean

      ######## determine ovelap cnv in the chrss region ######

      olaps <- findOverlaps(cnv_ori_range, chrss_range)

      cnv_ovl<- as.data.frame(pintersect(cnv_ori_range[queryHits(olaps)], chrss_range[subjectHits(olaps)]))



      if(nrow(cnv_ovl)>0){
        cnv_ovl$size=cnv_ovl$width
        del_size_mean=mean(cnv_ovl[cnv_ovl$copy_loss=="deletion",]$size)

        del_size_mean=ifelse(is.na(del_size_mean),0,del_size_mean)

        del_mean_por=del_size_mean/sum(chrss_cluster$range_size)

        del_por=sum(cnv_ovl[cnv_ovl$copy_loss=="deletion",]$size)/sum(chrss_cluster$range_size)

        del_chr_por=sum(cnv_ovl[cnv_ovl$copy_loss=="deletion",]$size)/sum(chrss_cluster$chr_size)


        del_number_10Mb=nrow(cnv_ovl[cnv_ovl$copy_loss=="deletion",])/(sum(chrss_cluster$range_size)/10000000)

        del_number_length=nrow(cnv_ovl[cnv_ovl$copy_loss=="deletion",])/(sum(chrss_cluster$range_size))


        dup_size_mean=mean(cnv_ovl[cnv_ovl$copy_loss=="gain",]$size)

        dup_size_mean=ifelse(is.na(dup_size_mean),0,dup_size_mean)

        dup_mean_por=dup_size_mean/sum(chrss_cluster$range_size)

        dup_por=sum(cnv_ovl[cnv_ovl$copy_loss=="gain",]$size)/sum(chrss_cluster$range_size)


        dup_chr_por=sum(cnv_ovl[cnv_ovl$copy_loss=="gain",]$size)/sum(chrss_cluster$chr_size)

        dup_number_10Mb=nrow(cnv_ovl[cnv_ovl$copy_loss=="gain",])/(sum(chrss_cluster$range_size)/10000000)

        dup_number_length=nrow(cnv_ovl[cnv_ovl$copy_loss=="gain",])/(sum(chrss_cluster$range_size))

        ###### cnv state per chr ####
        chr_uni=unique(cnv_ovl$seqnames)
        cnv_state=vector()
        for (k in 1:length(chr_uni)){


          cnv_uni=cnv_ovl[cnv_ovl$seqnames==chr_uni[k],]
          cnv_state_uni=length(unique(cnv_uni$total_cn))
          cnv_state=c(cnv_state,cnv_state_uni)
        }
        cnv_state_max=max(cnv_state)
        cnv_state_mean=mean(cnv_state)
        cnv_max=max(as.numeric(cnv_ovl$total_cn))
        cluster$CN_state[i] <- cnv_state_max
        cluster$cnv_state_mean[i] <- cnv_state_mean
        cluster$max_CNV[i] <- cnv_max
        cluster$del_size[i] <- del_size_mean
        cluster$Loss_size_percentage[i]=del_por
        cluster$Loss_chr_percentage[i]=del_chr_por

        cluster$dup_size[i] <- dup_size_mean
        cluster$Gain_size_percentage[i]=dup_por
        cluster$Gain_chr_percentage[i]=dup_chr_por

        cluster$Loss_number_10Mb[i]=del_number_10Mb
        cluster$Gain_number_10Mb[i]=dup_number_10Mb

        cluster$Loss_size_mean_ratio[i]=del_mean_por
        cluster$Gain_size_mean_ratio[i]=dup_mean_por
        cluster$Loss_number_length[i]=del_number_length
        cluster$Gain_number_length[i]=dup_number_length

      }

      cluster$telo_vs_chrss_del_mean_ratio[i]=(telo_deletion_size_mean)/(cluster$del_size[i]+1)
      cluster$telo_vs_chrss_del_max_ratio[i]=(telo_deletion_size_max)/(cluster$del_size[i]+1)
      ######################  sv ##############
      chrss_chr=chrss_cluster$chr
      sample=unique(chrss_cluster$sample)
      sv=complex_sv[complex_sv$sample==sample&(complex_sv$chrom1 %in% chrss_chr|complex_sv$chrom2 %in% chrss_chr),]

      if (nrow(sv)>1){
        brk_n=nrow(sv[sv$chrom1 %in% chrss_chr,])+nrow(sv[sv$chrom2 %in% chrss_chr,])
        brk10=brk_n/(sum(chrss_cluster$range_size)/10000000)
        brk1=brk_n/(sum(chrss_cluster$range_size)/1000000)
        ### size ###
        sv_avg=mean(sv[sv$svtype!="TRA",]$size)
        sv_avg=ifelse(is.na(sv_avg),0,sv_avg)

        sv$min_range=""

        ### get all break ####
        all_dispersion=vector()
        all_brk_dispersion=vector()
        all_sd=vector()
        all_brk_sd=vector()
        all_interval_total=vector()
        all_rd=vector()
        for (p in 1:length(chrss_cluster$chr)){
          sv_chr=chrss_cluster$chr[p]
          all_brk=sort(c(sv[sv$chrom1==sv_chr,]$start1,sv[sv$chrom2==sv_chr,]$start2))


          all_interval=diff(all_brk)
          all_interval_total=c(all_interval_total,all_interval)

          ##### mean absolute deviation with normalize ########

          dispersion_mean=mean(abs(as.numeric(all_interval-mean(all_interval))))/mean(as.numeric(all_interval))
          dispersion_mean=ifelse(is.na(dispersion_mean),0,dispersion_mean)
          all_dispersion=c(all_dispersion,dispersion_mean)

          sd_mean=sd(all_interval)/mean(as.numeric(all_interval))
          sd_mean=ifelse(is.na(sd_mean),0,sd_mean)
          all_sd=c(all_sd,sd_mean)

          brk_dispersion_mean=mean(abs(as.numeric(all_brk-mean(all_brk))))/mean(as.numeric(all_brk))
          brk_dispersion_mean=ifelse(is.na(brk_dispersion_mean),0,brk_dispersion_mean)
          all_brk_dispersion=c(all_brk_dispersion,brk_dispersion_mean)


        }
        dispersion_mean_total=mean(abs(as.numeric(all_interval_total-mean(all_interval_total))))/mean(as.numeric(all_interval_total))

        dispersion_mean_total=ifelse(is.na(dispersion_mean_total),0,dispersion_mean_total)
        all_dispersion_mean=mean(all_dispersion)

        cluster$Brk_dispersion_MAD[i] <- all_dispersion_mean

        cluster$brk_10mb[i] <-  brk10
        cluster$brk_1mb[i] <-  brk1
        cluster$sv_avg_size[i] <- sv_avg
        cluster$Brk_dispersion_MAD_mean_total[i] <- dispersion_mean_total

      }
    }
  }
print("CGR feature computing is done!")
cluster=cluster[c("sample","cluster_id","max_CNV","Loss_size_percentage","Gain_size_percentage","max_telo_loss_percentage","Brk_dispersion_MAD_mean_total")]
 # chrss_cluster_feature=list("cluster_feature"=cluster,"cnv_baseline"=cnvtotal)
 chrss_cluster_feature=list("cluster_feature"=cluster)
 filename=ifelse(prefix=="","CGR_feature_matrix.csv",paste0(prefix,"_CGR_feature_matrix.csv"))
 write.csv(cluster,filename,row.names = F)
 return (chrss_cluster_feature)
}








