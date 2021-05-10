#' starfish_sig
#'
#' This function loads the CGR feature matrix and perform clustering and classification.
#'
#'
#' @param cluster_feature the CGR feature matrix, output from starfish_feature
#' @param prefix prefix for intermediate files, default is none
#' @param class_method "class" based on a pre-constructed classifier or "cluster" based on de-novo unsupervised clustering, default is "class"
#' @return The signature classification table and plot if "class" is selected.The clustering results are stored under "CGR_cluster" folder if "cluster" is selected. $cluster_id contains CGR event IDs, and  $CGR_signature shows the signature names used in the PCAWG paper.
#' @export
starfish_sig=function(cluster_feature,prefix="",class_method="class",pcawg_feature=pcawg_chrss_raw,class_model=nn_model){

  chrss_select=cluster_feature

  chrss_select$log2_CN_state=ifelse(chrss_select$CN_state!=0,log(chrss_select$CN_state,2),0)
  chrss_select$log_max_CN=ifelse(chrss_select$max_CNV!=0,log(chrss_select$max_CNV,10),0)

  chrss_select$Loss_number_10Mb=ifelse(chrss_select$Loss_number_10Mb!=0,chrss_select$Loss_number_10Mb,min(chrss_select$Loss_number_10Mb[chrss_select$Loss_number_10Mb!=min(chrss_select$Loss_number_10Mb)])/10)
  chrss_select$Gain_number_10Mb=ifelse(chrss_select$Gain_number_10Mb!=0,chrss_select$Gain_number_10Mb,min(chrss_select$Gain_number_10Mb[chrss_select$Gain_number_10Mb!=min(chrss_select$Gain_number_10Mb)])/10)

  chrss_select$Loss_number_10Mb_log=log(chrss_select$Loss_number_10Mb,10)
  chrss_select$Gain_number_10Mb_log=log(chrss_select$Gain_number_10Mb,10)

  chrss_select$dataset="non-pcawg"

  pcawg_feature$dataset="pcawg"
  pcawg_col=colnames(pcawg_feature)
  chrss_select=chrss_select[,colnames(chrss_select) %in% c(pcawg_col)]
  select_col=colnames(chrss_select)
  pcawg_feature=pcawg_feature[,colnames(pcawg_feature) %in% c(select_col)]

  chrss_select_all=data.frame()

  ######### scale with pcawg feature #######

  for (n in 1:nrow(chrss_select)){
    chrss_select_i=chrss_select[n,]
    chrss_select_i_b=rbind(chrss_select_i,pcawg_feature)
    chrss_select_i_b[,c("Brk_dispersion_MAD_mean_total","Loss_size_percentage","Gain_size_percentage","log_max_CN","max_telo_loss_percentage")]=apply(chrss_select_i_b[,c("Brk_dispersion_MAD_mean_total","Loss_size_percentage","Gain_size_percentage","log_max_CN","max_telo_loss_percentage")],2,scale)
    chrss_select_i_n=chrss_select_i_b[chrss_select_i_b$dataset=="non-pcawg",]
    chrss_select_all=rbind(chrss_select_all,chrss_select_i_n)
}
  chrss_select_scale=chrss_select_all[chrss_select_all$dataset=="non-pcawg",]

  ######### classifier ###########
  if (class_method=="class") {




    feature_vector=c("Brk_dispersion_MAD_mean_total","Loss_size_percentage","Gain_size_percentage","log_max_CN","max_telo_loss_percentage")

    chrss_class_select=chrss_select_scale[c("cluster_id",feature_vector)]
    x=chrss_class_select[,-1]
    nn_train <- neuralnet::compute(nn_model, x)
    train_result=as.data.frame(nn_train$net.result)
    colnames(train_result)=c("hc_cluster1","hc_cluster2","hc_cluster3","hc_cluster4","hc_cluster5","hc_cluster6")
    chrss_class_select=cbind(chrss_class_select,train_result)
    chrss_class_select$CGR_signature=colnames(chrss_class_select[,c("hc_cluster1","hc_cluster2","hc_cluster3","hc_cluster4","hc_cluster5","hc_cluster6")])[apply(chrss_class_select[,c("hc_cluster1","hc_cluster2","hc_cluster3","hc_cluster4","hc_cluster5","hc_cluster6")],1,which.max)]

    chrss_class_select$CGR_signature=sapply(chrss_class_select$CGR_signature,switch,
                                    "hc_cluster1"="6 Hourglass",
                                    "hc_cluster2"="5 Large gain",
                                    "hc_cluster3"="1 ecDNA/double minutes",
                                    "hc_cluster4"="4 Micronuclei",
                                    "hc_cluster5"="3 Large loss",
                                    "hc_cluster6"="2 BFB cycles/chromatin bridge")


    chrss_class_select=chrss_class_select[order(chrss_class_select$CGR_signature,chrss_class_select$cluster_id),]

    sampleorder=chrss_class_select$cluster_id
    chrss_class_select$cluster_id=factor(chrss_class_select$cluster_id,levels = sampleorder)

    chrss_class_select$CGR_signature=factor(chrss_class_select$CGR_signature,levels = c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge", "3 Large loss", "4 Micronuclei","5 Large gain","6 Hourglass"))

    chrss_class_select[,c(2:6)]=apply(chrss_class_select[,c(2:6)],2,function(x) x - min(x))

    p4<-ggplot(chrss_class_select,aes(cluster_id,y=1,fill=CGR_signature))+geom_tile()+scale_fill_manual("",values = c("1 ecDNA/double minutes"="#e78ac3","2 BFB cycles/chromatin bridge"="#fc8d62", "3 Large loss"="#a6d854", "4 Micronuclei"="#66c2a5","5 Large gain"="#ffd92f","6 Hourglass"="#8da0cb"))+theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),legend.position="right",legend.text=element_text(size=15),legend.title = element_blank(),legend.spacing=unit(5,'cm'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border=element_blank(),axis.line.x=element_blank(), axis.line.y=element_blank())

    stack2=chrss_class_select[c("cluster_id","Brk_dispersion_MAD_mean_total","Loss_size_percentage","Gain_size_percentage","log_max_CN","max_telo_loss_percentage")]

    stack2=melt(stack2)

    stack2$variable=sapply(stack2$variable,switch,
                           "Brk_dispersion_MAD_mean_total"="Breakpoint dispersion score",
                           "Loss_size_percentage"="Loss percentage",
                           "Gain_size_percentage"="Gain percentage",
                           "log_max_CN"="log10(max CN)",
                           "max_telo_loss_percentage"="Telomere loss percentage")

    stack2$variable=factor(stack2$variable,levels=c("Breakpoint dispersion score","Loss percentage","Gain percentage","Telomere loss percentage","log10(max CN)"))

    stack2$cluster_id=factor(stack2$cluster_id,levels=sampleorder)

    p3=ggplot(stack2,aes(x=cluster_id,y=value),fill="darkgrey")+geom_bar(stat = "identity",width=1.1)+theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.x = element_blank(),strip.text.y = element_text(size=20,angle = 0,hjust = 0), strip.background = element_blank(), axis.ticks.y = element_blank())+xlab("")+ylab("")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),legend.text=element_text(size=10),panel.border=element_blank(),legend.position="right",axis.line.x=element_blank(), axis.line.y=element_blank())+ guides(fill = guide_legend(nrow = 1))+guides(fill=guide_legend(nrow=6,byrow=TRUE))+facet_grid(variable~.,scale="free")

    leg3=get_legend(p3)

    leg4=get_legend(p4)

    blank_p <- patchwork::plot_spacer() + theme_void()




    tree_total=cowplot::plot_grid(p4+ theme(legend.position = "none"),p3+ theme(legend.position = "none"), align = "v", ncol = 1, rel_heights = c(1,6),axis=("lr"))


    tree_final=ggdraw(tree_total) + draw_plot(leg4, x = 0.35, y = 0.41, scale = 0.5)

    cluster_pdf=ifelse(prefix=="","pcawg_6signatures_class.pdf",paste0(prefix,"_pcawg_6signatures_class.pdf"))
    pdf(cluster_pdf,height = 9,width = 12)
    print(tree_final)
    dev.off()
    print("Signature classification is done!")

    filename=ifelse(prefix=="","pcawg_6signatures_class.csv",paste0(prefix,"_pcawg_6signatures_class.csv"))
    write.csv(chrss_class_select,filename,row.names = F)
    return(chrss_class_select)
  } else if (class_method=="cluster") {

    feature_vector=c("Brk_dispersion_MAD_mean_total","Loss_size_percentage","Gain_size_percentage","log_max_CN","max_telo_loss_percentage")
    chrss_cluster=chrss_select_scale[c("cluster_id",feature_vector)]
    chrss_cluster2=setNames(data.frame(t(chrss_cluster[,-1])), chrss_cluster[,1])
    row.names(chrss_cluster2)=colnames(chrss_cluster[,-1])
    ssnmfw2=as.matrix(as.data.frame(chrss_cluster2))

    if(ncol(ssnmfw2)<10){
      print("Too few samples to run clustering, please try with at least 10 samples or run classification with 'class'.")
    } else {
    ########### ConsensusCluster ########
    cluster_results = ConsensusClusterPlus(ssnmfw2,maxK=10,reps=5000,pItem=0.9,pFeature=1,clusterAlg="pam",distance="euclidean",innerLinkage = "ward.D2",finalLinkage = "ward.D2",plot="pdf",writeTable=TRUE,title="CGR_cluster")
    print("Clustering is done!")
}

  } else if (class_method!="class"|class_method!="cluster") {

    print("Wrong method, please select from 'class' or 'cluster'.")

  }

}
