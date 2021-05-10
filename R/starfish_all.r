#'  starfish_all
#'
#' This function loads SV, CNV and gender files, identifies "connected" complex genome rearrangement (CGR) regions and complex SVs, then constructs a CGR event vs. feature matrix and infers CGR signatures based on clustering and classification.
#'
#'
#' @param sv_file a SV dataframe with 8 columns: "chrom1","end1", "chrom2","end2","svtype" (DEL,DUP,h2hINV,t2tINV,TRA),"strand1" (+/-) and "strand2" (+/-),"sample". Other svtypes like INV, INS, BND are not accepted
#' @param cnv_file a CNV dataframe with 5 columns: "chromosome","start","end","total_cn", and "sample". "total_cn" should contain absolute copy numbers.
#' @param gender_file a sample table with 2 columns: "sample" and "gender" ("Female, "female","F","f","Male","male","M","m"). If gender is unknown, any other characters could be put here like "unknown".
#' @param prefix_total the prefix for all intermediate files, default is none
#' @param genome_version which genome assembly was used to call SV and CNV. It should be "hg19" or "hg38", default is "hg19"
#' @param plot the logical value of plotting "connected" CGRs, default is TRUE
#' @param smethod "class" based on a pre-constructed classifier from PCAWG dataset or "cluster" based on de-novo unsupervised clustering, default is "class"
#' @return a list of files: $complex_sv contains complex SVs, $starfish_call contains connected CGRs, $feature_matrix contains CGR feature matrix, and $chrss_signature shows the signatures computed and inferred from PCAWG dataset. The signature classification table and plot will be generated if "class" is selected, otherwise clustering matrices and plots will be stored under "CGR_cluster" folder with K from 2 to 10.
#' @export

starfish_all=function(sv_file,cnv_file,gender_file,prefix_total="",genome_version="hg19",plot=TRUE,smethod="class"){

  starfish_link=starfish_link(sv_file,prefix = prefix_total)


  starfish_feature=starfish_feature(starfish_link$starfish_call,starfish_link$interleave_tra_complex_sv,cnv_file,gender_file,prefix=prefix_total,genome_v=genome_version)


  starfish_sig=starfish_sig(starfish_feature$cluster_feature,prefix=prefix_total,class_method=smethod,pcawg_feature=pcawg_chrss_raw,class_model=nn_model)

  if(smethod=="class"){

    starfish_output=list("starfish_call"=starfish_link$starfish_call,"complex_sv"=starfish_link$interleave_tra_complex_sv,"feature_matrix"=starfish_feature$cluster_feature,"chrss_signature"=starfish_sig)


  } else if (smethod=="cluster"){

    starfish_output=list("starfish_call"=starfish_link$starfish_call,"complex_sv"=starfish_link$interleave_tra_complex_sv,"feature_matrix"=starfish_feature$cluster_feature)

  }
  if (plot==TRUE) {

    starfish_plot(sv_file,cnv_file,starfish_link$starfish_call,genome_v=genome_version)

  } else {

  }

  return(starfish_output)
}
