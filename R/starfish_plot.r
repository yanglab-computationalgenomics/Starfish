#' starfish_plot
#'
#' This function loads a SV dataframe, a CNV dataframe and "connected" CGR regions reported by starfish_link, to draw the connected regions in a linear setting.
#'
#'
#' @param sv_file a SV dataframe with "chrom1","end1", "chrom2","end2","svtype","strand1","strand2","sample"
#' @param cnv_file a CNV dataframe with "chromosome","start","end","total_cn","sample"
#' @param cgr "connected" CGR regions, output from starfish_link
#' @param genome_v genome version of "hg19" or "hg38", default is "hg19"
#' @return plots of "connected" CGR regions.
#' @export


starfish_plot=function(sv_file,cnv_file,cgr,genome_v="hg19"){

`%!in%` = Negate(`%in%`)
chrlist=c(as.character(1:22),"X")

cnv_file$chromosome=as.character(cnv_file$chromosome)
cnv_file$chromosome=gsub("Chr|chr","",cnv_file$chromosome)
sv_file$chrom1=as.character(sv_file$chrom1)
sv_file$chrom1=gsub("Chr|chr","",sv_file$chrom1)
sv_file$chrom2=as.character(sv_file$chrom2)
sv_file$chrom2=gsub("Chr|chr","",sv_file$chrom2)
cnv_total=cnv_file[cnv_file$chromosome %in% chrlist,]
cnv_total=cnv_total[order(cnv_total$sample,cnv_total$start,cnv_total$total_cn),]


chrss=cgr

if(nrow(cnv_total[!cnv_total$sample %in% chrss$sample, ])>0){
  cnv_no_sample=cnv_total[!cnv_total$sample %in% chrss$sample, ]
  print(paste0("There is no CNV file for sample ",unique(cnv_no_sample$sample)))
}

if(nrow(sv_file[!sv_file$sample %in% chrss$sample, ])>0){
  sv_no_sample=sv_file[!sv_file$sample %in% chrss$sample, ]
  print(paste0("There is no SV file for sample ",unique(sv_no_sample$sample)))
}



cluster_list=unique(chrss$cluster_id)

hg19_chr=unique(hg19_chr_length[c("chrom","chr_size")])
hg38_chr=unique(hg38_chr_length[c("chrom","chr_size")])
hg19_chr=hg19_chr[order(as.numeric(hg19_chr$chrom)),]
hg38_chr=hg38_chr[order(as.numeric(hg38_chr$chrom)),]
hg19_chr$chrom=paste0("chr",hg19_chr$chrom)
hg38_chr$chrom=paste0("chr",hg38_chr$chrom)

if(genome_v=="hg19"){
  genome <- GenomeInfoDb::Seqinfo(
    seqnames = hg19_chr$chrom,
    seqlengths = hg19_chr$chr_size,
    genome = "hg19"
  )
  hg_genome=hg19_chr
  } else if(genome_v=="hg38"){
    genome <- GenomeInfoDb::Seqinfo(
      seqnames = hg38_chr$chrom,
      seqlengths = hg38_chr$chr_size,
      genome = "hg38"
    )
    hg_genome=hg38_chr
    } else {
    print("Wrong genome version!")
  }


GenomeInfoDb::seqlevelsStyle(genome) <- "NCBI"
chromosomes = as.character(c(seq(22), "X"))
genome = genome[chromosomes]


genome.ranges = GRanges(
  seqnames = seqnames(genome),
  strand = "*",
  ranges = IRanges(start = 1, width = seqlengths(genome)),
  seqinfo = genome
)


for (j in 1:length(cluster_list)){
  chrss_j=chrss[chrss$cluster_id==cluster_list[j],]

  seq=unlist(strsplit(unique(chrss_j$link_chromosome), split="_"))
  seq=gsub("X","23",seq)

  seq=as.numeric(seq)

  cnv=cnv_total[cnv_total$sample==unique(chrss_j$sample),]
  cnv_raw=cnv
  cnv_backup=cnv[c("chromosome","start","end","total_cn")]
  cnv_backup=cnv_backup[cnv_backup$chromosome %in% chrlist,]
  cnv$chromosome=gsub("X","23",cnv$chromosome)


  sv.data.raw <- sv_file[sv_file$sample==unique(chrss_j$sample),]
  sv.data=sv.data.raw

  sv.data.raw=sv.data.raw[sv.data.raw$chrom1 %in% chrlist & sv.data.raw$chrom2 %in% chrlist, ]
  # colnames(sv.data)[11] <- "svtype"
  sv.data$chrom1=gsub("X","23",sv.data$chrom1)
  sv.data$chrom2=gsub("X","23",sv.data$chrom2)


  sv.data=sv.data[sv.data$chrom1 %in% seq | sv.data$chrom2 %in% seq,]
  sv.data$pos1=ifelse(sv.data$chrom1 %!in% seq,sv.data$end2-2,sv.data$end1)
  sv.data$pos2=ifelse(sv.data$chrom2 %!in% seq,sv.data$end1+2,sv.data$end2)
  sv.data$chrom2=ifelse(sv.data$chrom2 %!in% seq,sv.data$chrom1,sv.data$chrom2)
  sv.data$chrom1=ifelse(sv.data$chrom1 %!in% seq,sv.data$chrom2,sv.data$chrom1)

  sv.data=sv.data[!duplicated(sv.data[,c('chrom1','pos1','chrom2','pos2','svtype')]),]


  sv.data$chrom1=gsub("23","X",sv.data$chrom1)
  sv.data$chrom2=gsub("23","X",sv.data$chrom2)


  if (length(sv.data$chrom1)>0){
    cnv=cnv[cnv$chromosome %in% seq,]
    region=chrss[chrss$cluster_id==cluster_list[j],]
    region$chr=gsub("X","23",region$chr)



    ##### Get the SV chromosomes...
    chr1 <- paste0("chr", sv.data$chrom1)
    chr2 <- paste0("chr", sv.data$chrom2)

    ##### ...positions
    pos1 <- sv.data$pos1
    pos2 <- sv.data$pos2

    ##### ... and the events type
    type <- sv.data$svtype

    genome.length <- sum(hg_genome$chr_size[seq])
    chrs_fake_starts <- vector("list", length(seq))
    chrs_fake_starts  <- setNames(chrs_fake_starts,  hg_genome$chrom[seq] )
    chrs_fake_starts[[hg_genome$chrom[seq[1]]]] <- 0
    length_sum <- 0
    if (length(seq)>1){
      for ( i in 2:length(seq) ) {


        length_sum <- length_sum + as.numeric(hg_genome$chr_size[[seq[i-1]]])+50000000
        chrs_fake_starts[[hg_genome$chrom[seq[i]]]] <- length_sum
      }
    }

    # assign 50k-500k to cnv and region to make the segement visible

    if (length(seq)<3){

      for (i in 1:length(seq)) {

        cnv[cnv$chromosome==seq[i],]$start=cnv[cnv$chromosome==seq[i],]$start+chrs_fake_starts[[i]]
        cnv[cnv$chromosome==seq[i],]$end=cnv[cnv$chromosome==seq[i],]$end+chrs_fake_starts[[i]]+50000

      }

      for (i in 1:length(seq)) {

        region[region$chr==seq[i],]$start=region[region$chr==seq[i],]$start+chrs_fake_starts[[i]]
        region[region$chr==seq[i],]$end=region[region$chr==seq[i],]$end+chrs_fake_starts[[i]]+50000

      }
    } else {

      for (i in 1:length(seq)) {

        cnv[cnv$chromosome==seq[i],]$start=cnv[cnv$chromosome==seq[i],]$start+chrs_fake_starts[[i]]
        cnv[cnv$chromosome==seq[i],]$end=cnv[cnv$chromosome==seq[i],]$end+chrs_fake_starts[[i]]+500000

      }


      for (i in 1:length(seq)) {

        region[region$chr==seq[i],]$start=region[region$chr==seq[i],]$start+chrs_fake_starts[[i]]
        region[region$chr==seq[i],]$end=region[region$chr==seq[i],]$end+chrs_fake_starts[[i]]+500000

      }
    }



    chrs_fake_label.pos <- vector("list", length(seq))
    chrs_fake_label.pos <- setNames(chrs_fake_label.pos,  hg_genome$chrom[seq] )

    for ( i in 1:length(chrs_fake_starts) ) {

      chrs_fake_label.pos[[hg_genome$chrom[seq[i]]]] <- hg_genome$chr_size[[seq[i]]]/2 + chrs_fake_starts[[hg_genome$chrom[[seq[i]]]]]


    }


    chrs_boundary.pos <- vector("list", length(seq))
    chrs_boundary.pos <- setNames(chrs_fake_label.pos,  hg_genome$chrom[seq] )

    for ( i in 1:length(chrs_fake_starts) ) {

      chrs_boundary.pos[[hg_genome$chrom[seq[i]]]] <- hg_genome$chr_size[[seq[i]]] + chrs_fake_starts[[hg_genome$chrom[[seq[i]]]]]


    }



    pos1_fake <- vector("list", nrow(sv.data))
    pos2_fake <- vector("list", nrow(sv.data))

    for ( i in 1:nrow(sv.data) ) {


      pos1_fake[[i]] <- chrs_fake_starts[[chr1[i]]] + pos1[i]
      pos2_fake[[i]] <- chrs_fake_starts[[chr2[i]]] + pos2[i]
    }


    beziers.height <- runif(nrow(sv.data), 0.5, 1)
    beziers.mid <- unlist(pos1_fake)+(unlist(pos2_fake)-unlist(pos1_fake))/2

    beziers <- data.frame(
      x = c(rbind( unlist(pos1_fake), beziers.mid, unlist(pos2_fake) )),
      y = c(rbind( 0.2, beziers.height, 0.2 ) ),
      group = rep( paste( chr1, pos1, chr2, pos2,type, sep="_" ), each=3),
      svtype = rep( type, each=3)

    )


    beziers_intra=beziers[beziers$svtype!='TRA',]
    beziers_inter=beziers[beziers$svtype=='TRA',]
    beziers_inter$y=beziers_inter$y+0.5

    cnv=na.omit(cnv)
    cnv$chromosome=gsub("23","X",cnv$chromosome)



    boundary=as.data.frame(c(unlist(chrs_boundary.pos),unlist(chrs_fake_starts)))
    colnames(boundary)="position"
    ymax=max(cnv$total_cn)

    title_position=tail(unlist(chrs_boundary.pos),n=1)/2

    plot_range=ifelse(length(seq)>1,(tail(unlist(chrs_boundary.pos),n=1)+50000000),(tail(unlist(chrs_boundary.pos),n=1)+5000000))

    p2=ggplot(cnv)+ggplot2::geom_segment(aes(x=(start), y=total_cn, xend=(end), yend=total_cn),size=1)+ggplot2::geom_segment(data=boundary,aes(x = position , xend = position, y = 0, yend = as.numeric(ymax)),linetype=2, colour = 'black', size = 0.2)+ scale_x_continuous(expand = c(0, 0),breaks = NULL,limits=c(0,plot_range)) +scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +ylab("CNV")+xlab("") +theme_bw()+theme(axis.title.x=element_blank(), axis.title.y=element_text(size=6),axis.text.x= element_blank(), axis.ticks.x=element_blank(),axis.text.y= element_text(size=3),panel.background=element_blank(),plot.margin = margin(0, 0.1, 0.1, 0.1, "cm"))

    label_size=3-log2(length(seq))/2
    title_size=1.5

    region$chr=gsub("23","X",region$chr)


    link_region=region[region$chrss_status=="link",]
    chrss_region=region[region$chrss_status=="chrss",]


    title_text=paste0(unique(chrss_j$histology_abbreviation)," ", unique(chrss_j$sample))

    if (length(link_region$chr)==0) {
      pcolor=ggplot() + geom_bezier(aes(x= x, y = y, group = group, color = svtype ), data = beziers_inter, show.legend = TRUE, size = 0.1)+geom_bezier(aes(x= x, y = y, group = group, color = svtype ), data = beziers_intra, show.legend = TRUE, size = 0.1) + scale_x_continuous(expand = c(0, 0),limits=c(0,plot_range))+xlab("")+

        ##### Remove default axes labels and grey backgroud
        theme(axis.title.x=element_text(size=title_size), axis.text.x= element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=7), axis.text.y= element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank(),panel.grid = element_blank(),
              ##### ...and the grey backgroud
              panel.background = element_rect(fill = NA),
              ##### ...change the legend parameters
              legend.title=element_text(size=6), legend.background=element_blank(),legend.text=element_text(size=6), legend.key.size = unit(0.3,"line"), legend.key= element_blank(), legend.position = c(0.5,0.75),legend.direction="horizontal",plot.margin = margin(-0.1, 0.1, -0.2, 0.5, "cm")) +

        ##### Set the axes limits

        scale_y_continuous(limits = c(-0.2, 1.8)) +

        ##### Add chromosomes boundaries

        ggplot2::geom_segment(aes(x = c(unlist(chrs_boundary.pos),unlist(chrs_fake_starts)) , xend = c(unlist(chrs_boundary.pos),unlist(chrs_fake_starts)), y = 0, yend = 0.5),linetype=2, colour = 'black', size = 0.2)+ggplot2::geom_segment(data=chrss_region,aes(x=start, y=0.1, xend=end, yend=0.1),color="red",size=0.5)+

        labs( color = "svtype") +scale_color_manual("SV type ",values = c("DEL"="#428BCA","DUP"="#d9534f","h2hINV"="#FFC425","t2tINV"="#5cb85c","TRA"="#aa96dc"))+annotate(geom = 'text', label = title_text, x = title_position, y = 1.7, size = title_size)+annotate(geom = 'text', label = names(chrs_fake_label.pos), x = unlist(chrs_fake_label.pos), y = -0.1, size = label_size)+ylab("SV")




    } else if (length(chrss_region$chr)==0){
      pcolor=ggplot() + geom_bezier(aes(x= x, y = y, group = group, color = svtype ), data = beziers_inter, show.legend = TRUE, size = 0.1)+geom_bezier(aes(x= x, y = y, group = group, color = svtype ), data = beziers_intra, show.legend = TRUE, size = 0.1) + scale_x_continuous(expand = c(0, 0),limits=c(0,plot_range))+xlab("")+

        ##### Remove default axes labels and grey backgroud
        theme(axis.title.x=element_text(size=title_size), axis.text.x= element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=7), axis.text.y= element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank(),panel.grid = element_blank(),panel.background = element_rect(fill = NA),legend.title=element_text(size=6), legend.background=element_blank(),legend.text=element_text(size=6), legend.key.size = unit(0.3,"line"), legend.key= element_blank(), legend.position = c(0.5,0.75),legend.direction="horizontal",plot.margin = margin(-0.1, 0.1, -0.2, 0.5, "cm")) +

        ##### Set the axes limits

        scale_y_continuous(limits = c(-0.2, 1.8)) +

        ##### Add chromosomes boundaries

        ggplot2::geom_segment(aes(x = c(unlist(chrs_boundary.pos),unlist(chrs_fake_starts)) , xend = c(unlist(chrs_boundary.pos),unlist(chrs_fake_starts)), y = 0, yend = 0.5),linetype=2, colour = 'black', size = 0.2)+ggplot2::geom_segment(data=link_region,aes(x=start, y=0.1, xend=end, yend=0.1),color="blue",size=0.5)+

        labs( color = "svtype") +

        ##### Add chromosomes labels
        scale_color_manual("SV type ",values = c("DEL"="#428BCA","DUP"="#d9534f","h2hINV"="#FFC425","t2tINV"="#5cb85c","TRA"="#aa96dc"))+annotate(geom = 'text', label = title_text, x = title_position, y = 1.7, size = title_size)+annotate(geom = 'text', label = names(chrs_fake_label.pos), x = unlist(chrs_fake_label.pos), y = -0.1, size = label_size)+ylab("SV")



    } else {
      pcolor=ggplot() + geom_bezier(aes(x= x, y = y, group = group, color = svtype ), data = beziers_inter, show.legend = TRUE, size = 0.1)+geom_bezier(aes(x= x, y = y, group = group, color = svtype ), data = beziers_intra, show.legend = TRUE, size = 0.1) + scale_x_continuous(expand = c(0, 0),limits=c(0,plot_range))+xlab("")+

        ##### Remove default axes labels and grey backgroud
        theme(axis.title.x=element_text(size=title_size), axis.text.x= element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=7), axis.text.y= element_blank(), axis.ticks.y=element_blank(),panel.border = element_blank(),panel.grid = element_blank(),
              ##### ...and the grey backgroud
              panel.background = element_rect(fill = NA),
              ##### ...change the legend parameters
              legend.title=element_text(size=6), legend.background=element_blank(),legend.text=element_text(size=6), legend.key.size = unit(0.3,"line"), legend.key= element_blank(), legend.position = c(0.5,0.75),legend.direction="horizontal",plot.margin = margin(-0.1, 0.1, -0.2, 0.5, "cm")) +

        ##### Set the axes limits

        scale_y_continuous(limits = c(-0.2, 1.8)) +

        ##### Add chromosomes boundaries

        ggplot2::geom_segment(aes(x = c(unlist(chrs_boundary.pos),unlist(chrs_fake_starts)) , xend = c(unlist(chrs_boundary.pos),unlist(chrs_fake_starts)), y = 0, yend = 0.5),linetype=2, colour = 'black', size = 0.2)+ggplot2::geom_segment(data=link_region,aes(x=start, y=0.1, xend=end, yend=0.1),size=0.5,color="blue")+ggplot2::geom_segment(data=chrss_region,aes(x=start, y=0.1, xend=end, yend=0.1),size=0.5,color="red")+

        labs( color = "svtype") +

        ##### Add chromosomes labels
        scale_color_manual("SV type ",values = c("DEL"="#428BCA","DUP"="#d9534f","h2hINV"="#FFC425","t2tINV"="#5cb85c","TRA"="#aa96dc"))+annotate(geom = 'text', label = title_text, x = title_position, y = 1.7, size = title_size)+annotate(geom = 'text', label = names(chrs_fake_label.pos), x = unlist(chrs_fake_label.pos), y = -0.1, size = label_size)+ylab("SV")


    }

    p3=plot_grid(pcolor,p2,ncol=1, align="v",rel_heights =c(2,1))



    pdfname=paste0(unique(chrss_j$code4),cluster_list[j],"_chrss.pdf")
    ggbio::ggsave(pdfname,plot=p3,width=3,height=3,device = cairo_pdf)

  }
}

print("Plotting is done!")
}




