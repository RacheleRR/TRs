#?TSS
 #?libraires 
        library(dplyr)
        library(GenomicRanges)
        library(ggplot2)
        library(tidyr)

 #?LOAD
        refflat <- read.delim("/home/rachele/Downloads/refFlat.txt", stringsAsFactors = F, header = F)
        known.expansion <- read.delim("/home/rachele/Downloads/bret_stuff/UCSC_simple_repeats_hg38.period_lte20.txt", stringsAsFactors = F, header = F)
        # known.expansion <- read.delim("~/Downloads/Human Full Genome_repeats.bed", header=FALSE)
        # known.expansion <- fread("/home/rachele/Downloads/MSDBG000001_perf.tsv")
        #!NOT SURE WHAT TO USE AS KNOWN EXPANSION TRF or MSDBdatabase or UCSC databse ??(NOT USE THE USCS DATABASE)


        detected.expansion <- read.delim("/home/rachele/EHDN_DBSCAN_new_CALLS_OCT/Result/output_DBSCAN_Oct_24/DBSCAN/EHdn_for_DBSCAN_combined_counts_compared.bed", header = F)
        outliers_1_case <- read.delim("~/outliers_1_case_no_split.tsv")
        outliers_case <- read.delim("~/outliers_case_over1.tsv")
        outliers_control_mixed <- read.delim("~/outliers_control_mixed.tsv")


 #?CLEAN UP 
  #CHANGE+ UNIFY  CHROMOSome names 
        known.expansion <- known.expansion[known.expansion$V1 %in% paste0("chr", c(1:22)), ]
        detected.expansion <- detected.expansion[detected.expansion$V1 %in% paste0("chr", c(1:22)), ]
        outliers_1_case <- outliers_1_case[outliers_1_case$contig %in% paste0("chr", c(1:22)), ]
        outliers_case <- outliers_case[outliers_case$contig %in% paste0("chr", c(1:22)), ]
        outliers_control_mixed <- outliers_control_mixed[outliers_control_mixed$contig %in% paste0("chr", c(1:22)), ]
                
  #CHANGE NAMES 
        common.expansion <- known.expansion
        colnames(common.expansion) <-c("contig", "start","end")
        colnames(detected.expansion) <-c("contig", "start","end")
        #for when common.expansion is hugh knwonstr file
        # colnames(common.expansion)[1] <- "contig"
        #  colnames(common.expansion)[2] <- "start"
        # colnames(common.expansion)[3] <- "end"

  #ADJUST SIZES
        refflat$tis100kbstart <- ifelse(refflat$V4 == "+", refflat$V5-10000, refflat$V6-10000)
        refflat$tis100kbend <- ifelse(refflat$V4 == "+", refflat$V5+10000, refflat$V6+10000)
        refflat.g <- GRanges(refflat$V3, IRanges(refflat$tis100kbstart, refflat$tis100kbend))
 #? GRANGES
 #GRANGES (WITH REDUCE)
        common.expansion <- data.frame(reduce(GRanges(common.expansion$contig, IRanges(common.expansion$start, common.expansion$end), "*")))
        detected.expansion <- data.frame(reduce(GRanges(detected.expansion$contig, IRanges(detected.expansion$start, detected.expansion$end), "*")))
        outliers_1_case <- data.frame(reduce(GRanges(outliers_1_case$contig, IRanges(outliers_1_case$start, outliers_1_case$end), "*")))
        outliers_case <- data.frame(reduce(GRanges(outliers_case$contig, IRanges(outliers_case$start, outliers_case$end), "*")))
        outliers_control_mixed <- data.frame(reduce(GRanges(outliers_control_mixed$contig, IRanges(outliers_control_mixed$start, outliers_control_mixed$end), "*")))


            
 #GRANGES WITHOUT REDUCE
            common.expansion.g <- GRanges(common.expansion$seqnames, IRanges(common.expansion$start, common.expansion$end), "*")
            detected.expansion.g <- GRanges(detected.expansion$seqnames, IRanges(detected.expansion$start, detected.expansion$end), "*")
            outliers_1_case.g <- GRanges(outliers_1_case$seqnames, IRanges(outliers_1_case$start, outliers_1_case$end), "*")
            outliers_case.g <- GRanges(outliers_case$seqnames, IRanges(outliers_case$start, outliers_case$end), "*")
            outliers_control_mixed.g <- GRanges(outliers_control_mixed$seqnames, IRanges(outliers_control_mixed$start, outliers_control_mixed$end), "*")

 #? FUNTION
 #CREATE FUNCTION
            getdistance <- function(refflat.g, Outliers_2.g, refflat, Outliers_2){
                    rare.olap <- data.frame(findOverlaps(Outliers_2.g, refflat.g))
                    rare.olap$strain <- refflat$V4[rare.olap$subjectHits]
                    rare.olap$tss <- ifelse(rare.olap$strain == "+", 
                                            refflat$tis100kbstart[rare.olap$subjectHits] + 10000, 
                                            refflat$tis100kbend[rare.olap$subjectHits] - 10000)
                    
                    rare.olap$expansion.start <- Outliers_2$start[rare.olap$queryHits]
                    rare.olap$expansion.end <- Outliers_2$end[rare.olap$queryHits]
                    rare.olap$mid.point <- (rare.olap$expansion.start + rare.olap$expansion.end) / 2
                        
                    rare.olap$distance <- rare.olap$tss - rare.olap$mid.point
                    rare.olap$distance <- ifelse(rare.olap$strain == "+", rare.olap$distance, -rare.olap$distance)
                    rare.olap <- rare.olap[order(abs(rare.olap$distance)), ]
                    rare.olap <- rare.olap[!duplicated(rare.olap$queryHits), ]
                    return(rare.olap)
                    }

 #APPLYFUNCTION +create rare column
            outliers_1_case.distance <- getdistance(refflat.g, outliers_1_case.g, refflat, outliers_1_case)
            outliers_1_case.distance$rarity <- "outliers_1_case"

            outliers_case.distance <- getdistance(refflat.g, outliers_case.g, refflat, outliers_case)
            outliers_case.distance$rarity <- "outliers_case"

            outliers_control_mixed.distance <- getdistance(refflat.g, outliers_control_mixed.g, refflat, outliers_control_mixed)
            outliers_control_mixed.distance$rarity <- "outliers_control_mixed"

            common.distance <- getdistance(refflat.g, common.expansion.g, refflat, common.expansion)
            common.distance$rarity <- "known STRs"

            detected.distance <- getdistance(refflat.g, detected.expansion.g, refflat, detected.expansion)
            detected.distance$rarity <- "all expansions"

 #DO THE statistical test
            wilcox.test(abs(outliers_1_case.distance$distance), abs(outliers_case.distance$distance), alternative = "less")$p.value
            wilcox.test(abs(outliers_1_case.distance$distance), abs(outliers_control_mixed.distance$distance), alternative = "less")$p.value
            wilcox.test(abs(outliers_1_case.distance$distance), abs(common.distance$distance), alternative = "less")$p.value
            wilcox.test(abs(outliers_1_case.distance$distance), abs(detected.distance$distance), alternative = "less")$p.value
 #? GRAPH FULL
 # bind them all togther
            distance <- rbind(
            outliers_1_case.distance,
            outliers_case.distance,
            outliers_control_mixed.distance,
            common.distance,
            detected.distance
            )


 #PLOT
                   TSS <- ggplot(distance, aes(x = distance, y = after_stat(density), color = rarity)) +
                    geom_density(alpha = .15, adjust = 1/10) +
                    geom_vline(xintercept = c(0), lty = 2) + xlab("Distance from TSS") + theme_classic() +
                    theme(panel.border = element_rect(fill = NA), 
                            legend.position = "bottom", 
                            axis.title.x = element_text(size = 14),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14), 
                            axis.title.y = element_text(size = 14)) +
                    geom_hline(yintercept = 0, size = 1, color = "white") +
                    coord_cartesian(xlim = c(-5000, 5000)) + scale_x_continuous(breaks = c(-5000, 0, 5000)) +
                    scale_y_continuous(breaks = c(0, 0.0002)) +
                    scale_color_manual(values = c("#4393C3", "#338833", "#D6604D", "#FF7F00", "#984EA3", "#6A3D9A", "#CC0033" ))
  print(TSS)


 #? GRAPH LESS
 # bind them all togther
            distance <- rbind(
            outliers_1_case.distance,
            outliers_case.distance,
            outliers_control_mixed.distance,
            common.distance,
            detected.distance
            )


 #plot 
                   TSS_less <- ggplot(distance, aes(x = distance, y = after_stat(density), color = rarity)) +
                    geom_density(alpha = .15, adjust = 1/10) +
                    geom_vline(xintercept = c(0), lty = 2) + xlab("Distance from TSS") + theme_classic() +
                    theme(panel.border = element_rect(fill = NA), 
                            legend.position = "bottom", 
                            axis.title.x = element_text(size = 14),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14), 
                            axis.title.y = element_text(size = 14)) +
                    geom_hline(yintercept = 0, size = 1, color = "white") +
                    coord_cartesian(xlim = c(-5000, 5000)) + scale_x_continuous(breaks = c(-5000, 0, 5000)) +
                    scale_y_continuous(breaks = c(0, 0.0002)) +
                    scale_color_manual(values = c("#4393C3", "#338833", "#D6604D", "#FF7F00",  "#984EA3", "#CC0033"))


 print(TSS_less)
 #? SAVE GRAPHS 
    ggsave("TSS_less_plot.png", TSS_less, width = 10, height = 6, dpi = 300)
    ggsave("TSS_full_plot.png", TSS, width = 10, height = 6, dpi = 300)
#



#? SPLICING JUNTION
 #?load lib
        library(dplyr)
        library(GenomicRanges)
        library(ggplot2)
        library(tidyr)
 #?load dataframes
        refflat <- read.delim("/home/rachele/Downloads/refFlat.txt", stringsAsFactors = F, header = F)
        known.expansion <- read.delim("/home/rachele/Downloads/bret_stuff/UCSC_simple_repeats_hg38.period_lte20.txt", stringsAsFactors = F, header = F)
        # known.expansion <- read.delim("~/Downloads/Human Full Genome_repeats.bed", header=FALSE)
        # known.expansion <- fread("/home/rachele/Downloads/MSDBG000001_perf.tsv")
        #NOT SURE WHGAT TO USE AS KNOWN EXPANSION TRF or MSDBdatabase or UCSC databse ??(NOT USE THE USCS DATABASE)
        introns <- read.delim("/home/rachele/Downloads/bret_stuff/hg38_intron_refFlat 1.txt", stringsAsFactors = F)
        appris <- rbind(read.delim("/home/rachele/Downloads/bret_stuff/appris_data.principal.refseq108.hg38.txt", stringsAsFactors = F, header = F),
                        read.delim("/home/rachele/Downloads/bret_stuff/appris_data.principal.refseq109.hg38.txt", stringsAsFactors = F, header = F))
        exons <- read.delim("/home/rachele/Downloads/bret_stuff/hg38_exon_refFlat.txt", stringsAsFactors = F)
        
        detected.expansion <- read.delim("/home/rachele/EHDN_DBSCAN_new_CALLS_OCT/Result/output_DBSCAN_Oct_24/DBSCAN/EHdn_for_DBSCAN_combined_counts_compared.bed", header = F)
        outliers_1_case <- read.delim("~/outliers_1_case_no_split.tsv")
        outliers_case <- read.delim("~/outliers_case_over1.tsv")
        outliers_control_mixed <- read.delim("~/outliers_control_mixed.tsv")


        
 #?CLEAN
  #POSITONS
        appris$V3 <- sapply(sapply(appris$V3, strsplit, "\\."), "[", 1)
        introns <- introns[introns$isoform %in% appris$V3, ]


  # CHANGE+ UNIFY  CHROMOSome names 
        common.expansion <- known.expansion[known.expansion$V1 %in% paste0("chr", c(1:22)), ]
        outliers_1_case <- outliers_1_case[outliers_1_case$contig %in% paste0("chr", c(1:22)), ]
        outliers_case <- outliers_case[outliers_case$contig %in% paste0("chr", c(1:22)), ]
        outliers_control_mixed <- outliers_control_mixed[outliers_control_mixed$contig %in% paste0("chr", c(1:22)), ]
        detected.expansion <- detected.expansion[detected.expansion$V1 %in% paste0("chr", c(1:22)), ]
        common.expansion <- known.expansion
        colnames(common.expansion) <-c("contig", "start","end")
        colnames(detected.expansion) <-c("contig", "start","end")

  #ADJUST SIZES
        outliers_1_case$size <- outliers_1_case$end - outliers_1_case$start
        outliers_1_case <- outliers_1_case[outliers_1_case$size < 10000, ]
        outliers_case$size <- outliers_case$end - outliers_case$start
        outliers_case <- outliers_case[outliers_case$size < 10000, ]
        outliers_control_mixed$size <- outliers_control_mixed$end - outliers_control_mixed$start
        outliers_control_mixed <- outliers_control_mixed[outliers_control_mixed$size < 10000, ]
        
        detected.expansion$size <- detected.expansion$end - detected.expansion$start
        detected.expansion <- detected.expansion[detected.expansion$size < 10000, ]
        common.expansion <- common.expansion[common.expansion$end-common.expansion$start < 10000, ]

  #FIND UNIQUES
        outliers_1_case <- unique(outliers_1_case[, c("contig", "start", "end")])
        outliers_case <- unique(outliers_case[, c("contig", "start", "end")])
        outliers_control_mixed <- unique(outliers_control_mixed[, c("contig", "start", "end")])
        common.expansion <- unique(common.expansion[, c("contig", "start", "end")])
        detected.expansion <- unique(detected.expansion[, c("contig", "start", "end")])

  #CHANGE NAMES 
        colnames(outliers_1_case)[1] <- "chr"
        colnames(outliers_case)[1] <- "chr"
        colnames(outliers_control_mixed)[1] <- "chr"
        colnames(detected.expansion)[1] <- "chr"
        names(common.expansion) <- names(outliers_1_case)


  #gRANGES
        rare.expansion.g <- GRanges(outliers_1_case$chr, IRanges(outliers_1_case$start, outliers_1_case$end), "*")
        common.expansion.g <- GRanges(common.expansion$chr, IRanges(common.expansion$start, common.expansion$end), "*")
        detected.expansion.g <- GRanges(detected.expansion$chr, IRanges(detected.expansion$start, detected.expansion$end), "*")
        outliers_case.g <- GRanges(outliers_case$chr, IRanges(outliers_case$start, outliers_case$end), "*")
        outliers_control_mixed.g <- GRanges(outliers_control_mixed$chr, IRanges(outliers_control_mixed$start, outliers_control_mixed$end), "*")
        introns.g <- GRanges(introns$chr, IRanges(introns$start, introns$end), "*")
        common.expansion.g <- reduce(common.expansion.g)
        common.expansion <- data.frame(common.expansion.g)

  #CHANGE NAMES 
        names(common.expansion) <- c("chr", "start", "end", "width", "strand")
                
 #? FUNTION
  #DEFINE THE FUCNTION
        getsplicedistance <- function(introns.g, rare.expansion.g, introns, rare.expansion, exons){
        introns.real.g <- GRanges(introns$chr, IRanges(introns$start, introns$end), "*")
        rare.olap <- data.frame(findOverlaps(rare.expansion.g, introns.real.g))
        rare.olap$isoform <- introns$isoform[rare.olap$subjectHits]
        dup <- rare.olap$queryHits[which(duplicated(rare.olap[, c("queryHits", "isoform")]))]
        
        
        exons.g <- GRanges(exons$chr, IRanges(exons$start, exons$end), "*")
        rare.exon <-  data.frame(findOverlaps(rare.expansion.g, exons.g))
        rare.exon$sizeOlap <- width(pintersect(rare.expansion.g[rare.exon$queryHits], 
                                                exons.g[rare.exon$subjectHits]))
        rare.exon$sizeExon <- width(rare.expansion.g[rare.exon$queryHits])
        rare.exon <- rare.exon[rare.exon$sizeOlap == rare.exon$sizeExon, ]
        
        dup <- union(dup, rare.exon$queryHits)
        
        rare.olap <- data.frame(findOverlaps(rare.expansion.g, introns.g))
        rare.olap <- rare.olap[!rare.olap$queryHits %in% rare.olap$queryHits[dup], ]
        rare.olap$strand <- introns$strand[rare.olap$subjectHits]
        rare.olap$chr <- rare.expansion$chr[rare.olap$queryHits]
        rare.olap$start.expansion <- rare.expansion$start[rare.olap$queryHits]
        rare.olap$end.expansion <- rare.expansion$end[rare.olap$queryHits]
        rare.olap$start.intron <- introns$start[rare.olap$subjectHits]
        rare.olap$end.intron <- introns$end[rare.olap$subjectHits]
        rare.olap$mid.point <- round((rare.olap$start.expansion + rare.olap$end.expansion) / 2)
        
        rare.olap$donor.distance <- ifelse(rare.olap$strand == "+", rare.olap$mid.point - rare.olap$start.intron, 
                                            rare.olap$end.intron - rare.olap$mid.point)
        rare.olap$acceptor.distance <- ifelse(rare.olap$strand == "+", rare.olap$end.intron - rare.olap$mid.point, 
                                                rare.olap$mid.point - rare.olap$start.intron)
        
        
        rare.olap$nearest.site <- ifelse(abs(rare.olap$donor.distance) < abs(rare.olap$acceptor.distance), "Donor", "Acceptor")
        rare.olap$distance <- ifelse(rare.olap$nearest.site == "Donor", rare.olap$donor.distance, rare.olap$acceptor.distance) %>% abs()
        
        rare.olap <- rare.olap[order(abs(rare.olap$distance)), ]
        rare.olap <- rare.olap[!duplicated(rare.olap$queryHits), ]
        
        rare.olap <- rare.olap[abs(rare.olap$distance) < 10000, ]

        rare.olap$distance[which(rare.olap$nearest.site == "Acceptor")] <-
            -rare.olap$distance[which(rare.olap$nearest.site == "Acceptor")]
        return(rare.olap)
        }

  #apply the fucntion
        rare.distance <- getsplicedistance(introns.g, rare.expansion.g, introns, outliers_1_case, exons)
        detected.distance <- getsplicedistance(introns.g, detected.expansion.g, introns, detected.expansion, exons)
        common.distance <- getsplicedistance(introns.g, common.expansion.g, introns, common.expansion, exons)
        outliers_case.distance <- getsplicedistance(introns.g, outliers_case.g, introns, outliers_case, exons)
        outliers_control_mixed.distance <- getsplicedistance(introns.g, outliers_control_mixed.g, introns, outliers_control_mixed, exons)

  #create rare column
        outliers_case.distance$rarity <- "outliers_case"
        rare.distance$rarity <- "outliers_1_case"
        common.distance$rarity <- "known STRs"
        detected.distance$rarity <- "all expansions"
        outliers_control_mixed.distance$rarity <- "outliers_control_mixed"

  #DO STATISTICAL TEST
        wilcox.test(abs(rare.distance$distance),  abs(detected.distance$distance), alternative = "less")$p.value
        wilcox.test(abs(rare.distance$distance),  abs(common.distance$distance), alternative = "less")$p.value
        wilcox.test(abs(rare.distance$distance), abs(outliers_case.distance$distance), alternative = "less")$p.value
        wilcox.test(abs(rare.distance$distance), abs(outliers_control_mixed.distance$distance), alternative = "less")$p.value

 #?GRAPH FULL
  #BIND THEM ALL TOGETHER
        distance <- rbind( 
                    rare.distance,
                    outliers_case.distance,
                    common.distance,
                    outliers_control_mixed.distance,
                    detected.distance
                    )


  #PLOT
        junction <- ggplot(distance, aes(x = distance, y = ..density.., color = rarity)) +
        geom_density(alpha = .15, adjust = 1/4) +
        geom_vline(xintercept = c(0), lty = 2) + xlab("Distance from splice junction") + theme_classic() +
        theme(panel.border = element_rect(fill = NA), 
                axis.title.x = element_text(size = 14),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14), 
                axis.title.y = element_text(size = 14),
                legend.position = "bottom") +  
        geom_hline(yintercept = 0, linewidth = 1, color = "white") + 
        scale_x_continuous(breaks = c(-5000, 0, 5000)) +
        scale_y_continuous(breaks = c(0, 0.0002)) +
        scale_color_manual(values = c("#4393C3", "#338833", "#D6604D", "#FF7F00", "#984EA3", "#6A3D9A", "#CC0033" ))

 #?GRAPH LESS
  #BIND A SUBGROUP TOGETHER
        distance <- rbind( 
                    rare.distance,
                    outliers_case.distance,
                    common.distance,
                    outliers_control_mixed.distance,
                    detected.distance
                    )

  #PLOT
        junction_less<- ggplot(distance, aes(x = distance, y = ..density.., color = rarity)) +
        geom_density(alpha = .15, adjust = 1/4) +
        geom_vline(xintercept = c(0), lty = 2) + xlab("Distance from splice junction") + theme_classic() +
        theme(panel.border = element_rect(fill = NA), 
                axis.title.x = element_text(size = 14),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14), 
                axis.title.y = element_text(size = 14),
        legend.position = "bottom") +  
        geom_hline(yintercept = 0, size = 1, color = "white") + 
        scale_x_continuous(breaks = c(-5000, 0, 5000)) +
        scale_y_continuous(breaks = c(0, 0.0002)) +
        scale_color_manual(values = c("#4393C3", "#338833", "#D6604D", "#FF7F00","#6A3D9A", "#CC0033" ))


 #? SAVE GRAPHS
        ggsave("junction_full_plot.png", junction, width = 10, height = 6, dpi = 300)
        ggsave("junction_less_plot.png", junction_less, width = 10, height = 6, dpi = 300)








#

