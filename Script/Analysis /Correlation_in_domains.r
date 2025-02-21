#! 10) Correlation in domains 
#TODO: To understand the distribution of the tandem repeat and get insight into possibile funcional importance ,by doing a logistical regression 
#? library 
    library(data.table)
    library(GenomicRanges)
    library(ggforce)
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    library(magrittr)
    library(phastCons100way.UCSC.hg38)
    library(ordinal)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(Biostrings)

#? Prepare data 
    setwd("/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_without_QC_NOUHR/AFTER_DBSCAN/")
    
    gs <- phastCons100way.UCSC.hg38

    bin.dt <- read.delim("/home/rachele/Downloads/bret_stuff/genome.bin.1k.tsv", stringsAsFactors = F)
    bin.dt <- bin.dt[bin.dt$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]
    bin.g <- GRanges(bin.dt$seqnames, IRanges(bin.dt$start, bin.dt$end), "*")

    tmp <- score(gs, bin.g)
    bin.dt$PhastCons <- tmp

    expansion <- read.delim("outliers_1_case_no_split.tsv", stringsAsFactors = F)
    expansion.g <- GRanges(expansion$contig, IRanges(expansion$start, expansion$end), "*")
    olap <- data.frame(findOverlaps(bin.g, expansion.g))
    olap <- aggregate(subjectHits ~ queryHits, olap, length)
    bin.dt$expan <- 0
    bin.dt$expan[olap$queryHits] <- olap$subjectHits

    phylop <- read.delim("/home/rachele/Downloads/bret_stuff/subset_phylop_1kb_bins.tsv", header = T)
    phylop.g <- GRanges(phylop$seqnames, IRanges(phylop$start, phylop$end), "*")
    olap <- data.frame(findOverlaps(bin.g, phylop.g))
    olap$score <- phylop$phylop.max[olap$subjectHits]
    olap <- aggregate(score ~ queryHits, olap, mean)
    bin.dt$phylop <- NA
    bin.dt$phylop[olap$queryHits] <- olap$score


    fragile <- read.delim("/home/rachele/Downloads/fagile_sites_hg38.tsv",  stringsAsFactors = F)
    fragile$V1 <- sub("chry", "chrY", fragile$V1)
    fragile$V1 <- sub("chrx", "chrX", fragile$V1)
    fragile <- fragile[fragile$V2 != "-", ]
    fragile.g <- GRanges(fragile$V1, IRanges(as.numeric(fragile$V2), as.numeric(fragile$V3)), "*")
    olap <- data.frame(findOverlaps(bin.g, fragile.g))
    bin.dt$fragile <- 0
    bin.dt$fragile[unique(olap$queryHits)] <- 1

    known.str <- read.delim("/home/rachele/Downloads/bret_stuff/UCSC_simple_repeats_hg38.period_lte20.txt", stringsAsFactors = F, header = F)
    known.str <- data.frame(known.str) ### 1031708
    known.str.g <- GRanges(known.str$V1, IRanges(known.str$V2, known.str$V3), "*")
    olap <- data.frame(findOverlaps(bin.g, known.str.g))
    bin.dt$str <- 0
    bin.dt$str[unique(olap$queryHits)] <- 1

    bin.dt <- bin.dt[bin.dt$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]
    features <- c("fragile", "GC", "str", "PhastCons", "phylop")
    bin.dt[is.na(bin.dt)] <- 0
    bin.dt$EHdn <- factor(bin.dt$expan > 0)
    bin.dt$knownSTR <- factor(bin.dt$str > 0)
    dt.out <- data.frame()
    p <- list()

#? Logistical regression 
    for(type in c("EHdn", "knownSTR"))
    for(f in features){
    bin.dt$feat <- bin.dt[, f]#scale(bin.dt[, f])#
    
    lm <- glm(sprintf("%s ~ feat", type), data = bin.dt, family = binomial(link = "logit"))
    pvalue <- summary(lm)$coefficients["feat", "Pr(>|z|)"]
    
    conf <- confint.default(lm)
    dt.out <- rbind(dt.out, data.frame("feature" = f, "Odds ratio" = exp(lm$coefficients["feat"]),
                                        "OR.lower" = exp(conf["feat", 1]),
                                        "OR.upper" = exp(conf["feat", 2]),
                                        "type" = type,
                                        "pvalue" = pvalue, stringsAsFactors = F))
    
    p[[length(p) + 1]] <- ggplot(bin.dt, aes(x = feat, y = expan)) + geom_point(shape = 1) + ylab("expansions") + xlab(f) +
        geom_abline(intercept = lm$coefficients[1], slope = lm$coefficients[2], lty = 2, color = "red")
    }

    dt.out$feature <- c("fragile sites", "GC content", "known STRs", 
                        "PhastCons", "phyloP")
    dt.out$feature <- factor(dt.out$feature,
                            levels = c("GC content", "phyloP",  "PhastCons", "fragile sites", "known STRs"))
#? Graph 
# Create the plot with the legend
domains<-ggplot(dt.out[dt.out$feature != "known STRs", ], 
       aes(x = feature, y = Odds.ratio, fill = type)) + 
    geom_bar(stat = "identity", width = .5, position = position_dodge(width = .5), color = "black") +
    geom_errorbar(aes(ymin = OR.lower, ymax = OR.upper), 
                  position = position_dodge(width = .5), linewidth = .4, width = .2) +
    geom_hline(yintercept = 1, lty = 2) +
    theme_classic() + ylab("Odds ratio") + xlab("") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
          axis.text.y = element_text(size = 14), 
          axis.title.y = element_text(size = 14), 
          panel.border = element_rect(fill = NA)) +
    scale_fill_manual(values = c("#D6604D", "#4393C3")) +
    labs(title = "Correlation of functional domains")+
    scale_y_continuous(breaks = c(0, 1)) +
    guides(fill = guide_legend(title = "Type"))

ggsave("correlation_domain.png", domains, width = 10, height = 6, dpi = 300)