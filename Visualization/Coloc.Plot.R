library(eQTpLot)
library(data.table)
library(stringr)

###########################################################################################

combine_data <- function(GWAS.df, QTL.df, Genes.df, gene, trait, sigpvalue_GWAS, sigpvalue_QTL, 
                         rangebp , gbuild){
  
  QTL.data <- QTL.df[which(QTL.df$Gene.Symbol == gene & QTL.df$P.Value < sigpvalue_QTL & !(is.na(QTL.df$NES)) & !(is.na(QTL.df$P.Value))), ]
  if (dim(QTL.data)[1] == 0){ 
    stop("Sorry, QTL.df does not have any data for the gene ", paste(gene), " meeting your sigpvalue_QTL threshold")
  }
  startpos <- min(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == gbuild), ] %>% dplyr::select(Start)) - rangebp
  stoppos <- max(Genes.df[which(Genes.df$Gene == gene & Genes.df$Build == gbuild), ] %>% dplyr::select(Stop)) + rangebp
  chromosome <- Genes.df$CHR[Genes.df$Gene == gene & Genes.df$Build == gbuild]
  if("PHE" %in% names(GWAS.df)){
    GWAS.df <- GWAS.df[GWAS.df$PHE == trait,]
  }
  
  gwas.data <- GWAS.df[which(GWAS.df$CHR == chromosome  & GWAS.df$BP >= startpos & GWAS.df$BP <= stoppos & !(is.na(GWAS.df$P)) & !(is.na(GWAS.df$BETA))), ]
  
  NoFisher <- FALSE
  if (dim(gwas.data[which(gwas.data$P < sigpvalue_GWAS), ])[1] == 0) {
    NoFisher <- TRUE
    print(paste(sep = "", "WARNING: GWAS.df does not contain any SNPs with p-value < sigpvalue_GWAS within the ", 
                range, "kb flanking the gene ", gene, " for the trait ", trait, ". QTL Enrcihment Plot statistics will not be calculated"))
  }
  
  QTL.data <- dplyr::ungroup(QTL.data)
  gwas.data$SNP <- as.factor(gwas.data$SNP)
  QTL.data$SNP.Id <- as.factor(QTL.data$SNP.Id)
  combinedSNPS <- sort(union(levels(gwas.data$SNP), levels(QTL.data$SNP.Id)))
  
  Combined.Data <- dplyr::left_join(dplyr::mutate(gwas.data, SNP = factor(SNP, levels = combinedSNPS)), 
                                    dplyr::mutate(QTL.data, SNP.Id = factor(SNP.Id, levels = combinedSNPS)) %>% 
                                      dplyr::rename(SNP = SNP.Id), by = "SNP")
  Combined.Data$DirectionOfEffect_GWAS <- ifelse(Combined.Data$BETA < 0, "Negative", ifelse(Combined.Data$BETA > 0, "Positive", NA))
  Combined.Data$DirectionOfEffect_QTL <- ifelse(Combined.Data$NES < 0, "DOWN", ifelse(Combined.Data$NES > 0, "UP", NA))
  Combined.Data$Congruence <- (Combined.Data$BETA * Combined.Data$NES)
  Combined.Data$Congruence <- ifelse(Combined.Data$Congruence < 0, "Incongruent-QTL", ifelse(Combined.Data$Congruence > 0, "Congruent-QTL", NA))
  Combined.Data$NeglogQTLpValue <- -(log10(Combined.Data$P.Value))
  Combined.Data$Neglog10pvalue_GWAS <- -(log10(Combined.Data$P))
  Combined.Data <- Combined.Data[which(!(is.na(Combined.Data$P))), ]
  Combined.Data$significance <- ifelse(Combined.Data$P >= sigpvalue_GWAS, "Non-significant", "Significant")
  Combined.Data$Congruence[is.na(Combined.Data$Congruence)] <- "Non-QTL"
  Combined.Data <- Combined.Data %>% dplyr::mutate(Congruence = factor(Congruence, levels = c("Non-QTL","Congruent-QTL", "Incongruent-QTL"), ordered = TRUE))
  Combined.Data$isQTL <- Combined.Data$Congruence
  Combined.Data$isQTL <- ifelse(Combined.Data$isQTL == "Non-QTL", "Non-QTL", "QTL")  
  
  return(list(Combined.Data=Combined.Data,Fisher=NoFisher))
}
##########################################################################################################

plot.coloc <- function(GWAS.df, QTL.df, Genes.df, gene, trait, sigpvalue_GWAS, sigpvalue_QTL, 
                       rangebp , gbuild, congruence, GeneList){
  if(!GeneList){
    if(length(gene) >= 2){
      stop("When you set GeneList to FALSE you can't generate plots for more than one gene! Please set GeneList to TRUE or enter only one gene.")
    }
    temp <- combine_data(GWAS.df, QTL.df, Genes.df, gene, trait, sigpvalue_GWAS, sigpvalue_QTL, 
                         rangebp, gbuild)
    Combined.QTL.GWAS.Data <- temp$Combined.Data
    NoFisher <- temp$Fisher
    
    minpos <- min(Combined.QTL.GWAS.Data$BP, na.rm = TRUE)
    maxpos <- max(Combined.QTL.GWAS.Data$BP, na.rm = TRUE)
    
    if (dim(Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL"), ])[1] == 0) {
      Congruentdata <- FALSE
    }else {
      Congruentdata <- TRUE
    }
    if (dim(Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL"), ])[1] == 0) {
      Incongruentdata <- FALSE
    }else {
      Incongruentdata <- TRUE
    }
    if (gbuild == "hg19") {
      hostname <- "https://grch37.ensembl.org"
    }
    if (gbuild == "hg38") {
      hostname <- "https://apr2020.archive.ensembl.org"
    }
    bm <- biomaRt::useMart(host = hostname, biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    biomTrack <- Gviz::BiomartGeneRegionTrack(genome = gbuild, chromosome = median(Combined.QTL.GWAS.Data$CHR, na.rm = TRUE), 
                                              start = minpos, end = maxpos, filter = list(with_refseq_mrna = TRUE), 
                                              name = "ENSEMBL", background.panel = "gray95", biomart = bm,margin = c(-3, -3))
    gtrack <- Gviz::GenomeAxisTrack(fontcolor = "#000000", fontsize = 14, margin = c(-3, -3))
    itrack <- Gviz::IdeogramTrack(genome = gbuild, chromosome = paste0("chr",Genes.df$CHR[Genes.df$Gene == gene]),margin = c(-3, -3))
    genetracks <- patchwork::wrap_elements(panel = (grid::grid.grabExpr(Gviz::plotTracks(list(itrack,biomTrack, gtrack), collapseTranscripts = "meta", transcriptAnnotation = "symbol", 
                                                                                         chromosome = median(Combined.QTL.GWAS.Data$CHR, na.rm = TRUE), 
                                                                                         from = min(Combined.QTL.GWAS.Data$BP, na.rm = TRUE), to = max(Combined.QTL.GWAS.Data$BP, na.rm = TRUE), 
                                                                                         showTitle = FALSE, labelPos = "below", distFromAxis = 10, innermargin = 0, maxHeight = (2 * 10), 
                                                                                         minHeight = (2 * 10), sizes = c(1,2, 1), margin = c(-3, -3)))))
    
    p1 <- ggplot2::ggplot(data = Combined.QTL.GWAS.Data, aes(stroke = 0)) + 
      ggplot2::coord_cartesian(xlim = c(minpos, maxpos), expand = FALSE) + 
      ggplot2::geom_point(data = subset(Combined.QTL.GWAS.Data, Congruence == "Non-QTL"), shape = 15, color = "black", alpha = 0.2, aes(x = BP, y = Neglog10pvalue_GWAS)) + 
      ggplot2::xlab("") + ggplot2::ylab(bquote(-log10(P[GWAS]))) + 
      ggplot2::ggtitle(paste("GWAS of ", trait, ", colored by QTL data for ", gene, "\n(Significance thresholds: GWAS, ", sigpvalue_GWAS, "; QTL, ", sigpvalue_QTL, ")", sep = "")) + 
      ggplot2::scale_shape_manual("GWAS Direction of Effect", values = c(Negative = 25, Positive = 24), na.value = 22) + 
      ggplot2::guides(alpha = FALSE, size = guide_legend("QTL Normalized Effect Size", override.aes = list(shape = 24, color = "black", fill = "grey"), title.position = "top", order = 2), 
                      shape = guide_legend(title.position = "top", direction = "vertical", order = 1, override.aes = list(size = 3, fill = "grey"))) + 
      ggplot2::theme(legend.direction = "horizontal", legend.key = element_rect(fill = NA, colour = NA, size = 0.25)) + 
      ggplot2::geom_vline(xintercept = Genes.df$Start[Genes.df$Gene == gene], linetype = "twodash", color = "blue", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = -log10(sigpvalue_GWAS), linetype = "solid", color = "red", linewidth = 0.5) +
      ggplot2::geom_text(aes(x=Genes.df$Start[Genes.df$Gene == gene], label=paste0("\n",gene), y=2.5), colour="blue", angle=90) +
      ggplot2::theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
      ggplot2::theme(plot.margin = unit(c(0, 1, -0.8, 0), "cm")) 
    
    if (congruence == TRUE) {
      if (Congruentdata == TRUE) {
        p1 <- p1 + ggplot2::geom_point(data = subset(Combined.QTL.GWAS.Data, Congruence == "Congruent-QTL"), alpha = 1, 
                                       aes(x = BP, y = Neglog10pvalue_GWAS, fill = NeglogQTLpValue, 
                                           alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(NES))) + 
          ggplot2::scale_fill_gradient(bquote(atop(-log10(P[QTL]), paste("Congruous SNPs"))), low = "#000099", 
                                       high = "#33FFFF", guide = guide_colorbar(title.position = "top"), 
                                       limits = c(min(Combined.QTL.GWAS.Data %>% 
                                                        dplyr::select(NeglogQTLpValue), na.rm = TRUE), 
                                                  max(Combined.QTL.GWAS.Data %>% dplyr::select(NeglogQTLpValue), na.rm = TRUE)))
      }
      if (Congruentdata == TRUE & Incongruentdata == TRUE) {
        p1 <- p1 + ggnewscale::new_scale_fill()
      }
      if (Incongruentdata == TRUE) {
        p1 <- p1 + ggplot2::geom_point(data = subset(Combined.QTL.GWAS.Data, Congruence == "Incongruent-QTL"), 
                                       alpha = 1, aes(x = BP,y = Neglog10pvalue_GWAS, fill = NeglogQTLpValue, 
                                                      alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(NES))) + 
          ggplot2::scale_fill_gradient(bquote(atop(-log10(P[QTL]), paste("Incongruous SNPs"))), low = "#990000", 
                                       high = "#FFCC33", guide = guide_colorbar(title.position = "top"), 
                                       limits = c(min(Combined.QTL.GWAS.Data %>% 
                                                        dplyr::select(NeglogQTLpValue), na.rm = TRUE), 
                                                  max(Combined.QTL.GWAS.Data %>% dplyr::select(NeglogQTLpValue),na.rm = TRUE)))
      }
    }else{
      p1 <- p1 + ggplot2::geom_point(data = subset(Combined.QTL.GWAS.Data, Congruence == "Congruent-QTL"), alpha = 1, 
                                     aes(x = BP, y = Neglog10pvalue_GWAS, fill = NeglogQTLpValue, alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(NES))) + 
        ggplot2::scale_fill_viridis_c((bquote(-log10(P[QTL]))), option = "C", guide = guide_colorbar(title.position = "top"), 
                                      limits = c(min(Combined.QTL.GWAS.Data %>% dplyr::select(NeglogQTLpValue), 
                                                     na.rm = TRUE), max(Combined.QTL.GWAS.Data %>% dplyr::select(NeglogQTLpValue), na.rm = TRUE)))
    }
    
    if (nrow(as.table(table(Combined.QTL.GWAS.Data$isQTL, Combined.QTL.GWAS.Data$significance))) < 2 | 
        ncol(as.table(table(Combined.QTL.GWAS.Data$isQTL, Combined.QTL.GWAS.Data$significance))) < 2) {
      NoFisher <- TRUE
      print("Not enough data to compute enrichment significance for Enrichment Plot")
    }
    if (NoFisher == FALSE) {
      fisher <- fisher.test(table(Combined.QTL.GWAS.Data$isQTL, Combined.QTL.GWAS.Data$significance))
      fpvalue <- fisher$p.value
    }
    
    if(congruence){
      p2 <- ggplot2::ggplot(Combined.QTL.GWAS.Data) + 
        ggplot2::aes(x = significance, y = 1, fill = Congruence) + 
        ggplot2::theme(legend.title = element_blank())+
        ggplot2::geom_bar(stat = "identity", position = "fill") + 
        ggplot2::ggtitle(paste("Enrichment of QTLs among\nGWAS-significant SNPs")) + 
        ggplot2::ylab("Proportion of SNPs\nthat are QTLs") + 
        ggplot2::xlab(paste("GWAS significance\n(threshold p <", sigpvalue_GWAS, ")")) + ggplot2::ylim(0, 1.2) + 
        if (NoFisher == FALSE) {
          ggpubr::geom_signif(y_position = c(1.1, 1.1), xmin = c("Non-significant"), xmax = c("Significant"), 
                              annotation = (paste("p =", formatC(fpvalue, format = "e", digits = 2))), tip_length = 0.05)
        }
    }else{
      p2 <- ggplot2::ggplot(Combined.QTL.GWAS.Data) + 
        ggplot2::aes(x = significance, y = 1, fill = isQTL) + 
        ggplot2::theme(legend.title = element_blank())+
        ggplot2::geom_bar(stat = "identity", position = "fill") + 
        ggplot2::ggtitle(paste("Enrichment of QTLs among\nGWAS-significant SNPs")) + 
        ggplot2::ylab("Proportion of SNPs\nthat are QTLs") + 
        ggplot2::xlab(paste("GWAS significance\n(threshold p <", sigpvalue_GWAS, ")")) + ggplot2::ylim(0, 1.2) + 
        if (NoFisher == FALSE) {
          ggpubr::geom_signif(y_position = c(1.1, 1.1), xmin = c("Non-significant"), xmax = c("Significant"), 
                              annotation = (paste("p =", formatC(fpvalue, format = "e", digits = 2))), tip_length = 0.05)
        }
    }             
    if(congruence){
      if(Congruentdata & Incongruentdata){
        
        df1 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL"), ]
        df2 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL"), ]
        if(nrow(df1) >= 3){
          pearson.congruent <- cor.test(df1$NeglogQTLpValue,df1$Neglog10pvalue_GWAS,method = "pearson")
        }else{
          print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in Congruent QTLs"))
          pearson.congruent <- list(estimate = 0,p.value=1)
        }
        if(nrow(df2) >= 3){
          pearson.Incongruent <- cor.test(df2$NeglogQTLpValue,df2$Neglog10pvalue_GWAS,method = "pearson")
        }else{
          print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in Incongruent QTLs"))
          pearson.Incongruent <- list(estimate = 0,p.value=1)
        }
        
        p3 <- ggplot(data =df1 , aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), ) + ggplot2::geom_point(color="#000099") +
          ggplot2::geom_smooth( data = df1, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, colour="Congruent-QTL"), 
                                formula = y ~ x, method = "lm") +
          ggplot2::geom_point(data = df2,aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), color = "#990000") +
          ggplot2::geom_smooth(data = df2, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, colour="Incongruent-QTL"), 
                               formula = y ~ x, method = "lm") +
          scale_colour_manual(name="Direction of effect", values=c("#000099","#990000")) + 
          ggplot2::xlab((bquote(-log10(P[QTL])))) + 
          ggplot2::ylab((bquote(-log10(P[GWAS])))) +
          ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.congruent$estimate, 3), "\np = ", 
                                                             formatC(pearson.congruent$p.value, format = "e", digits = 2)), 
                             color = "#000099", hjust = 1, vjust = 1) + 
          ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.Incongruent$estimate, 3), "\np = ", 
                                                             formatC(pearson.Incongruent$p.value, format = "e", digits = 2)), 
                             color = "#990000", hjust = 1, vjust = 2.5) 
      }
      if(Congruentdata & !Incongruentdata){
        
        df1 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL"), ]
        
        if(nrow(df1) >= 3){
          pearson.congruent <- cor.test(df1$NeglogQTLpValue,df1$Neglog10pvalue_GWAS,method = "pearson")
        }else{
          print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in Congruent QTLs"))
          pearson.congruent <- list(estimate = 0,p.value=1)
        }
        
        p3 <- ggplot(data =df1 , aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), ) + ggplot2::geom_point(color="#000099") +
          ggplot2::geom_smooth( data = df1, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, colour="Congruent-QTL"), 
                                formula = y ~ x, method = "lm") +
          scale_colour_manual(name="Direction of effect", values=c("#000099")) + 
          ggplot2::xlab((bquote(-log10(P[QTL])))) + 
          ggplot2::ylab((bquote(-log10(P[GWAS])))) +
          ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.congruent$estimate, 3), "\np = ", 
                                                             formatC(pearson.congruent$p.value, format = "e", digits = 2)), 
                             color = "#000099", hjust = 1, vjust = 1)
      }
      if(!Congruentdata & Incongruentdata){
        
        df2 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL"), ]
        
        if(nrow(df2) >= 3){
          pearson.Incongruent <- cor.test(df2$NeglogQTLpValue,df2$Neglog10pvalue_GWAS,method = "pearson")
        }else{
          print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in Incongruent QTLs"))
          pearson.Incongruent <- list(estimate = 0,p.value=1)
        }
        
        p3 <- ggplot(data =df2 , aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), ) + ggplot2::geom_point(color="#990000") +
          ggplot2::geom_point(data = df2,aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), color = "#990000") +
          ggplot2::geom_smooth(data = df2, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, colour="Incongruent-QTL"), 
                               formula = y ~ x, method = "lm") +
          scale_colour_manual(name="Direction of effect", values=c("#990000")) + 
          ggplot2::xlab((bquote(-log10(P[QTL])))) + 
          ggplot2::ylab((bquote(-log10(P[GWAS])))) +
          ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.Incongruent$estimate, 3), "\np = ", 
                                                             formatC(pearson.Incongruent$p.value, format = "e", digits = 2)), 
                             color = "#990000", hjust = 1, vjust = 1) 
      }
    }else{
      df <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & 
                                           (Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL" | Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL")), ]
      
      if(nrow(df) >= 3){
        pearson.all <- cor.test(df$NeglogQTLpValue,df$Neglog10pvalue_GWAS, method = "pearson")
      }else{
        print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in all QTLs"))
        pearson.all <- list(estimate = 0,p.value=1)
      }
      
      p3 <- ggplot(data =df , aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue), ) + ggplot2::geom_point(color = "#360f70") +
        ggplot2::geom_smooth(data = df, aes(y = Neglog10pvalue_GWAS, x = NeglogQTLpValue, color = Congruence), method = "lm", formula = (y ~ x), color = "#ffee00") +
        scale_colour_manual(name="Direction of effect", values=c("#ffee00")) + 
        ggplot2::geom_text(x = Inf, y = Inf, label = paste(sep = "", "r = ", round(pearson.all$estimate, 3), "\np = ", formatC(pearson.all$p.value, format = "e", digits = 2)), 
                           color = "#360f70", hjust = 1, vjust = 1) 
    }
    
    p4 <- p1 + genetracks + patchwork::plot_spacer() + (p2 + p3 + patchwork::plot_layout(ncol = 2, widths = c(2, 3))) + 
      patchwork::plot_layout(ncol = 1, height = c(4, 2, 0.1, 2)) + 
      patchwork::plot_annotation(tag_levels = "A", tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12))
    
    p1 <- p1 + genetracks + patchwork::plot_spacer() + patchwork::plot_layout(ncol = 1, height = c(4, 2, 0.1))
    return(list(plot.merdeg = p4, plot.coloc=p1,plot.cor=p3,plot.enrich=p2)) 
  }else{
    cor.result <- as.data.frame(matrix(data = NA,nrow = length(gene), ncol = 7))
    names(cor.result) <- c("gene","congruent.cor","Congruent.Pval","Incongruent.cor","Incongruent.Pval","All.cor","All.Pval")
    for(i in 1:length(gene)){
      Combined.QTL.GWAS.Data <- combine_data(GWAS.df, QTL.df, Genes.df, gene[i], trait, sigpvalue_GWAS, sigpvalue_QTL, 
                                             rangebp, gbuild)$Combined.Data
      cor.result$gene[i] <- gene[i]
      
      df1 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL"), ]
      df2 <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL"), ]
      df <- Combined.QTL.GWAS.Data[which(!is.na(Combined.QTL.GWAS.Data$NES) & 
                                           (Combined.QTL.GWAS.Data$Congruence == "Congruent-QTL" | Combined.QTL.GWAS.Data$Congruence == "Incongruent-QTL")), ]
      if(nrow(df1) >= 3){
        pearson.congruent <- cor.test(df1$NeglogQTLpValue,df1$Neglog10pvalue_GWAS,method = "pearson")
        cor.result$congruent.cor[i] <- pearson.congruent$estimate
        cor.result$Congruent.Pval[i] <- pearson.congruent$p.value
      }else{
        print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in Congruent QTLs"))
      }
      if(nrow(df2) >= 3){
        pearson.Incongruent <- cor.test(df2$NeglogQTLpValue,df2$Neglog10pvalue_GWAS,method = "pearson")
        cor.result$Incongruent.cor[i] <- pearson.Incongruent$estimate
        cor.result$Incongruent.Pval[i] <- pearson.Incongruent$p.value
      }else{
        print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in Incongruent QTLs"))
      }
      if(nrow(df) >= 3){
        pearson.all <- cor.test(df$NeglogQTLpValue,df$Neglog10pvalue_GWAS, method = "pearson")
        cor.result$All.cor[i] <- pearson.all$estimate
        cor.result$All.Pval[i] <- pearson.all$p.value
      }else{
        print(paste("Not enough data to compute pearson correlation of gene",gene[i],"in all QTLs"))
      }
      
    }
    
    return(cor.result)
  }
  
}

#################################################################################
args <- commandArgs(T)

qtl_file=args[1]
cpg_file=args[2]

gwas_file=args[3]
gwas_trait = args[4]

out_prefix=args[5]
cpg.interest = args[6]
pvalue_GWAS  = args[7]
pvalue_QTL = args[8]
range = args[9]
gbuild = args[10]
congruence = args[11]



qtl <- fread(file = qtl_file,stringsAsFactors = F,header = T)

gwas <- fread(gwas_file,header = T,stringsAsFactors = F)

cpg <- fread(cpg_file,header = T,stringsAsFactors = F)
cpg = cpg[cpg$Gene %in% qtl$Gene.Symbol,]

cpg.interest <- str_split_1(cpg.interest, pattern = ",")

if(length(cpg.interest) > 1){
  result1 <- plot.coloc(GWAS.df = gwas[gwas$P<pvalue_GWAS,],QTL.df = qtl,Genes.df = cpg,gene =cpg.interest ,trait = gwas_trait,
                        sigpvalue_GWAS = pvalue_GWAS,sigpvalue_QTL = pvalue_QTL,rangebp = range,gbuild = gbuild,
                        congruence = congruence,GeneList = T)
  write.csv(result1, file = paste0(out_prefix,".",gwas_trait,".coloc.plot.csv"),row.names = F )
  result1 <- result1[!is.na(result1$All.Pval),]
  cpg.interest <- result1$gene[result1$All.Pval==min(result1$All.Pval)]
  
  result2 <- plot.coloc(GWAS.df = gwas,QTL.df = qtl,Genes.df = cpg,gene =cpg.interest ,trait = gwas_trait,
                        sigpvalue_GWAS = pvalue_GWAS,sigpvalue_QTL = pvalue_QTL,rangebp = range,gbuild = gbuild,
                        congruence = congruence,GeneList = F)
  tiff(filename = paste0(out_prefix,".",gwas_trait,".",cpg.interest,".coloc.plot.tif"),units = "in",height = 12,width = 14,res = 500)
  print(result2$plot.coloc)
  graphics.off()
}else{
  result2 <- plot.coloc(GWAS.df = gwas,QTL.df = qtl,Genes.df = cpg,gene =cpg.interest ,trait = gwas_trait,
                        sigpvalue_GWAS = pvalue_GWAS,sigpvalue_QTL = pvalue_QTL,rangebp = range,gbuild = gbuild,
                        congruence = congruence,GeneList = F)
  tiff(filename = paste0(out_prefix,".",gwas_trait,".",cpg.interest,".coloc.plot.tif"),units = "in",height = 12,width = 14,res = 500)
  print(result2$plot.coloc)
  graphics.off()
}

