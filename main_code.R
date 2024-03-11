library(data.table)
library(MendelianRandomization)
library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)
library(gsmr)
library(MRPRESSO)
library(tidyverse)
library(vcfR)
library(dplyr)
library(cause)
library(stringr)
library(ggplot2)
library(ggsci)

options(scipen = 200)
options(stringsAsFactors = F)

###functioin
harmonise_ld_dat <- function(x, ld){
  snpnames <- do.call(rbind, strsplit(rownames(ld), split="_"))
  rownames(snpnames) <- snpnames[,1]
  rownames(ld) <- colnames(ld) <- snpnames[,1]
  snp <- intersect(snpnames[,1], x$SNP)
  snpnames <- snpnames[snp,]
  ld <- ld[snp,snp]
  rownames(x) <- x$SNP
  x <- x[snp,]
  stopifnot(all(snpnames[,1]== x$SNP))
  x$effect_allele.exposure <- as.character(x$effect_allele.exposure)
  x$other_allele.exposure <- as.character(x$other_allele.exposure)
  # Set1 x and ld alleles match
  snpnames <- data.frame(snpnames, stringsAsFactors=FALSE)
  snpnames <- merge(subset(x, select=c(SNP, effect_allele.exposure, other_allele.exposure)), snpnames, by.x="SNP", by.y="X1")
  snpnames <- snpnames[match(x$SNP, snpnames$SNP),]
  snpnames$keep <- (snpnames$X2 == snpnames$effect_allele.exposure & snpnames$X3 == snpnames$other_allele.exposure) |
    (snpnames$X3 == snpnames$effect_allele.exposure & snpnames$X2 == snpnames$other_allele.exposure)
  
  # What happens if everything is gone?
  if(nrow(x) == 0)
  {
    message(" - none of the SNPs could be aligned to the LD reference panel")
    return(NULL)
  }
  
  if(any(!snpnames$keep))
  {
    message(" - the following SNPs could not be aligned to the LD reference panel: \n- ", paste(subset(snpnames, !keep)$SNP, collapse="\n - "))
  }
  
  
  snpnames$flip1 <- snpnames$X2 != snpnames$effect_allele.exposure
  x <- subset(x, SNP %in% snpnames$SNP)
  temp1 <- x$effect_allele.exposure[snpnames$flip1]
  temp2 <- x$other_allele.exposure[snpnames$flip1]
  x$beta.exposure[snpnames$flip1] <- x$beta.exposure[snpnames$flip1] * -1
  x$beta.outcome[snpnames$flip1] <- x$beta.outcome[snpnames$flip1] * -1
  x$effect_allele.exposure[snpnames$flip1] <- temp2
  x$other_allele.exposure[snpnames$flip1] <- temp1
  
  rownames(ld) <- snpnames$SNP
  colnames(ld) <- snpnames$SNP
  
  if(any(!snpnames$keep))
  {
    message("Removing ", sum(!snpnames$keep), " variants due to harmonisation issues")
    ld <- ld[snpnames$keep, snpnames$keep]
    x <- x[snpnames$keep, ]
    
  }
  return(list(x=x, ld=ld))
}

format_dat_to_MRInput <- function (dat, get_correlations = FALSE,plink_file,bfile_dir) {
  out <- plyr::dlply(dat, c("exposure", "outcome"), function(x) {
    x <- plyr::mutate(x)
    message("Converting:")
    message(" - exposure: ", x$exposure[1])
    message(" - outcome: ", x$outcome[1])
    if (get_correlations) {
      message(" - obtaining LD matrix")
      if (length(x$SNP)>500) {
        ld <- ld_matrix(x$SNP, plink_bin = plink_file, bfile = bfile_dir)
      }else{
        ld <- ld_matrix(x$SNP, pop = "EUR")
      }
      # print(ld)
      out <- harmonise_ld_dat(x, ld)
      if (is.null(out)) {
        return(NULL)
      }
      x <- out$x
      ld <- out$ld
      MendelianRandomization::mr_input(bx = x$beta.exposure, 
                                       bxse = x$se.exposure, by = x$beta.outcome, byse = x$se.outcome, 
                                       exposure = x$exposure[1], outcome = x$outcome[1], 
                                       snps = x$SNP, effect_allele = x$effect_allele.exposure, 
                                       other_allele = x$other_allele.exposure, eaf = x$eaf.exposure, 
                                       correlation = ld)
    }
    else {
      MendelianRandomization::mr_input(bx = x$beta.exposure, 
                                       bxse = x$se.exposure, by = x$beta.outcome, byse = x$se.outcome, 
                                       exposure = x$exposure[1], outcome = x$outcome[1], 
                                       snps = x$SNP, effect_allele = x$effect_allele.exposure, 
                                       other_allele = x$other_allele.exposure, eaf = x$eaf.exposure)
    }
  })
  return(out)
}

mr_forest_plot <- function(singlesnp_results, exponentiate=FALSE)
{
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
    levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
    am <- grep("All", d$SNP, value=TRUE)
    d$up <- d$b + 1.96 * d$se
    d$lo <- d$b - 1.96 * d$se
    d$tot <- 0.01
    d$tot[d$SNP %in% am] <- 1
    d$SNP <- as.character(d$SNP)
    nom <- d$SNP[! d$SNP %in% am]
    nom <- nom[order(d$b)]
    d <- rbind(d, d[nrow(d),])
    d$SNP[nrow(d)-1] <- ""
    d$b[nrow(d)-1] <- NA
    d$up[nrow(d)-1] <- NA
    d$lo[nrow(d)-1] <- NA
    d$SNP <- ordered(d$SNP, levels=c(am, "", nom))
    
    xint <- 0
    if(exponentiate)
    {
      d$b <- exp(d$b)
      d$up <- exp(d$up)
      d$lo <- exp(d$lo)
      xint <- 1
    }
    
    ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
      ggplot2::geom_vline(xintercept=xint, linetype="dashed") +
      # ggplot2::geom_errorbarh(ggplot2::aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="gray",linetype = "dashed") +
      # ggsci::scale_color_nejm() +
      ggplot2::scale_colour_manual(values=c("#374E55FF", "#B24745FF")) +
      ggplot2::scale_size_manual(values=c(0.3, 1)) +
      # xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
      ggplot2:: theme_bw() + theme(
        legend.position="none", 
        axis.text.y=ggplot2::element_text(size=8), 
        axis.ticks.y=ggplot2::element_line(size=0),
        axis.title.x=ggplot2::element_text(size=14)) +
      ggplot2::labs(y="", x=paste0("MR effect size for ", d$exposure[1], " on ", d$outcome[1], ""))
  })
  res
}

mr_leaveoneout_plot <- function(leaveoneout_results)
{
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    # Need to have at least 3 SNPs because IVW etc methods can't be performed with fewer than 2 SNPs
    if(sum(!grepl("All", d$SNP)) < 3) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    d$up <- d$b + 1.96 * d$se
    d$lo <- d$b - 1.96 * d$se
    d$tot <- 1
    d$tot[d$SNP != "All"] <- 0.01
    d$SNP <- as.character(d$SNP)
    nom <- d$SNP[d$SNP != "All"]
    nom <- nom[order(d$b)]
    d <- rbind(d, d[nrow(d),])
    d$SNP[nrow(d)-1] <- ""
    d$b[nrow(d)-1] <- NA
    d$up[nrow(d)-1] <- NA
    d$lo[nrow(d)-1] <- NA
    d$SNP <- ordered(d$SNP, levels=c("All", "", nom))
    
    ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
      ggplot2::geom_vline(xintercept=0, linetype="dashed") +
      # ggplot2::geom_errorbarh(ggplot2::aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="gray",linetype = "dashed") +
      ggplot2::scale_colour_manual(values=c("#374E55FF", "#B24745FF")) +
      ggplot2::scale_size_manual(values=c(0.3, 1)) +
      # xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
      ggplot2::theme_bw() + theme(
        legend.position="none", 
        axis.text.y=ggplot2::element_text(size=12,family = "serif"),
        # axis.text.y=ggplot2::element_text(size=6,family = "serif"),
        axis.ticks.y=ggplot2::element_line(size=0),
        axis.title.x=ggplot2::element_text(size=16)) +
      ggplot2::labs(y="", x=paste0("MR leave-one-out sensitivity analysis for\n", d$exposure[1], " on ", d$outcome[1], ""))
  })
  res
}

mr_funnel_plot <- function(singlesnp_results)
{
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    am <- grep("All", d$SNP, value=TRUE)
    d$SNP <- gsub("All - ", "", d$SNP)
    am <- gsub("All - ", "", am)
    ggplot2::ggplot(subset(d, ! SNP %in% am), ggplot2::aes(y = 1/se, x=b)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(data = subset(d, SNP %in% am),  ggplot2::aes(xintercept=b, color = SNP, linetype = SNP), key_glyph = "smooth") +
      # ggsci::scale_color_nejm() +
      ggplot2::scale_colour_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF",
                                              "#6F99ADFF","#FFDC91FF","#EE4C97FF","black")) +
      ggplot2::labs(y=expression(paste("Inverse of the standard error of each IV","(",1/SE[IV],")",sep="")),
                    x=expression(paste("Effect size of each IV","(",b[IV],")",sep=""))) +
      ggplot2::theme_bw() + theme(#legend.justification=c(1,1), legend.position=c(0.98,0.98), 
        legend.position="top",
        legend.direction="horizontal",
        legend.box = "horizontal",
        legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        # axis.text.y=ggplot2::element_text(size=12),
        # axis.text.x=ggplot2::element_text(size=12),
        axis.text.y=ggplot2::element_text(size=12),
        axis.text.x=ggplot2::element_text(size=12), 
        axis.title.x=ggplot2::element_text(size=16),
        axis.title.y=ggplot2::element_text(size=16),
        #legend.direction="vertical",
        # legend.background = element_rect(color = "gray"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
      ggplot2::guides(colour=ggplot2::guide_legend(nrow=1)) +
      ggplot2::labs(linetype ="MR Method", color = "MR Method")
  })
  res
}

MR_cor_fun <- function(inputds,presso = FALSE){
  
  if (nrow(inputds)==1) {
    b = inputds$beta.outcome/inputds$beta.exposure
    se = inputds$se.outcome/abs(inputds$beta.exposure)
    pval = pnorm(abs(b)/se, lower.tail = FALSE) * 2
    rlt_wald <- data.frame(method = "Wald ratio", b, se, pval, nsnp = 1)
    
    rlt <- rlt_wald
    rlt$exposure <- inputds$exposure
    rlt$outcome <- inputds$outcome
    return(rlt)
  }else{
    mrds <- format_dat_to_MRInput(dat = inputds, get_correlation=TRUE, 
                                  plink_file = "/Users/wangxiang/Library/CloudStorage/OneDrive-共享的库-Onedrive/ref_code/plink_mac/plink",
                                  bfile_dir = "/Users/wangxiang/Library/CloudStorage/OneDrive-共享的库-Onedrive/我的项目/MR/1kg/EUR")
    
    res_ivw <- try(MendelianRandomization::mr_ivw(mrds[[1]], correl=TRUE))
    rlt_ivw <- data.frame(method = "MR-IVW")
    if (nrow(res_ivw@Correlation)==2) {
      if (class(res_ivw) != 'try-error') {
        rlt_ivw <- data.frame(method = "MR-IVW", b = res_ivw@Estimate, 
                              se = res_ivw@StdError, pval = res_ivw@Pvalue, 
                              nsnp = res_ivw@SNPs)
        rlt <- rlt_ivw
        rlt$exposure <- mrds[[1]]@exposure
        rlt$outcome <- mrds[[1]]@outcome
        return(rlt)
      }else{
        rlt <- data.frame(method = "NA", b = NA, se = NA, pval =NA, nsnp = NA, exposure = mrds[[1]]@exposure, outcome = mrds[[1]]@outcome)
      }
    }else{
      
      if (class(res_ivw) != 'try-error') {
        rlt_ivw <- data.frame(method = "MR-IVW", b = res_ivw@Estimate, 
                              se = res_ivw@StdError, pval = res_ivw@Pvalue, 
                              nsnp = res_ivw@SNPs)
      } 
      
      res_egger <- try(MendelianRandomization::mr_egger(mrds[[1]], correl=TRUE))
      rlt_egger <- data.frame(method = "MR-Egger", b = NA, se = NA, pval =NA, nsnp = NA)
      if (class(res_egger) != 'try-error') {
        rlt_egger <- data.frame(method = c("MR-Egger","MR-Egger Intercept"),
                                b = c(res_egger@Estimate,res_egger@Intercept),
                                se = c(res_egger@StdError.Est,res_egger@StdError.Int),
                                pval = c(res_egger@Pvalue.Est,res_egger@Pvalue.Int), 
                                nsnp = res_egger@SNPs)
      }
      
      res_median <- try(MendelianRandomization::mr_median(mrds[[1]]))
      rlt_median <- data.frame(method = "MR-median", b = NA, se = NA, pval =NA, nsnp = NA)
      if (class(res_median) != 'try-error') {
        rlt_median <- data.frame(method = "MR-median", b = res_median@Estimate, 
                                 se = res_median@StdError, pval = res_median@Pvalue, 
                                 nsnp = res_median@SNPs)
      }
      
      res_mbe <- try(MendelianRandomization::mr_mbe(mrds[[1]]))
      rlt_mbe <- data.frame(method = "Mode-based method", b = NA, se = NA, pval =NA, nsnp = NA)
      if (class(res_mbe) != 'try-error') {
        rlt_mbe <- data.frame(method = "Mode-based method", b = res_mbe@Estimate,
                              se = res_mbe@StdError, pval = res_mbe@Pvalue,
                              nsnp = res_mbe@SNPs)
      }
      
      res_maxlik <- try(MendelianRandomization::mr_maxlik(mrds[[1]]))
      rlt_maxlik <- data.frame(method = "Maximum-likelihood method", 
                               b = NA, se = NA, pval =NA, nsnp = NA)
      if (class(res_maxlik) != 'try-error') {
        rlt_maxlik <- data.frame(method = "Maximum-likelihood method", b = res_maxlik@Estimate, 
                                 se = res_maxlik@StdError, pval = res_maxlik@Pvalue, 
                                 nsnp = res_maxlik@SNPs)
      }
      
      # res_raps <- try(TwoSampleMR::mr(inputds, method_list = "mr_raps"))
      
      res_raps<-try(mr.raps(inputds$beta.exposure, inputds$beta.outcome, inputds$se.exposure, inputds$se.outcome))
      
      rlt_raps <- data.frame(method = "Robust adjusted profile score (RAPS)",
                             b = NA, se = NA, pval =NA, nsnp = NA)
      if (class(res_raps) != 'try-error') {
        
        # rlt_raps <- res_raps[,c(5,7,8,9,6)]
        
        rlt_raps<-data.frame(method = "Robust adjusted profile score (RAPS)", b = res_raps$beta.hat, 
                             se = res_raps$beta.se, pval = res_raps$beta.p.value, 
                             nsnp = nrow(inputds))
      }
      
      res_radial <- try(TwoSampleMR::mr(inputds, method_list = "mr_ivw_radial"))
      rlt_radial <- data.frame(method = "IVW Radial", b = NA, se = NA, pval =NA, nsnp = NA)
      if (class(res_radial) != 'try-error') {
        rlt_radial <- data.frame(method = "IVW Radial", b = res_radial$b,
                                 se = res_radial$se, pval = res_radial$pval,
                                 nsnp = res_radial$nsnp)
      }
      
      
      inputds_p <- inputds[which(inputds$SNP %in% colnames(mrds[[1]]@correlation)),]
      
      res_gsmr <- try(gsmr::gsmr(bzx = inputds_p$beta.exposure, bzx_se = inputds_p$se.exposure, bzx_pval = inputds_p$pval.exposure,
                                 bzy = inputds_p$beta.outcome, bzy_se = inputds_p$se.outcome, bzy_pval = inputds_p$pval.outcome,
                                 ldrho = mrds[[1]]@correlation[inputds_p$SNP,inputds_p$SNP], snpid = inputds_p$SNP, heidi_outlier_flag = T, n_ref = 504,
                                 gwas_thresh=5e-8, nsnps_thresh=1, ld_r2_thresh=1, ld_fdr_thresh=1))
      rlt_gsmr <- data.frame(method = "GSMR", b = NA, se = NA, pval =NA, nsnp = NA)
      if (class(res_gsmr) != 'try-error') {
        rlt_gsmr <- data.frame(method = "GSMR", b = res_gsmr$bxy,
                               se = res_gsmr$bxy_se, pval = res_gsmr$bxy_pval,
                               nsnp = length(res_gsmr$used_index))
      }
      if (presso == TRUE) {
        res_presso <- try(MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",BetaExposure = "beta.exposure", 
                                              SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                                              data = inputds, OUTLIERtest = T,
                                              DISTORTIONtest = T,
                                              SignifThreshold = 0.05,
                                              seed = 12345,
                                              NbDistribution = nrow(inputds)/0.05))
        rlt_presso <- data.frame(method = "MR-PRESSO", b = NA, se = NA, pval =NA, nsnp = NA)
        if (class(res_presso) != 'try-error') {
          presso_rlt <- data.frame(res_presso$`Main MR results`)
          if (is.na(presso_rlt[2,3])) {
            rlt_presso <- data.frame(method = "MR-PRESSO", b = presso_rlt[1,3],
                                     se = presso_rlt[1,4], pval = presso_rlt[1,6],
                                     nsnp = nrow(inputds))
          }else{
            outlier_nsnp <- length(res_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
            nsnp <- nrow(inputds)-outlier_nsnp
            rlt_presso <- data.frame(method = "MR-PRESSO", b = presso_rlt[2,3],
                                     se = presso_rlt[2,4], pval = presso_rlt[2,6],
                                     nsnp = nsnp)
          }
        }
        rlt <- rbind(rlt_egger,rlt_ivw,rlt_median,rlt_mbe, 
                     rlt_maxlik,rlt_raps,rlt_radial,rlt_gsmr,rlt_presso)
      }else{
        rlt <- rbind(rlt_egger,rlt_ivw,rlt_median,rlt_mbe, rlt_maxlik,rlt_raps,rlt_radial,rlt_gsmr)
      }
      rlt$exposure <- mrds[[1]]@exposure
      rlt$outcome <- mrds[[1]]@outcome
      return(rlt)
      
    }
  }
}

mr_scatter_plot <- function(mr_results, dat, xlab, ylab)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0,) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE, key_glyph = "smooth") +
      ggsci::scale_color_nejm() +
      # scale_x_continuous(limits = c(min(xmin),max(xmax)),n.breaks = 5) +
      # scale_y_continuous(limits =c(min(ymin),max(ymax)), n.breaks = 5) +
      ggplot2::scale_colour_manual(values = c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF",
                                              "#6F99ADFF","#FFDC91FF","#EE4C97FF","black")) +
      ggplot2::labs(colour="MR Method", x = xlab, y = ylab) +
      ggplot2::theme_bw() + theme(#legend.justification=c(0,0), legend.position=c(0.02,0.02), 
        legend.position="top",
        # legend.direction="horizontal",
        # legend.box = "horizontal",
        # legend.background = element_rect(color = "gray"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=20),
        axis.text.y=ggplot2::element_text(size=16),
        axis.text.x=ggplot2::element_text(size=16), 
        axis.title.x=ggplot2::element_text(size=20),
        axis.title.y=ggplot2::element_text(size=20),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) + 
      ggplot2::guides(colour=ggplot2::guide_legend(nrow=2)) 
  })
  mrres
}

plot_func <- function(rltds,harmds,exp_trait,out_trait){
  
  rltds$b <- as.numeric(rltds$b)
  
  rltds$id.exposure <- rep(unique(harmds$id.exposure),nrow(rltds))
  
  rltds$id.outcome <- rep(unique(harmds$id.out),nrow(rltds))
  
  rltds$method[2:4] <- c("MR-Median","MR-RAPS","MR-Radial")
  
  plotds <- rltds[rltds$method!="MR-Egger Intercept",]
  
  plotds$a <- ifelse(plotds$method=="MR-Egger",rltds$b[rltds$method=="MR-Egger Intercept"],0)
  
  # plotds <- plotds[plotds$method!="MR-Egger",]
  p1 <- mr_scatter_plot(plotds, harmds, xlab = paste0("SNPs effect on ",exp_trait),
                        ylab = paste0("SNPs effect on ",out_trait))
  
  
  ggsave(p1[[1]], file=paste0("./result_data/figure/",exp_trait,"_",out_trait,"_MR.pdf"), width=10, height=10)
  
  res_single <- mr_singlesnp(harmds,all_method = c("mr_ivw"))
  
  res_single_mr <- plotds[,c(1:4,6:9)]
  
  res_single_mr$method <- paste("All - ",res_single_mr$method,sep = "")
  
  colnames(res_single_mr)[c(1,4)] <- c("SNP","p")
  
  res_single <- rbind(res_single[-nrow(res_single),c(1:4,6:9)],res_single_mr)
  
  p2 <- mr_funnel_plot(res_single)
  
  ggsave(p2[[1]], file=paste0("./result_data/figure/",exp_trait,"_",out_trait,"-funnel",".pdf"), width=12, height=9)
  
  p3 <- mr_forest_plot(res_single)
  
  ggsave(p3[[1]], file=paste0("./result_data/figure/",exp_trait,"_",out_trait,"-forest",".pdf"), width=9, height=12)
  
  res_loo <- mr_leaveoneout(harmds)
  
  p4 <- mr_leaveoneout_plot(res_loo)
  
  ggsave(p4[[1]], file=paste0("./result_data/figure/",exp_trait,"_",out_trait,"-LOO",".pdf"), width=9, height=12)
  
}

setwd("/Users/wangxiang/Library/CloudStorage/OneDrive-个人/我的项目/兼职/银屑病与多发性硬化")

###data source
#L12_PSORIASIS Psoriasis XII Diseases of the skin and subcutaneous tissue (L12_) 8075 330975 
#https://storage.googleapis.com/finngen-public-data-r8/summary_stats/finngen_R8_L12_PSORIASIS.gz

PSORIASIS <- fread("raw_data/finngen_R8_L12_PSORIASIS.gz")

PSORIASIS$chr_pos <- paste0(PSORIASIS$`#chrom`,":",PSORIASIS$pos)

PSORIASIS_nona <- PSORIASIS[-which(PSORIASIS$rsids==""),]

###hg38 to hg19 by liftover
lift = data.frame(v1 = paste0("chr",PSORIASIS_nona$`#chrom`),v2 = PSORIASIS_nona$pos,v3 = PSORIASIS_nona$pos+1,v4 = PSORIASIS_nona$chr_pos)

write.table(lift,"anal_data/PSORIASIS.ucsc.bed",row.names = F,col.names = F,quote=F,sep= "\t")

lift_rename <- fread("anal_data/PSORIASIS_rename.ucsc.bed")

lift_rename$hg19_chr_pos <- paste0(lift_rename$V1,":",lift_rename$V2)
  
lift_rename$hg19_chr_pos <- str_replace(lift_rename$hg19_chr_pos,"chr","")

lift_rename_nodup <-lift_rename[-which(duplicated(lift_rename$hg19_chr_pos)),]

PSORIASIS_nodup <- PSORIASIS_nona[-which(duplicated(PSORIASIS_nona$chr_pos)),]

PSORIASIS_h19 <- left_join(PSORIASIS_nodup,lift_rename_nodup[,c("V4","hg19_chr_pos")],by=c("chr_pos"="V4"))

PSORIASIS_h19_nona <- PSORIASIS_h19[-which(is.na(PSORIASIS_h19$hg19_chr_pos)),]

save(PSORIASIS_h19_nona, file="anal_data/PSORIASIS_h19_nona.RData")

load("anal_data/PSORIASIS_h19_nona.RData")


###data source
#https://gwas.mrcieu.ac.uk/datasets/ieu-b-18/
Sclerosis <- vcfR::read.vcfR("raw_data/ieu-b-18.vcf.gz")

Sclerosis_meta <- data.frame(Sclerosis@meta)

Sclerosis_fix <- data.frame(Sclerosis@fix)

Sclerosis_gt <- data.frame(Sclerosis@gt)

Sclerosis_fix$ES <- sapply(str_split(Sclerosis_gt$ieu.b.18, ":"),'[',1)

Sclerosis_fix$SE <- sapply(str_split(Sclerosis_gt$ieu.b.18, ":"),'[',2)

Sclerosis_fix$LP <- sapply(str_split(Sclerosis_gt$ieu.b.18, ":"),'[',3)

Sclerosis_fix$LP <- as.numeric(Sclerosis_fix$LP)

Sclerosis_fix$P <- 10^(-Sclerosis_fix$LP)

Sclerosis_fix$chr_pos <- paste0(Sclerosis_fix$CHROM,":",Sclerosis_fix$POS)

Sclerosis_fix <- Sclerosis_fix[,c("chr_pos","ID","CHROM","POS","ALT","REF","ES","SE","P")]

colnames(Sclerosis_fix) <- c("chr_pos","rs_number","chr","bp","A1","A2","beta","se","pval")

save(Sclerosis_fix, file="anal_data/Sclerosis_fix.RData")


load("anal_data/Sclerosis_fix.RData")

Sclerosis_sig <- Sclerosis_fix[Sclerosis_fix$pval<=5e-8,]

###select snp
Sclerosis_PSORIASIS_sig <- Sclerosis_sig[Sclerosis_sig$chr_pos %in% PSORIASIS_h19_nona$hg19_chr_pos,]

###LD clumping
Sclerosis_PSORIASIS_sig$pheno <- "Sclerosis"

Sclerosis_PSORIASIS_exp <- format_data(Sclerosis_PSORIASIS_sig,type = "exposure",snp_col = "rs_number", phenotype_col = "pheno",
                           beta_col = "beta",se_col = "se", eaf_col = "MAF",
                           effect_allele_col = "A1",other_allele_col = "A2", 
                           pval_col = "pval",samplesize_col = "N",info_col = "chr_pos")

Sclerosis_clump <- clump_data(Sclerosis_PSORIASIS_exp, clump_kb = 10000, clump_r2 = 0.1, clump_p1 = 1, clump_p2 = 1, pop = "EUR")

PSORIASIS_clump <- merge(PSORIASIS_h19_nona,Sclerosis_clump[,c("info.exposure","SNP")],by.x = "hg19_chr_pos",by.y="info.exposure")

PSORIASIS_clump$pheno <- "PSORIASIS"

PSORIASIS_clump <- data.frame(PSORIASIS_clump)

PSORIASIS_out <- format_data(PSORIASIS_clump,type = "outcome",snp_col = "rsids", phenotype_col = "pheno",
                        beta_col = "beta",se_col = "sebeta", eaf_col = "af_alt",
                        effect_allele_col = "alt",other_allele_col = "ref", 
                        pval_col = "pval",samplesize_col = "N")

###In order to perform MR the effect of a SNP on an outcome and exposure must be harmonised to be relative to the same allele
Sclerosis_PSORIASIS_harm <- harmonise_data(Sclerosis_clump, PSORIASIS_out, action = 2)

Sclerosis_PSORIASIS_harm <- Sclerosis_PSORIASIS_harm[Sclerosis_PSORIASIS_harm$mr_keep,]

##MR result1 for Heterogeneity test using Cochran's Q-statistic to identify outliers
res_radial <- RadialMR::ivw_radial(Sclerosis_PSORIASIS_harm)

if(length(res_radial$outliers)>1){
  
  Sclerosis_PSORIASIS_harm_nopleio <- Sclerosis_PSORIASIS_harm[!Sclerosis_PSORIASIS_harm$SNP%in%res_radial$outliers$SNP,]
  
}else if(res_radial$outliers == "No significant outliers") {
  
  Sclerosis_PSORIASIS_harm_nopleio <- Sclerosis_PSORIASIS_harm
  
}else{
  
  Sclerosis_PSORIASIS_harm_nopleio <- Sclerosis_PSORIASIS_harm[!Sclerosis_PSORIASIS_harm$SNP%in%res_radial$outliers$SNP,]
  
}

mr_pleiotropy_test(Sclerosis_PSORIASIS_harm_nopleio)

##MR result1
Sclerosis_PSORIASIS_rlt <- MR_cor_fun(Sclerosis_PSORIASIS_harm_nopleio,presso = T)

Sclerosis_PSORIASIS_rlt_s <- rbind(Sclerosis_PSORIASIS_rlt,
                                 c("CAUSE",-0.07,0.03316327,0.045,64,"Sclerosis","PSORIASIS"))

Sclerosis_PSORIASIS_rlt_s$b <- as.numeric(Sclerosis_PSORIASIS_rlt_s$b)

Sclerosis_PSORIASIS_rlt_s$se <- as.numeric(Sclerosis_PSORIASIS_rlt_s$se)

#######
Sclerosis_PSORIASIS_harm_nopleio2 <- Sclerosis_PSORIASIS_harm_nopleio[-which(Sclerosis_PSORIASIS_harm_nopleio$SNP %in% c("rs1112718", "rs4676756", "rs9275602")),]

Sclerosis_PSORIASIS_rlt2 <- MR_cor_fun(Sclerosis_PSORIASIS_harm_nopleio2,presso = T)

Sclerosis_PSORIASIS_rlt2_s <- Sclerosis_PSORIASIS_rlt2[-which(Sclerosis_PSORIASIS_rlt2$method %in% c("Mode-based method","GSMR","MR-Egger Intercept")),]

Sclerosis_PSORIASIS_rlt2_s <- rbind(Sclerosis_PSORIASIS_rlt2_s,
                                   c("CAUSE",-0.07,0.03316327,0.045,61,"Sclerosis","PSORIASIS"))

Sclerosis_PSORIASIS_rlt2_s$b <- as.numeric(Sclerosis_PSORIASIS_rlt2_s$b)

Sclerosis_PSORIASIS_rlt2_s$se <- as.numeric(Sclerosis_PSORIASIS_rlt2_s$se)

mr_res <- Sclerosis_PSORIASIS_rlt2_s

mr_res$upper <- sprintf("%.2f",exp(mr_res$b+1.96*mr_res$se))

mr_res$lower <- sprintf("%.2f",exp(mr_res$b-1.96*mr_res$se))

mr_res$OR <- sprintf("%.2f",exp(mr_res$b))

mr_res$ci <- paste0(mr_res$OR, " (", mr_res$lower, ", ", mr_res$upper,")")

mr_res <- mr_res[c(8,1,2,3,4,5,6,7),c("method","ci","pval")]

write.csv(mr_res,"result_data/Sclerosis_PSORIASIS_mr_res2.csv")

save(Sclerosis_PSORIASIS_harm_nopleio,Sclerosis_PSORIASIS_rlt_s, Sclerosis_PSORIASIS_harm_nopleio2,Sclerosis_PSORIASIS_rlt2_s, file="result_data/Sclerosis_PSORIASIS_MR.RData")

###select snp
PSORIASIS_sig <- PSORIASIS_h19_nona[which(PSORIASIS_h19_nona$pval<=5e-8),]

PSORIASIS_Sclerosis_sig <- PSORIASIS_sig[PSORIASIS_sig$hg19_chr_pos %in% Sclerosis_fix$chr_pos,]

##LD clumping
PSORIASIS_Sclerosis_sig$pheno <- "PSORIASIS"

PSORIASIS_Sclerosis_exp <- format_data(PSORIASIS_Sclerosis_sig,type = "exposure",snp_col = "rsids", phenotype_col = "pheno",
                                       beta_col = "beta",se_col = "sebeta", eaf_col = "af_alt",
                                       effect_allele_col = "alt",other_allele_col = "ref", 
                                       pval_col = "pval",samplesize_col = "N",info_col = "hg19_chr_pos")

PSORIASIS_clump <- clump_data(PSORIASIS_Sclerosis_exp, clump_kb = 10000, clump_r2 = 0.1, clump_p1 = 1, clump_p2 = 1, pop = "EUR")

Sclerosis_clump <- merge(Sclerosis_fix,PSORIASIS_clump[,c("info.exposure","SNP")],by.x = "chr_pos",by.y="info.exposure")

Sclerosis_clump$pheno <- "Sclerosis"

Sclerosis_out <- format_data(Sclerosis_clump,type = "outcome",snp_col = "rs_number", phenotype_col = "pheno",
                             beta_col = "beta",se_col = "se", eaf_col = "MAF",
                             effect_allele_col = "A1",other_allele_col = "A2", 
                             pval_col = "pval",samplesize_col = "N")

###In order to perform MR the effect of a SNP on an outcome and exposure must be harmonised to be relative to the same allele
PSORIASIS_Sclerosis_harm <- harmonise_data(PSORIASIS_clump, Sclerosis_out, action = 2)

PSORIASIS_Sclerosis_harm <- PSORIASIS_Sclerosis_harm[PSORIASIS_Sclerosis_harm$mr_keep,]

###MR result1 for Heterogeneity test using Cochran's Q-statistic to identify outliers
res_radial <- RadialMR::ivw_radial(PSORIASIS_Sclerosis_harm)

if(length(res_radial$outliers)>1){
  
  PSORIASIS_Sclerosis_harm_nopleio <- PSORIASIS_Sclerosis_harm[!PSORIASIS_Sclerosis_harm$SNP%in%res_radial$outliers$SNP,]
  
}else if(res_radial$outliers == "No significant outliers") {
  
  PSORIASIS_Sclerosis_harm_nopleio <- PSORIASIS_Sclerosis_harm
  
}else{
  
  PSORIASIS_Sclerosis_harm_nopleio <- PSORIASIS_Sclerosis_harm[!PSORIASIS_Sclerosis_harm$SNP%in%res_radial$outliers$SNP,]
  
}

mr_pleiotropy_test(PSORIASIS_Sclerosis_harm_nopleio)#多效性检验

###MR result1
PSORIASIS_Sclerosis_rlt <- MR_cor_fun(PSORIASIS_Sclerosis_harm_nopleio,presso = T)

PSORIASIS_Sclerosis_rlt_s <- rbind(PSORIASIS_Sclerosis_rlt,
                                   c("CAUSE",-0.33,0.1530612,0.0014,29,"PSORIASIS","Sclerosis"))

PSORIASIS_Sclerosis_rlt_s$b <- as.numeric(PSORIASIS_Sclerosis_rlt_s$b)

PSORIASIS_Sclerosis_rlt_s$se <- as.numeric(PSORIASIS_Sclerosis_rlt_s$se)

###
PSORIASIS_Sclerosis_harm_nopleio2 <- PSORIASIS_Sclerosis_harm_nopleio[-which(PSORIASIS_Sclerosis_harm_nopleio$SNP %in% c("rs12663590", "rs2395471")),]

PSORIASIS_Sclerosis_rlt2 <- MR_cor_fun(PSORIASIS_Sclerosis_harm_nopleio2,presso = T)

PSORIASIS_Sclerosis_rlt2_s <- PSORIASIS_Sclerosis_rlt2[-which(PSORIASIS_Sclerosis_rlt2$method %in% c("Mode-based method","GSMR","MR-Egger Intercept")),]

PSORIASIS_Sclerosis_rlt2_s <- rbind(PSORIASIS_Sclerosis_rlt2_s,
                                    c("CAUSE",-0.30,0.1632653,0.006,27,"Sclerosis","PSORIASIS"))

PSORIASIS_Sclerosis_rlt2_s$b <- as.numeric(PSORIASIS_Sclerosis_rlt2_s$b)

PSORIASIS_Sclerosis_rlt2_s$se <- as.numeric(PSORIASIS_Sclerosis_rlt2_s$se)

mr_res <- PSORIASIS_Sclerosis_rlt2_s

mr_res$upper <- sprintf("%.2f",exp(mr_res$b+1.96*mr_res$se))

mr_res$lower <- sprintf("%.2f",exp(mr_res$b-1.96*mr_res$se))

mr_res$OR <- sprintf("%.2f",exp(mr_res$b))

mr_res$ci <- paste0(mr_res$OR, " (", mr_res$lower, ", ", mr_res$upper,")")

mr_res <- mr_res[c(8,1,2,3,4,5,6,7),c("method","ci","pval")]

write.csv(mr_res,"result_data/PSORIASIS_Sclerosis_mr_res2.csv")

save(PSORIASIS_Sclerosis_harm_nopleio, PSORIASIS_Sclerosis_rlt_s, 
     PSORIASIS_Sclerosis_harm_nopleio2, PSORIASIS_Sclerosis_rlt2_s, file="result_data/PSORIASIS_Sclerosis_MR.RData")


###plot
load("result_data/Sclerosis_PSORIASIS_MR.RData")

# Sclerosis_PSORIASIS_rlt2_s <- Sclerosis_PSORIASIS_rlt2_s[-which(Sclerosis_PSORIASIS_rlt2_s$method %in% c("Mode-based method","GSMR","MR-Egger Intercept")),]

Sclerosis_PSORIASIS_rlt2_s <- Sclerosis_PSORIASIS_rlt2_s[-which(Sclerosis_PSORIASIS_rlt2_s$method %in% c("MR-Egger","Maximum-likelihood method")),]

plot_func(rltds = Sclerosis_PSORIASIS_rlt2_s, harmds = Sclerosis_PSORIASIS_harm_nopleio, exp_trait = "Sclerosis", out_trait = "PSORIASIS")

mr_heterogeneity(Sclerosis_PSORIASIS_harm_nopleio)

mr_pleiotropy_test(Sclerosis_PSORIASIS_harm_nopleio)

mr_egger(mr_input(bx = Sclerosis_PSORIASIS_harm_nopleio$beta.exposure, bxse = Sclerosis_PSORIASIS_harm_nopleio$se.exposure, 
                  by = Sclerosis_PSORIASIS_harm_nopleio$beta.outcome, byse = Sclerosis_PSORIASIS_harm_nopleio$se.outcome))

mr_res <- Sclerosis_PSORIASIS_rlt2_s

mr_res$upper <- sprintf("%.2f",exp(mr_res$b+1.96*mr_res$se))

mr_res$lower <- sprintf("%.2f",exp(mr_res$b-1.96*mr_res$se))

mr_res$OR <- sprintf("%.2f",exp(mr_res$b))

mr_res$ci <- paste0(mr_res$OR, " (", mr_res$lower, ", ", mr_res$upper,")")

mr_res <- mr_res[c(6,2,1,3,4,5),c("method","ci","pval")]

Sclerosis_PSORIASIS_SNP <- Sclerosis_PSORIASIS_harm_nopleio[,c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                                                               "beta.exposure", "se.exposure", "pval.exposure")]

colnames(Sclerosis_PSORIASIS_SNP) <- c("SNP", "Alt", "Ref","Beta", "SE", "Pval")

Sclerosis_PSORIASIS_SNP$F<-(Sclerosis_PSORIASIS_SNP$Beta^2)/(Sclerosis_PSORIASIS_SNP$SE^2)

write.csv(Sclerosis_PSORIASIS_SNP,"result_data/Sclerosis_PSORIASIS_SNP.csv",row.names = F)

####plot
load("result_data/PSORIASIS_Sclerosis_MR.RData")

PSORIASIS_Sclerosis_rlt2_s <- rbind(PSORIASIS_Sclerosis_rlt2_s[-8,], PSORIASIS_Sclerosis_rlt_s[11,])

PSORIASIS_Sclerosis_rlt2_s <- PSORIASIS_Sclerosis_rlt2_s[-which(PSORIASIS_Sclerosis_rlt2_s$method %in% c("MR-Egger","Maximum-likelihood method")),]

plot_func(rltds = PSORIASIS_Sclerosis_rlt2_s, harmds = PSORIASIS_Sclerosis_harm_nopleio, exp_trait = "PSORIASIS", out_trait = "Sclerosis")

mr_heterogeneity(PSORIASIS_Sclerosis_harm_nopleio)

mr_pleiotropy_test(PSORIASIS_Sclerosis_harm_nopleio)

mr_res <- PSORIASIS_Sclerosis_rlt2_s

mr_res$upper <- sprintf("%.2f",exp(mr_res$b+1.96*mr_res$se))

mr_res$lower <- sprintf("%.2f",exp(mr_res$b-1.96*mr_res$se))

mr_res$OR <- sprintf("%.2f",exp(mr_res$b))

mr_res$ci <- paste0(mr_res$OR, " (", mr_res$lower, ", ", mr_res$upper,")")

mr_res <- mr_res[c(6,2,1,3,4,5),c("method","ci","pval")]

mr_egger(mr_input(bx = PSORIASIS_Sclerosis_harm_nopleio$beta.exposure, bxse = PSORIASIS_Sclerosis_harm_nopleio$se.exposure, 
                  by = PSORIASIS_Sclerosis_harm_nopleio$beta.outcome, byse = PSORIASIS_Sclerosis_harm_nopleio$se.outcome))

PSORIASIS_Sclerosis_SNP <- PSORIASIS_Sclerosis_harm_nopleio[,c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                                                               "beta.exposure", "se.exposure", "pval.exposure")]

colnames(PSORIASIS_Sclerosis_SNP) <- c("SNP", "Alt", "Ref","Beta", "SE", "Pval")

PSORIASIS_Sclerosis_SNP$F<-(PSORIASIS_Sclerosis_SNP$Beta^2)/(PSORIASIS_Sclerosis_SNP$SE^2)

write.csv(PSORIASIS_Sclerosis_SNP,"result_data/PSORIASIS_Sclerosis_SNP.csv",row.names = F)


library(cause)

####Sclerosis_PSORIASIS
#Step 1: Format Data for CAUSE

Sclerosis_PSORIASIS_all <- inner_join(Sclerosis_fix[,c("chr_pos","rs_number","A1","A2","beta","se","pval")],
                                      PSORIASIS_h19_nona[,c("hg19_chr_pos","rsids","alt","ref","beta","sebeta","pval")],by=c("chr_pos"="hg19_chr_pos"))

save(Sclerosis_PSORIASIS_all, file="anal_data/Sclerosis_PSORIASIS_all.RData")

load("anal_data/Sclerosis_PSORIASIS_all.RData")

colnames(Sclerosis_PSORIASIS_all) <- c("chr_pos", "rs_number", "A1", "A2", "beta_hat_1", "seb1", "p1", "rsids", "alt", "ref", "beta_hat_2",   
                                       "seb2", "p2")

Sclerosis_PSORIASIS_all_keep <- Sclerosis_PSORIASIS_all[-which(Sclerosis_PSORIASIS_all$A1!=Sclerosis_PSORIASIS_all$alt | Sclerosis_PSORIASIS_all$A2!=Sclerosis_PSORIASIS_all$ref),]

Sclerosis_PSORIASIS_all_pan <- Sclerosis_PSORIASIS_all[which(Sclerosis_PSORIASIS_all$A2==Sclerosis_PSORIASIS_all$alt & Sclerosis_PSORIASIS_all$A1==Sclerosis_PSORIASIS_all$ref),]

Sclerosis_PSORIASIS_all_pan$beta_hat_2 <- (-Sclerosis_PSORIASIS_all_pan$beta_hat_2)

Sclerosis_PSORIASIS_all2 <- rbind(Sclerosis_PSORIASIS_all_keep, Sclerosis_PSORIASIS_all_pan)

Sclerosis_PSORIASIS_all2 <- Sclerosis_PSORIASIS_all2[,c("rs_number","beta_hat_1","seb1", "p1","beta_hat_2","seb2", "p2")]

colnames(Sclerosis_PSORIASIS_all2) <- c("snp","beta_hat_1", "seb1", "p1","beta_hat_2","seb2", "p2")

Sclerosis_PSORIASIS_all2$beta_hat_1 <- as.numeric(Sclerosis_PSORIASIS_all2$beta_hat_1)

Sclerosis_PSORIASIS_all2$seb1 <- as.numeric(Sclerosis_PSORIASIS_all2$seb1)

#Step 2: Calculate nuisance parameters
set.seed(100)

varlist <- with(Sclerosis_PSORIASIS_all2, sample(snp, size=1000000, replace=FALSE))

params <- est_cause_params(Sclerosis_PSORIASIS_all2, varlist)

res <- cause(X=Sclerosis_PSORIASIS_all2, variants = Sclerosis_PSORIASIS_harm_nopleio$SNP, param_ests = params)

res$elpd

summary(res, ci_size=0.95)

res2 <- cause(X=Sclerosis_PSORIASIS_all2, variants = Sclerosis_PSORIASIS_harm_nopleio2$SNP, param_ests = params)

res$elpd

summary(res, ci_size=0.95)

save(res, res2, file="result_data/Sclerosis_PSORIASIS_cause.RData")

load("result_data/Sclerosis_PSORIASIS_cause.RData")

####PSORIASIS_Sclerosis
PSORIASIS_Sclerosis_all <- Sclerosis_PSORIASIS_all2[,c("snp","beta_hat_2", "seb2", "p2","beta_hat_1","seb1", "p1")]

colnames(PSORIASIS_Sclerosis_all) <- c("snp","beta_hat_1", "seb1", "p1","beta_hat_2","seb2", "p2")

PSORIASIS_Sclerosis_all$beta_hat_1 <- as.numeric(PSORIASIS_Sclerosis_all$beta_hat_1)

PSORIASIS_Sclerosis_all$beta_hat_1 <- as.numeric(PSORIASIS_Sclerosis_all$beta_hat_1)

PSORIASIS_Sclerosis_all$seb1 <- as.numeric(PSORIASIS_Sclerosis_all$seb1)

#Step 2: Calculate nuisance parameters
set.seed(100)

varlist <- with(PSORIASIS_Sclerosis_all, sample(snp, size=100000, replace=FALSE))

params <- est_cause_params(PSORIASIS_Sclerosis_all, varlist)

res <- cause(X=PSORIASIS_Sclerosis_all, variants = PSORIASIS_Sclerosis_harm_nopleio$SNP, param_ests = params)

res$elpd

summary(res, ci_size=0.95)

res2 <- cause(X=PSORIASIS_Sclerosis_all, variants = PSORIASIS_Sclerosis_harm_nopleio2$SNP, param_ests = params)

res2$elpd

summary(res2, ci_size=0.95)

save(res, res2, file="result_data/PSORIASIS_Sclerosis_cause.RData")