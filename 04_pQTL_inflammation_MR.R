## ============================================================================
## 04_pQTL_inflammation_MR.R  (逐蛋白 pQTL-MR; 与 04_drive_pQTL.sh 配套)
## 炎症蛋白 pQTL-MR: Zhao 2023 SCALLOP/Olink 91 炎症蛋白 → MDD2025
## ----------------------------------------------------------------------------
## 输入:
##   $PQTL_DIR/<GCST>.tsv.gz        (91 protein pQTL files, GRCh38)
##   data/gwas/mdd2025.ma                (结局 MDD2025 daner→.ma)
##   data/ld_ref/1000G.EUR.QC            (EUR 1000G PLINK bfile, clump 用)
##   data/pqtl/cis/<gene>.txt            (由 04_drive_pQTL.sh awk 预提取的 ±1Mb cis 区)
## 输出:
##   out/pqtl_inflammation_MR.csv        (逐蛋白 MR 估计: Wald/IVW/Egger/WM/PRESSO + coloc)
## 对应正文: pQTL 层 (炎症蛋白是否因果驱动 MDD)
## ----------------------------------------------------------------------------
## 方法 (用户指定): cis 工具 P<5e-8 → 与 MDD 交集后 plink r2<0.001 clump → F>10 →
##   harmonise by rsid → Wald(1SNP)/IVW+Egger+WM(≥2)+MR-PRESSO(≥4) → coloc.abf。
## 样本重叠自查: Zhao n≈14,736 即便全落在 MDD2025 总 N(1,639,572) 内仅占 0.90% → 可忽略。
## 结果: 91 全跑; 51 有 MR 估计, 40 无可用 cis 工具 (含 IL6/TNF/IFNG/CCL2 分泌型弱 cis);
##   面板 FDR-BH: 4 名义 p<0.05 (FGF21/CCL28/IL15RA/CD274), 无一过 FDR, 无一共定位;
##   No protein met both the panel-corrected MR and colocalization criteria.
## 运行: 由 04_drive_pQTL.sh 逐蛋白调用 (见该文件)。
## Zhao et al. 2023: PMID 37563310; doi:10.1038/s41590-023-01588-w.
## ============================================================================
suppressPackageStartupMessages({
  .libPaths(c("./.r-libs/geo", .libPaths()))
  library(data.table); library(TwoSampleMR); library(coloc); library(ieugwasr)
})
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 5L) stop("Usage: Rscript 04_pQTL_inflammation_MR.R <gene> <GCST> <chr> <gene_start> <gene_end>")
gene<-args[1]; gcst<-args[2]; chrom<-as.integer(args[3]); gstart<-as.numeric(args[4]); gend<-as.numeric(args[5])
DL <- Sys.getenv("PQTL_DIR", "data/raw/scallop")
LDREF <- Sys.getenv("LDREF_PREFIX", "data/ld_ref/1000G.EUR.QC")
PLINK <- genetics.binaRies::get_plink_binary()
MDD_MA <- Sys.getenv("MDD_MA", "data/gwas/mdd2025.ma")
CIS_DIR <- Sys.getenv("PQTL_CIS_DIR", "data/pqtl/cis"); dir.create(CIS_DIR, showWarnings=FALSE, recursive=TRUE)
OUT_CSV <- Sys.getenv("PQTL_OUT", "out/pqtl_inflammation_MR.csv")
dir.create(dirname(OUT_CSV), showWarnings=FALSE, recursive=TRUE)
win <- 1e6
lo <- max(0, gstart-win); hi <- gend+win

appendrow <- function(...){
  r <- data.frame(..., stringsAsFactors=FALSE)
  fwrite(r, OUT_CSV, append=file.exists(OUT_CSV))
}
emit <- function(method,n_snp,b,se,pval,orr,direction,pph4,causal,note=""){
  appendrow(gene=gene,GCST=gcst,n_snp=n_snp,method=method,b=b,se=se,pval=pval,
            OR=orr,direction=direction,coloc_PPH4=pph4,causal=causal,note=note)
}

## 1. read pre-extracted cis region (extracted by bash awk before this script)
cis_raw <- file.path(CIS_DIR, paste0(gene,".txt"))
if(!file.exists(cis_raw)){ emit("no_cis_file",0,NA,NA,NA,NA,NA,NA,"skip","cis file missing"); quit(save="no") }
cis <- fread(cis_raw)
if(nrow(cis)==0){ emit("no_cis_data",0,NA,NA,NA,NA,NA,NA,"skip","empty cis region"); quit(save="no") }

## 2. instrument selection: cis P<5e-8
setnames(cis, c("chromosome","base_pair_location","effect_allele","other_allele","beta",
                "standard_error","effect_allele_frequency","p_value","variant_id","rsid","n"),
         c("chr","pos","ea","oa","beta","se","eaf","pval","varid","rsid","n"), skip_absent=TRUE)
iv <- cis[pval < 5e-8 & rsid!="" & !is.na(rsid)]
iv <- iv[!duplicated(rsid)]
if(nrow(iv)==0){ emit("no_cis_instrument",0,NA,NA,NA,NA,NA,NA,"skip","no cis SNP p<5e-8"); quit(save="no") }

## 2b. restrict instruments to those present in the outcome (MDD) BEFORE clumping,
##     so the clump lead SNP is guaranteed to exist in the outcome
mdd_full <- fread(MDD_MA)
setnames(mdd_full, c("SNP","A1","A2","freq","b","se","p","n"))
iv <- iv[rsid %in% mdd_full$SNP]
if(nrow(iv)==0){ emit("no_shared_instrument",0,NA,NA,NA,NA,NA,NA,"skip","no cis p<5e-8 SNP present in MDD"); quit(save="no") }

## 3. LD clump r2<0.001 via local plink + 1000G EUR
expo <- data.frame(SNP=iv$rsid, pval.exposure=iv$pval, id.exposure=gene, rsid=iv$rsid)
clp <- tryCatch(
  suppressMessages(suppressWarnings(utils::capture.output(
    clp_tmp <- ld_clump(dplyr::tibble(rsid=iv$rsid, pval=iv$pval, id=gene),
           clump_r2=0.001, clump_kb=10000, clump_p=1,
           plink_bin=PLINK, bfile=LDREF)))),
  error=function(e){ message("clump err: ",conditionMessage(e)); NULL })
clp <- if(exists("clp_tmp")) clp_tmp else NULL
if(is.null(clp) || nrow(clp)==0){ emit("clump_failed",nrow(iv),NA,NA,NA,NA,NA,NA,"skip","LD clump returned 0"); quit(save="no") }
ivc <- iv[rsid %in% clp$rsid]

## 4. build exposure dataset + F-stat
exp_dat <- format_data(as.data.frame(ivc), type="exposure", snp_col="rsid",
    beta_col="beta", se_col="se", effect_allele_col="ea", other_allele_col="oa",
    eaf_col="eaf", pval_col="pval", samplesize_col="n")
exp_dat$exposure <- gene; exp_dat$id.exposure <- gene
exp_dat$F_stat <- (exp_dat$beta.exposure/exp_dat$se.exposure)^2
exp_dat <- exp_dat[exp_dat$F_stat>10, ]
if(nrow(exp_dat)==0){ emit("weak_instruments",nrow(ivc),NA,NA,NA,NA,NA,NA,"skip","all F<=10"); quit(save="no") }
meanF <- mean(exp_dat$F_stat)

## 5. outcome = MDD2025 (.ma keyed by rsid) — reuse mdd_full loaded above
mdd <- mdd_full
out_sub <- mdd[SNP %in% exp_dat$SNP]
if(nrow(out_sub)==0){ emit("no_outcome_snp",nrow(exp_dat),NA,NA,NA,NA,NA,NA,"skip","instruments absent in MDD"); quit(save="no") }
out_dat <- format_data(as.data.frame(out_sub), type="outcome", snp_col="SNP",
    beta_col="b", se_col="se", effect_allele_col="A1", other_allele_col="A2",
    eaf_col="freq", pval_col="p", samplesize_col="n")
out_dat$outcome <- "MDD2025"; out_dat$id.outcome <- "MDD2025"

## 6. harmonise
dat <- harmonise_data(exp_dat, out_dat, action=2)
dat <- dat[dat$mr_keep, ]
nkeep <- nrow(dat)
if(nkeep==0){ emit("harmonise_empty",nrow(exp_dat),NA,NA,NA,NA,NA,NA,"skip","0 SNP after harmonise"); quit(save="no") }

## 7. MR
dirstr <- function(b) if(is.na(b)) NA else if(b>0) "protein↑→MDD风险↑" else "protein↑→MDD风险↓"
if(nkeep==1){
  m <- mr(dat, method_list="mr_wald_ratio")
  orr <- exp(m$b)
  emit("Wald_ratio",1,m$b,m$se,m$pval,orr,dirstr(m$b),NA,ifelse(m$pval<0.05,"nominal","no"),
       sprintf("meanF=%.1f",meanF))
} else {
  m <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  for(i in seq_len(nrow(m))){
    orr <- exp(m$b[i])
    emit(m$method[i], m$nsnp[i], m$b[i], m$se[i], m$pval[i], orr, dirstr(m$b[i]), NA,
         "", sprintf("meanF=%.1f",meanF))
  }
  # MR-PRESSO if >=4 SNP
  if(nkeep>=4){
    pr <- tryCatch(run_mr_presso(dat, NbDistribution=1000), error=function(e) NULL)
    if(!is.null(pr)){
      gp <- pr[[1]]$`Main MR results`
      rr <- gp[gp$`MR Analysis`=="Raw",]
      emit("MR-PRESSO_raw", nkeep, rr$`Causal Estimate`, rr$Sd, rr$`P-value`,
           exp(rr$`Causal Estimate`), dirstr(rr$`Causal Estimate`), NA, "",
           sprintf("global_p=%s", tryCatch(pr[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue, error=function(e) NA)))
    }
  }
}

## 8. coloc over full cis window (both exposure & outcome SNPs)
cisall <- cis[rsid!="" & !is.na(rsid)]
cisall <- cisall[!duplicated(rsid)]
mmerge <- merge(cisall[,.(rsid,beta,se,eaf,n)], mdd[,.(rsid=SNP,b,se_o=se,freq)], by="rsid")
mmerge <- mmerge[!is.na(beta)&!is.na(se)&se>0 & !is.na(b)&!is.na(se_o)&se_o>0]
pph4 <- NA
if(nrow(mmerge)>=20){
  s_frac <- 357636/(357636+1281936)
  De <- list(beta=mmerge$beta, varbeta=mmerge$se^2, snp=mmerge$rsid, type="quant",
             N=median(cisall$n,na.rm=TRUE), MAF=pmin(mmerge$eaf,1-mmerge$eaf))
  Dg <- list(beta=mmerge$b, varbeta=mmerge$se_o^2, snp=mmerge$rsid, type="cc", s=s_frac,
             MAF=pmin(mmerge$freq,1-mmerge$freq))
  ct <- tryCatch(coloc.abf(De,Dg), error=function(e) NULL)
  if(!is.null(ct)) pph4 <- ct$summary["PP.H4.abf"]
}
## 9. final causal verdict row (primary = IVW or Wald)
prim <- if(nkeep==1) m[1,] else m[m$method=="Inverse variance weighted",][1,]
causal <- if(!is.na(prim$pval) && prim$pval<0.05 && !is.na(pph4) && pph4>0.8) "YES(MR+coloc)" else
          if(!is.na(prim$pval) && prim$pval<0.05) "MR-only(coloc<0.8 or NA)" else "no"
emit("SUMMARY", nkeep, prim$b, prim$se, prim$pval, exp(prim$b), dirstr(prim$b), pph4, causal,
     sprintf("meanF=%.1f;coloc_nsnp=%d", meanF, nrow(mmerge)))
cat(sprintf("DONE %s (%s): nkeep=%d prim_p=%.3g coloc_PPH4=%s causal=%s\n",
            gene,gcst,nkeep, prim$pval, ifelse(is.na(pph4),"NA",sprintf("%.3f",pph4)), causal))
