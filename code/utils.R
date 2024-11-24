library(Matrix)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggsci)
library(viridis)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(uwot)
library(ggrepel)
library(uwot)
library(ggsci)
library(ggseqlogo)
library(lme4)
library(broom)
library(broom.mixed)
library(metafor)

cmem_col_scheme <- function(){
  return(c("olivedrab4", "turquoise4", "aquamarine", "dodgerblue3", "navyblue", "cadetblue1", "black", "orangered3"))
}

convert_to_quantile <- function (vals, step=0.1){
  ptl = quantile(vals[!(is.na(vals))], probs=seq(step, 1, step))
  ret = sapply(vals, function(x) names(ptl)[min(which(ptl>=x))])
  ret = as.numeric(as.character(gsub("%", "", ret)))
  return(ret)
}

get_metabeta <- function(betas, ses){
  we = 1 / (ses)^2
  return(sum(betas * we) / sum(we))
}

get_metase <- function(ses){
  we = 1 / (ses)^2
  return(sqrt(1/sum(we)))
}

get_metap <- function(betas, ses){
  we = 1 / (ses)^2
  metabeta = sum(betas * we) / sum(we)
  metase = sqrt(1/sum(we))
  metaz = metabeta / metase
  metap = 2 * pnorm(-abs(metaz))
  return(metap)
}

my_umap_theme <- function(g, cont = FALSE, show.legend=TRUE, aes=36){
  g = g + geom_point_rast(size=0.00001, show.legend=show.legend, shape=".") + xlab("UMAP1") + ylab("UMAP2")
  g = g + theme_classic(base_size=20) + theme(plot.title.position = "plot", plot.title=element_text(size=100, face="bold"))
  if (cont==FALSE){
    return(g + guides(color = guide_legend(override.aes = list(size=aes))))
  }
  return(g + theme(legend.key.size = unit(1.5, 'cm')))
}

my_cluster_umap <- function(md, cl_name, palette="Dark2", index = NULL, col_scheme=NULL){
  md$cl = md[,which(colnames(md)==cl_name)]
  nc = length(unique(md$cl))
  g = ggplot()
  g = g + geom_point_rast(aes(md$UMAP1, md$UMAP2, color=md$cl), show.legend=FALSE, size=0.001)
  ctds = md %>% group_by(cl) %>% dplyr::summarise(mUMAP1 = median(UMAP1), mUMAP2=median(UMAP2))
  g = g + geom_label_repel(aes(x=ctds$mUMAP1, y=ctds$mUMAP2, color=ctds$cl, label=ctds$cl), show.legend=FALSE, max.overlaps = Inf, size=8)
  if (!(is.null(col_scheme))){
    g = g + scale_color_manual(values=col_scheme)
  } else if (nc<=8){
    g = g + scale_color_brewer(palette=palette)
  } else {
    getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
    g = g + scale_color_manual(values = getPalette(nc))
  }
  if (!(is.null(index))) { g = g + ggtitle(index) }
  return(my_umap_theme(g, show.legend=FALSE))
}

cluster_boxplot <- function(df, clustname, yval){
  df$score = as.numeric(as.character(df[,which(colnames(df)==yval)]))
  df = df[!(is.na(df$score)),]
  df$cluster = df[,which(colnames(df)==clustname)]
  gr = df %>% group_by(cluster) %>% summarise(median_score = median(score))
  gr = gr[order(gr$median_score),]
  df$cluster = factor(df$cluster, levels=unique(gr$cluster))
  g = ggplot(df, aes(cluster, y=score))
  g = g + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title=element_text(hjust=0.5))
  return(g)
}

create_seqlogo <- function(betas, poss, ymin, ymax){
  nposs = as.numeric(as.character(sapply(poss, function(x) strsplit(x, "p")[[1]][2])))
  nposs = sapply(1:length(poss), function(x) ifelse(grepl("cdr3", poss[x]), 103+nposs[x], nposs[x]))
  df = data.frame(poss, nposs)
  nposs = sort(as.numeric(as.character(nposs)))
  aminos <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  mat = matrix(nrow=length(aminos), ncol=length(nposs))
  colnames(mat) = nposs
  rownames(mat) = aminos
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      b = betas[betas$val==aminos[i] & betas$model==df$poss[df$nposs==nposs[j]],]
      if (nrow(b)==1){
        mat[i,j] = b$estimate[1]
      } else if (nrow(b)==0){
        mat[i,j] = 0
      } else {
        print("too many rows")
        print(b)
        stopifnot(FALSE)
      }
    }
  }
  g = ggseqlogo(mat, method='custom', seq_type='aa', col_scheme=amino_color_scheme()) + ylim(c(ymin, ymax))
  g = g + scale_x_continuous(breaks=1:length(nposs), labels=sort(nposs))
  g = g + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title=element_text(hjust=0.5))
  g = g + xlab("IMGT position") + ylab("R")
  return(g)
}

amino_color_scheme <- function(){
  return(make_col_scheme(chars=aminos, groups=c('neutral', 'hydrophobic', 'acidic', 'acidic', 'very hydrophobic', 'neutral', 'basic', 'hydrophobic', 'basic', 'very hydrophobic', 'hydrophobic', 'neutral', 'neutral', 'neutral', 'basic', 'neutral', 'neutral', 'neutral', 'very hydrophobic', 'very hydrophobic'),
                         
                         cols=c('black', 'darkgoldenrod', 'navyblue', 'navyblue', 'darkgoldenrod2', 'black', 'coral2', 'darkgoldenrod', 'coral2', 'darkgoldenrod2', 'darkgoldenrod', 'black', 'black', 'black', 'coral2', 'black', 'black', 'black', 'darkgoldenrod2', 'darkgoldenrod2')))
  
}

label_invariants <- function(md){
  colnames(md) = gsub("TCR_", "", colnames(md))
  md$my_aNKT = (md$v_gene_TRA %in% c("TRAV10*01", "TRAV10*02"))  & md$j_gene_TRA=="TRAJ18*01"
  md$my_abNKT = md$my_aNKT  & md$v_gene_TRB=="TRBV25-1*01"
  md$my_aMAIT = grepl("TRAV1-2", md$v_gene_TRA) & (grepl("TRAJ33", md$j_gene_TRA) | grepl("TRAJ20", md$j_gene_TRA) | grepl("TRAJ12", md$j_gene_TRA))
  md$my_abMAIT = md$my_aMAIT & (grepl("TRBV6-", md$v_gene_TRB) | grepl("TRBV20", md$v_gene_TRB))  
  return(md)
}

label_invariants2 <- function(md){
  md$my_aNKT = md$TCRA_vgene=="TRAV10" & md$TCRA_jgene=="TRAJ18"
  md$my_abNKT = md$my_aNKT  & md$TCRB_vgene=="TRBV25-1"
  md$my_aMAIT = grepl("TRAV1-2", md$TCRA_vgene) & (grepl("TRAJ33", md$TCRA_jgene) | grepl("TRAJ20", md$TCRA_jgene) | grepl("TRAJ12", md$TCRA_jgene))
  md$my_abMAIT = md$my_aMAIT & (grepl("TRBV6-", md$TCRB_vgene) | grepl("TRBV20", md$TCRB_vgene))  
  return(md)
}

process_igor_output <- function(file){
  igor.suo = readRDS(file)
  igor.suo$va_del = gsub("\\(", "", igor.suo$igorA_Deletion_V_gene_Three_prime_prio5_size21)
  igor.suo$va_del = as.numeric(gsub("\\)", "", igor.suo$va_del))
  
  igor.suo$ja_del = gsub("\\(", "", igor.suo$igorA_Deletion_J_gene_Five_prime_prio5_size23)
  igor.suo$ja_del = as.numeric(gsub("\\)", "", igor.suo$ja_del))
  
  igor.suo$a_del = igor.suo$va_del + igor.suo$ja_del
  
  igor.suo$a_ins = gsub("\\(", "", igor.suo$igorA_Insertion_VJ_gene_Undefined_side_prio4_size41)
  igor.suo$a_ins = as.numeric(gsub("\\)", "", igor.suo$a_ins))
  
  igor.suo$a_abschange = igor.suo$a_ins + igor.suo$a_del
  igor.suo$a_change = igor.suo$a_ins - igor.suo$a_del
  
  igor.suo$vb_del = gsub("\\(", "", igor.suo$igorB_Deletion_V_gene_Three_prime_prio5_size21)
  igor.suo$vb_del = as.numeric(gsub("\\)", "", igor.suo$vb_del))
  
  igor.suo$db_del1 = gsub("\\(", "", igor.suo$igorB_Deletion_D_gene_Five_prime_prio5_size21)
  igor.suo$db_del1 = as.numeric(gsub("\\)", "", igor.suo$db_del1))
  
  igor.suo$db_del2 = gsub("\\(", "", igor.suo$igorB_Deletion_D_gene_Three_prime_prio5_size21)
  igor.suo$db_del2 = as.numeric(gsub("\\)", "", igor.suo$db_del2))
  
  igor.suo$jb_del = gsub("\\(", "", igor.suo$igorB_Deletion_J_gene_Five_prime_prio5_size23)
  igor.suo$jb_del = as.numeric(gsub("\\)", "", igor.suo$jb_del))
  
  igor.suo$b_del = igor.suo$vb_del + igor.suo$jb_del + igor.suo$db_del1 + igor.suo$db_del2
  
  igor.suo$b_ins1 = gsub("\\(", "", igor.suo$igorB_Insertion_VD_genes_Undefined_side_prio4_size31)
  igor.suo$b_ins1 = as.numeric(gsub("\\)", "", igor.suo$b_ins1))
  
  igor.suo$b_ins2 = gsub("\\(", "", igor.suo$igorB_Insertion_DJ_gene_Undefined_side_prio2_size31)
  igor.suo$b_ins2 = as.numeric(gsub("\\)", "", igor.suo$b_ins2))
  
  igor.suo$b_ins = igor.suo$b_ins1 + igor.suo$b_ins2
  
  igor.suo$b_abschange = igor.suo$b_ins1 + igor.suo$b_ins2 + igor.suo$b_del
  igor.suo$b_change = (igor.suo$b_ins1 + igor.suo$b_ins2) - igor.suo$b_del
  
  igor.suo$tdt_indep = (igor.suo$a_ins + igor.suo$b_ins)==0
  igor.suo$ab_ins = igor.suo$a_ins + igor.suo$b_ins
  return(igor.suo)
}

get_mode <- function(vec){
  t = sort(table(vec))
  return(names(t)[length(t)])
}

twin_match_binomial_test <- function(twin_cells, res="cl0.5"){
  twin_cells$cl = twin_cells[,which(colnames(twin_cells)==res)]
  ntrials = length(unique(twin_cells$clonotype))
  twin_tcrs = twin_cells %>% group_by(COMBAT_ID, clonotype) %>% dplyr::summarise(ncells = length(unique(barcode_id)), nstates = length(unique(cl)), mode_state = get_mode(cl))
  null_freq = prop.table(table(twin_tcrs$mode_state))
  match_table = twin_tcrs %>% group_by(clonotype) %>% dplyr::summarise(match = length(unique(mode_state))==1)
  nmatches = nrow(match_table[match_table$match==TRUE,])
  prob_success = sum(sapply(null_freq, function(x) x^2))
  p = pbinom(nmatches, ntrials, prob_success, lower.tail=FALSE)
  exp = ntrials * prob_success
  return(list(p=p, nmatches=nmatches, ntrials = ntrials, expmatches = exp))
}
  
get_meta <- function(df){
  ses = df$std.error
  betas = df$estimate
  we = 1 / (ses)^2
  metabeta = sum(betas * we) / sum(we)
  metase = sqrt(1/sum(we))
  metaz = metabeta / metase
  metap = 2 * pnorm(-abs(metaz))
  lOR = exp(metabeta + qnorm(0.025)*metase)
  hOR = exp(metabeta + qnorm(0.975)*metase)
  tt = df$tt[1]
  split = df$split[1]
  return(data.frame(condition = c("meta"), estimate = c(metabeta), std.error = c(metase),
                    p.value = c(metap), OR = c(exp(metabeta)), tt = c(tt), split = c(split), lOR = c(lOR),
                    hOR = c(hOR), meta = c(TRUE), file=c("meta")))
}

aggregate_glmer <- function(dir, pd=FALSE){
  savepath = getwd()
  setwd(dir)
  files = list.files()
  files = files[files!="models_perdonor"]
  files = files[!(grepl("interaction", files))]
  if (pd){
    df = data.frame(matrix(nrow=0, ncol=14))
    colnames(df) = c("term", "estimate", "std.error", "statistic", "p.value",
                     "CV", "level", "dataset", "test", "target", "file", "npos", "nneg", "conv")
  } else {
  df = data.frame(matrix(nrow=0, ncol=11))
  colnames(df) = c("term", "estimate", "std.error", "statistic", "p.value",
                   "CV", "level", "dataset", "test", "target", "file")
  }
  for (i in 1:length(files)){
    if (pd) { 
      load(files[i])
      res = stats
      res$conv = conv
      res$npos = npos
      res$nneg = nneg
    } else {
      res = readRDS(files[i])
    }
    res$effect <- NULL
    res$group <- NULL
    info = strsplit(files[i], "_")[[1]]
    axis = info[1]
    res$CV = strsplit(axis, "\\.")[[1]][3]
    res$level = strsplit(axis, "\\.")[[1]][2]
    res$dataset = strsplit(axis, "\\.")[[1]][1]
    res$condition = info[4]
    res$file = files[i]
    res$split = info[3]
    df = rbind(df, res)
  }
  df = df[df$term=="CV",]
  setwd(savepath)
  return(df)
}

plot_forest <- function(df, CR, color, trainonly=FALSE, rev=FALSE, collapse_test = TRUE, splitCD48=TRUE, min=0.8, max=1.5, rm_labels=FALSE, de_exp = FALSE, base_size=NULL){
  if (splitCD48){
    df = df[df$split!="allcells",]
  } else {
    df = df[df$split=="allcells",]
  }
  df$p.value[df$p.value==0] = 1e-300
  if (rev){
    df$estimate = -df$estimate
  }
  df$tt = ifelse(grepl("train", df$file), "train", "test")
  df$tt = factor(df$tt, levels=c("train", "test"))
  if (collapse_test){
    df$condition = as.character(df$condition)
    df = df[df$tt=="train" | df$condition=="res.RData",]
    df$condition[df$tt=="test"] = "individuals held-out\nfrom training"
  }
  df$OR = exp(df$estimate)
  df$lOR = exp(df$estimate + qnorm(0.025)*df$std.error)
  df$hOR = exp(df$estimate + qnorm(0.975)*df$std.error)
  df$meta = FALSE
  df$condition = as.character(df$condition)
  df$condition = gsub("COMBAT", "Dataset 1", df$condition)
  df$condition = gsub("Ren et al.", "Dataset 2", df$condition)
  df$condition = gsub("none", "none of the above", df$condition)
  df$tt = as.character(df$tt)
  if (trainonly){
    df = df[df$tt=="train",]
  } 
  df$tt = factor(df$tt, levels=c("train", "test"))
  df = df[,c("condition", "estimate", "std.error", "p.value", "OR", "tt", "split", "lOR", "hOR", "meta", "file")]
  if (splitCD48){
    m1 = get_meta(df[df$tt=="train" & df$split=="CD4",])
    m2 = get_meta(df[df$tt=="train" & df$split=="CD8",])
    df = rbind(df, m1, m2)
  } else if (nrow(df)>1) {
    get_meta(df)
    m1 = get_meta(df[df$tt=="train",])
    df = rbind(df, m1)
    if (!trainonly & !collapse_test){
      m2 = get_meta(df[df$tt=="test",])
      df = rbind(df, m2)
    }
  }
  df$condition[df$condition=="meta"] = "meta-analysis"
  if (CR){
    if (!collapse_test){
      df$condition = factor(df$condition, levels=rev(c("meta-analysis", "COVID+; Dataset 1", "Sepsis; Dataset 1", "Influenza; Dataset 1", "none of the above; Dataset 1", "COVID+; Dataset 2", "COVID-; Dataset 2")))
    } else {
      df$condition = factor(df$condition, levels=rev(c("meta-analysis", "COVID+; Dataset 1", "Sepsis; Dataset 1", "Influenza; Dataset 1", "none of the above; Dataset 1", "COVID+; Dataset 2", "COVID-; Dataset 2", "individuals held-out\nfrom training")))
    }
  } else {
    df$condition = factor(df$condition, levels=rev(c("meta-analysis", "COVID+", "Sepsis", "Influenza", "none")))
  }
  df$col = df$meta==TRUE | df$condition=="test"
  print(df[df$col==TRUE,c("condition", "estimate", "std.error", "OR", "lOR", "hOR", "p.value", "file", "split", "tt")])
  xint = 0
  g = ggplot(df, aes(estimate, condition, shape=meta, color=col, size=meta))
  g = g + geom_vline(xintercept=xint, linetype="dashed") + scale_color_manual(values=c("black", color))
  if (!(is.null(base_size))){
    g = g + theme_bw(base_size=base_size) + geom_point(show.legend=FALSE, size=3) + geom_errorbar(aes(xmin=estimate + qnorm(0.025)*std.error, xmax=estimate + qnorm(0.975)*std.error), size=1, show.legend = FALSE)
  } else {
    g = g + theme_bw() + geom_point(show.legend=FALSE) + geom_errorbar(aes(xmin=estimate + qnorm(0.025)*std.error, xmax=estimate + qnorm(0.975)*std.error), size=0.5, show.legend = FALSE)
  }
  g = g + coord_cartesian(xlim = c(min, max)) + scale_shape_manual(values=c(16,18)) + scale_size_manual(values=c(1, 2)) + xlab("beta") + ylab("subset of individuals")
  if (splitCD48){
    if (!trainonly & !collapse_test){
      g = g + facet_grid(rows=vars(tt), cols=vars(split))
    } else {
      g = g + facet_wrap(~split, ncol=1)
    }
  } else if (!trainonly & !collapse_test) {
    g = g + facet_wrap(~tt, ncol=1)
  }
  if (rm_labels){
    g = g + theme(axis.text.y=element_blank()) + ylab("")
  }
  return(list(g=g, df=df))
}

aggregate_tcr_cors_simp <- function(dir){
  savepath = getwd()
  setwd(dir)
  files = list.files()
  filenames = files
  filenames = gsub("TRAcdr1p", "TRAcdr1_p", filenames)
  filenames = gsub("TRAcdr2p", "TRAcdr2_p", filenames)
  filenames = gsub("TRBcdr1p", "TRBcdr1_p", filenames)
  filenames = gsub("TRBcdr2p", "TRBcdr2_p", filenames)
  estimate = numeric(length(files))
  CV = character(length(files))
  p.value = numeric(length(files))
  model = character(length(files))
  val = numeric(length(files))
  cases = numeric(length(files))
  for (i in 1:length(files)){
    load(files[i])
    estimate[i] = test$estimate
    p.value[i] = test$p.value
    info = strsplit(filenames[i], "_")[[1]]
    model[i] = paste0(info[3:4], collapse="_")
    cases[i] = tab["TRUE"]
    val[i] = info[5]
    CV[i] = info[2]
  }
  df = data.frame(estimate, CV, p.value, model, val, cases)
  df$type = sapply(df$model, function(x) strsplit(x, "_")[[1]][1])
  setwd(savepath)
  return(df)
}

twogenecor_repplot <- function(combatfile, renfile, labeled, rev=FALSE){
  cgenes = readRDS(combatfile)
  rgenes = readRDS(renfile)
  if (rev){
    cgenes$R = -cgenes$R
    rgenes$R = -rgenes$R
  }
  colnames(cgenes)[2:ncol(cgenes)] = paste("COMBAT", colnames(cgenes)[2:ncol(cgenes)], sep="_")
  colnames(rgenes)[2:ncol(rgenes)] = paste("Ren", colnames(rgenes)[2:ncol(rgenes)], sep="_")
  cgenes = left_join(cgenes, rgenes)
  cgenes = cgenes[order(-cgenes$COMBAT_R),]
  print(head(cgenes[!(is.na(cgenes$COMBAT_R)),]))
  cgenes$lbl = ""
  cgenes$lbl[cgenes$gene %in% labeled] = as.character(cgenes$gene[cgenes$gene %in% labeled])
  g = ggplot(cgenes, aes(COMBAT_R, Ren_R, label=lbl))
  g = g + geom_smooth(method="lm", color="gray") + geom_point_rast(size=0.001) + theme_bw() + geom_text_repel(max.overlaps=Inf)
  return(g + stat_cor())
}

facet_seq_logos <- function(sequences, min=10, max=18){
  sequences = sequences[!(is.na(sequences))]
  df = data.frame(seq = sequences, len = sapply(sequences, function(x) nchar(x)))
  print(table(df$len))
  lens = as.numeric(as.character(unique(df$len)))
  lens = lens[lens>=min & lens<=max]
  lens = lens[!(is.na(lens))]
  lens = sort(lens)
  print("doing lens:")
  print(lens)
  plots = list()
  for (i in 1:length(lens)){
    plots[[i]] = ggseqlogo(as.character(df$seq[df$len==lens[i]]), seq_type='aa', col_scheme=amino_color_scheme(), method='prob') 
  }
  return(plots)
}

aminos <- function(){
return(c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))
}

amino_color_scheme <- function(){
return(make_col_scheme(chars=aminos(), groups=c('neutral', 'hydrophobic', 'acidic', 'acidic', 'very hydrophobic', 'neutral', 'basic', 'hydrophobic', 'basic', 'very hydrophobic', 'hydrophobic', 'neutral', 'neutral', 'neutral', 'basic', 'neutral', 'neutral', 'neutral', 'very hydrophobic', 'very hydrophobic'),
cols=c('black', 'darkgoldenrod', 'navyblue', 'navyblue', 'darkgoldenrod2', 'black', 'coral2', 'darkgoldenrod', 'coral2', 'darkgoldenrod2', 'darkgoldenrod', 'black', 'black', 'black', 'coral2', 'black', 'black', 'black', 'darkgoldenrod2', 'darkgoldenrod2')))
}

scale_variables <- function(data, mns=NULL, sds=NULL){
  if (is.null(mns)){
    mns = sapply(1:ncol(data), function(x) mean(data[,x], na.rm=TRUE))
    names(mns) = colnames(data)
    sds = sapply(1:ncol(data), function(x) sd(data[,x], na.rm=TRUE))
    names(sds) = colnames(data)
  }
  #mns = mns[!(is.na(mns)) & !(is.na(sds)) & sds!=0]
  #sds = sds[names(sds) %in% names(mns)]
  #data = data[,as.character(names(mns))]
  #print(table(colnames(data)==names(mns)))
  data = sweep(data, 2, mns)
  data = sweep(data, MARGIN=2, FUN="/", sds)
  data[is.na(data)] <- 0
  return(data)
}

add_adjacent_ints <- function(x, prefix){
  print("in adjacent")
  factors = paste("AF", seq(1, 5), sep="")
  factor_grid = expand.grid(factors, factors)
  print(factor_grid)
  for (i in 2:18){
    pos2 = paste(prefix, i, sep="")
    if (!(paste(pos2, "AF1", sep="_") %in% colnames(x))) { next }
    pos1 = paste(prefix, (i-1), sep="")
    for (j in 1:nrow(factor_grid)){
      ind2 = which(colnames(x)==paste(pos2, factor_grid$Var1[j], sep="_"))
      ind1 = which(colnames(x)==paste(pos1, factor_grid$Var2[j], sep="_"))
      new = x[,ind1] * x[,ind2]
      x = cbind(x, new)
      colnames(x)[ncol(x)] = paste(colnames(x)[ind1], colnames(x)[ind2], sep="_by_")
    }
  }
  return(x)
}

get_pos_entropy_grid <- function(sequences, min_len=11, max_len=18){
  sequences = sequences[!(is.na(sequences))]
  seq_len = sapply(sequences, function(x) nchar(x))
  lens = unique(seq_len)
  lens = lens[lens>=min_len & lens<=max_len]
  if (length(lens)==0) { return(NULL) }
  lens = lens[!(is.na(lens))]
  print(lens)
  m = matrix(nrow=length(lens), ncol=max(lens))
  for (i in 1:nrow(m)){
    for (j in 1:ncol(m)){
      if (j>lens[i]) { 
        m[i,j] = NA
      } else {
        m[i,j] = Entropy(table(sapply(sequences[seq_len==lens[i]], function(x) substr(x, j, j))))
      }
    }
  }
  rownames(m) = lens
  colnames(m) = seq(1, ncol(m))
  return(m)
}

get_pos_lk_grid <- function(sequences, min_len=11, max_len=18){
  seq_len = sapply(sequences, function(x) nchar(x))
  lens = unique(seq_len)
  lens = lens[lens>=min_len & lens<=max_len]
  if (length(lens)==0) { return(NULL) }
  lens = lens[!(is.na(lens))]
  lk = matrix(nrow=length(lens), ncol=max(lens))
  for (i in 1:nrow(lk)){
    for (j in 1:ncol(lk)){
      if (j>lens[i]) { 
        lk[i,j] = NA
      } else {
        lk[i,j] = log(length(unique(sapply(sequences[seq_len==lens[i]], function(x) substr(x, j, j)))), base=2)
      }
    }
  }
  rownames(lk) = lens
  colnames(lk) = seq(1, ncol(lk))
  return(lk)
}

get_long_format_entropy <- function(sequences, min_len=11, max_len=18, bg_sequences=NULL){
  sequences = gsub("\\.", "", sequences)
  grid = get_pos_entropy_grid(sequences, min_len, max_len)
  if (!(is.null(bg_sequences))){
    bg_sequences = gsub("\\.", "", bg_sequences)
    lk = get_pos_lk_grid(bg_sequences, min_len, max_len)
    grid = grid/lk
  }
  if (is.null(grid)){
    df = data.frame(matrix(nrow=0, ncol=2))
    colnames(df) = c("pos", "ent")
  } else {
    ent = numeric(nrow(grid)*ncol(grid))
    pos = numeric(nrow(grid)*ncol(grid))
    c = 1
    for (i in 1:nrow(grid)){
      for (j in 1:ncol(grid)){
        if (!(is.na(grid[i,j]))){
          ent[c] = grid[i,j]
          pos[c] = paste(rownames(grid)[i], colnames(grid)[j], sep="_")
          c = c + 1
        }
      }
    }
    df = data.frame(pos, ent)
  }
  return(df)
}

clean_antigen_labels <- function(tp){
  tp$meta = grepl("meta", tp$antigen)
  tp$lab = sapply(tp$antigen, function(x) ifelse(x=="meta-analysis", "meta-analysis", strsplit(x, "_")[[1]][2]))
  tp$lab = sapply(tp$antigen, function(x) ifelse(x=="meta-analysis", "meta-analysis", strsplit(x, "_")[[1]][2]))
  tp$hla = sapply(tp$antigen, function(x) ifelse(x=="meta-analysis", "meta-analysis", strsplit(x, "_")[[1]][1]))
  tp$source = "human"
  tp$source[grepl("HIV", tp$antigen)] = "HIV"
  tp$source[grepl("HPV", tp$antigen)] = "HPV"
  tp$source[grepl("CMV", tp$antigen)] = "CMV"
  tp$source[grepl("EBV", tp$antigen)] = "EBV"
  tp$source[grepl("HTLV", tp$antigen)] = "HTLV"
  tp$source[grepl("Kanamycin.B.dioxygenase", tp$antigen)] = "fungal"
  tp$source[grepl("Influenza", tp$antigen)] = "Influenza"
  tp$source[grepl("meta", tp$antigen)] = "meta"
  tp$lab = paste(tp$source, tp$lab, sep=", ")
  tp$lab[grepl("meta-analysis", tp$lab)] = "meta-analysis"
  ants = unique(tp$lab)
  tp$lab = factor(tp$lab, levels=c(ants[ants!="meta-analysis"], "meta-analysis"))
  return(tp)
}

antigen_forest_plot <- function(df, ylog10=FALSE, xlims = NULL, colorHLA = FALSE, colorself = FALSE, colorsource = FALSE, colorpower = FALSE, meta_hline=FALSE){
  betas = df$b
  antigens = df$antigen
  ses = df$se
  ps = df$p
  rma.res = rma(df$b, sei=df$se, measure="OR", method="ML")
  mb = as.numeric(as.character(rma.res$beta))
  p = as.numeric(as.character(rma.res$pval))
  se = as.numeric(as.character(rma.res$se))
  tp = data.frame(c(antigens, "meta-analysis"), c(betas, mb), c(ses, se), c(ps, p), c(df$power, NA))
  colnames(tp) = c("antigen", "beta", "se", "P", "power")
  tp = clean_antigen_labels(tp)
  pal = brewer.pal(8, "Dark2")
  if (colorHLA){
    g = ggplot(tp, aes(lab, beta, shape=meta, color=hla))
    g = g + scale_color_manual(values=c(pal, "black"))
  } else if (colorself){
    g = ggplot(tp, aes(lab, beta, shape=meta, color=source=="self"))
    g = g + scale_color_manual(values=c(pal, "black"))
  } else if (colorsource){
    g = ggplot(tp, aes(lab, beta, shape=meta, color=source))
    g = g + scale_color_manual(values=c(pal, "black"))
  } else if (colorpower){
    g = ggplot(tp, aes(lab, beta, shape=meta, color=power))
    g = g + scale_color_viridis(option="plasma")
  } else {
    g = ggplot(tp, aes(lab, beta, shape=meta, color=meta))
    g = g + scale_color_manual(values=c(pal, "black"))
  }
  g = g + geom_errorbar(aes(ymin=beta+qnorm(0.05)*se, ymax=beta+qnorm(0.95)*se), show.legend = FALSE, color="black") + geom_point(size=3) + geom_hline(yintercept=0)
  g = g + theme_classic() + scale_shape_manual(values=c(16, 18)) 
  if (!(is.null(xlims))){
    g = g + coord_cartesian(xlim=(c(xlims[1], xlims[2])))
  }
  g = g + theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1), title=element_text(hjust=0.5))
  if (meta_hline){
    g = g + geom_hline(yintercept=tp$beta[tp$meta==TRUE], color="#e7288a", linetype="dashed")
  }
  return(g + xlab("antigen"))
}

antigen_power_analysis_plot <- function(betas, meta_power){
  betas = betas[order(betas$b),]
  tp = data.frame(betas[,c("antigen", "power")])
  tp = rbind(tp, data.frame(antigen="meta-analysis", power=meta_power))
  tp = clean_antigen_labels(tp)
  g = ggplot(tp, aes(lab, power, fill=grepl("meta", antigen)))
  g = g + geom_bar(stat="identity", color="black") + theme_classic()
  g = g + xlab("antigen specificity") + scale_fill_manual(values=c("gray", "hotpink3"))
  g = g + theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1), title=element_text(hjust=0.5))
  return(g)
}

get_bin_stats <- function(data, cols, quantile_bin = NULL, lcat=NULL, randeff=TRUE){
  if (!(is.null(lcat))){
    data = data[data$lcat == lcat,]
  }
  vals = unique(data$bin)
  vals = sort(vals)
  if (is.null(quantile_bin)){
    ref = "0"
  } else {
    ref = vals[length(vals)/2]
  }
  data$bin = factor(data$bin, levels=c(ref, vals[vals!=ref]))
  if (randeff){
    fit = glmer(target ~ bin + (1|Donor), data=data, family="binomial")
  } else {
    fit = glm(target ~ bin, data=data, family="binomial")
  }
  st = tidy(fit)
  st$OR = exp(st$estimate)
  st$lowOR = exp(st$estimate+qnorm(0.025)*st$std.error)
  st$highOR = exp(st$estimate+qnorm(0.975)*st$std.error)
  st$bin = gsub("bin", "", st$term)
  if (!(is.null(lcat))){
    st$lcat = lcat
    st = st[,cols]
    st = rbind(st, c(ref, 1, 1, 1, "ref", lcat))
  } else {
    st = st[,cols]
    st = rbind(st, c(ref, 1, 1, 1, "ref"))
  }
  st$OR = as.numeric(as.character(st$OR))
  st$lowOR = as.numeric(as.character(st$lowOR))
  st$highOR = as.numeric(as.character(st$highOR))
  return(st)
}

get_odds <- function(p){
  return(p/(1-p))
}

plot_TCRscore_validation_bybin <- function(data, cv, lcats, minx=-3, maxx=4, plot_together = FALSE, fate, colors = NULL, xlab = NULL, ylab = NULL, quantile_bin = NULL, ymin=NULL, ymax=NULL, thresh =10, show.legend=TRUE, do_mm=FALSE, rescale=FALSE, barplot_only=FALSE){##}, lcat){
  if (is.null(colors)){
    colors = c("darkred", "cadetblue", "darkgoldenrod3", "pink3")
    key = c("blood", "lymph", "nonlymph", "thymus")
    colors = colors[key %in% lcats]
  } 
  data = data[data$lcat %in% lcats,]
  data$TCRscore = data[,which(colnames(data)==paste("X", cv, sep=""))]
  data$target = data[,which(colnames(data)==paste("target", cv, sep=""))]
  data = data[!(is.na(data$target)),]
  if (rescale){
    for (lc in 1:length(lcats)){
      data$TCRscore[data$lcat==lcats[lc]] = scale(data$TCRscore[data$lcat==lcats[lc]])[,1]
    }
  }
  if (is.null(quantile_bin)){
    data$bin = sapply(data$TCRscore, function(x) round(x))
  } else {
    data$qtl = convert_to_quantile(data$TCRscore, step=quantile_bin)
    map = data %>% group_by(qtl) %>% dplyr::summarise(bin = mean(TCRscore), minbin = min(TCRscore), maxbin = max(TCRscore))
    data = left_join(data, map)
  }
  tmp = data
  if (rescale){
    tmp$TCRscore = scale(tmp$TCRscore)[,1]
  }
  fit = glmer(target ~ TCRscore + (1|Donor), data=tmp, family="binomial")
  st = tidy(fit)
  print(st)
  b = numeric(length(lcats))
  se = numeric(length(lcats))
  for (lc in 1:length(lcats)){
    tmp = data[data$lcat==lcats[lc],]
    fit = glmer(target ~ TCRscore + (1|Donor), data=tmp, family="binomial")
    b[lc] = tidy(fit)$estimate[2]
    se[lc] = tidy(fit)$std.error[2]
  }
  stats = data.frame(b, se, lcat=lcats)
  if (plot_together){
    data$lcat = "foo"
    lcats = c("foo")
  }
  gr = data %>% group_by(bin, lcat) %>% dplyr::summarise(npos = length(target[target==TRUE]), nneg = length(target[target==FALSE]))
  gr$prop = gr$npos/(gr$npos+gr$nneg)
  if (is.null(xlab)){ xlab = "TCR score" }
  if (is.null(ylab)) { ylab = paste0(c("proportion of cells\n in", fate, "state"), collapse=" ") }
  g = ggplot(gr, aes(bin, prop, fill=lcat))
  #g = g + coord_cartesian(ylim=c(0.85,1))
  g = g + geom_bar(stat="identity", position="dodge", show.legend = show.legend, size=4) + theme_bw(base_size=15) + scale_fill_manual(values=colors)##c("darkred", "cadetblue", "darkgoldenrod3"))
  if (!(is.null(quantile_bin))){
    g = g + xlim(c(minx,maxx))
  }
  g1 = g + xlab(xlab) + ylab(ylab) 
  gr = gr[gr$nneg>=thresh,]
  gr = gr[gr$npos>=thresh,]
  data = data[data$bin %in% gr$bin,]
  st = data.frame(matrix(nrow=0, ncol=6))
  colnames(st) = c("bin", "OR", "lowOR", "highOR", "term", "lcat")
  t = table(gr$lcat)
  colors = colors[lcats %in% names(t)[t>1]]
  lcats = lcats[lcats %in% names(t)[t>1]]
  for (lc in 1:length(lcats)){
    res = get_bin_stats(data=data, cols = colnames(st), quantile_bin=quantile_bin, lcat=lcats[lc])
    colnames(res) = colnames(st)
    st = rbind(st, res)
    sub = gr[gr$lcat==lcats[lc],]
    ref = res$bin[res$term=="ref"]
    data$mm = "other"
    data$mm[data$bin==min(data$bin)] = "minbin"
    data$mm[data$bin==max(data$bin)] = "maxbin"
    data$ptl = convert_to_quantile(data$TCRscore, step=0.01)
    data$mmp = "other"
    data$mmp[data$ptl==min(data$ptl)] = "minptl"
    data$mmp[data$ptl==max(data$ptl)] = "maxptl"
    d = data[data$mm %in% c("minbin", "maxbin") & data$lcat==lcats[lc],]
    d$mm = factor(d$mm, levels=c("minbin", "maxbin"))
    tmp = data %>% group_by(mmp) %>% dplyr::summarise(prop.pos = length(target[target==TRUE])/length(target[!(is.na(target))]))
    print(tmp)
    if (do_mm){
      fit = glmer(target ~ mm + (1|Donor), data=d, family="binomial")
      print(tidy(fit))
    }
  }
  
  st = st[!(grepl("Intercept", st$term)),]
  st$bin = as.numeric(as.character(st$bin))
  ylab = paste0(c("OR for", fate, "fate"), collapse=" ") 
  g = ggplot() + geom_hline(yintercept=1, linetype="dashed", color="gray")
  g = g + geom_point(aes(x=st$bin, y=st$OR, color=st$lcat), show.legend = show.legend, size=2) + theme_classic(base_size=15) 
  g = g + geom_errorbar(aes(x=st$bin, y=st$OR, ymin=st$lowOR, ymax=st$highOR, color=st$lcat), show.legend = show.legend)
  g = g + scale_color_manual(values=colors) + xlim(c(minx,maxx))
  if (!(is.null(quantile_bin))){
    map$bin = as.numeric(as.character(map$bin))
    st = left_join(st, map)
    print(st)
    #g = g + geom_errorbar(aes(x=st$bin, y=st$OR, xmin=st$minbin, xmax=st$maxbin, color=st$lcat), show.legend = show.legend)
  } 
  if (cv==4){
    g = g + scale_y_continuous(trans="log10", limits=c(0.7, 1.6), breaks=c(0.8, 1, 1.2, 1.5))
  } else {
    g = g + scale_y_continuous(trans="log10")
  }
  g2 = g + xlab(xlab) + ylab(ylab) 
  if (barplot_only) { 
    return(g1)
  } else {
    return(g1 + g2 + plot_layout(ncol=1))
  }
}

call_antigens <- function(R, thresholds, data){
  bin = readRDS("data/10xG_CD8.3.0.2_aggregated_binarized_matrix_full0927.rds")
  ants = colnames(R)
  counts = bin[,19:62]
  hlas = sapply(ants, function(x) strsplit(x, "_")[[1]][1])
  match_donor_list = list()
  match_donor_list[[1]] = c("donor2")
  match_donor_list[[2]] = c("donor1", "donor2")
  match_donor_list[[3]] = c("donor4")
  match_donor_list[[4]] = c("donor1")
  match_donor_list[[5]] = c("donor3")
  match_donor_list[[6]] = c("donor4")
  match_donor_list[[7]] = c("donor2")
  match_donor_list[[8]] = c("donor1")
  names(match_donor_list) = sort(unique(hlas))
  call = data.frame(matrix(nrow=nrow(R), ncol=ncol(R), 0))
  c = 1
  for (i in 1:length(ants)){
    hla = strsplit(ants[i], "_")[[1]][1]
    peptide = strsplit(ants[i], "_")[[1]][2]
    drs = match_donor_list[[which(names(match_donor_list)==hla)]]
    for (j in 1:length(drs)){
      case = paste(ants[i], drs[j])
      if (case %in% names(thresholds)){
        th = thresholds[which(names(thresholds)==case)]
        call[,i][data$donor==drs[j] & R[,i]>th] = 1
        c = c + 1
      }
    }
  }
  colnames(call) = colnames(R)
  rownames(call) = rownames(R)
  print("number of calls per cell:")
  print(table(rowSums(call)))
  ##for cells that see one antigen, return assignment and dextramer staining
  conf = call[rowSums(call)==1,]
  inds = sapply(1:nrow(conf), function(x) which(conf[x,]==1))
  stain = R[rowSums(call)==1,]
  stain_vals = sapply(1:nrow(stain), function(x) stain[x,inds[x]])
  antigens = sapply(1:length(inds), function(x) colnames(R)[inds[x]])
  df = data.frame(cell = rownames(conf), antigen = antigens, stain = stain_vals)
  return(df)
}

testTCRscore_perantigen <- function(call, md10xg, cell_state="binary", thresh=10, stain_cov=FALSE, stain_only=FALSE, remove_CD4s = FALSE){
  ants_called = unique(call$antigen)
  call = suppressMessages(left_join(call, md10xg[,c("cell", "donor", "AE", "TCR.mem", "cmem")]))
  call = call[!(is.na(call$donor)),]
  call = call[!(is.na(call$cmem)),]
  if (remove_CD4s){
    call = call[!(grepl("CD4", call$cmem)),]
  }
  print(paste("analyzing", paste(nrow(call), "cells")))
  print(paste("analyzing", paste(length(unique(call$antigen)), "antigens prior to N/M thresh")))
  b = vector(mode="numeric")
  se = vector(mode="numeric")
  p = vector(mode="numeric")
  b.stain = vector(mode="numeric")
  se.stain = vector(mode="numeric")
  p.stain = vector(mode="numeric")
  antigen = vector(mode="character")
  n = vector(mode="numeric")
  c = 1
  for (i in 1:length(ants_called)){
    dat = call[call$antigen==ants_called[i],]
    dat = dat[!(is.na(dat$AE)),]
    dat = dat[!(is.na(dat$TCR.mem)),]
    dat$stain = scale(dat$stain)[,1]
    if (nrow(dat[dat$AE,])<1 | nrow(dat[!(dat$AE),])<thresh) { next }
    if (length(unique(dat$donor))>1){
      if (stain_cov){
        fit = glm(AE ~ TCR.mem + stain + donor, data=dat, family="binomial")
      } else if (stain_only){
        fit = glm(AE ~ stain + donor, data=dat, family="binomial")
      } else {
        fit = glm(AE ~ TCR.mem + donor, data=dat, family="binomial")
      }
    } else {
      if (stain_cov){
        fit = glm(AE ~ TCR.mem + stain, data=dat, family="binomial")
      } else if (stain_only){
        fit = glm(AE ~ stain, data=dat, family="binomial")
      } else {
        fit = glm(AE ~ TCR.mem, data=dat, family="binomial")
      }
    }
    st = tidy(fit)
    b[c] = st$estimate[2]
    se[c] = st$std.error[2]
    p[c] = st$p.value[2]
    if (stain_cov){
      b.stain[c] = st$estimate[3]
      se.stain[c] = st$std.error[3]
      p.stain[c] = st$p.value[3]
    }
    n[c] = nrow(dat)
    antigen[c] = ants_called[i]
    c = c +1
  }
  print(paste(length(antigen), "antigens retained"))
  print("positive beta?")
  print(table(b>0))
  print("(nominally significant):")
  print(table(b[p<0.05]>0))
  df = data.frame(antigen, b, se, p, n)
  if (stain_cov){
    df = cbind(df, data.frame(b.stain, se.stain, p.stain))
  }
  df = df[order(df$b),]
  return(df)
}

get_cells_included <- function(clust, covid_status="any", T_subset="allT"){
  if (T_subset=="CD4"){
    clust = clust[clust$X %in% readRDS("combatCD4s_QC1007_TH10505.rds"),]
  }
  if (T_subset=="CD8"){
    clust = clust[clust$X %in% readRDS("combatCD8s_QC1007_TH10505.rds"),]
  }
  if (covid_status=="neg"){
    clust = clust[clust$SARSCoV2PCR==0,]
  }
  if (covid_status=="pos"){
    clust = clust[clust$SARSCoV2PCR==0,]
  }
  return(as.character(clust$X))
}

twin_analysis <- function(md, covid_status="any", T_subset="allT", cl_res="cl0.5"){
  md$cl = as.character(md[,which(colnames(md)==cl_res)])
  cells_included = get_cells_included(md, covid_status, T_subset)
  possible_clusters = unique(as.character(md$cl[md$X %in% cells_included]))
  twin_tcrs = get_twin_TCRs(md[md$barcode_id %in% cells_included,])
  expobs = compute_expobs(twin_tcrs, possible_clusters)
  print(expobs$test$p.value)
  clustdf = expobs$clustdf
  clustdflong = clustdf %>% pivot_longer(cols=c("exp", "obs"), names_to="type", values_to="count")
  getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
  clustdflong$x = ifelse(clustdflong$type=="exp", "expected", "observed")
  g = ggplot(clustdflong, aes(x, count, color=clust, group=clust))
  g = g + geom_point() + geom_line() + theme_bw(base_size=15) + xlab("") + ylab("number of matches")
  lineplot = g + scale_color_manual(values = getPalette(length(unique(clustdflong$clust)))) + labs(color="transcriptional\n cluster")
  return(list(twin_tcrs = twin_tcrs, clustdf = clustdf, lineplot = lineplot))
}

p_via_combinations <- function(table){
  denom = choose(sum(table), 2)
  p = numeric(length(table))
  for (i in 1:length(table)){
    p[i] = choose(table[i], 2)/denom
  }
  names(p) = names(table)
  return(p)
}

compute_expobs <- function(tt, possible_clusters, exc=c("")){
  t = table(tt$clonotype)
  tt = tt[tt$clonotype %in% names(t)[t==2],]
  mt = tt %>% group_by(clonotype) %>% dplyr::summarise(match = length(unique(mode_state))==1)
  mt$clust = sapply(1:nrow(mt), function(x) ifelse(mt$match[x], tt$mode_state[tt$clonotype==mt$clonotype[x]][1], "mixed"))
  
  null_freq  = p_via_combinations(table(tt$mode_state))
  prob_success = sum(null_freq)
  exp_perclust = null_freq*(nrow(tt)/2)
  names(exp_perclust) = names(null_freq)
  clustdf = data.frame(clust = names(exp_perclust), exp = exp_perclust)
  msg = unique(as.character(possible_clusters[!(possible_clusters %in% clustdf$clust)]))
  if (length(msg)>0){
    msg = msg[!(msg %in% exc)]
    for (i in 1:length(msg)){
      clustdf = rbind(clustdf, c(msg[i], 0))
    }
  }
  clustdf$exp = as.numeric(as.character(clustdf$exp))
  clustdf$obs = sapply(clustdf$clust, function(x) nrow(mt[mt$match==TRUE & mt$clust==x,]))
  bt = binom.test(sum(clustdf$obs), nrow(tt)/2, prob_success, alternative="greater")
  print(bt)
  return(list(clustdf=clustdf, test = bt))
}

select_clonotypes_toviz <- function(twin_tcrs){
  match_table = twin_tcrs %>% group_by(clonotype) %>% dplyr::summarise(match = length(unique(mode_state))==1)
  unambig = twin_tcrs[twin_tcrs$ncells==1,]
  gr = unambig %>% group_by(clonotype) %>% dplyr::summarise(ndonor = length(unique(COMBAT_ID)))
  unambig = unambig[unambig$clonotype %in% gr$clonotype[gr$ndonor==2],]
  unambig = unambig[unambig$clonotype %in% match_table$clonotype[match_table$match==TRUE],]
  #7 states we can visualize (including cl8)
  set.seed(27)
  viz_clonotypes = character(7)
  states = unique(unambig$mode_state)
  for (i in 1:length(states)){
    viz_clonotypes[i] = sample(unambig$clonotype[unambig$mode_state==states[i]], 1)
  }
  return(viz_clonotypes)
}

highlight_twins <- function(md){
  md$cl = md$mRNA_cluster
  twin_tcrs = get_twin_TCRs(md)
  viz_clonotypes = select_clonotypes_toviz(twin_tcrs)
  g = ggplot() + ggtitle("d")
  g = g + geom_point_rast(aes(md$UMAP1[!(md$clonotype %in% viz_clonotypes)], md$UMAP2[!(md$clonotype %in% viz_clonotypes)]), size=0.001, color="black")
  g = g + geom_point_rast(aes(md$UMAP1[(md$clonotype %in% viz_clonotypes)], md$UMAP2[(md$clonotype %in% viz_clonotypes)], color=md$clonotype[(md$clonotype %in% viz_clonotypes)]), size=4, show.legend=FALSE)
  g = g + scale_color_brewer(palette = "Set2") + xlab("UMAP1") + ylab("UMAP2")
  g = g + theme_classic(base_size=20) + theme(plot.title.position = "plot", plot.title=element_text(size=30, face="bold"))
  return(g)
}

process_tcr_info <- function(singleAB = TRUE){
  tcr_info = readRDS("data/tcr_chain_information.rds")
  tcr_info$nbeta = 0
  tcr_info$nbeta[tcr_info$chain_composition %in% c("double_alpha_beta", "single_beta", "triple_2_alpha")] = 1
  tcr_info$nbeta[tcr_info$chain_composition %in% c("double_beta", "quad_2_alpha_2_beta", "triple_2_beta")] = 2
  
  tcr_info$nalpha = 0
  tcr_info$nalpha[tcr_info$chain_composition %in% c("double_alpha_beta", "single_alpha", "triple_2_beta")] = 1
  tcr_info$nalpha[tcr_info$chain_composition %in% c("double_alpha", "quad_2_alpha_2_beta", "triple_2_alpha")] = 2
  
  tcr_info = tcr_info[tcr_info$nbeta==1 & tcr_info$nalpha %in% c(1,2),]
  if (singleAB){
    tcr_info = tcr_info[tcr_info$nbeta==1 & tcr_info$nalpha==1,]
  }
  tcr_info$clonotypeA = paste(tcr_info$v_gene_TRA, tcr_info$cdr3_TRA)
  tcr_info$clonotypeA = paste(tcr_info$clonotypeA, tcr_info$j_gene_TRA)
  tcr_info$clonotypeA = paste(tcr_info$clonotypeA, tcr_info$v_gene_TRA2)
  tcr_info$clonotypeA = paste(tcr_info$clonotypeA, tcr_info$cdr3_TRA2)
  tcr_info$clonotypeA = paste(tcr_info$clonotypeA, tcr_info$j_gene_TRA2)
  tcr_info$clonotypeB = paste(tcr_info$v_gene_TRB, tcr_info$cdr3_TRB)
  tcr_info$clonotypeB = paste(tcr_info$clonotypeB, tcr_info$j_gene_TRB)
  tcr_info$clonotype = paste(tcr_info$clonotypeA, tcr_info$clonotypeB, sep="-")
  return(tcr_info)
}

get_twin_TCRs <- function(tcr_info){
  dr_counts = tcr_info %>% group_by(clonotype) %>% dplyr::summarise(nindiv = length(unique(COMBAT_ID)))
  
  twins = tcr_info[tcr_info$clonotype %in% dr_counts$clonotype[dr_counts$nindiv==2],]
  
  twin_tcrs = twins %>% group_by(COMBAT_ID, clonotype) %>% dplyr::summarise(ncells = length(unique(barcode_id)), nstates = length(unique(cl)), mode_state = get_mode(cl))
  
  return(twin_tcrs)
}

display_combat_reference <- function(sref, gene_markers, protein_markers, exp_file = "/Users/klagattu/Downloads/combat_expnorm_incNKT.rds", prot_file="/Users/klagattu/Downloads/combat_protnorm_incNKT.rds"){
  gene_names = gene_markers
  gene_markers[gene_names=="PLZF"] = "ZBTB16"
  protein_names = protein_markers
  protein_names[protein_markers=="AB_TCR_Va24_Ja18"] = "AB_MAITTCR"
  md = sref$meta_data
  ##done once, stored in object:
  #V = t(sref$Z_orig)
  #set.seed(27)
  #uw_orig = uwot::umap(V, n_neighbors = 30, learning_rate = 0.5, init = "laplacian", metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1, min_dist = 0.3, n_threads = 4, ret_model = TRUE)
  #md$UMAP1_orig = uw_orig$embedding[,1]
  #md$UMAP2_orig = uw_orig$embedding[,2]
  
  md$UMAP1 = sref$umap$embedding[,1]
  md$UMAP2 = sref$umap$embedding[,2]
  
  g = ggplot(md[sample(1:nrow(md), nrow(md), replace=FALSE),], aes(UMAP1_orig, UMAP2_orig, color=Pool_ID))
  g = g + scale_color_futurama() + ggtitle("a")
  pool1 = my_umap_theme(g) 
  
  g = ggplot(md[sample(1:nrow(md), nrow(md), replace=FALSE),], aes(UMAP1_orig, UMAP2_orig, color=Institute))
  g = g + ggtitle("c")
  inst1 = my_umap_theme(g)
  
  g = ggplot(md[sample(1:nrow(md), nrow(md), replace=FALSE),], aes(UMAP1, UMAP2, color=Pool_ID))
  g = g + scale_color_futurama() + ggtitle("b")
  pool2 = my_umap_theme(g, show.legend=FALSE)
  
  g = ggplot(md[sample(1:nrow(md), nrow(md), replace=FALSE),], aes(UMAP1, UMAP2, color=Institute))
  g = g + ggtitle("d")
  inst2 = my_umap_theme(g, show.legend=FALSE)
  
  g = ggplot(md, aes(UMAP1, UMAP2, color=QC_pct_mitochondrial))
  g = g + scale_color_viridis(option="cividis") + labs(color="% mito") + ggtitle("e")
  mito = my_umap_theme(g, cont=TRUE)
  
  g = ggplot(md, aes(UMAP1, UMAP2, color=QC_ngenes))
  g = g + scale_color_viridis(option="cividis") + labs(color="n genes") + ggtitle("f")
  ngenes = my_umap_theme(g, cont=TRUE) 
  
  g = ggplot(md, aes(UMAP1, UMAP2, color=QC_scrub_doublet_scores))
  g = g + scale_color_viridis(option="cividis") + labs(color="doublet score") + ggtitle("g")
  doublet = my_umap_theme(g, cont=TRUE)
  
  qc_plots = list(pool1, pool2, inst1, inst2, mito, ngenes, doublet)
  
  exp = readRDS(exp_file)
  exp = exp[rownames(exp) %in% gene_markers,]
  exp = t(exp)
  exp = exp[as.character(md$X),]
  gene_plots= list()
  for (i in 1:length(gene_markers)){
    md$gene = exp[,which(colnames(exp)==gene_markers[i])]
    g = ggplot(md, aes(UMAP1, UMAP2, color=gene))
    g = g + scale_color_viridis() + labs(color=gene_names[i]) + ggtitle(letters[length(qc_plots)+i])
    gene_plots[[i]] = my_umap_theme(g, cont=TRUE)
  }
  rm(exp)
  
  prot = readRDS(prot_file)
  prot = data.frame(t(as.matrix(prot)))
  prot = prot[as.character(md$X),]
  prot_plots= list()
  for (i in 1:length(protein_markers)){
    md$prot = prot[,which(colnames(prot)==protein_markers[i])]
    g = ggplot(md, aes(UMAP1, UMAP2, color=prot))
    g = g + scale_color_viridis(option="plasma") + labs(color=protein_names[i]) + ggtitle(letters[length(qc_plots) + length(gene_plots)+i])
    prot_plots[[i]] = my_umap_theme(g, cont=TRUE)
  }
  rm(prot)
  
  return(list(qc_plots = qc_plots, gene_plots = gene_plots, prot_plots = prot_plots))
}

check_different_TCRresolutions <- function(tcr_fres_dir){
  setwd(tcr_fres_dir)
  files = list.files()
  file = vector(mode="character")
  cor = vector(mode="numeric")
  CV = vector(mode="numeric")
  for (i in 1:length(files)){
    load(files[i])
    print(files[i])
    print(dim(res$loadings$X))
    CV = c(CV, seq(1,4))
    cor = c(cor, cor_test[1:4])
    file = c(file, rep(files[i], 4))
  }
  df = data.frame(file, CV, cor)
  
  df$resolution = "CDRs"
  df$resolution[grepl("cdr3_only", df$file)] = "CDR3 only"
  df$resolution[grepl("FRandCDR", df$file)] = "CDRs and FRs"
  df$resolution[grepl("ints52", df$file)] = "CDRs + adjacent interactions"
  df$resolution[grepl("ints86", df$file)] = "CDRs + adjacent cross-chain interactions"
  res = unique(df$resolution)
  df$resolution = factor(df$res, levels=c("CDR3 only", "CDRs and FRs", "CDRs", "CDRs + adjacent interactions"))
  df = df[order(df$resolution),]
  df$align = "alignment: middle"
  df$align[grepl("LR", df$file)] = "alignment: flanking"
  g = ggplot(df, aes(CV, cor, color=resolution)) 
  g = g + geom_point() + geom_line() + theme_bw() + facet_wrap(~align)
  g = g + ylab("Canonical Correlation") + labs(color="TCR feature resolution")
  return(g)
}

plot_clone_invariance <- function(dir){
  setwd(dir)
  files = list.files()
  print(load(files[1]))
  ldgs = res$loadings$X[,1:4]
  colnames(ldgs) = paste(paste("run", paste(1), sep=""), paste("CV", seq(1,4), sep=""), sep="_")
  for (i in 2:length(files)){
    print(load(files[i]))
    sim = res$loadings$X[,1:4]
    colnames(sim) = paste(paste("run", paste(i), sep=""), paste("CV", seq(1,4), sep=""), sep="_")
    ldgs = cbind(ldgs, sim)
  }
  cm = cor(ldgs)
  anndf = data.frame(CV = as.factor(rep(seq(1,4), 100)))
  rownames(anndf) = rownames(cm)
  return(pheatmap(cm, annotation_row = anndf, annotation_col=anndf, show_rownames = FALSE, show_colnames = FALSE))
}

aggregate_tcr_cors <- function(dir){
  setwd(dir)
  files = list.files()
  all = data.frame((matrix(nrow=0, ncol=6)))
  colnames(all) = c("term", "estimate", "std.error", "statistic", "p.value", "file")
  for (i in 1:length(files)){
    ex = readRDS(files[i])
    ex$file = files[i]
    all = rbind(all, ex)
  }
  all$model = ""
  all$CCA = ""
  all$CV = ""
  for (i in 1:nrow(all)){
    info = strsplit(all$file[i], "_")[[1]]
    if (info[2] %in% c("v", "j")){
      all$model[i] = paste0(info[2:4], collapse="_")
      all$CV[i] = info[5]
      all$CCA[i] = info[6]
    } else if (info[2] %in% c("TRA", "TRB", "TRAcdr3", "TRBcdr3")){
      all$model[i] = paste0(info[2:3], collapse="_")
      all$CV[i] = info[4]
      all$CCA[i] = info[5]
    } else {
      all$model[i] = info[2]
      all$CV[i] = info[3]
      all$CCA[i] = info[4]
    }
  }
  all$type = "other"
  all$type[all$model %in% c("CDR3alen", "CDR3alenQ", "CDR3blen", "CDR3blenQ")] = "length"
  all$type[grepl("gene", all$model)] = all$model[grepl("gene", all$model)]
  all$type[grepl("TRAcdr3", all$model)] = "TRAcdr3"
  all$type[grepl("TRBcdr3", all$model)] = "TRBcdr3"
  all$type[all$model %in% c("TRA_p27", "TRA_p28", "TRA_p29", "TRA_p30", "TRA_p36", "TRA_p37", "TRA_p38")] = "TRAcdr1"
  all$type[all$model %in% c("TRA_p56", "TRA_p57", "TRA_p58", "TRA_p63", "TRA_p64", "TRA_p65")] = "TRAcdr2"
  all$type[all$model %in% c("TRB_p27", "TRB_p28", "TRB_p29",  "TRB_p36", "TRB_p37", "TRB_p38")] = "TRBcdr1"
  all$type[all$model %in% c("TRB_p56", "TRB_p57", "TRB_p58", "TRB_p59", "TRB_p63", "TRB_p64", "TRB_p65")] = "TRBcdr2"
  
  all = all[!(grepl("Intercept", all$term)),]
  all = all[all$p.value < 0.05/nrow(all),]
  
  all$estimate[all$CV=="X1"] = -all$estimate[all$CV=="X1"]
  all$estimate[all$CV=="X3"] = -all$estimate[all$CV=="X3"]
  return(all)
}

train_ren_combat_separately <- function(combat_file, ren_file){
  sref = readRDS("sref_combat_full_authTplusNKT_20hPCs_tcrfilt0607_nvargenes200_sampTH1.instTH0.5.poolTH0.5_origuwot.rds")
  ctcs = t(sref$Z_corr)
  colnames(ctcs) = paste("hPC", seq(1, 20), sep="")
  rownames(ctcs) = sref$meta_data$X
  tcr = readRDS("CRtrtest_0602/CR_xtrain.rds")
  g1 = plot_cca_replication(ren_file, tcr, ctcs)
  mapped = readRDS("ren_mappedtocombat_full_authTplusNKT__500g_20hPCs_tcrfilt0607_nvargenes200_theta10505.rds")
  rtcs = t(mapped$Z)
  colnames(rtcs) = paste("hPC", seq(1, 20), sep="")
  rownames(rtcs) = mapped$meta_data$cellName
  g2 = plot_cca_replication(combat_file, tcr, rtcs)
  return(list(g1=g1, g2=g2))
}

plot_cca_replication <- function(cca_file, tcr, tcs){
  colors = c("orangered3", "dodgerblue3", "turquoise4", "darkmagenta")
  int = intersect(rownames(tcr), rownames(tcs))
  print(length(int))
  tcr = tcr[as.character(int),]
  tcs = tcs[as.character(int),]
  
  print(load(cca_file))
  cca = res
  tcs = scale_variables(tcs, mns_y, sds_y)
  tcr = scale_variables(tcr, mns_x, sds_x)
  
  x = as.matrix(tcr) %*% as.matrix(cca$loadings$X[,1:4])
  y = as.matrix(tcs) %*% as.matrix(cca$loadings$Y[,1:4])
  Rp = numeric(4)
  Rr = numeric(4)
  Rlr = numeric(4)
  Rhr = numeric(4)
  Dp = numeric(4)
  Dr = numeric(4)
  Dlr = numeric(4)
  Dhr = numeric(4)
  for (i in 1:4){
    Dtest = cor.test(cca$variates$X[,i], cca$variates$Y[,i])
    Dp[i] = Dtest$p.value
    Dr[i] = Dtest$estimate
    Dlr[i] = Dtest$conf.int[1]
    Dhr[i] = Dtest$conf.int[2]
    Rtest = cor.test(x[,i], y[,i])
    Rp[i] = Rtest$p.value
    Rr[i] = Rtest$estimate
    Rlr[i] = Rtest$conf.int[1]
    Rhr[i] = Rtest$conf.int[2]
    print(i)
  }
  tp = data.frame(Dp, Dr, Dlr, Dhr, Rp, Rr, Rlr, Rhr)
  tp$CV = seq(1,4)
  tp$pv = sapply(as.character(tp$Rp), function(x) strsplit(x, "e")[[1]][1])
  tp$lbl = sapply(1:nrow(tp), function(x) gsub(tp$pv[x], round(as.numeric(as.character(tp$pv[x])), 1), as.character(tp$Rp[x])))
  tp$lbl = paste("P =", tp$lbl)
  tp$lbl[tp$lbl=="P = 0"] = "P < 1e-300"
  g = ggplot(tp, aes(Dr, Rr, color=factor(CV), label=lbl))  + geom_hline(yintercept=0) + geom_vline(xintercept=0)
  g = g + geom_point(size=0.5) + geom_errorbar(aes(xmin=Dlr, xmax=Dhr), width=0.005)
  g = g + geom_errorbar(aes(ymin=Rlr, ymax=Rhr), width=0.005) + geom_text_repel() + scale_color_manual(values=colors)
  g = g + theme_bw() + xlim(c(0,0.65)) ##+ ylim(c(0,0.))
  return(g)
}

viz_all_posTCRcors <- function(betas, ymin=-8, ymax=8, case_thresh=NULL, rev=FALSE){
  if (rev){
    betas$estimate = -betas$estimate
  }
  if (!(is.null(case_thresh))){
    betas = betas[betas$cases>=case_thresh,]
  }
  cdr1a = create_seqlogo(betas, sort(unique(betas$model[betas$type=="TRAcdr1"])), ymin, ymax)
  cdr2a = create_seqlogo(betas, sort(unique(betas$model[betas$type=="TRAcdr2"])), ymin, ymax)
  #^mce
  cdr1b = create_seqlogo(betas, sort(unique(betas$model[betas$type=="TRBcdr1"])), ymin, ymax)
  cdr2b = create_seqlogo(betas, sort(unique(betas$model[betas$type=="TRBcdr2"])), ymin, ymax)
  possCDR3a = paste("TRAcdr3_p", seq(2, 17), sep="")
  cdr3a = create_seqlogo(betas, possCDR3a, ymin, ymax)
  possCDR3b = paste("TRBcdr3_p", seq(2, 18), sep="")
  cdr3b = create_seqlogo(betas, possCDR3b, ymin, ymax)
  g = (cdr1a + cdr2a + cdr3a)/(cdr1b + cdr2b + cdr3b)
  return(g)
}

arrange_score_comps <- function(target, cdrscore_dir, fullscore_dir){
  savepath = getwd()
  setwd(cdrscore_dir)
  targets$target = targets[,which(colnames(targets)==target)]
  load(paste(target, "CDR1a_LR_TCRscore.RData", sep="_"))
  cdr1a = preds_min
  load(paste(target, "CDR2a_LR_TCRscore.RData", sep="_"))
  cdr2a = preds_min
  load(paste(target, "CDR1b_LR_TCRscore.RData", sep="_"))
  cdr1b = preds_min
  load(paste(target, "CDR2b_LR_TCRscore.RData", sep="_"))
  cdr2b = preds_min
  load(paste(target, "CDR3a_LR_TCRscore.RData", sep="_"))
  cdr3a = preds_min
  load(paste(target, "CDR3b_LR_TCRscore.RData", sep="_"))
  cdr3b = preds_min
  setwd("..")
  setwd(fullscore_dir)
  load(paste(target, "all_alphabeta_LR_TCRscore.RData", sep="_"))
  total = preds_min
  df = data.frame(cdr1a, cdr2a, cdr1b, cdr2b, cdr3a, cdr3b, total)
  colnames(df) = c("cdr1a", "cdr2a", "cdr1b", "cdr2b", "cdr3a", "cdr3b", "total")
  setwd(savepath)
  return(df)
}

plot_var_comp <- function(vc){
  tp = data.frame(val = 100*vc, chain = rep(c("alpha", "beta"), 3), loop = c("CDR1", "CDR1", "CDR2", "CDR2", "CDR3", "CDR3"))
  tp$region = tp$loop
  tp$region[tp$loop %in% c("CDR1", "CDR2")] = "CDR1-2"
  tp = tp %>% group_by(region, chain) %>% dplyr::summarise(val = sum(val))
  tp$chain = factor(tp$chain, levels=rev(levels(factor(tp$chain))))
  tp$region = factor(tp$region, levels=rev(levels(factor(tp$region))))
  g = ggplot(tp, aes(val, chain, fill=region))
  g = g + geom_bar(stat="identity", color="black") + theme_classic() + scale_fill_manual(values=c("powderblue", "salmon4"))
  return(g + xlab("% contribution") + xlim(c(0,75)))
}

get_var_comps_LR <- function(data, order=c("cdr1a", "cdr1b", "cdr2a", "cdr2b", "cdr3a", "cdr3b")){
  r2 = numeric(length(order))
  form = as.formula(paste("total ~", order[1]))
  r2[1] = summary(lm(form, data=data))$r.squared
  for (i in 2:length(order)){
    form = as.formula(paste("total ~", paste0(order[1:i], collapse=" + ")))
    r2[i] = summary(lm(form, data=data))$r.squared - sum(r2[1:(i-1)])
  }
  return(r2)
}

get_ccascore_umap <- function(combat_md, ren_mapped_file, cca_res, CV, side="Y", mp, trunc=NULL, order=FALSE, rev=FALSE){  ##clean up later
  md = combat_md
  md$dataset = "combat"
  md$cell = md$X
  md = md[,c("cell", "UMAP1", "UMAP2", "dataset")]
  
  mapped = readRDS(ren_mapped_file)
  renmd = data.frame(cell = mapped$meta_data$cellName, UMAP1 = mapped$umap[,1], UMAP2 = mapped$umap[,2], dataset="ren")
  crmd = rbind(md, renmd)
  
  fullCRres = cca_res
  yscores = data.frame(fullCRres$variates$Y)
  colnames(yscores) = paste("Y", seq(1:ncol(yscores)), sep="")
  yscores$cell = rownames(yscores)
  crmd = left_join(crmd, yscores)
  xscores = data.frame(fullCRres$variates$X)
  colnames(xscores) = paste("X", seq(1:ncol(xscores)), sep="")
  xscores$cell = rownames(xscores)
  crmd = left_join(crmd, xscores)
  tp = crmd[!(is.na(crmd$Y1)),]
  tp$color = tp[,which(colnames(tp)==paste(side, CV, sep=""))]
  colors = c("orangered3", "dodgerblue3", "darkmagenta", "turquoise4", "purple4")
  clr = colors[CV]
  if (rev) { tp$color = -tp$color }
  if (!(is.null(trunc))){ tp$color= sapply(tp$color, function(x) ifelse(x>=trunc, trunc, x))}
  if (order) { tp = tp[order(tp$color),] }
  g = ggplot(tp, aes(UMAP1, UMAP2, color=color))
  g = g + geom_point_rast(size=0.001) + scale_color_gradient2(low="palegoldenrod", mid="palegoldenrod", high=clr, midpoint=mp) + theme_classic() + labs(color=paste("CV", CV, sep=""))
  g = g +  theme(legend.text = element_text(size=15), legend.title = element_text(size=25))
  return(g)
}

geneprot_rankplot <- function(R, cv, labeled, rev=FALSE){
  df = data.frame(marker = colnames(R), R = R[cv,])
  ##NAs are due to 0 variance genes
  df = df[!(is.na(df$R)),]
  if (rev) { df$R = -df$R }
  df$lbl = df$marker
  df$lbl[!(df$marker %in% labeled)] = ""
  df = df[order(-df$R),]
  df$rank = seq(1, nrow(df))
  df$gp = ifelse(grepl("AB_", df$marker), "protein", "gene")
  print(head(df))
  print(tail(df))
  g = ggplot(df, aes(rank, R, color=gp, label=lbl))
  g = g + geom_point_rast(size=0.001) + scale_color_manual(values=c("darkslateblue", "deeppink2"))
  g = g + geom_text_repel(max.overlaps = Inf) + theme_classic()
  return(g)
}






