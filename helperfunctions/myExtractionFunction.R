myExtractionFunction = function(data,SH,set){
  # data = copy(gx_assoc_pre_filt)
  # SH = mySH
  # set = "females_pre"
  
  mytab1 = data[[SH]]
  mytab2 = copy(mytab1$assoc)
  myNames = c("PROBE_ID","symbol_INGENUITY","description_INGENUITY","location_INGENUITY",
              "types_INGENUITY","biomarkers_INGENUITY","drugs_INGENUITY", "ilmn_chr_hg19",
              "ilmn_start_hg19","ilmn_end_hg19","ilmn_genelength_hg19","ilmn_entrezID_INGENUITY",
              "logFC","CI.L","CI.R","AveExpr","t","P.Value","adj.P.Val","var","FC","n_chips")
  colsOut<-setdiff(colnames(mytab2),myNames)
  mytab2[,get("colsOut"):=NULL]
  
  setnames(mytab2,"logFC","beta")
  setnames(mytab2,"n_chips","n_samples")
  setnames(mytab2,"var","phenotype")
  setnames(mytab2,"FC","FoldChange")
  setnames(mytab2,"t","t_statistic")
  
  mytab2[,se := beta / t_statistic]
  mytab2[,CI.L := NULL]
  mytab2[,CI.R := NULL]
  mytab2[,setting := set]
  
  return(mytab2)
}
