# FOR EXECUTING UNDER HPC
if (length(commandArgs(T)) > 0 & commandArgs(T)[1]=="HPC")
{
  print(paste(Sys.time(),"[Starting]"))

  library(plyr)
  library(reshape2)
  library(ggplot2)
  library(coda)

  load('index.Rdata')
  load('design.Rdata')
  load('parameters.Rdata')

  files = list.files(path="samplestats",pattern="^[0-9]+\\.Rdata")


  test_samples <- levels(design$Sample)[levels(design$Sample) != tolower(levels(design$Sample))]
  test_samples <- test_samples[2:length(test_samples)]

  samps <- mdply(files, .id=NULL, function(f) {
    load(paste0("samplestats/",f))
    samps <- data.frame(t(colMeans(s.Sol[,colnames(s.Sol) %in% paste0('Sample', levels(design$Sample)),drop=F])))
    colnames(samps) <- sub('Sample', '', colnames(samps))
    samps$ProteinID <- factor(as.integer(gsub("\\.Rdata","",f)))
    #samps$itt <- seq(1,nrow(samps))
    samps
  })


  stats <- mdply(files, .id=NULL, function(f) {
    load(paste0("samplestats/",f))
    stats$ProteinID <- as.integer(gsub("\\.Rdata","",f))
    stats
  })
  stats$Sample <- factor(stats$Sample)

  results <- merge(data.index, stats)

  res <- unique(results[,c("ProteinID","N","Protein","Peptides","Spectra","Sample","mean")])

  res <- dcast(res, ProteinID + N + Protein + Peptides + Spectra ~ Sample)
  res <- res[with(res, order(ProteinID)),]

  write.csv(res,paste0(parameters$Value[parameters$Key=="id"],"_Samples.csv"), row.names=F)


  print(paste(Sys.time(),"[Finished]"))
}
