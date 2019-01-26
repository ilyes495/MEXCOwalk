cancers = c('BLCA','BRCA','CRC', 'GBM', 'HNSC', 'KIRC', 'LAML', 'LUAD', 'LUSC', 'OV', 'UCEC')
for (cancer in cancers){
  
  CancerProject <- sprintf("TCGA-%s", cancer)
  DataDirectory <- paste0("../data/TCGA/",gsub("-","_",CancerProject))
  FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")
  
  
  query <- GDCquery(project = CancerProject,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  samplesDown <- getResults(query,cols=c("cases"))
  if (cancer == 'LAML') {
    smTP = 'TB'
  } else{
    smTP = 'TP'
    }
  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                    typesample = smTP)
  
  
  
  queryDown <- GDCquery(project = CancerProject, 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts", 
                        barcode = c(dataSmTP))
  GDCdownload(query = queryDown,
              directory = DataDirectory)
  
  
  dataPrep <- GDCprepare(query = queryDown, 
                         save = TRUE, 
                         directory =  DataDirectory,
                         save.filename = FileNameData)
  
  
  GBM<- assay(dataPrep) 
}
