meta <- data.frame(
  Title = "Spatial allelic expression counts for fly cross embryo",
  Description = "Allelic expression counts of spatial slices of a fly embryo from Combs & Fraser (2018), a D melanogaster x D simulans reciprocal cross",
  BiocVersion = "3.14",
  Genome = "dm6", 
  SourceType = "TXT",
  SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102233",
  SourceVersion = "v1",
  Species = "Drosophila melanogaster",
  TaxonomyId = 7227,
  Coordinate_1_based = TRUE,
  DataProvider = "Fraser Lab, Stanford",
  Maintainer = "Michael Love <michaelisaiahlove@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "Rda",
  RDataPath = "spatialDmelxsim/v1/spatialDmelxsim.rda",
  Tags = "allelic:ASE:spatial:embryo:patterning",
  Notes = "")

write.csv(meta, file="metadata.csv", row.names=FALSE)
