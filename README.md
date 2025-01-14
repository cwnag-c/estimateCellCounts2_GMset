# estimateCellCounts
The estimateCellCounts2 in FlowSorted.Blood.EPIC was modified to support GenomicMethylSet (output from minfi::readGEORawFile for GEO datasets). This simplified version only supports blood cell types, enforces preprocessQuantile normalization, and uses probeSelect = "IDOL" by default.
