# estimateCellCounts3
The estimateCellCounts2 in FlowSorted.Blood.EPIC was modified to support GenomicMethylSet (output from minfi::readGEORawFile for GEO datasets). This simplified version only supports blood cell types, enforces preprocessQuantile normalization, and uses probeSelect = "IDOL" by default.

**Only need to input a GenomicMethylSet in minfi format.**

This code depends on `FlowSorted.Blood.EPIC`. Please install this package and cite the **original publication** of `FlowSorted.Blood.EPIC` when using this code: [https://github.com/immunomethylomics/FlowSorted.Blood.EPIC](url).

```
cell_P<-estimateCellCounts3(GenomicMethylSet)
cell_P<-cell_p$prop
```
