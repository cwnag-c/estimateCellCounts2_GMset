# estimateCellCounts2_GMset
The function ```estimateCellCounts2_GMset``` extends ```estimateCellCounts2``` from **FlowSorted.Blood.EPIC** with specialized support for **GenomicMethylSet** objects , particularly addressing the output format from minfi::readGEORawFile. (Only for blood sample ,  preprocessQuantile normalization and utilizes probeSelect = "IDOL".)

This code depends on `FlowSorted.Blood.EPIC`. Please install this package and cite the **original publication** of `FlowSorted.Blood.EPIC` when using this code: [https://github.com/immunomethylomics/FlowSorted.Blood.EPIC](url).

```
#Step1:Download ref and load
library(ExperimentHub)
library(ExperimentHub)
eh <- ExperimentHub()
query_results <- query(eh, "FlowSorted.Blood.EPIC")
reference_data <-eh[[query_results$ah_id]]
save(reference_data, file = "FlowSorted.Blood.EPIC.rda")
referenceset<-load("./FlowSorted.Blood.EPIC.rda")##Load reference data locally
```

```
#Step2: read GenomicMethylSet from intensity
GenomicMethylSet$Sex<-targets$Gender ##Ensure that the gender information is coded as "F" and "M" (F represents Female, M represents Male)
cellp <-estimateCellCounts2_GMset(GenomicMethylSet,referenceset=referenceset)
cellp <-cellp2$prop
```

```
#Step3: eastimate
cellp <-estimateCellCounts2_GMset(GenomicMethylSet,referenceset=referenceset)
cellp <-cellp2$prop
```
