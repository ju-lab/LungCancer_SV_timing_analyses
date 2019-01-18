
This script generates nonamer from SNV of indel calls.

###Step1

```
python 01_SNPannotate_ju_9merProtein.py <VCF>
```

###Step2

```
python 02_make_fsa.py <annotated VCF> <sampleid> <column number of Step1 annotation>
```

