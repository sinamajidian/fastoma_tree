
# FastOMA

FastOMA is a novel and fast pipeline for species tree inference using placement in Hierarchical Orthologous Groups (HOGs).


# Pipeline

Consider you have the proteome of two species as two fasta files `WoodlandKingfisher.fa` and `AtlanticCanary.fa` which are provided in the `test_data` folder.


The first step is to find the placement of each protein of each species in the HOG dataset.  To do so , we use OMAmer. Please install [OMAmer](https://github.com/DessimozLab/omamer/blob/master/README.md) using `pip install omamer`. We need to download the OMAmer dataset. Then, we can run OMAmer:

```
wget https://zenodo.org/record/4593702/files/LUCA.h5

omamer search --db  LUCA.h5 --query WoodlandKingfisher.fa --nthreads 1  > WoodlandKingfisher.hogmap
omamer search --db  LUCA.h5 --query AtlanticCanary.fa --nthreads 1 > AtlanticCanary.hogmap
```
After few minutes, you will have the output.
