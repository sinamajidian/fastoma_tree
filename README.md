
# FastOMA

FastOMA is a novel and fast pipeline for species tree inference using placement in Hierarchical Orthologous Groups (HOGs).


# prerequisites

Please install these

- [OMAmer](https://github.com/DessimozLab/omamer/blob/master/README.md) using `pip install omamer`.
- [Pyoma](https://github.com/DessimozLab/pyoma) using `pip install pyoma`
- [zoo wrapper] ??
- [mafft](https://mafft.cbrc.jp/)


and packages datetime, numpy, ast, pickle .

Download these  

- [OMAmer dataset](https://zenodo.org/record/4593702/files/LUCA.h5)  (10Gb)  ??
- [HOH-OG-map.dic] >> my zenodo   (100Mb)
-  [OMA dataset](https://omabrowser.org/All/OmaServer.h5) and [its index file](https://omabrowser.org/All/OmaServer.h5.idx) (95+60Gb).

```
cd test_data; mkdir dataset; cd dataset
wget https://zenodo.org/record/4593702/files/LUCA.h5     ??
wget https://omabrowser.org/All/OmaServer.h5
wget https://omabrowser.org/All/OmaServer.h5.idx
```
Be careful about different OMA releases. The OMAmer dataset and OMA dataset should be from the same release since HOG id changes in every release.

# Pipeline

Consider you have the proteome of two species as two fasta files `WoodlandKingfisher.fa` (Halcyon senegalensis, HALSEN) and `AtlanticCanary.fa` (Serinus canaria, SERCAN) which are provided in the `test_data` folder. Note that the fasta file shouldn't contain U, you could replace it with X with `sed s/U/X/g` otherwise.


The first step is to find the placement of each protein of each species in the HOG dataset.  To do so , we use OMAmer.  We need to download the OMAmer dataset. Then, we can run OMAmer:


```
cd ..; mkdir hogmap

number_threads=20

omamer search --db  dataset/LUCA.h5 --query proteome/WoodlandKingfisher.fa --nthreads ${number_threads}  --out hogmap/WoodlandKingfisher.hogmap
omamer search --db  dataset/LUCA.h5 --query proteome/AtlanticCanary.fa --nthreads ${number_threads}  --out hogmap/AtlanticCanary.hogmap
```

After few minutes, the HOG maps of species are ready in the `hogmap` folder.

```
number_threads=20
python fastoma-tree.py  ./test_data/dataset/   ./test_data/  ${number_threads}
```
The second argument should contain two folders `proteome` and `hogmap`. The output will be the multiple sequence alignment supermatrix.


Finally, you could run a tree inference method
```
fasttree _20_msa_concatanated.txt

(CHLSB:0.007967770,PAPAN:0.009104705,(HUMAN:0.012101579,((AtlanticCanary:0.327152404,(MOUSE:0.121086106,CHILA:0.087558595)0.971:0.015068557)0.167:0.000847169,WoodlandKingfisher:0.353362694)1.000:0.041260314)1.000:0.010528717);
```





# bird data set

We used following thresholds for bird dataset:
```
keep_og_treshold_species_query =   360
threshold_least_query_sepecies_in_OG = 320
kept_oma_species_num = 20
ogs_keep_number = 100
```






done!
