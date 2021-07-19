# LongPass



## get started

```python
python longpass/LongPass.py -o outcluster/out_cluster.txt --clustering paraclu --normalization raw  --params longpass/params.txt --replicate --bamfile sequinB_rep1.trimmed.sort.bam sequinB_rep2.trimmed.sort.bam
```



You can adjust parameters in params.txt file:

example of parameter file:

```shell
powerlaw:
    value_range = [10,1000]
distclu:
    maxdist = 20
paraclu:
    minStability = 1
    maxLength = 500
    removeSingletons = True
    keepSingletonAbove = 15
    reducetoNoneoverlap = True
```



## requirement

```
```



## input file

### types:

```text
bam
lrtsp
range
```



### examples of lrtsp:

1-base file.

#### tbsfile:

first column: chromosome

second column: position

third column: strand

fourth column: TSS count at this position

fifth column: TES count at this position

```
chrIS   77834   +       29      0
chrIS   90912   +       12      0
chrIS   94351   +       114     0
```



#### tssfile/tesfile:

first column: chromosome

second column: position

third column: strand

fourth column: TSS/TES count at this position

```
chrIS   77834   +       29
chrIS   90912   +       12
chrIS   94351   +       114
```



## OUTPUT FILE

example of output file:

0-base file

```
chrIS   +       2298169 2298173 2       2298169 190.0   200.0   3.33333 0.07033 tss     SRR13057605.7787,SRR13057605.1473
chrIS   +       4771512 4771514 2       4771512 16.0    29.0    13.0    11.23858        tss     SRR13057605.525671,SRR13057605.1473
chrIS   +       4778059 4778062 3       4778059 90.0    217.0   63.5    29.96   tss     SRR13057605.412450,SRR13057605.1473
chrIS   +       3998391 3998466 9       3998465 61.0    198.0   1.6875  0.13665 tss     SRR13057605.436442,SRR13057605.1473
chrIS   +       4777903 4777908 4       4777907 109.0   249.0   35.0    30.07874        tss     SRR13057605.1332445,SRR13057605.91742,SRR13057605.1473
```

first column: chromosome

second column: gene strand

third column: cluster start site

fourth column: cluster end site

fifth column: peak number in cluster

sixth column: dominant site in cluster

seventh column: number of reads at dominant site

eighth column: number of total reads in the cluster

ninth column: max density of the cluster

tenth column: min density of the cluster

11st column: tss/tes

12nd column: reads ids assign to this cluster

