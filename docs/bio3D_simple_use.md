# Bio3D
Bio3D is an R package containing utilities for the analysis of protein structure, sequence and trajectory data.We can use bio3D to perform some simple tasks。
It is best to have a basic understanding of pdb files before using subsequent code.
## Calculate distance

We can certainly use Chimerax or Pymol to obtain the residues at the binding interface, but calculating the distance through bio3D provides another strategy for determining the residues at the binding interface.
``` r
## basic example code
library(wlab.block)

## caculate interchain distance
interchain_distance <- bio3D_minimum_interchain_distances(
    input_file = "path/to/pdb",
    chain_query = "A",
    chain_target = "B")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis <- bio3D_chain_hetatm_distances(
    input_file = "path/to/pdb",
    chain_query = "A",
    hetatm_target = "GNP"
)
```
Using distance information can easily determine the position of residues. Here is a simple example.
``` r
library(data.table)

interchain_distance[scHAmin_ligand <5 , type := 'interface']
id<-stringr::str_c(interchain_distance[type=='interface',Pos], collapse = ",")
```
## B factor
B-factors in PDB files commonly are seen as a measure of (local) mobility in the (macro)molecule.By using biopython, we can convert bfactor into residue depth or Solvent Accessible Surface Area (SASA), and then use bio3D to obtain relevant data.
``` r

res_depth <- bio3D_temperature_factor(
    input_file = "path/to/res_depth.pdb",
    chain = "A"
)

sasa <- bio3D_temperature_factor(
    input_file = "path/to/sasa.pdb",
    chain = "A"
)
```
## Other use
Here, some simple applications have been added, such as obtaining secondary structures.
``` r
ss8 <- bio3D_secondary_structure(
    input_file = "path/to/pdb",
    chain = "A"
)
```

