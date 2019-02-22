## Gene_connect

Script designed for visually represent syntenic regions based on orthology relations between
genes in a single reference scaffold of one species and orthologous genes of another species
which may appear in multiple scaffolds.

## Usage

Usage:

```
perl draw_colinearity.pl orthology_list.txt refSpecies.geneIDs querySpecies.geneIds suffix
```
The main input file of the program is a conection file with two columns separated by a space or tab.
The first column contains the ids of the genes of the reference scaffold. The second column contains
orthologous gene ids of the query species. Each line represents a single orthology relation.

The format for gene names needs to be:

"ScaffoldName__StartPosition__EndPosition__GeneName"

Note the double-underscore characters, there should not be another double-underscore in the ScaffoldName or geneName.
We call this string the geneID of a gene.

If a gene does not have an ortholog, a "-" may appear in the other column.
Alternatively, two additional files may be provided that contain the complete lists of geneIDs in the respective
scaffolds of each of the two species (or just all geneIDs of all scaffolds of each species) in one column.
The orthology list as well as the additional geneID lists will be sorted by the script, i.e. they
do not need to be sorted initially.

The "suffix" will be the last part of the output image filename before the ".png" extension. If omitted the image
files will have the name of the reference scaffold (and extension ".png"), otherwise "refScafName.suffix.png".

The input may be either one ortholgy file only, or one orthology file plus two geneID files plus suffix string.
If only an orthology file is provided, only the genes contained therein will be drawn.

Formatting example for an orthology list including dashes for missing orthologs:
 "0020.scaffold00069__1__4111__Bv_000500_kues.t1	-
 0020.scaffold00069__32216__33048__Bv_000510_ffys.t1	-
 0020.scaffold00069__34996__59920__Bv_000520_xuag.t1	scaffold1170__32936__58560__mar_g9410.t1
 0020.scaffold00069__69976__71790__Bv_000530_rsxf.t1	scaffold1170__63346__68792__mar_g9411.t1
 0020.scaffold00069__71766__75423__Bv_000540_hoty.t1	-
 \-	scaffold1170__63346__68792__mar_g9411.t1
 0020.scaffold00069__82096__92059__Bv_000550_zstx.t1	scaffold1170__73726__83471__mar_g9412.t1
 \-	scaffold1170__83974__85481__mar_g9413.t1
 0020.scaffold00069__95506__96127__Bv_000570_pmqw.t1	scaffold1170__87426__88060__mar_g9414.t1"

Formatting example for an orthology list showing orthologs only:
 0020.scaffold00069__34996__59920__Bv_000520_xuag.t1	scaffold1170__32936__58560__mar_g9410.t1
 0020.scaffold00069__69976__71790__Bv_000530_rsxf.t1	scaffold1170__63346__68792__mar_g9411.t1
 0020.scaffold00069__82096__92059__Bv_000550_zstx.t1	scaffold1170__73726__83471__mar_g9412.t1
 0020.scaffold00069__95506__96127__Bv_000570_pmqw.t1	scaffold1170__87426__88060__mar_g9414.t1

Example for gene lists to show all genes of the scaffold (species 1):
 0020.scaffold00069__1__4111__Bv_000500_kues.t1
 0020.scaffold00069__32216__33048__Bv_000510_ffys.t1
 0020.scaffold00069__34996__59920__Bv_000520_xuag.t1
 0020.scaffold00069__69976__71790__Bv_000530_rsxf.t1
 0020.scaffold00069__71766__75423__Bv_000540_hoty.t1
 0020.scaffold00069__82096__92059__Bv_000550_zstx.t1
 0020.scaffold00069__95506__96127__Bv_000570_pmqw.t1
 
Example for a gene list of species 2:
 scaffold1170__32936__58560__mar_g9410.t1
 scaffold1170__63346__68792__mar_g9411.t1
 scaffold1170__63346__68792__mar_g9411.t1
 scaffold1170__73726__83471__mar_g9412.t1
 scaffold1170__83974__85481__mar_g9413.t1
 scaffold1170__87426__88060__mar_g9414.t1

# paralelize -> juliane command

## Output


