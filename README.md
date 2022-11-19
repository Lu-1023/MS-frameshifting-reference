**<font color="grey"><font size=10>Make frameshifting reference </font></font>**
<font size=5><font color="steelblue"><p align="right">2022.11.17</p></font></font>
# <font color="steelblue">(Pipe for making reference of MS after first database search) </font>



[TOC]

***
##  <font size=6>1   Make frameshifting peptide</font>
If the proteins were expressed in the data, we selected these proteins to make references for frameshift peptides of plus and minus.

```shell
python make_ref_shift0.py ${cds_file} ${pep_file} ${shift_file} ${proteinGroup_file} ${outfile} 
#cds_file is the DNA sequence of the coding sequence
#pep_file is the Protein sequence of the coding sequence
#shift_file is the CRIF-Finder result of codon repeat sites which can significantly induce frameshifting
#proteinGroup_file is the first database search results of maxquant
#outfile is the output file  
```
[make_ref_shift0.py](https://github.com/Lu-1023/MS-frameshifting-reference/blob/main/mass_ref_shift0.py)

##  <font size=6>2   Add original protein as background</font>
In order to screen out the real frameshift peptide, the expressed original protein was added to reference as background.
```shell
python origin_protein.py ${pep_file} ${proteinGroup_file} ${outfile}
#pep_file is the Protein sequence of the coding sequence
#proteinGroup_file is the first database search results of maxquant
#outfile is the output file 
```
[origin_protein.py](https://github.com/Lu-1023/MS-frameshifting-reference/blob/main/origin_protein.py)


##  <font size=6>3 Maxquant autorun   </font>
Because of the large number of raw data samples, it needs to be run 56 times, so the python script is used to run it automatically.

```shell
python autorun.py

```
[autorun.py](https://github.com/Lu-1023/MS-frameshifting-reference/blob/main/autorun.py)
<!--/TOC-->
