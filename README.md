

# COCOS: Codon Consequence Scanner 
COCOS is a plugin for the Ensembl Variant Effect Predictor (VEP) plugin for annotating reading frame changes.
The plugin captures Amino Acid sequence alterations stemming from variants that produce an altered reading frame, e.g. stop-lost variants and small genetic Insertion and Deletions (InDels).  As a result, the stop codon of the transcript may be changed and loses its functionality. The coding sequence is terminated when a new subsequent stop-codon is reached, terminates the translation process. 


## Install VEP on your system
find installation instructions here (<a href=http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html>link</a>)

## Example commandline to run the COCOS plugin with VEP

```
./variant_effect_predictor.pl -i /path/your/input.vcf --cache --offline --force_overwrite --dir_plugins /path/to/cocos/ --plugin cocos
```

## Testing COCOS

COCOS can be tested by using the example .vcf file in the 'example' folder.


## Testing COCOS using the Ensembl virtual machine via VirtualBox

Ensembl provides a virtual machine (ftp.ensembl.org/pub/current_virtual_machine) for <a href="https://www.virtualbox.org/">VirtualBox</a>. On the virtual machine, VEP is pre-installed and can be used immediately. The user can download and install the virtual machine machine with 5 simple steps outline on this <a href=http://www.ensembl.org/info/data/virtual_machine.html>Ensembl website</a>.

Note: the virtual machine image is about 2.6GB in size. The following step are tested on a ubuntu 14.04 linux environment. Once you have VirtualBox installed and the Ensembl virtual machine started, you may proceed with the follwoing steps:

<li>start</li>
double click the icon 'VEP terminal' on the Ensembl ubuntu desktop. It will lead you to the path: /home/ensembl/ensembl-api-folder-86/ensembl-tools/scripts/variant_effect_predictor

<li>install vep cache file for homosapiens</li>
```
perl INSTALL.pl 
```
skip installing the API; download the VEP cache file for homosapiens `homo_sapiens_vep_86_GRCh37.tar.gz` (choice: 44)

<li>get cocos plugin</li>
```
git clone git://github.com/butkiem/COCOS.git cocos-github
```
This will create a clone cocos github repository in the folder: ./cocos-github

<li>run VEP with cocos</li>
```
./variant_effect_predictor.pl -i ./cocos-github/example/example.vcf --cache --offline --force_overwrite --dir_plugins ./cocos-github/ --plugin cocos
```
