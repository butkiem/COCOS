

# COCOS: Codon Consequence Scanner 
COCOS is a plugin for the Ensembl Variant Effect Predictor (VEP) plugin for annotating reading frame changes.
The plugin captures Amino Acid sequence alterations stemming from variants that produce an altered reading frame, e.g. stop-lost variants and small genetic Insertion and Deletions (InDels).  As a result, the stop codon of the transcript may be changed and loses its functionality. The coding sequence is terminated when a new subsequent stop-codon is reached, terminates the translation process. 


## Install VEP on your system
find installation instructions here (<a href=http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html>link</a>)

## Example commandline to run the COCOS plugin with VEP

```
./variant_effect_predictor.pl -i /path/your/input.vcf --cache --offline --plugin cocos
```

## Testing COCOS

COCOS can be tested by using the example .vcf file in the 'example' folder.


## Testing COCOS using a virtual machine

Ensembl provides a virtual machine (ftp.ensembl.org/pub/current_virtual_machine) for virtual box (free). VEP is pre-installed and can be used immediately. The user can download and install the virtual machine machine with 5 simple steps outline on this <a href=http://www.ensembl.org/info/data/virtual_machine.html>Ensembl website</a>.

Note: the virtual machine image is about 2.6GB in size.

```
<li>ftp://ftp.ensembl.org/pub/current_virtual_machine</li>
```


