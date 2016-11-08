

# COCOS: Codon Consequence Scanner 
COCOS is a plugin for the Ensembl Variant Effect Prefictor (VEP) plugin for annotating reading frame changes.
The COCOS plugin captures Amino Acid sequence alterations stemming from variants that produce an altered reading frame, e.g. stop-lost variants and small genetic Insertion and Deletions (InDels).  As a result, the stop codon of the transcript may be changed and loses its functionality. The coding sequence is terminated when a new subsequent stop-codon is reached, terminates the translation process. 


## Install VEP on your system
find installation instructions here (<a href=http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html>link</a>)

## Example commandline to run the COCOS plugin with VEP

```
./variant_effect_predictor.pl -i /path/your/input.vcf --cache --offline --plugin cocos
```

## Testing COCOS

COCOS can be tested by using the exa


## Testing using a virtual machine

Ensembl provides a <a href="ftp://ftp.ensembl.org/pub/current_virtual_machine" class="external-link" rel="nofollow">virtual_machine</a> for virtual box (free). VEP is pre-installed and can be readily used.
The user can install the virtual machine machine with 5 simple steps outline on the Ensembl website - <a href=http://www.ensembl.org/info/data/virtual_machine.html>link</a>.
<li>Obtain VirtualBox</li>
<li>Download and import the virtual machine</li>
<li>Create shared folders</li>
<li>Start and verify the Ensembl installation</li>
<li>(Optional) Resize virtual disk</li> 

