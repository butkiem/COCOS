

# COCOS: Codon Consequence Scanner 
COCOS is a plugin for the Ensembl Variant Effect Predictor (VEP) plugin for annotating reading frame changes.
The plugin captures Amino Acid sequence alterations stemming from variants that produce an altered reading frame, e.g. stop-lost variants and small genetic Insertion and Deletions (InDels).  As a result, the stop codon of the transcript may be changed and loses its functionality. The coding sequence is terminated when a new subsequent stop-codon is reached, terminates the translation process. 


## Install Ensembl VEP on your system
find installation instructions here (<a href=http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html>link</a>)

## Quickstart example to run the COCOS plugin with VEP

Once VEP is installed (example is given for a linux environment), change to the directory where VEP is installed.
It contains the executable variant_effect_predictor.pl perl script. The '/path/to' placeholder have to be appropriately replaced.
```
cd /path/to/VEP
```

Install the COCOS plugin from github
```
git clone git://github.com/butkiem/COCOS.git /path/to/cocos/cocos-github
```

Run VEP with the COCOS plugin
```
./variant_effect_predictor.pl -i /path/to/cocos-github/example/example.vcf --cache --offline --force_overwrite --dir_plugins /path/to/cocos-github --plugin cocos,"./path/to/my_cocos.fasta"
```

COCOS can be tested by using the example .vcf file in the 'example' folder as you input. The COCOS output file can be specified as an argument to the `--plugin` flag as shown in the example commandline.


## Testing COCOS using the Ensembl virtual machine via VirtualBox
In case you do not have VEP installed but would still like to try the COCOS plugin, you can use the Ensembl virtual machine using VirtualBox.

### installing VirtualBox
VirtualBox is a freely availble general-purpose full virtualizer, that allows user to run virtual machines on their systems.
Note: Make sure that Intel Virtualization Technology (also known as Intel VT) or AMD-V depending on the brand of the processor is activated in your BIOS.

VirtualBox can be download from  <a href="https://www.virtualbox.org/">here</a>. The 'Downloads' section contains information on how to install the software on you particular system.

### install and run the Ensembl virtual machine 

Ensembl provides a virtual machine (ftp.ensembl.org/pub/current_virtual_machine) for <a href="https://www.virtualbox.org/">VirtualBox</a>.  The user can download and install the virtual machine machine with 5 simple steps outline on this <a href=http://www.ensembl.org/info/data/virtual_machine.html>Ensembl website</a>.

Note: the virtual machine image is about 2.6GB in size (version 86). 

### download and install COCOS on the Ensembl virtual machine 

Start the Ensembl virtual machine. VEP is pre-installed and can be used immediately.
The following steps are tested on a ubuntu 14.04 linux host environment.

#### Start
Once you started the virtual machine, double click the icon 'VEP terminal' on the Ensembl ubuntu desktop. It will lead you to the `path: /home/ensembl/ensembl-api-folder-86/ensembl-tools/scripts/variant_effect_predictor` where VEP is installed.

#### Install vep cache file for homosapiens
VEP provides the option to install local caches spific for different species. 

To install the VEP cache for homosapiens, you can run:
```
perl INSTALL.pl 
```
skip installing the API; download the VEP cache file for homosapiens `homo_sapiens_vep_86_GRCh37.tar.gz` (choice: 44)

Note: The download process might take some time.

####Download and install COCOS plugin

You can clone the COCOS github repository and save it in the directory `./cocos-github`
```
git clone git://github.com/butkiem/COCOS.git ./cocos-github
```

####Run VEP with COCOS

The COCOS github repository has an example vcf file that you can use for testing purposes.
Now you can run VEP with the COCOS plugin:

```
./variant_effect_predictor.pl -i ./cocos-github/example/example.vcf --cache --offline --force_overwrite --dir_plugins ./cocos-github/ --plugin cocos,"./cocos.fasta"
```

When successful, COCOS will generate a fasta file with the results called `./cocos.fasta`.




