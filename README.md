# About

Similator that genarate artifitial SNV and CNV at the same time at given tumor content rtaio.

# Licence

 Apache 2


# Install

- Download karkinosSim.jar from web sites.
- Install Java Runtime 1.5 or later.

# Required files

- Bam files for normal reads Alignment.
- 2.bit reference file
- CaptureTarget Information (bed format,optional)
- Low depth region (bed format,optional)

# Run

usage: karkinosSim <command> options

possible commands are; devideBAM ; assignSNVIndel ; similate ; checkAnswer 

1. Devide normal bam files in to 3 bam files.
   karkinosSim need 3 normal bamfiles from same sample.
   devideBAM can devide one bam file to 3.
   If you can provide 3 normal bamfiles, skip this process.


```
usage: karkinosSim.jar devideBAM -bam <arg> -out <arg>

 -bam,--bam <arg>   input bam file
 -out,--out <arg>   output directory

```

2. Create positions Vcf for random artifitial SNV and Indel.


```
usage: karkinosSim.jar assignSNVIndel -target <arg> -r <arg> -numSnv <arg>
       -numIndel <arg> -indelfrac <arg> -lowdepth <arg> -outvcf <arg>

 -target,--target bed <arg>     target bed file
 -r,--ref <arg>                 2bit reference file
 -numSnv,--numSnv <arg>         number of snv to similate
 -numIndel,--numIndel <arg>     number of Indel to similate
 -lowdepth,--lowdepth <arg>     low depth vcf
 -outvcf,--outvcf <arg>         output vcf file

```

3. Genarate similated BAM and FastQ files.

```   
usage: karkinosSim.jar similate -id <arg> -target <arg> -r <arg> -SNP
       <arg> -SNVgen <arg> -SNVgenSub <arg> -CNVgen <arg> -tc <arg>
       -normalbam1 <arg> -normalbam2 <arg> -out <arg>
 -id,--id <arg>                   id for this similation
 -target,--target bed <arg>       target bed file
 -r,--ref <arg>                   2bit reference file
 -SNP,--SNPVcf <arg>              SNP list exclude SNP positions
 -SNVgen,--SNVgen <arg>           SNV to generate, vcf
 -SNVgenSub,--SNVgenSub <arg>     subpolulation SNV to generate, vcf
 -CNVgen,--CNVgen <arg>           CNV to generate, vcf
 -tc,--tumorContent <arg>         tumor content ration for similation
 -normalbam1,--normalbam1 <arg>   normalbam1
 -normalbam2,--normalbam2 <arg>   normalbam2
 -out,--out <arg>                 output directory

```

4. Check Answer,Summary Output to  System.out


```
usage: karkinosSim.jar checkAnswer -answer <arg> -check <arg>
 -answer,--answer <arg>   answer vcf file
 -check,--check <arg>     vcf file to check

```


