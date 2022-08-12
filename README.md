# Postassm
This pipeline offers a variety of evaluations for genome assembly outcomes, such as accuracy, continuity, and completeness.

## Quick start 

  1. change the config.yaml according to you environment. 

       
           vim config.yaml 

2. run with snakemake    

          
          nohup snakemake -s postassm.smk -j 10 -k --ri >sublog 2>&1 &
          
3. run on a cluster

  
          nohup snakemake -s postassm.smk -j 10 -k --ri --cluster "qsub -l nodes=1:ppn=20 -l walltime=999:00:00" >sublog 2>&1 & 
  
 ## Configuration
  * project_name: project  name
  * dir_work: work directory
  * software: software path (absolute path)
      * busco:
      * merqury:   
      * samtools: 
      * Rscript: 
      * xf_stat: 
      * seqtk: 
      * ...

  * lib_dir: database directory 
      * busco: path of BUSCO lib, for example /path/to/busco/datasets/mammalia_odb10
      * busco_conf: config of BUSCO lib, for example /path/to/busco/conf/busco-master/config/config.ini
 
  * samples:
      * sample1:
        * assm: sample1.fasta
        * meryl_lab: /path/to/meryl_db/sample1.meryl
      * sample2:
        * assm: /path/to/sample2.fa
        * meryl_lab: /path/to/meryl_db/sample2.meryl
      * ...
  * threads:
    * busco: 48
    * default: 2
    * ...
    
 ## Support tools 
 * BUSCO 
 * Merquery
 ## Contribution 
 If you want to apply other tools to evaluate the genome, we encourage you to pull a request or email us. 
 
 ## Contact  
 * Peng Jia (pengjia@stu.xjtu.edu.cn)

