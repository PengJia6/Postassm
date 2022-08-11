# postassm
This pipeline offers a variety of evaluations for genome assembly outcomes, such as accuracy, continuity, and completeness.

## How to run? 

1. change the config.yam according to you environment. 

``` vim config.yaml ```

2. run with snakemake    

  ``` 
  nohup snakemake -s postassm.smk -j 10 -k --ri >sublog 2>&1 &
  ```
3. run on a cluster

  ```
  nohup snakemake -s postassm.smk -j 10 -k --ri --cluster "qsub -l nodes=1:ppn=20 -l walltime=999:00:00" >sublog 2>&1 & 
  ```
 ## Support tools 
 * BUSCO 
 * Merquery
 ## contribution 
 If you want to apply other tools to evaluate the genome, we encourage you pull a request or email us 
 
 ## Contact  
 * Peng Jia (pengjia@stu.xjtu.edu.cn)

