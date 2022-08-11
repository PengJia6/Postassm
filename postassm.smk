# ======================================================================================================================
# Project:  PostAssm
# Script :  postassm.smk
# Author : Peng Jia
# Date   : 2021.07.06
# Email  : pengjia@stu.xjtu.edu.cn
# Description: post analysis of assembly
# ======================================================================================================================
configfile: "configs.yaml"
import pandas as pd
import os


# software and work path

proect_name =  config["project_name"] if "project_name" in  config else "my_postassm"
dir_work = "postassm_workdir/" if "dir_work" not in config else config["dir_work"]
dir_work = dir_work + "/" if dir_work[-1] != "/" else dir_work
busco_conf = config["lib_dir"]["busco_conf"]
busco = config["software"]["busco"]
merquery =  config["software"]["merqury"]

Rscript = config["software"]["Rscript"] if "Rscript" in config["software"] else "Rscript"
xf_stat = config["software"]["xf_stat"] if "xf_stat" in config["software"] else "xf_stat"
samtools = config["software"]["samtools"] if "samtools" in config["software"] else "samtools"
seqtk = config["software"]["seqtk"] if "seqtk" in config["software"] else "seqtk"
samples = [sample for sample, info in config["samples"].items()]

wildcard_constraints:
    sample="|".join(samples)


################################################## final output ########################################################
#
rule all:
    input:
        dir_work + f"{proect_name}.continuity.csv",
        dir_work + f"{proect_name}.completeness.csv",
        dir_work + f"{proect_name}.accuracy.csv",
        dir_work+ f"{proect_name}.stat.csv",


def get_fa(wildcards):
    return config["samples"][wildcards.sample]["assm"]


rule copy_fa:
    input:
        get_fa
    output:
        dir_work + "fasta/{sample}.fa"
    run:
        shell("cp {input} {output}")


################################################## completeness ########################################################

def get_all_busco(wildcards):
    all_stats = []
    for sample, info in config["samples"].items():
        all_stats.append(dir_work + f"busco/{sample}.busco")
    return all_stats


rule busco:
    input:
        ref=dir_work + "fasta/{sample}.fa",
        busco_lib=config["lib_dir"]["busco"]
    output:
        dir_work + "busco/{sample}.busco"
    log: dir_work + "logs/{sample}.busco"
    threads: config["threads"]["busco"] if "busco" in config["threads"] else config["default"]
    run:
        import os

        out_dir = "/".join(str(output).split("/")[:-1])
        out_tmp = str(output).split("/")[-1] + "_dir"
        busco_db = f"{input.busco_lib}"
        if os.path.exists(str(output) + "_dir"):
            shell("rm -rf {output}_dir")
        busco_prefix = "/".join(f"{busco}".split("/")[:-1])

        shell(
            "cd {out_dir} && "
            " export BUSCO_CONFIG_FILE={busco_conf}&&"
            " export PATH={busco_prefix}:$PATH &&"
            " {busco} -o {out_tmp} -i {input.ref} -l {busco_db} -m genome -c {threads} 2>{log}")
        shell("cp {output}_dir/short_summary.specific.mammalia_odb10.{out_tmp}.txt {output} ")

rule merge_all_busco:
    input:
        busco=get_all_busco,
    output:
        dir_work + f"{proect_name}.busco.csv",
    run:
        output_data = pd.DataFrame()
        for item in input.busco:
            name = item.split("/")[-1][:-6]
            busco_info = [i.lstrip() for i in open(item).readlines() if i.lstrip().startswith("C")][0]
            busco_C = busco_info.split(":")[1].split("%")[0]
            busco_S = busco_info.split(":")[2].split("%")[0]
            busco_D = busco_info.split(":")[3].split("%")[0]
            busco_F = busco_info.split(":")[4].split("%")[0]
            busco_M = busco_info.split(":")[5].split("%")[0]
            output_data.loc[name, "busco_C"] = busco_C
            output_data.loc[name, "busco_S"] = busco_S
            output_data.loc[name, "busco_D"] = busco_D
            output_data.loc[name, "busco_F"] = busco_F
            output_data.loc[name, "busco_M"] = busco_M
        output_data.to_csv(str(output),sep="\t")


def get_all_merqury(wildcards):
    all_stats = []
    for sample, info in config["samples"].items():
        all_stats.append(dir_work + f"merqury/{sample}.merqury")
    return all_stats


def get_merdb(wildcards):
    return config["samples"][wildcards.sample]["meryl_lab"]


rule merqury:
    input:
        fa=dir_work + "fasta/{sample}.fa",
        db=get_merdb
    output:
        dir_work + "merqury/{sample}.merqury"
    threads: config["threads"]["merqury"] if "merqury" in config["threads"] else config["threads"]["default"]
    priority: 50
    run:
        output_dir = "/".join(str(output).split("/")[:-1])
        output_pre = str(output).split("/")[-1][:-3]
        software_pre = merquery[:-11]
        if not os.path.exists(str(input) + "sta"):
            shell("ln -fs {input.fa} {output_dir}/assembly.fasta")
        shell("export PATH=/data/home/pengjia/mybin:$PATH && export MERQURY={software_pre} && "
              "cd {output_dir} && {merquery} {input.db} assembly.fasta {output_pre}")
        # print("Run End")
        # print("export PATH=/data/home/pengjia/mybin:$PATH && export MERQURY={software_pre} && "
        #       "cd {output_dir} && {merquery} {input.db} assembly.fasta {output_pre}")
        # print("software_pre",software_pre)
        # print("output_dir ",output_dir)
        # print("input.db",input.db)
        # print("output_pre",output_pre)
        print(input.fa)

        shell("touch {output}")

rule merge_all_merqury:
    input:
        merquery=get_all_merqury,
    output:
        dir_work + f"{proect_name}.merqury.csv",
    run:
        import pandas as pd

        output_data = pd.DataFrame()
        for item in input.merquery:
            name = item.split("/")[-1][:-11]
            output_data.loc[name, "merquery_qv"] = open(item[:-2] + "qv").read().split("\t")[3]
            output_data.loc[name, "merquery_C"] = open(item[:-2] + "completeness.stats").read().split("\t")[4][:-1]
        output_data.to_csv(str(output),sep="\t")

rule merged_completenessy:
    input:
        merqury=dir_work + "{proect_name}.merqury.csv",
        busco=dir_work + "{proect_name}.busco.csv",
    output:
        dir_work + "{proect_name}.completeness.csv",
    run:
        merqury_csv = pd.read_csv(f"{input.merqury}",index_col=0)
        busco_csv = pd.read_csv(f"{input.busco}",index_col=0)
        merqury_csv = merqury_csv["merquery_C"]
        data = pd.concat([busco_csv, merqury_csv])
        data.to_csv(f"{output}")


################################################## continuity ########################################################

def get_all_gaps(wildcards):
    all_stats = []
    for sample, info in config["samples"].items():
        all_stats.append(dir_work + f"gaps/{sample}.gaps.csv")
    return all_stats


rule merge_gaps:
    input:
        get_all_gaps
    output:
        dir_work + f"{proect_name}.gaps.csv",
    run:
        import pandas as pd

        output_data = pd.DataFrame()
        for item in input:
            name = item.split("/")[-1][:-9]
            output_data.loc[name, "gaps_num"] = len(open(item).readlines())
        output_data.to_csv(str(output))

rule get_gaps:
    input:
        fa=dir_work + "fasta/{sample}.fa",
    output:
        dir_work + "gaps/{sample}.gaps.csv",
    run:
        shell("{seqtk} cutN -n1 -g {input.fa} >{output}")


def get_all_NG50_stat(wildcards):
    all_stats = []
    for sample, info in config["samples"].items():
        all_stats.append(dir_work + f"N50/{sample}.N50.csv")
    return all_stats


rule merge_all_NG50:
    input:
        stat=get_all_NG50_stat,
    output:
        dir_work + "{proect_name}.N50.csv",

    priority: 50
    run:
        import pandas as pd

        output_data = pd.DataFrame()
        for item in input.stat:
            name = item.split("/")[-1][:-8]
            this_data = pd.read_table(item)
            print(this_data)
            output_data.loc[name, "Assembly_length"] = str(int(this_data.loc[0, "length"]))
            output_data.loc[name, "Scaffold_number"] = str(int(this_data.loc[0, "rank"]))
            output_data.loc[name, "Scaffold_N50"] = str(int(this_data.loc[5, "length"]))
            output_data.loc[name, "Scaffold_N90"] = str(int(this_data.loc[9, "length"]))
            output_data.loc[name, "Contig_number"] = str(int(this_data.loc[10, "rank"]))
            output_data.loc[name, "Config_N50"] = str(int(this_data.loc[15, "length"]))
            output_data.loc[name, "Config_N90"] = str(int(this_data.loc[19, "length"]))
        output_data.to_csv(str(output))

rule get_scaffoldN50:
    input:
        fa=dir_work + "fasta/{sample}.fa",
        fai=dir_work + "fasta/{sample}.fa",
    output:
        dir_assm=dir_work + "N50/{sample}.N50.csv"
    run:
        tmp = str(output.dir_assm).rstrip(".stat") + "_tmp"
        shell("mkdir -p {tmp}")
        shell("{seqtk} cutN -n 10 {input.fa} >{tmp}/contigs.fa")
        shell("{samtools} faidx  {tmp}/contigs.fa")
        shell("cp {input.fa} {tmp}/scaffolds.fa")
        shell("{samtools} faidx {tmp}/scaffolds.fa")
        shell("cd {tmp} && {Rscript} {xf_stat} scaffolds.fa.fai contigs.fa.fai")
        shell("mv {tmp}/assembly_Statistics.stat {output.dir_assm}")

rule merged_continuity:
    input:
        gaps=dir_work + "{proect_name}.N50.csv",
        N50=dir_work + "{proect_name}.gaps.csv",
    output:
        dir_work + "{proect_name}.continuity.csv",
    run:
        gaps_csv = pd.read_csv(f"{input.gaps}",index_col=0)
        N50_csv = pd.read_csv(f"{input.N50}",index_col=0)
        data = pd.concat([gaps_csv, N50_csv],axis=1)

        data.to_csv(f"{output}")

################################################## Accuracy  ########################################################

rule merged_accuracy:
    input:
        merqury=dir_work + "{proect_name}.merqury.csv",
    output:
        dir_work + "{proect_name}.accuracy.csv",
    run:
        merqury_csv = pd.read_csv(f"{input.merqury}",index_col=0)
        merqury_csv = merqury_csv["merquery_qv"]
        merqury_csv.to_csv(f"{output}")

################################################## Other rules ########################################################
rule merge_all_stats:
    input:
        dir_work + f"{proect_name}.continuity.csv",
        dir_work + f"{proect_name}.completeness.csv",
        dir_work + f"{proect_name}.accuracy.csv",
    output:
        dir_work + f"{proect_name}.stat.csv",
    run:
        data = pd.DataFrame()
        for item in input:
            sub_data = pd.read_csv(item,index_col=0)
            data = pd.concat[data, sub_data]
        data.to_csv(f"{output}")


rule fasta_index:
    input:
        "{prefix}.fa"
    output:
        "{prefix}.fa.fai"
    run:
        shell("{samtools} faidx {input}")

rule fasta_index2:
    input:
        "{prefix}.fasta"
    output:
        "{prefix}.fasta.fai"
    run:
        shell("{samtools} faidx {input}")
