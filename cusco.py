#!/usr/bin/env python3
# ver.0.1.1

import argparse
from Bio import Phylo
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from collections import defaultdict
import copy
import csv
from ete3 import Tree
import glob
import itertools
from io import StringIO
import multiprocessing
import os
import pandas as pd
from pathlib import Path
import re
import subprocess
import sys


'''prepare input file for cureatscg'''
def wgs2gff(threads, ref_cds, w, outsthreshold, input_dir):
    if not os.path.exists(f"{input_dir}/gff"):
        os.makedirs(f"{input_dir}/gff")
    if not os.path.exists(f"{input_dir}/paf"):
        os.makedirs(f"{input_dir}/paf")

    wgs_files = [f for f in sorted(os.listdir(f'{input_dir}/species')) if f.endswith(('.fasta', '.fas', '.fa', '.fna', '.seq', '.aln', '.maf', '.nt'))] # List of WGS files

    for wgs_file in wgs_files:
        wgs_path = os.path.join(f"{input_dir}/species", wgs_file)
        sample_name = os.path.splitext(os.path.basename(wgs_path))[0]
        index_file = f"{wgs_path}.fai"

        if not os.path.exists(index_file): # If the index file does not exist, run samtools faidx
            subprocess.run(['samtools', 'faidx', wgs_path], check=True)

        if w == True: # when polyploid or WGD sample in ./species, extract number of ploidy form threshold file
            threshold_file = None
            for fname in os.listdir(f'{input_dir}'):
                if fname.lower() in ['threshold.tsv', 'threshold.csv']:
                    threshold_file = os.path.join(input_dir, fname)
                    break
            if threshold_file is None:
                raise FileNotFoundError("There is no threshold file.")
            ext = os.path.splitext(threshold_file)[1].lower()
            delimiter = '\t' if ext == '.tsv' else ','
            # read threshold.csv or threshold.tsv
            threshold_df = pd.read_csv(threshold_file, sep=delimiter, header=None, index_col=0)
            threshold_df = threshold_df[1].astype(int)

            if sample_name in threshold_df: # Annotate with low float only ploidy genome
                if threshold_df[sample_name] == 1:
                    cmd_wgs2gff = f"miniprot --gff --outs=0.97 --no-cs -Iut {threads} {wgs_path} {input_dir}/{ref_cds} > {input_dir}/{sample_name}.out"
                    subprocess.run(cmd_wgs2gff, shell=True, check=True)  # Run miniprot to GFF and PAF

                else:
                    cmd_wgs2gff = f"miniprot --gff --outs={outsthreshold} --no-cs -Iut {threads} {wgs_path} {input_dir}/{ref_cds} > {input_dir}/{sample_name}.out"
                    subprocess.run(cmd_wgs2gff, shell=True, check=True)  # Run miniprot to GFF and PAF

            else:
                print(f"{sample_name} not found in threshold list.")

        else:
            cmd_wgs2gff = f"miniprot --gff --outs=0.97 --no-cs -Iut {threads} {wgs_path} {input_dir}/{ref_cds} > {input_dir}/{sample_name}.out"
            subprocess.run(cmd_wgs2gff, shell=True, check=True)  # Run miniprot to GFF and PAF

        paf_file, gff_file = f"{input_dir}/{sample_name}.paf", f"{input_dir}/{sample_name}.gff"

        with open(f"{input_dir}/{sample_name}.out", "r", encoding="utf-8") as infile, \
             open(paf_file, "w", encoding="utf-8") as file_paf, \
             open(gff_file, "w", encoding="utf-8") as file_gff:

            for line in infile:
                if line.startswith("##PAF\t"):
                    file_paf.write(line[6:])
                else:
                    file_gff.write(line)

        cmd_mv1 = f"mv {input_dir}/{sample_name}.gff {input_dir}/gff"
        subprocess.run(cmd_mv1, shell=True, check=True)

        cmd_mv2 = f"mv {input_dir}/{sample_name}.paf {input_dir}/paf"
        subprocess.run(cmd_mv2, shell=True, check=True)

        cmd_rm = f"rm {input_dir}/{sample_name}.out"
        subprocess.run(cmd_rm, shell=True, check=True)

def insert_id(args):
    sample_name, file_path = args
    cmd_iid = f"sed -i -E 's/^>(MP[0-9]+)/>{sample_name}_\\1/g' {file_path}"
    subprocess.run(cmd_iid, shell=True, check=True)  # Run sed

def process_gff2fna(args_na):
    gff_file, input_dir = args_na
    gff_path = os.path.join(f"{input_dir}/gff", gff_file)
    sample_name = gff_file.split('.gff')[0]

    extensions = ['.fasta', '.fas', '.fa', '.fna', '.seq', '.aln', '.maf', '.nt']
    for ext in extensions:
        fasta_files = glob.glob(f"{input_dir}/species/{sample_name}{ext}")
        if fasta_files:
            fasta_file = fasta_files[0]
            break
    else:
        raise FileNotFoundError(f"No FASTA file found for {sample_name} in {input_dir}/species/")

    cmd_fna = f"gffread -x {sample_name}.fna -g {fasta_file} {gff_path}"
    subprocess.run(cmd_fna, shell=True, check=True)  # Run gffread

    file_path = sample_name + ".fna"
    insert_id((sample_name, file_path))

    cmd_mv = f"mv {sample_name}.fna {input_dir}/input"
    subprocess.run(cmd_mv, shell=True, check=True)

def process_gff2faa(args_aa):
    gff_file, input_dir = args_aa
    gff_path = os.path.join(f"{input_dir}/gff", gff_file)
    sample_name = gff_file.split('.gff')[0]

    extensions = ['.fasta', '.fas', '.fa', '.fna', '.seq', '.aln', '.maf', '.nt']
    for ext in extensions:
        fasta_files = glob.glob(f"{input_dir}/species/{sample_name}{ext}")
        if fasta_files:
            fasta_file = fasta_files[0]
            break
    else:
        raise FileNotFoundError(f"No FASTA file found for {sample_name} in {input_dir}/species/")

    cmd_faa = f"gffread -y {sample_name}.faa -g {fasta_file} {gff_path}"
    subprocess.run(cmd_faa, shell=True, check=True)  # Run gffread

    file_path = sample_name + ".faa"
    insert_id((sample_name, file_path))

    cmd_mv = f"mv {sample_name}.faa {input_dir}/input"
    subprocess.run(cmd_mv, shell=True, check=True)

def gff2fasta(threads, p, input_dir):
    if not os.path.exists(f"{input_dir}/input"):
        os.makedirs(f"{input_dir}/input")
    gff_files = [f for f in sorted(os.listdir(f'{input_dir}/gff')) if f.endswith('.gff')] # List of GFF files

    if p == True:
        with multiprocessing.Pool(threads) as pool_aa: # processing that convert gff to fna
            args_aa = [(gff_file, input_dir) for gff_file in gff_files]
            pool_aa.map(process_gff2faa, args_aa)
    else:
        with multiprocessing.Pool(threads) as pool_na: # processing that convert gff to fna
            args_na = [(gff_file, input_dir) for gff_file in gff_files]
            pool_na.map(process_gff2fna, args_na)

def paf2tsv(input_dir): # make CNVtable from paf
    cmd_gfa = f"./pangene-1.1-bin/bin_x64-linux/pangene {input_dir}/paf/*paf > {input_dir}/graph.gfa"
    subprocess.run(cmd_gfa, shell=True, check=True)  # Run pangene

    cmd_tsv = f"./pangene-1.1-bin/bin_x64-linux/k8 ./pangene-1.1-bin/scripts/pangene.js gfa2matrix -c {input_dir}/graph.gfa > {input_dir}/pre_geneCNV.Rtab"
    subprocess.run(cmd_tsv, shell=True, check=True)  # Run pangene

    with open(f"{input_dir}/pre_geneCNV.Rtab", "r", encoding="utf-8") as infile, \
         open(f"{input_dir}/geneCNV.Rtab", "w", encoding="utf-8", newline="") as outfile:

        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        for row_index, row in enumerate(reader):
            if row_index == 0:
                row[1:] = [col[:-2] if col.endswith("#0") else col for col in row[1:]]
            writer.writerow(row)

    cmd_rm = f"rm {input_dir}/pre_geneCNV.Rtab"
    subprocess.run(cmd_rm, shell=True, check=True)


'''cureate single copy genes'''
'''1st step of cureate single copy genes'''
def tsv2scg(w, input_dir): # first step : extract single copy genes from Rtab table
    ori_Rtab = pd.read_csv(f"{input_dir}/geneCNV.Rtab", sep='\t')

    if w == True: # when polyploid or WGD sample in ./species, extract multiple copy genes for the sample
        threshold_file = None
        for fname in os.listdir(f'{input_dir}'):
            if fname.lower() in ['threshold.tsv', 'threshold.csv']:
                threshold_file = os.path.join(input_dir, fname)
                break
        if threshold_file is None:
            raise FileNotFoundError("There is no threshold file.")
        ext = os.path.splitext(threshold_file)[1].lower()
        delimiter = '\t' if ext == '.tsv' else ','

        # read threshold and cast to int
        #threshold_df = pd.read_csv(threshold_file, sep=delimiter, index_col=0).squeeze().astype(int) # read threshold.csv or threshold.tsv
        threshold_df = pd.read_csv(threshold_file, sep=delimiter, header=None, index_col=0)
        threshold_df = threshold_df[1].astype(int)
        #print(threshold_df)

        # gene copy numbers only
        data_only = ori_Rtab.iloc[:, 1:].copy()

        # check if all threshold species are in Rtab
        if not all(species in data_only.columns for species in threshold_df.index):
            missing = [s for s in threshold_df.index if s not in data_only.columns]
            raise ValueError(f"Missing species in geneCNV.Rtab: {missing}")

        # align both in same order and cast data_only to int
        data_filtered = data_only[threshold_df.index].astype(int)
        match_df = data_filtered.eq(threshold_df.values, axis=1)

        filtered_Rtab = ori_Rtab[match_df.all(axis=1)]
        pre_scg_step1 = filtered_Rtab.iloc[:, 0].to_list()

        #print(f"Filtered {len(pre_scg_step1)} genes based on threshold file.")
        #print(filtered_Rtab[threshold_df.index].head())  # <-- デバッグ表示用
        #filtered_Rtab.to_csv('test.csv')

    else: # extract all 1 copy genes
        filtered_Rtab = ori_Rtab[ori_Rtab.iloc[:, 1:].eq(1).all(axis=1)]
        pre_scg_step1 = filtered_Rtab.iloc[:, 0].to_list()
        print(f"Filtered {len(pre_scg_step1)} genes with all copy number = 1.")

    return pre_scg_step1

'''2nd step of cureate single copy genes'''
def list_to_df_and_find_duplicates(input_list):
    find_dup_df = pd.DataFrame(input_list, columns=["Column1"])
    duplicates_list = find_dup_df[find_dup_df.duplicated(subset="Column1", keep=False)]["Column1"].unique().tolist()

    return duplicates_list

def process_gff_file(args1):
    gff_file, pre_scg, input_dir = args1
    results = []
    file_name = gff_file.split('.gff')[0]  # extract sample name
    gff_path = f'{input_dir}/gff/' + gff_file

    with open(gff_path, 'r') as gff_input:
        for line in gff_input:
            if not line.strip():  # skip empty row
                continue

            columns = re.split(r'\s+', line.strip())

            if len(columns) > 8 and columns[2] == "mRNA":  # mRNA in third column
                col_9 = columns[8]
                target_match = re.search(r"Target=([^;]+)", col_9)  # check original gene name
                target_value = target_match.group(1) if target_match else ""

                if target_value in pre_scg:  # check whether the original gene name is a single copy gene
                    results.append(f"{target_value}")

    return results

def gff2scg(pre_scg, threads, input_dir): #second step : trimming single copy genes of first step using GFF files
    gff_files = [f for f in sorted(os.listdir(f'{input_dir}/gff')) if f.endswith('.gff')] # List of GFF files

    with multiprocessing.Pool(threads) as pool: # processing that find duplicate genes with multicore
        args1 = [(gff_file, pre_scg, input_dir) for gff_file in gff_files]
        results = pool.map(process_gff_file, args1)

    duplist1, duplist2 = [], []

    for result in results:
        X = list_to_df_and_find_duplicates(result)
        duplist1 = duplist1 + X
        duplist2 = duplist1 + result

    duplist3 = list_to_df_and_find_duplicates(duplist2)
    scg_step2 = [i for i in pre_scg if i not in duplist3]
    scg_step2.sort()

    return scg_step2


def process_gff_file_df(args2):
    gff_file, pre_scg2, scg_orthologs1, input_dir = args2
    results2 = []
    file_name = gff_file.split('.gff')[0]  # extract sample name
    gff_path = f'{input_dir}/gff/' + gff_file

    with open(gff_path, 'r') as gff_input:
        for line in gff_input:
            if not line.strip():  # skip empty row
                continue

            columns = re.split(r'\s+', line.strip())

            if len(columns) > 8 and columns[2] == "mRNA":  # mRNA in third column
                col_9 = columns[8]
                target_match = re.search(r"Target=([^;]+)", col_9)  # check original gene name
                target_value = target_match.group(1) if target_match else ""
                id_match = re.search(r"ID=([^;]+)", col_9) # check orthologs name
                id_value = id_match.group(1) if id_match else ""

                if target_value in pre_scg2:  # check whether the original gene name is a single copy gene
                    if target_value and id_value:
                        results2.append(f"{target_value}\t{id_value}")

        ori_ortho = [line.split("\t") for line in results2] # list to df
        ori_ortho_df = pd.DataFrame(ori_ortho, columns=["reference_gene", "file_name"])

    merged_ori_ortho_df = pd.merge(scg_orthologs1, ori_ortho_df, on="reference_gene", how='right') # merge df
    merged_ori_ortho_df = merged_ori_ortho_df[["file_name"]].fillna('')
    merged_ori_ortho_df.columns = [file_name]

    return merged_ori_ortho_df

def gff2scg_df(pre_scg2, threads, input_dir): # make orthologs table
    gff_files = [f for f in sorted(os.listdir(f'{input_dir}/gff')) if f.endswith('.gff')] # List of GFF files
    scg_orthologs1 = pd.DataFrame(pre_scg2, columns=["reference_gene"]) # main df

    with multiprocessing.Pool(threads) as pool2: # Processing that create merge df with multicore
        args2 = [(gff_file, pre_scg2, scg_orthologs1, input_dir) for gff_file in gff_files]
        result_dfs = pool2.map(process_gff_file_df, args2)

        final_df = scg_orthologs1.copy() # merge all df
        for df in result_dfs:
            final_df = pd.concat([final_df, df], axis=1)

        #final_df = final_df.sort_index(axis=1) # sort df

    return final_df


def process_pep_file(args3):
    ref_gene, scg_orthologs1, input_dir, output_dir = args3
    pep_files = [f for f in sorted(os.listdir(f'{input_dir}/input')) if f.endswith(('.fasta', '.fas', '.fa', '.fna', '.faa', '.pep', '.seq', '.aln', '.maf', '.nt', '.aa'))]
    temp_files = []

    if ref_gene in scg_orthologs1["reference_gene"].values:
        target_row = scg_orthologs1[scg_orthologs1["reference_gene"] == ref_gene]

        for col in scg_orthologs1.columns[1:]:
            ortho_gene = target_row[col].values[0] # extract orthologs name
            file_name = col # extract sample name

            pep_files = glob.glob(os.path.join(f"{input_dir}/input", f"{file_name}.*"))

            if pep_files:
                for pep_path in pep_files:
                    file_extension = os.path.splitext(pep_path)[1] # get extension
                    if os.path.exists(pep_path):
                        output_pep = f"{ref_gene}_{file_name}{file_extension}"
                        temp_files.append(output_pep)

                        cmd = [
                            "seqkit", "grep", "-n", "-r", "-p", ortho_gene, pep_path, "-o", output_pep
                        ]
                        subprocess.run(cmd, check=True) # Run the seqkit grep
            else:
                print(f"No files found = {pep_files} in {input_dir}/input")

        combined_pep = os.path.join(f"{output_dir}/single_copy_genes", f"{ref_gene}.fa")
        with open(combined_pep, 'w') as out_f: # merge temporary files
            for temp_file in temp_files:
                with open(temp_file, 'r') as tf:
                    out_f.write(tf.read())

        for temp_file in temp_files: # remove temporary files
            os.remove(temp_file)

def fasta2scg(pre_scg2, scg_orthologs1, threads, input_dir, output_dir): # Create single copy gene FASTA
    if not os.path.exists(f"{output_dir}/single_copy_genes"): # create output directory
        os.makedirs(f"{output_dri}/single_copy_genes")
    ref_genes = scg_orthologs1["reference_gene"].tolist() # reference genes

    with multiprocessing.Pool(threads) as pool3: # Processing that create singel copy gene fasta with multicore
        args3 = [(ref_gene, scg_orthologs1, input_dir, output_dir) for ref_gene in ref_genes]
        pool3.map(process_pep_file, args3)


### The following steps for polyploid or whole geneme dupulication ###
def load_expected_counts(input_dir):  # Read expected count from threshold.(tsv|csv)
    candidates = [(f"{input_dir}/threshold.tsv", '\t'), (f"{input_dir}/threshold.csv", ',')]

    for filepath, sep in candidates:
        if os.path.exists(filepath):
            df = pd.read_csv(filepath, sep=sep, header=None, names=["sample", "expected_count"])
            return dict(zip(df["sample"], df["expected_count"]))

    raise FileNotFoundError("Neither 'threshold.tsv' nor 'threshold.csv' was found in the current directory.")

def extract_strict_scg(results_dict, expected_counts_dict): # Extract only those whose SCG occurrence count matches the expected value
    sample_valid_scg_sets = []

    for sample, scg_list in results_dict.items():
        expected_count = expected_counts_dict.get(sample, None)
        if expected_count is None:
            continue

        scg_counter = Counter(scg_list)
        valid_scg = {scg for scg, count in scg_counter.items() if count == expected_count}
        sample_valid_scg_sets.append(valid_scg)

    if sample_valid_scg_sets:
        strict_scg = set.intersection(*sample_valid_scg_sets)
    else:
        strict_scg = set()

    return sorted(strict_scg)

def gff2scg_w(pre_scg, threads, input_dir):
    gff_files = [f for f in sorted(os.listdir(f'{input_dir}/gff')) if f.endswith('.gff')]  # List of GFF files

    with multiprocessing.Pool(threads) as pool:
        args1 = [(gff_file, pre_scg, input_dir) for gff_file in gff_files]
        results_list = pool.map(process_gff_file, args1)

    # Combine the GFF file name (and extract sample names from it) and the result list to convert it into a dictionary
    sample_names = [f.split('.gff')[0] for f in gff_files]
    results_dict = dict(zip(sample_names, results_list))

    return results_dict


def process_gff_file_df_w(args2):
    gff_file, pre_scg2, scg_orthologs1, input_dir = args2
    results2 = []
    file_name = gff_file.split('.gff')[0]
    gff_path = f'{input_dir}/gff/' + gff_file

    with open(gff_path, 'r') as gff_input:
        for line in gff_input:
            if not line.strip():
                continue

            columns = re.split(r'\s+', line.strip())

            if len(columns) > 8 and columns[2] == "mRNA":
                col_9 = columns[8]
                target_match = re.search(r"Target=([^;]+)", col_9)
                id_match = re.search(r"ID=([^;]+)", col_9)
                rank_match = re.search(r"Rank=([^;]+)", col_9)

                target_value = target_match.group(1) if target_match else ""
                id_value = id_match.group(1) if id_match else ""
                rank_value = rank_match.group(1) if rank_match else ""

                if target_value in pre_scg2 and id_value and rank_value:
                    column_name = f"{file_name}_{rank_value}"
                    results2.append([target_value, column_name, id_value])

    if not results2:
        return pd.DataFrame()  # If there is no data, it returns an empty DataFrame

    df = pd.DataFrame(results2, columns=["reference_gene", "column_name", "ortholog"])
    pivot_df = df.pivot(index="reference_gene", columns="column_name", values="ortholog")
    pivot_df = pivot_df.reset_index()
    merged_df = pd.merge(scg_orthologs1, pivot_df, on="reference_gene", how="right")

    return merged_df

def gff2scg_df_w(strict_scg, threads, input_dir):
    gff_files = [f for f in sorted(os.listdir(f'{input_dir}/gff')) if f.endswith('.gff')]
    scg_orthologs1 = pd.DataFrame(strict_scg, columns=["reference_gene"])

    with multiprocessing.Pool(threads) as pool2:
        args2 = [(gff_file, strict_scg, scg_orthologs1, input_dir) for gff_file in gff_files]
        result_dfs = pool2.map(process_gff_file_df_w, args2)

    final_df = scg_orthologs1.copy()
    for df in result_dfs:
        if not df.empty:
            final_df = pd.merge(final_df, df, on="reference_gene", how="outer")

    final_df = final_df.fillna('')  # Fill gaps with nulls (if necessary)

    return final_df


def resolve_sample_names(columns):
    prefixes = [re.sub(r'_\d+$', '', col) for col in columns]
    prefix_counts = Counter(prefixes)
    col_to_sample = {}
    for col in columns:
        prefix = re.sub(r'_\d+$', '', col)
        if prefix_counts[prefix] == 1:
            col_to_sample[col] = prefix
        else:
            col_to_sample[col] = prefix  # Use prefix for search
    return col_to_sample

def process_pep_file_w(args3):
    ref_gene, scg_orthologs1, input_dir, output_dir = args3
    temp_files = []

    if ref_gene not in scg_orthologs1["reference_gene"].values:
        return

    target_row = scg_orthologs1[scg_orthologs1["reference_gene"] == ref_gene]
    col_to_sample_prefix = resolve_sample_names(scg_orthologs1.columns[1:])

    # Record the column number for each sample to keep track of the indexes of samples with the same name
    sample_col_indices = {}
    for idx, col in enumerate(scg_orthologs1.columns[1:], start=1):
        prefix = col_to_sample_prefix[col]
        if prefix not in sample_col_indices:
            sample_col_indices[prefix] = []
        sample_col_indices[prefix].append(col)

    for col in scg_orthologs1.columns[1:]:
        ortho_gene = target_row[col].values[0]
        sample_prefix = col_to_sample_prefix[col]  # e.g., "samplename"
        sample_fullname = col                     # e.g., "samplename_2"

        pep_files = glob.glob(os.path.join(f"{input_dir}/input", f"{sample_prefix}.*"))

        if pep_files:
            # Index assignment: Only when there are multiple identical prefixes
            suffix_index = ""
            if len(sample_col_indices[sample_prefix]) > 1:
                col_index = sample_col_indices[sample_prefix].index(col) + 1
                suffix_index = f"_{col_index}"

            for pep_path in pep_files:
                file_extension = os.path.splitext(pep_path)[1]
                output_pep = f"{sample_fullname}_{ortho_gene}{file_extension}"
                temp_files.append(output_pep)

                # Filtering to change FASTA IDs
                cmd = ["seqkit", "grep", "-n", "-r", "-p", ortho_gene, pep_path]
                result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True)

                modified_output = []
                for line in result.stdout.splitlines():
                    if line.startswith(">"):
                        new_id = f">{sample_prefix}{suffix_index}_{ortho_gene}"
                        modified_output.append(new_id)
                    else:
                        modified_output.append(line)

                with open(output_pep, 'w') as f:
                    f.write("\n".join(modified_output) + "\n")
        else:
            print(f"[Warning] No files found for sample '{sample_prefix}' in {input_dir}/input")

    if temp_files:
        combined_pep = os.path.join(f"{output_dir}/single_copy_genes", f"{ref_gene}.fa")
        with open(combined_pep, 'w') as out_f:
            for temp_file in temp_files:
                with open(temp_file, 'r') as tf:
                    out_f.write(tf.read())

        for temp_file in temp_files:
            os.remove(temp_file)

def fasta2scg_w(pre_scg2, scg_orthologs1, threads, input_dir, output_dir): # Create single copy gene FASTA
    if not os.path.exists(f"{output_dir}/single_copy_genes"): # create output directory
        os.makedirs(f"{output_dir}/single_copy_genes")
    ref_genes = scg_orthologs1["reference_gene"].tolist() # reference genes

    with multiprocessing.Pool(threads) as pool3: # Processing that create singel copy gene fasta with multicore
        args3 = [(ref_gene, scg_orthologs1, input_dir, output_dir) for ref_gene in ref_genes]
        pool3.map(process_pep_file_w, args3)
### up to here ###

'''3rd step of cureate single copy genes'''
def count_sequences_in_fasta(file_path):
    try:
        result = subprocess.run(
            ["seqkit", "fx2tab", "-n", file_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        # Count lines in the output
        return len(result.stdout.strip().split("\n"))
    except subprocess.CalledProcessError as e:
        print(f"Error processing {file_path}: {e.stderr}")
        return None

def last_trimming(pre_scg2, scg_orthologs1, output_dir): # third step : trimming single copy genes of second step using FASTA
    scg_orthologs1 = scg_orthologs1.set_index("reference_gene")
    sample_num = len(scg_orthologs1.columns)
    fasta_files = [f for f in os.listdir(f"{output_dir}/single_copy_genes") if f.endswith(".fa")]

    for fasta_file in fasta_files:
        fasta_path = os.path.join(f"{output_dir}/single_copy_genes", fasta_file)
        seq_count = count_sequences_in_fasta(fasta_path) # Count the number of sequences in a FASTA file using seqkit
        if seq_count is not None and seq_count != sample_num:
            print(f"{fasta_file} has {seq_count} sequences (not {sample_num})")
            gene_name = fasta_file.split('.fa')[0]
            scg_orthologs1 = scg_orthologs1.drop(index=gene_name)
            os.remove(fasta_path)

    return scg_orthologs1

'''formatting fasta'''
def eject_id(alignment_path):
    cmd_eid = f"sed -i -E 's/^(>.*)_([^_]+)$/\\1/' {alignment_path}"
    subprocess.run(cmd_eid, shell=True, check=True)  # Run sed

def process_alignment(args_mafft): # Run mafft
    fasta_file, pip, output_dir =  args_mafft
    fasta_path = os.path.join(f"{output_dir}/single_copy_genes", fasta_file)
    alignment_path = fasta_path.replace(".fa", ".fasta")

    cmd_mafft = f"mafft --quiet {fasta_path} > {alignment_path}"
    subprocess.run(cmd_mafft, shell=True, check=True)  # Run mafft

    if pip:
        eject_id(alignment_path)

    cmd_rm = f"rm {fasta_path}"
    subprocess.run(cmd_rm, shell=True, check=True)

def alignment_fasta(threads, pip, output_dir):
    fasta_files = [f for f in os.listdir(f"{output_dir}/single_copy_genes") if f.endswith(".fa")]

    with multiprocessing.Pool(threads) as poolmafft: # Processing that create singel copy gene fasta with multicore
        args_mafft = [(fasta_file, pip, output_dir) for fasta_file in fasta_files]
        poolmafft.map(process_alignment, args_mafft)


'''make phylogenetic tree'''
def process_genetree(args_gt): # Run iqtree
    fasta_file, sub_model, bootstrap, input_dir, output_dir =  args_gt

    fasta_path = os.path.join(f"{output_dir}/single_copy_genes", fasta_file)
    gene_name = fasta_file.split('.fasta')[0]

    if bootstrap == 0:
        cmd1 = f"iqtree2 -s {fasta_path} -T AUTO --quiet -m {sub_model} --prefix {gene_name}"
    else:
        cmd1 = f"iqtree2 -s {fasta_path} -T AUTO --quiet -m {sub_model} --ufboot {bootstrap} --bnni --prefix {gene_name}"
    subprocess.run(cmd1, shell=True, check=True)  # Run iqtree

    cmd2 = f"mv {gene_name}.* {input_dir}/gene_tree"
    subprocess.run(cmd2, shell=True, check=True)

def make_gene_tree(threads, sub_model, bootstrap, input_dir, output_dir): # Make gene tree for species tree
    if not os.path.exists(f"{input_dir}/gene_tree"):
        os.makedirs(f"{input_dir}/gene_tree")

    fasta_files = [f for f in os.listdir(f"{output_dir}/single_copy_genes") if f.endswith(".fasta")]
    fasta_files.sort()

    with multiprocessing.Pool(threads) as poolgt: # Processing that create singel copy gene fasta with multicore
        args_gt = [(fasta_file, sub_model, bootstrap, input_dir, output_dir) for fasta_file in fasta_files]
        poolgt.map(process_genetree, args_gt)

def make_species_tree(threads, bootstrap, w, pip, input_dir, output_dir): # Make species tree
    if threads == 1:
        javaheapspace = "500M"
    else:
        int(threads * 0.8)  # Approximately 80% of the threads for java memoly
        javaheapspace = str(int(threads*0.8)) + "G"

    if w == True: # Make multi label tree and extend newick file
        if pip and bootstrap == 0:
            cmd_cat = f'cat {input_dir}/gene_tree/*.treefile > {input_dir}/ploidy.new'
            subprocess.run(cmd_cat, shell=True, check=True)
        elif pip and bootstrap != 0:
            cmd_cat = f'cat {input_dir}/gene_tree/*.contree > {input_dir}/ploidy.new'
            subprocess.run(cmd_cat, shell=True, check=True)
        else:
            cmd_cat = f'cat {input_dir}/gene_tree/*.contree > {input_dir}/ploidy.new'
            subprocess.run(cmd_cat, shell=True, check=True)
            if cmd_cat:
                pass
            else:
                cmd_cat = f'cat {input_dir}/gene_tree/*.treefile > {input_dir}/ploidy.new'
                subprocess.run(cmd_cat, shell=True, check=True)

        species_count = {}
        with open(f"{input_dir}/threshold.tsv") as f:
            for line in f:
                if line.strip():
                    species, count = line.strip().split()
                    species_count[species] = int(count)

        with open(f"{input_dir}/ploidy.new") as f: # Prepare input newick file for polyphest
            lines = f.readlines()

        processed_lines = []
        for line in lines:
            for species, count in species_count.items():
                if count >= 2: # Replac species_number to species
                    pattern = rf"{species}_[0-9]+"
                    line = re.sub(pattern, species, line)
            processed_lines.append(line)

        with open(f"{input_dir}/gene_tree.new", "w") as f:
            f.writelines(processed_lines)

        cmd_rm = f'rm {input_dir}/ploidy.new'
        subprocess.run(cmd_rm, shell=True, check=True)

        df = pd.read_csv(f"{input_dir}/threshold.tsv", sep='\t', header=None, names=["sample", "ploidy"])
        expanded_list = df.loc[df.index.repeat(df["ploidy"]), "sample"].tolist() # Repeat the sample name as many times as the ploidy
        with open(f"{input_dir}/gene_tree.txt", "w") as f:
            for sample in expanded_list:
                f.write(sample + "\n")

        cmd_polyphest = f'python Polyphest/Polyphest.py --gene_tree_file {input_dir}/gene_tree.new --consensus_multiset_file {input_dir}/gene_tree.txt --filter_strategy percentile --percentile 75 --use_near_isomorphic True --isomorphic_threshold 0.2'
        subprocess.run(cmd_polyphest, shell=True, check=True)

        cmd_exsptre = f'head -n 1 {input_dir}/gene_tree-polyphest.txt > {output_dir}species_tree.tree && sed -i "s@multree:\ @@g" {output_dir}/species_tree.tree'
        subprocess.run(cmd_exsptre, shell=True, check=True)

        cmd_exsptere = f'head -n 2 {input_dir}/gene_tree-polyphest.txt | tail -n 1 > {output_dir}/species_tree.enewick && sed -i "s@network:\ @@g" {output_dir}/species_tree.enewick'
        subprocess.run(cmd_exsptere, shell=True, check=True)

    else:
        if bootstrap == 0:
            cmd_cat = f'cat {input_dir}/gene_tree/*.treefile > {input_dir}/gene_tree.treefile'
            subprocess.run(cmd_cat, shell=True, check=True)

            cmd_iq = f'java -Xmx{javaheapspace} -jar ./Astral/astral.5.7.8.jar -i {input_dir}/gene_tree.treefiile -o {output_dir}/species_tree.tree 2> {input_dir}/species_nobs_tree.log'
            subprocess.run(cmd_iq, shell=True, check=True)

        else:
            cmd_cat = f'cat {input_dir}/gene_tree/*.contree > {input_dir}/gene_tree.contree'
            subprocess.run(cmd_cat, shell=True, check=True)

            cmd_find = f'find {input_dir}/gene_tree/*.ufboot > {input_dir}/gene_tree.ufboot'
            subprocess.run(cmd_find, shell=True, check=True)

            cmd_iq = f'java -Xmx{javaheapspace} -jar ./Astral/astral.5.7.8.jar -i {input_dir}/gene_tree.contree -b {input_dir}/gene_tree.ufboot -r {bootstrap} -o {input_dir}/species_bs{bootstrap}_tree.contree 2> {input_dir}/species_bs{bootstrap}_tree.log'
            subprocess.run(cmd_iq, shell=True, check=True)

            cmd_tree = f'tail -n 1 {input_dir}/species_bs{bootstrap}_tree.contree > {output_dir}/species_tree.tree'
            subprocess.run(cmd_tree, shell=True, check=True)


'''find potential marker gene'''
def rf_calculation(pip, bootstrap, species_tree, input_dir, output_dir):
    if pip:
        genetree_files = [f for f in os.listdir(f"{input_dir}/gene_tree") if f.endswith(".treefile" if bootstrap == 0 else ".contree")]
    else:
        genetree_files = [f for f in os.listdir(f"{input_dir}/gene_tree") if f.endswith(".contree")]
        if not genetree_files:
            genetree_files = [f for f in os.listdir(f"{input_dir}/gene_tree") if f.endswith(".treefile")]

    normalized_rf_dict = {}
    species_tree_path = os.path.join(f"{output_dir}/", species_tree)

    species_tree = Tree(species_tree_path, format=1)
    species_tree.set_outgroup(species_tree.get_midpoint_outgroup())

    for genetree_file in genetree_files:
        genetree_path = os.path.join(f"{input_dir}/gene_tree", genetree_file)

        gene_tree = Tree(genetree_path, format=1)
        #species_tree = Tree(species_tree_path, format=1)
        gene_tree.set_outgroup(gene_tree.get_midpoint_outgroup())
        #species_tree.set_outgroup(species_tree.get_midpoint_outgroup())

        try:
            rf_values = gene_tree.robinson_foulds(species_tree) # Get RF distance and maximum RF distance from dictionary
            rf_distance = rf_values[0]
            max_rf = rf_values[1]

            normalized_rf = rf_distance / max_rf
            gene_name = os.path.splitext(os.path.basename(genetree_file))[0]
            normalized_rf_dict[gene_name] = normalized_rf

        except Exception as e:
            print(f"error: {e}")
            normalized_rf_dict[gene_name] = None

    sorted_rf = sorted([(key, value) for key, value in normalized_rf_dict.items() if value is not None], key=lambda x: x[1])

    output_list = f"{output_dir}/marker.list"
    with open(output_list, "w") as f:
        f.write("gene\tRF\n")
        for key, value in sorted_rf:
            f.write(f"{key}\t{value:.2f}\n")

def parse_newick_labels(tree_str): # Parse the tree string and extract the labels in order
    tree = Phylo.read(StringIO(tree_str), "newick")
    labels = [term.name for term in tree.get_terminals()]
    return labels

def find_duplicates(labels):
    count = defaultdict(int)
    for label in labels:
        count[label] += 1
    return {label: n for label, n in count.items() if n > 1}

def generate_labelings(dup_species_counts):
    permutations_dict = {}
    for species, count in dup_species_counts.items():
        permutations_dict[species] = list(itertools.permutations(range(1, count+1)))
    return permutations_dict

def generate_all_combinations(permutations_dict):
    keys = list(permutations_dict.keys())
    all_combinations = list(itertools.product(*[permutations_dict[k] for k in keys]))
    return keys, all_combinations

def relabel_tree(tree_str, keys, combo):
    tree = Phylo.read(StringIO(tree_str), "newick")
    label_indices = defaultdict(int)
    new_tree = copy.deepcopy(tree)

    for clade in new_tree.get_terminals(): # Cycle through each terminal node and label them in order
        if clade.name in keys:
            species = clade.name
            idx = combo[keys.index(species)][label_indices[species]]
            clade.name = f"{species}_{idx}"
            label_indices[species] += 1
    return new_tree

def make_candidates(input_dir, output_dir):
    if not os.path.exists(f"{input_dir}./candidates_sptree"):
        os.makedirs(f"{input_dir}/candidates_sptree")

    with open(f"{output_dir}/species_tree.tree") as f:
        tree_str = f.read().strip()

    labels = parse_newick_labels(tree_str)
    duplicates = find_duplicates(labels)

    if not duplicates:
        print("No overlap species were found.")
        return

    perms = generate_labelings(duplicates)
    keys, all_combos = generate_all_combinations(perms)

    for i, combo in enumerate(all_combos, 1):
        new_tree = relabel_tree(tree_str, keys, combo)
        with open(f"{input_dir}/candidate{i}.tree", "w") as f:
            Phylo.write(new_tree, f, "newick")

    cmd_mv = f"mv {input_dir}/candidate*.tree {input_dir}/candidates_sptree"
    subprocess.run(cmd_mv, shell=True, check=True)

def rf_calculation_multiple_species(pip, bootstrap, candidates_dir, input_dir, output_dir):
    # Selection of gene phylogenetic trees
    if pip:
        genetree_files = [f for f in os.listdir(f"{input_dir}/gene_tree") if f.endswith(".treefile" if bootstrap == 0 else ".contree")]
    else:
        genetree_files = [f for f in os.listdir(f"{input_dir}/gene_tree") if f.endswith(".contree")]
        if not genetree_files:
            genetree_files = [f for f in os.listdir(f"{input_dir}/gene_tree") if f.endswith(".treefile")]

    # candidates_dirからspecies_treeファイルを取得
    species_tree_files = [f for f in os.listdir(candidates_dir) if f.endswith(".tree")]

    # Dictionary for storing results： { gene_name: { species_tree_name: nRF, ... }, ... }
    result_dict = {}

    for species_file in species_tree_files:
        species_path = os.path.join(candidates_dir, species_file)
        species_tree = Tree(species_path, format=1)
        species_tree.set_outgroup(species_tree.get_midpoint_outgroup())

        for genetree_file in genetree_files:
            gene_name = os.path.splitext(genetree_file)[0]
            genetree_path = os.path.join(f"{input_dir}/gene_tree", genetree_file)

            gene_tree = Tree(genetree_path, format=1)
            gene_tree.set_outgroup(gene_tree.get_midpoint_outgroup())

            try:
                rf_distance, max_rf = gene_tree.robinson_foulds(species_tree)[:2]
                normalized_rf = rf_distance / max_rf

            except Exception as e:
                print(f"Error for gene {gene_name} with species tree {species_file}: {e}")
                normalized_rf = None

            if gene_name not in result_dict:
                result_dict[gene_name] = {}
            result_dict[gene_name][species_file] = normalized_rf

    df = pd.DataFrame.from_dict(result_dict, orient='index')

    # Extract the column with the minimum values ​​for marker.list
    min_rf = df.min(axis=1, skipna=True)

    with open(f"{output_dri}/marker.list", "w") as f:
        f.write("gene\tRF\n")
        for gene, rf_val in min_rf.dropna().sort_values().items():
            f.write(f"{gene}\t{rf_val:.2f}\n")

    return df


'''make multi fasta for primer design'''
def split_fasta_by_segments(input_fasta, output_file):
    with open(input_fasta, "r") as fasta_file:
        record = next(SeqIO.parse(fasta_file, "fasta"))

    fasta_header = record.description
    fasta_sequence = str(record.seq)

    match = re.search(r'segs:([\d,-]+)', fasta_header)
    if not match:
        raise ValueError("Invalid header format. 'segs' information missing.")

    segments = match.group(1).split(',')
    segment_ranges = [tuple(map(int, seg.split('-'))) for seg in segments]

    records = []
    for i, (start, end) in enumerate(segment_ranges, 1):
        subseq = fasta_sequence[start-1:end]  # Convert 1-based to 0-based index
        record_part = SeqRecord(Seq(subseq), id=f"{record.id}_{start}-{end}", description="")
        records.append(record_part)

    SeqIO.write(records, output_file, "fasta")
    #print(f"Multi-FASTA file '{output_file}' created successfully.")

def extract_region_from_fasta(reference_fasta, query_header, output_fasta):
    match = re.search(r'loc:([\w\.]+)\(([-+])\)(\d+)-(\d+)', query_header)
    if not match:
        raise ValueError("Invalid header format. 'loc' information missing.")

    contig, strand, start, end = match.groups()
    start, end = int(start), int(end)
    start = start - 30
    end = end + 30

    with open(reference_fasta, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id == contig:
                subseq = record.seq[start-1:end] if strand == '+' else record.seq[start-1:end].reverse_complement()
                new_record = SeqRecord(subseq, id=f"{contig}_{start}-{end}", description=f"Extracted from {contig} ({strand}) {start}-{end}")
                SeqIO.write(new_record, output_fasta, "fasta")
                #print(f"Extracted sequence saved to '{output_fasta}'")
                return

    raise ValueError(f"Contig '{contig}' not found in reference FASTA.")

def extract_value_from_tsv_w(input_gene, sp, ploidy_value, output_dir):
    with open(f'{output_dir}/single_copy_genes.tsv', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row[reader.fieldnames[0]] == input_gene:
                sp_name_tsv = [f"{sp}_{i}" for i in range(1, ploidy_value + 1)]
                return [row.get(name, None) for name in sp_name_tsv]
    return None

def extract_value_from_tsv(input_gene, sp, output_dir):
    with open(f'{output_dir}/single_copy_genes.tsv', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row[reader.fieldnames[0]] == input_gene:
                return row.get(sp, None)
    return None

def primer(w, target_gene, input_species, input_dir, output_dir):
    if not os.path.exists(f"{output_dir}/primer_design"):
        os.makedirs(f"{output_dir}/primer_design")

    if w == True:
        threshold_file = None
        for fname in os.listdir(f'{input_dir}/'):
            if fname.lower() in ['threshold.tsv', 'threshold.csv']:
                threshold_file = fname
                break
        if threshold_file is None:
            raise FileNotFoundError(f"There is no {input_dir}/threshold file.")
        ext = os.path.splitext(threshold_file)[1].lower()
        delimiter = '\t' if ext == '.tsv' else ','

        # read threshold and cast to int
        threshold_df = pd.read_csv(threshold_file, sep=delimiter, header=None, index_col=0)
        threshold_df = threshold_df[1].astype(int)

    extensions = ['.fasta', '.fas', '.fa', '.fna', '.seq', '.aln', '.maf', '.nt']
    for num, sp in enumerate(input_species):
        for ext in extensions:
            wgs_files = glob.glob(f"{input_dir}/species/{sp}{ext}")
            if wgs_files:
                wgs_file = wgs_files[0]
                break
        else:
            raise FileNotFoundError(f"No FASTA file found for {sp} in {input_dir}/species/")

        if w == True:
            ploidy_value = threshold_df.loc[sp].item()
            geneMP = extract_value_from_tsv_w(target_gene, sp, ploidy_value, output_dir)
        else:
            geneMP = extract_value_from_tsv(target_gene, sp, output_dir)

        cmd_sp = f"gffread -W -x {sp}.tmp.fna -g {wgs_file} {input_dir}/gff/{sp}.gff"
        subprocess.run(cmd_sp, shell=True, check=True)

        if w == True:
            for i, gene in enumerate(geneMP, 1):
                if gene:  # Skip 'None' or 'empty'
                    outname = f"{sp}.tmp.{target_gene}_{i}.fna" if len(geneMP) > 1 else f"{sp}.tmp.{target_gene}.fna"
                    header_file = f"{sp}_header_{i}"
                    region_file = f"{sp}_region_{i}.fasta"

                    cmd_sp_g = f"seqkit grep -p {gene} {sp}.tmp.fna > {outname}"
                    subprocess.run(cmd_sp_g, shell=True, check=True)

                    if num == 0 and i == 1:
                        split_fasta_by_segments(outname, f'{sp}_{target_gene}_split.fasta')

                    ex_header = f"seqkit fx2tab -n {outname} | cat > {header_file}"
                    subprocess.run(ex_header, shell=True, check=True)

                    with open(header_file, "r") as file:
                        query_header = file.read().strip()

                    extract_region_from_fasta(wgs_file, query_header, region_file)

                    rmtmp = f"rm {outname} {header_file}"
                    subprocess.run(rmtmp, shell=True, check=True)

            rmtmp = f"rm {sp}.tmp.fna"
            subprocess.run(rmtmp, shell=True, check=True)

        else:
            if geneMP: # Skip 'None' or 'empty'
                outname = f"{sp}.tmp.{target_gene}.fna"
                cmd_sp_g = f"seqkit grep -p {geneMP} {sp}.tmp.fna > {outname}"
                subprocess.run(cmd_sp_g, shell=True, check=True)

                if num == 0: # first sample make exon fasta for primaer design
                    split_fasta_by_segments(outname, f'{sp}_{target_gene}_split.fasta')

            ex_header = f"seqkit fx2tab -n {outname} | cat > {sp}_header"
            subprocess.run(ex_header, shell=True, check=True)

            with open(sp + '_header', "r") as file:
                query_header = file.read().strip()

            extract_region_from_fasta(wgs_file, query_header, sp + '_region.fasta')

            rmtmp = f"rm {sp}.tmp.fna {outname} {sp}_header"
            subprocess.run(rmtmp, shell=True, check=True)

    catfasta = f"cat *fasta > {target_gene}.fasta"
    subprocess.run(catfasta, shell=True, check=True)

    mafftfasta = f"mafft --quiet {target_gene}.fasta > {target_gene}_primerdesign.fasta"
    subprocess.run(mafftfasta, shell=True, check=True)

    mvfasta = f"mv {target_gene}_primerdesign.fasta {output_dir}/primer_design/ && rm *fasta"
    subprocess.run(mvfasta, shell=True, check=True)


'''edid'''
def update_fasta(args):
    fasta_name, fasta_path, tab_data = args

    if fasta_name in tab_data: # Get the row corresponding to gene
        gene_data = tab_data[fasta_name]
    else:
        print(f"Warning: {fasta_name} not found in tab file.")
        return

    updated_records = []

    for record in SeqIO.parse(fasta_path, "fasta"): # Read FASTA file and update IDs
        sample = record.id
        if sample in gene_data:
            record.id = gene_data[sample]
            record.description = ""  # Reset description if necessary
        else:
            print(f"Warning: {sample} not found in tab file for {fasta_name}.")
        updated_records.append(record)

    with open(fasta_path, "w") as output_handle: # Save the updated FASTA file
        SeqIO.write(updated_records, output_handle, "fasta")

def parse_tab_file(tab_file_path): # Parses TSV and returns the data in dictionary format
    tab_data = {}
    with open(tab_file_path, "r") as file:
        lines = file.readlines()
        header = lines[0].strip().split("\t")[1:]  # Get sample name from first line
        for line in lines[1:]:
            parts = line.strip().split("\t")
            gene = parts[0]
            values = parts[1:]
            tab_data[gene] = {header[i]: values[i] for i in range(len(header))}

    return tab_data


def parse_species_name(file_path):
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
    if not first_line.startswith(">"):
        return None
    parts = first_line[1:].split()
    if len(parts) < 3:
        return None
    genus = parts[1]
    species = parts[2]
    return f"{genus[0].upper()}{species}"

def collect_species_counts(file_paths):
    species_counts = defaultdict(int)
    file_to_species = {}

    for file_path in file_paths:
        species = parse_species_name(file_path)
        if species:
            species_counts[species] += 1
            file_to_species[file_path] = species
        else:
            print(f"Skip (invalid header): {file_path}")

    return species_counts, file_to_species

def rename_all_wgs(file_paths):
    species_counts, file_to_species = collect_species_counts(file_paths)
    species_index = defaultdict(int)

    for file_path in sorted(file_paths):  # Sorting will result in a predictable order
        if file_path not in file_to_species:
            continue
        species = file_to_species[file_path]
        count = species_counts[species]
        directory = os.path.dirname(file_path)

        if count == 1:
            new_name = f"{species}.fna"
        else:
            species_index[species] += 1
            new_name = f"{species}_{species_index[species]}.fna"

        new_path = os.path.join(directory, new_name)
        os.rename(file_path, new_path)
        print(f"Renamed: {os.path.basename(file_path)} → {new_name}")


'''main commands'''
def prepare(threads, ref_cds, p, w, outsthreshold, input_dir):
    wgs2gff(threads, ref_cds, w, outsthreshold, input_dir)
    paf2tsv(input_dir)
    gff2fasta(threads, p, input_dir)

def curate(threads, w, input_dir, output_dir, pip):
    pre_scg1 = tsv2scg(w, input_dir)
    if w:
        results_dict = gff2scg_w(pre_scg1, threads, input_dir)
        expected_counts = load_expected_counts(input_dir)
        pre_scg2 = extract_strict_scg(results_dict, expected_counts)
        scg_orthologs1 = gff2scg_df_w(pre_scg2, threads, input_dir)
        fasta2scg_w(pre_scg2, scg_orthologs1, threads, input_dir, output_dir)
    else:
        pre_scg2 = gff2scg(pre_scg1, threads, input_dir)
        scg_orthologs1 = gff2scg_df(pre_scg2, threads, input_dir)
        fasta2scg(pre_scg2, scg_orthologs1, threads, input_dir, output_dir)
    scg_orthologs2 = last_trimming(pre_scg2, scg_orthologs1, output_dir)
    print(scg_orthologs2)
    scg_orthologs2.to_csv(f"{output_dir}/single_copy_genes.tsv", sep='\t')
    alignment_fasta(threads, pip, output_dir)

def phylo(threads, sub_model, bootstrap, w, input_dir, output_dir, species_tree, pip):
    make_gene_tree(threads, sub_model, bootstrap, input_dir, output_dir)
    if species_tree == "species_tree.tree":
        make_species_tree(threads, bootstrap, w, pip, input_dir, output_dir)
    else:
        pass

def marker(w, input_dir, output_dir, pip, bootstrap, species_tree):
    print(f"Calculate nRF using {species_tree}")
    if w == True:
        make_candidates(input_dir, output_dir)
        test_df = rf_calculation_multiple_species(pip, bootstrap, candidates_dir=f"{input_dir}/candidates_sptree", input_dir, output_dir)
        print(test_df)
        test_df.to_csv("test_df.csv")
    else:
        rf_calculation(pip, bootstrap, species_tree, input_dir, output_dir)

def pipeline(threads, ref_cds, p, w, outsthreshold, sub_model, bootstrap, input_dir, output_dir, species_tree):
    prepare(threads, ref_cds, p, w, outsthreshold, input_dir)
    pip = True
    curate(threads, w, input_dir, output_dir, pip)
    phylo(threads, sub_model, bootstrap, w, input_dir, output_dir, species_tree, pip)
    marker(w, input_dir, output_dir, pip, bootstrap, species_tree)

def edid(threads, species, insert, eject, restore, input_dir, output_dir):
    if species:
        file_paths = [os.path.join(f"{input_dir}/species", f) for f in os.listdir(f"{input_dir}/species")
              if f.endswith(('.fasta', '.fas', '.fa', '.fna', '.faa', '.pep', '.seq', '.aln', '.maf', '.nt', '.aa'))]
        rename_all_wgs(file_paths)

    if insert:
        fasta_files = [f for f in os.listdir(f"{input_dir}/input") if f.endswith(('.fasta', '.fas', '.fa', '.fna', '.faa', '.pep', '.seq', '.aln', '.maf', '.nt', '.aa'))]
        with multiprocessing.Pool(threads) as pool_in:
            arg_in = [(os.path.splitext(fasta_file)[0], os.path.join(f"{input_dir}/input", fasta_file)) for fasta_file in fasta_files]
            pool_in.map(insert_id, arg_in)

    if eject:
        fasta_files = [f for f in os.listdir(f"{output_dir}/single_copy_genes") if f.endswith(".fasta")]
        with multiprocessing.Pool(threads) as pool_ej:
            arg_ej = [os.path.join(f"{output_dir}/single_copy_genes", fasta_file) for fasta_file in fasta_files]
            pool_ej.map(eject_id, arg_ej)

    if restore:
        fasta_files = [f for f in os.listdir(f"{output_dir}/single_copy_genes") if f.endswith(".fasta")]
        tab_data = parse_tab_file(f"{output_dir}/single_copy_genes.tsv")
        with multiprocessing.Pool(threads) as pool_re:
            args_re = [(os.path.splitext(fasta_file)[0], os.path.join(f"{output_dir}/single_copy_genes", fasta_file), tab_data) for fasta_file in fasta_files]
            pool_re.map(update_fasta, args_re)
        cmd_rm = f"rm {output_dir}/single_copy_genes/sed*"
        subprocess.run(cmd_rm, shell=True, check=True)



class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def add_arguments(self, actions):
        # Override the `add_arguments` method to display nothing
        pass

def main():
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(dest="command")

    parser_1 = subparsers.add_parser("prepare", help="Making input from whole genome seq and proteins fasta.")
    parser_1.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")
    parser_1.add_argument("-r", "--ref_cds", required=True, metavar="STR", help="reference CDS fasta.")
    parser_1.add_argument("-p", action="store_true", help="Get amino acid seq of single copy genes (default: nucleotide seq).")
    parser_1.add_argument("-w", action="store_true", help="Providing information of ploidy as input.")
    parser_1.add_argument("-f", "--outsthreshold", type=float, default=0.50, help="Threshold for miniprot --outs when ploidy is detcted (default: 0.50).")
    parser_1.add_argument("-i", "--input_dir", type=str, default="./input", help="Change the input directory destination (default: ./input).")

    parser_2 = subparsers.add_parser("curate", help="Generating single copy gene fasta and table from input.")
    parser_2.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")
    parser_2.add_argument("-w", action="store_true", help="Providing information of ploidy as input.")
    parser_2.add_argument("-i", "--input_dir", type=str, default="./input", help="Change the input directory destination (default: ./input).")
    parser_2.add_argument("-o", "--output_dir", type=str, default="./output", help="Change the output directory destination (default: ./output).")

    parser_3 = subparsers.add_parser("phylo", help="Generating phylogenetic tree from single copy gene fasta.")
    parser_3.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")
    parser_3.add_argument("-m", "--sub_model", type=str, default="MFP", help="Change substitution model of iq-tree (default: MFP).")
    parser_3.add_argument("-bs", "--bootstrap", type=int, default=1000, help="Number of bootstrap (default: 1000).")
    parser_3.add_argument("-w", action="store_true", help="Providing information of ploidy as input.")
    parser_3.add_argument("-i", "--input_dir", type=str, default="./input", help="Change the input directory destination (default: ./input).")
    parser_3.add_argument("-o", "--output_dir", type=str, default="./output", help="Change the output directory destination (default: ./output).")

    parser_4 = subparsers.add_parser("marker", help="Generating list of potential marker resion from tree file.")
    parser_4.add_argument("-sp", "--species_tree", type=str, default="species_tree.tree", help="reference species tree (default: species_tree.tree).")
    parser_4.add_argument("-w", action="store_true", help="Providing information of ploidy as input.")
    parser_4.add_argument("-i", "--input_dir", type=str, default="./input", help="Change the input directory destination (default: ./input).")
    parser_4.add_argument("-o", "--output_dir", type=str, default="./output", help="Change the output directory destination (default: ./output).")

    parser_5 = subparsers.add_parser("primer", help="Generate multi fasta for primer design.")
    parser_5.add_argument("-w", action="store_true", help="Providing information of ploidy as input.")
    parser_5.add_argument("-g", "--target_gene", required=True, metavar="STR", help="Name of target gene.")
    parser_5.add_argument("-s", "--input_species", type=str, nargs='+', help="Name of species use for desgin primaer.")
    parser_5.add_argument("-i", "--input_dir", type=str, default="./input", help="Change the input directory destination (default: ./input).")
    parser_5.add_argument("-o", "--output_dir", type=str, default="./output", help="Change the output directory destination (default: ./output).")

    parser_6 = subparsers.add_parser("pipeline", help="Run pipleline prepare, curate, phylo and marker.")
    parser_6.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")
    parser_6.add_argument("-r", "--ref_cds", required=True, metavar="STR", help="reference CDS fasta.")
    parser_6.add_argument("-p", action="store_true", help="Get amino acid seq of single copy genes (default: nucleotide seq).")
    parser_6.add_argument("-w", action="store_true", help="Providing information of ploidy as input.")
    parser_6.add_argument("-f", "--outsthreshold", type=float, default=0.50, help="Threshold for miniprot --outs when ploidy is detcted (default: 0.50).")
    parser_6.add_argument("-m", "--sub_model", type=str, default="MFP", help="Change substitution model of iq-tree (default: MFP).")
    parser_6.add_argument("-bs", "--bootstrap", type=int, default=1000, help="Number of bootstrap (default: 1000).")
    parser_6.add_argument("-sp", "--species_tree", type=str, default="species_tree.tree", help="reference species tree (default: species_tree.tree).")
    parser_6.add_argument("-i", "--input_dir", type=str, default="./input", help="Change the input directory destination (default: ./input).")
    parser_6.add_argument("-o", "--output_dir", type=str, default="./output", help="Change the output directory destination (default: ./output).")

    parser_7 = subparsers.add_parser("edid", help="Edit fasta id.")
    parser_7.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")
    parser_7.add_argument("-species", action="store_true", help="Change the name of whole genome fasta to species_name.fna.")
    parser_7.add_argument("-insert", action="store_true", help="Insert sample_name to fasta_id.")
    parser_7.add_argument("-eject", action="store_true", help="Eject gene_id form fasta_id.")
    parser_7.add_argument("-restore", action="store_true", help="Restore sample_name_gene_id using tsv file.")
    parser_7.add_argument("-i", "--input_dir", type=str, default="./input", help="Change the input directory destination (default: ./input).")
    parser_7.add_argument("-o", "--output_dir", type=str, default="./output", help="Change the output directory destination (default: ./output).")

    if len(sys.argv) < 2 or sys.argv[1] in ["-h", "--help"]:
        print("cusco  Version=0.1.0\n\n\
Usage:\n\
    {} [command] option\n\n\
Command list:\n\
    prepare       Generate input files.\n\
    curate        Generate multi fasta and table of single copy genes.\n\
    phylo         Generate phylogenetic tree using single copy genes.\n\
    marker        Generate list of potential marker resion.\n\
    primer        Generate multi fasta for primer design.\n\
    pipeline      All processing is done automatically.\n\
    edid          Edit fasta id.\n".format(sys.argv[0]))
        sys.exit(0)

    args = parser.parse_args()

    if args.command == "prepare":
        prepare(args.threads, args.ref_cds, args.p, args.w, args.outsthreshold, args.input_dir)
    elif args.command == "curate":
        curate(args.threads, args.w, args.input_dir, args.output_dir, pip=False)
    elif args.command == "phylo":
        phylo(args.threads, args.sub_model, args.bootstrap, args.w, args.input_dir, args.output_dir, species_tree="species_tree.tree", pip=False)
    elif args.command == "marker":
        marker(args.w, args.input_dir, args.output_dir, pip=False, bootstrap=1000, species_tree=args.species_tree)
    elif args.command =="primer":
        primer(args.w, args.target_gene, args.input_species, args.input_dir, args.output_dir)
    elif args.command == "pipeline":
        pipeline(args.threads, args.ref_cds, args.p, args.w, args.outsthreshold, args.sub_model, args.bootstrap, args.input_dir, args.output_dir, species_tree=args.species_tree)
    elif args.command == "edid":
        edid(args.threads, args.species, args.insert, args.eject, args.restore, args.input_dir, args.output_dir)
    else:
        print("Error: {} [command] \n\n    [prepare][curate][phylo][marker][primer][pipeline][edid]\n".format(sys.argv[0]))
        sys.exit(1)

if __name__ == "__main__":
    main()
