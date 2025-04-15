#!/usr/bin/env python3
# ver.0.1.0

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from ete3 import Tree
import glob
import multiprocessing
import os
import pandas as pd
from pathlib import Path
import re
import subprocess
import sys


'''prepare input file for cureatscg'''
def wgs2gff(threads, ref_cds):
    if not os.path.exists("./gff"):
        os.makedirs("./gff")
    if not os.path.exists("./paf"):
        os.makedirs("./paf")

    wgs_files = [f for f in sorted(os.listdir('./species')) if f.endswith(('.fasta', '.fas', '.fa', '.fna', '.seq', '.aln', '.maf', '.nt'))] # List of WGS files

    for wgs_file in wgs_files:
        wgs_path = os.path.join("./species", wgs_file)
        sample_name = os.path.splitext(os.path.basename(wgs_path))[0]
        index_file = f"{wgs_path}.fai"

        if not os.path.exists(index_file): # If the index file does not exist, run samtools faidx
            subprocess.run(['samtools', 'faidx', wgs_path], check=True)

        cmd_wgs2gff = f"miniprot --gff --outs=0.97 --no-cs -Iut {threads} {wgs_path} {ref_cds} > {sample_name}.out"
        subprocess.run(cmd_wgs2gff, shell=True, check=True)  # Run miniprot to GFF and PAF

        paf_file, gff_file = f"{sample_name}.paf", f"{sample_name}.gff"

        with open(f"{sample_name}.out", "r", encoding="utf-8") as infile, \
             open(paf_file, "w", encoding="utf-8") as file_paf, \
             open(gff_file, "w", encoding="utf-8") as file_gff:

            for line in infile:
                if line.startswith("##PAF\t"):
                    file_paf.write(line[6:])
                else:
                    file_gff.write(line)

        cmd_mv1 = f"mv {sample_name}.gff ./gff"
        subprocess.run(cmd_mv1, shell=True, check=True)

        cmd_mv2 = f"mv {sample_name}.paf ./paf"
        subprocess.run(cmd_mv2, shell=True, check=True)

        cmd_rm = f"rm {sample_name}.out"
        subprocess.run(cmd_rm, shell=True, check=True)

def insert_id(args):
    sample_name, file_path = args
    cmd_iid = f"sed -i -E 's/^>(MP[0-9]+)/>{sample_name}_\\1/g' {file_path}"
    subprocess.run(cmd_iid, shell=True, check=True)  # Run sed

def process_gff2fna(arg_na):
    gff_path = os.path.join("./gff", arg_na)
    sample_name = arg_na.split('.gff')[0]

    extensions = ['.fasta', '.fas', '.fa', '.fna', '.seq', '.aln', '.maf', '.nt']
    for ext in extensions:
        fasta_files = glob.glob(f"./species/{sample_name}{ext}")
        if fasta_files:
            fasta_file = fasta_files[0]
            break
    else:
        raise FileNotFoundError(f"No FASTA file found for {sample_name} in ./species/")

    cmd_fna = f"gffread -x {sample_name}.fna -g {fasta_file} {gff_path}"
    subprocess.run(cmd_fna, shell=True, check=True)  # Run gffread

    file_path = sample_name + ".fna"
    insert_id((sample_name, file_path))

    cmd_mv = f"mv {sample_name}.fna ./input"
    subprocess.run(cmd_mv, shell=True, check=True)

def process_gff2faa(arg_aa):
    gff_path = os.path.join("./gff", arg_aa)
    sample_name = arg_aa.split('.gff')[0]

    extensions = ['.fasta', '.fas', '.fa', '.fna', '.seq', '.aln', '.maf', '.nt']
    for ext in extensions:
        fasta_files = glob.glob(f"./species/{sample_name}{ext}")
        if fasta_files:
            fasta_file = fasta_files[0]
            break
    else:
        raise FileNotFoundError(f"No FASTA file found for {sample_name} in ./species/")

    cmd_faa = f"gffread -y {sample_name}.faa -g {fasta_file} {gff_path}"
    subprocess.run(cmd_faa, shell=True, check=True)  # Run gffread

    file_path = sample_name + ".faa"
    insert_id((sample_name, file_path))

    cmd_mv = f"mv {sample_name}.faa ./input"
    subprocess.run(cmd_mv, shell=True, check=True)

def gff2fasta(threads, p):
    if not os.path.exists("./input"):
        os.makedirs("./input")
    gff_files = [f for f in sorted(os.listdir('./gff')) if f.endswith('.gff')] # List of GFF files

    if p == True:
        with multiprocessing.Pool(threads) as pool_aa: # processing that convert gff to fna
            arg_aa = [gff_file for gff_file in gff_files]
            pool_aa.map(process_gff2faa, arg_aa)
    else:
        with multiprocessing.Pool(threads) as pool_na: # processing that convert gff to fna
            arg_na = [gff_file for gff_file in gff_files]
            pool_na.map(process_gff2fna, arg_na)

def paf2tsv(): # make CNVtable from paf
    cmd_gfa = f"./pangene-1.1-bin/bin_x64-linux/pangene ./paf/*paf > graph.gfa"
    subprocess.run(cmd_gfa, shell=True, check=True)  # Run pangene

    cmd_tsv = f"./pangene-1.1-bin/bin_x64-linux/k8 ./pangene-1.1-bin/scripts/pangene.js gfa2matrix -c graph.gfa > pre_geneCNV.Rtab"
    subprocess.run(cmd_tsv, shell=True, check=True)  # Run pangene

    with open("pre_geneCNV.Rtab", "r", encoding="utf-8") as infile, \
         open("geneCNV.Rtab", "w", encoding="utf-8", newline="") as outfile:

        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        for row_index, row in enumerate(reader):
            if row_index == 0:
                row[1:] = [col[:-2] if col.endswith("#0") else col for col in row[1:]]
            writer.writerow(row)

    cmd_rm = "rm pre_geneCNV.Rtab"
    subprocess.run(cmd_rm, shell=True, check=True)


'''cureate single copy genes'''
'''1st step of cureate single copy genes'''
def tsv2scg(): # first step : extract single copy genes from Rtab table
    ori_Rtab = pd.read_csv("geneCNV.Rtab", sep='\t')
    filtered_Rtab = ori_Rtab[ori_Rtab.iloc[:, 1:].eq(1).all(axis=1)]
    pre_scg_step1 = filtered_Rtab.iloc[:, 0].to_list()

    return pre_scg_step1

'''2nd step of cureate single copy genes'''
def list_to_df_and_find_duplicates(input_list):
    find_dup_df = pd.DataFrame(input_list, columns=["Column1"])
    duplicates_list = find_dup_df[find_dup_df.duplicated(subset="Column1", keep=False)]["Column1"].unique().tolist()

    return duplicates_list

def list_to_df_and_find_unique(input_list):
    find_unique_df = pd.DataFrame(input_list, columns=["Column1"])
    unique_list = find_unique_df[~find_unique_df.duplicated(subset="Column1", keep=False)]["Column1"].tolist()

    return unique_list

def process_gff_file(args1):
    gff_file, pre_scg = args1
    results = []
    file_name = gff_file.split('.gff')[0]  # extract sample name
    gff_path = './gff/' + gff_file

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

def gff2scg(pre_scg, threads): #second step : trimming single copy genes of first step using GFF files
    gff_files = [f for f in sorted(os.listdir('./gff')) if f.endswith('.gff')] # List of GFF files

    with multiprocessing.Pool(threads) as pool: # processing that find duplicate genes with multicore
        args1 = [(gff_file, pre_scg) for gff_file in gff_files]
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
    gff_file, pre_scg2, scg_orthologs1 = args2
    results2 = []
    file_name = gff_file.split('.gff')[0]  # extract sample name
    gff_path = './gff/' + gff_file

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

def gff2scg_df(pre_scg2, threads): # make orthologs table
    gff_files = [f for f in sorted(os.listdir('./gff')) if f.endswith('.gff')] # List of GFF files
    scg_orthologs1 = pd.DataFrame(pre_scg2, columns=["reference_gene"]) # main df

    with multiprocessing.Pool(threads) as pool2: # Processing that create merge df with multicore
        args2 = [(gff_file, pre_scg2, scg_orthologs1) for gff_file in gff_files]
        result_dfs = pool2.map(process_gff_file_df, args2)

        final_df = scg_orthologs1.copy() # merge all df
        for df in result_dfs:
            final_df = pd.concat([final_df, df], axis=1)

        #final_df = final_df.sort_index(axis=1) # sort df

    return final_df


def process_pep_file(args3):
    ref_gene, scg_orthologs1 = args3
    pep_files = [f for f in sorted(os.listdir('./input')) if f.endswith(('.fasta', '.fas', '.fa', '.fna', '.faa', '.pep', '.seq', '.aln', '.maf', '.nt', '.aa'))]
    temp_files = []

    if ref_gene in scg_orthologs1["reference_gene"].values:
        target_row = scg_orthologs1[scg_orthologs1["reference_gene"] == ref_gene]

        for col in scg_orthologs1.columns[1:]:
            ortho_gene = target_row[col].values[0] # extract orthologs name
            file_name = col # extract sample name

            pep_files = glob.glob(os.path.join("./input", f"{file_name}.*"))

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
                print(f"No files found = {pep_files} in ./input")

        combined_pep = os.path.join("./single_copy_genes", f"{ref_gene}.fa")
        with open(combined_pep, 'w') as out_f: # merge temporary files
            for temp_file in temp_files:
                with open(temp_file, 'r') as tf:
                    out_f.write(tf.read())

        for temp_file in temp_files: # remove temporary files
            os.remove(temp_file)

def fasta2scg(pre_scg2, scg_orthologs1, threads): # Create single copy gene FASTA
    if not os.path.exists("./single_copy_genes"): # create output directory
        os.makedirs("./single_copy_genes")
    ref_genes = scg_orthologs1["reference_gene"].tolist() # reference genes

    with multiprocessing.Pool(threads) as pool3: # Processing that create singel copy gene fasta with multicore
        args3 = [(ref_gene, scg_orthologs1) for ref_gene in ref_genes]
        pool3.map(process_pep_file, args3)

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

def last_trimming(pre_scg2, scg_orthologs1): # third step : trimming single copy genes of second step using FASTA
    scg_orthologs1 = scg_orthologs1.set_index("reference_gene")
    sample_num = len(scg_orthologs1.columns)
    fasta_files = [f for f in os.listdir("./single_copy_genes") if f.endswith(".fa")]

    for fasta_file in fasta_files:
        fasta_path = os.path.join("./single_copy_genes", fasta_file)
        seq_count = count_sequences_in_fasta(fasta_path) # Count the number of sequences in a FASTA file using seqkit
        if seq_count is not None and seq_count != sample_num:
            print(f"{fasta_file} has {seq_count} sequences (not {sample_num})")
            gene_name = fasta_file.split('.fa')[0]
            scg_orthologs1 = scg_orthologs1.drop(index=gene_name)
            os.remove(fasta_path)

    return scg_orthologs1

'''formatting fasta'''
def eject_id(alignment_path):
    cmd_eid = f"sed -i -E 's/^>([^_]+)_.*/>\\1/g' {alignment_path}"
    subprocess.run(cmd_eid, shell=True, check=True)  # Run sed

def process_alignment(args_mafft): # Run mafft
    fasta_file, pip =  args_mafft
    fasta_path = os.path.join("./single_copy_genes", fasta_file)
    alignment_path = fasta_path.replace(".fa", ".fasta")

    cmd_mafft = f"mafft --quiet {fasta_path} > {alignment_path}"
    subprocess.run(cmd_mafft, shell=True, check=True)  # Run mafft

    if pip:
        eject_id(alignment_path)

    cmd_rm = f"rm {fasta_path}"
    subprocess.run(cmd_rm, shell=True, check=True)

def alignment_fasta(threads, pip):
    fasta_files = [f for f in os.listdir("./single_copy_genes") if f.endswith(".fa")]

    with multiprocessing.Pool(threads) as poolmafft: # Processing that create singel copy gene fasta with multicore
        args_mafft = [(fasta_file, pip) for fasta_file in fasta_files]
        poolmafft.map(process_alignment, args_mafft)


'''make phylogenetic tree'''

def process_genetree(args_gt): # Run iqtree
    fasta_file, bootstrap =  args_gt

    fasta_path = os.path.join("./single_copy_genes", fasta_file)
    gene_name = fasta_file.split('.fasta')[0]

    if bootstrap == 0:
        cmd1 = f"iqtree2 -s {fasta_path} -T AUTO --quiet -m MFP --prefix {gene_name}"
    else:
        cmd1 = f"iqtree2 -s {fasta_path} -T AUTO --quiet -m MFP --ufboot {bootstrap} --bnni --prefix {gene_name}"
    subprocess.run(cmd1, shell=True, check=True)  # Run iqtree

    cmd2 = f"mv {gene_name}.* ./gene_tree"
    subprocess.run(cmd2, shell=True, check=True)

def make_gene_tree(threads, bootstrap): # make gene tree for species tree
    if not os.path.exists("./gene_tree"):
        os.makedirs("./gene_tree")

    fasta_files = [f for f in os.listdir("./single_copy_genes") if f.endswith(".fasta")]
    fasta_files.sort()

    with multiprocessing.Pool(threads) as poolgt: # Processing that create singel copy gene fasta with multicore
        args_gt = [(fasta_file, bootstrap) for fasta_file in fasta_files]
        poolgt.map(process_genetree, args_gt)

def make_species_tree(threads, bootstrap): # make species tree
    if threads == 1:
        javaheapspace = "500M"
    else:
        int(threads * 0.8)  # Approximately 80% of the threads for java memoly
        javaheapspace = str(int(threads*0.8)) + "G"

    if bootstrap == 0:
        cmd_cat = 'cat ./gene_tree/*.treefile > gene_tree.treefile'
        subprocess.run(cmd_cat, shell=True, check=True)

        cmd_iq = f'java -Xmx{javaheapspace} -jar ./Astral/astral.5.7.8.jar -i gene_tree.treefiile -o species_tree.tree 2> species_nobs_tree.log'
        subprocess.run(cmd_iq, shell=True, check=True)

    else:
        cmd_cat = 'cat ./gene_tree/*.contree > gene_tree.contree'
        subprocess.run(cmd_cat, shell=True, check=True)

        cmd_find = 'find ./gene_tree/*.ufboot > gene_tree.ufboot'
        subprocess.run(cmd_find, shell=True, check=True)

        cmd_iq = f'java -Xmx{javaheapspace} -jar ./Astral/astral.5.7.8.jar -i gene_tree.contree -b gene_tree.ufboot -r {bootstrap} -o species_bs{bootstrap}_tree.contree 2> species_bs{bootstrap}_tree.log'
        subprocess.run(cmd_iq, shell=True, check=True)

        cmd_tree = f'tail -n 1 species_bs{bootstrap}_tree.contree > species_tree.tree'
        subprocess.run(cmd_tree, shell=True, check=True)


'''find potential marker gene'''
def rf_calculation(pip, bootstrap, species_tree):
    if pip and bootstrap == 0:
        genetree_files = [f for f in os.listdir("./gene_tree") if f.endswith(".treefile")]
    elif pip and bootstrap != 0:
        genetree_files = [f for f in os.listdir("./gene_tree") if f.endswith(".contree")]
    else:
        genetree_files = [f for f in os.listdir("./gene_tree") if f.endswith(".contree")]
        if genetree_files:
            pass
        else:
            genetree_files = [f for f in os.listdir("./gene_tree") if f.endswith(".treefile")]

    normalized_rf_dict = {}
    species_tree_path = os.path.join("./", species_tree)

    for genetree_file in genetree_files:
        genetree_path = os.path.join("./gene_tree", genetree_file)

        gene_tree = Tree(genetree_path, format=1)
        species_tree = Tree(species_tree_path, format=1)
        gene_tree.set_outgroup(gene_tree.get_midpoint_outgroup())
        species_tree.set_outgroup(species_tree.get_midpoint_outgroup())

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

    output_list = "marker.list"
    with open(output_list, "w") as f:
        f.write("gene\tRF\n")
        for key, value in sorted_rf:
            f.write(f"{key}\t{value:.2f}\n")


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

def extract_value_from_tsv(input_gene, sp):
    with open('single_copy_genes.tsv', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row[reader.fieldnames[0]] == input_gene:
                return row.get(sp, None)
    return None

def primer(target_gene, input_species):
    if not os.path.exists("./primer_design"):
        os.makedirs("./primer_design")

    extensions = ['.fasta', '.fas', '.fa', '.fna', '.seq', '.aln', '.maf', '.nt']
    for num, sp in enumerate(input_species):
        for ext in extensions:
            wgs_files = glob.glob(f"./species/{sp}{ext}")
            if wgs_files:
                wgs_file = wgs_files[0]
                break
        else:
            raise FileNotFoundError(f"No FASTA file found for {sp} in ./species/")

        geneMP = extract_value_from_tsv(target_gene, sp)

        cmd_sp = f"gffread -W -x {sp}.tmp.fna -g {wgs_file} ./gff/{sp}.gff"
        subprocess.run(cmd_sp, shell=True, check=True)

        cmd_sp_g = f"seqkit grep -p {geneMP} {sp}.tmp.fna > {sp}.tmp.{target_gene}.fna"
        subprocess.run(cmd_sp_g, shell=True, check=True)

        if num == 0: # first sample make exon fasta for primaer design
            split_fasta_by_segments(sp + '.tmp.' + target_gene + '.fna', sp + '_' + target_gene + '_split.fasta')

        ex_header = f"seqkit fx2tab -n {sp}.tmp.{target_gene}.fna | cat > {sp}_header"
        subprocess.run(ex_header, shell=True, check=True)

        with open(sp + '_header', "r") as file:
            query_header = file.read().strip()

        extract_region_from_fasta(wgs_file, query_header, sp + '_region.fasta')

        rmtmp = f"rm {sp}.tmp.fna {sp}.tmp.{target_gene}.fna {sp}_header"
        subprocess.run(rmtmp, shell=True, check=True)

    catfasta = f"cat *fasta > {target_gene}.fasta"
    subprocess.run(catfasta, shell=True, check=True)

    mafftfasta = f"mafft --quiet {target_gene}.fasta > {target_gene}_primerdesign.fasta"
    subprocess.run(mafftfasta, shell=True, check=True)

    mvfasta = f"mv {target_gene}_primerdesign.fasta ./primer_design/ && rm *fasta"
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


'''main commands'''
def prepare(threads, ref_cds, p):
    wgs2gff(threads, ref_cds)
    paf2tsv()
    gff2fasta(threads, p)

def curate(threads, pip):
    pre_scg1 = tsv2scg()
    pre_scg2 = gff2scg(pre_scg1, threads)
    scg_orthologs1 = gff2scg_df(pre_scg2, threads)
    fasta2scg(pre_scg2, scg_orthologs1, threads)
    scg_orthologs2 = last_trimming(pre_scg2, scg_orthologs1)
    print(scg_orthologs2)
    scg_orthologs2.to_csv("single_copy_genes.tsv", sep='\t')
    alignment_fasta(threads, pip)

def phylo(threads, bootstrap, species_tree):
    make_gene_tree(threads, bootstrap)
    if species_tree == "species_tree.tree":
        make_species_tree(threads, bootstrap)
    else:
        pass

def marker(pip, bootstrap, species_tree):
    print(species_tree)
    rf_calculation(pip, bootstrap, species_tree)

def pipeline(threads, ref_cds, p, bootstrap, species_tree):
    prepare(threads, ref_cds, p)
    pip = True
    curate(threads, pip)
    phylo(threads, bootstrap, species_tree)
    marker(pip, bootstrap, species_tree)

def edid(threads, insert, eject, restore):
    if insert:
        fasta_files = [f for f in os.listdir("./input") if f.endswith(('.fasta', '.fas', '.fa', '.fna', '.faa', '.pep', '.seq', '.aln', '.maf', '.nt', '.aa'))]
        with multiprocessing.Pool(threads) as pool_in:
            arg_in = [(os.path.splitext(fasta_file)[0], os.path.join("./input", fasta_file)) for fasta_file in fasta_files]
            pool_in.map(insert_id, arg_in)

    fasta_files = [f for f in os.listdir("./single_copy_genes") if f.endswith(".fasta")]
    if eject:
        with multiprocessing.Pool(threads) as pool_ej:
            arg_ej = [os.path.join("./single_copy_genes", fasta_file) for fasta_file in fasta_files]
            pool_ej.map(eject_id, arg_ej)
    if restore:
        tab_data = parse_tab_file("./single_copy_genes.tsv")
        with multiprocessing.Pool(threads) as pool_re:
            args_re = [(os.path.splitext(fasta_file)[0], os.path.join("./single_copy_genes", fasta_file), tab_data) for fasta_file in fasta_files]
            pool_re.map(update_fasta, args_re)
        cmd_rm = "rm ./single_copy_genes/sed*"
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

    parser_2 = subparsers.add_parser("curate", help="Generating single copy gene fasta and table from input.")
    parser_2.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")

    parser_3 = subparsers.add_parser("phylo", help="Generating phylogenetic tree from single copy gene fasta.")
    parser_3.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")
    parser_3.add_argument("-bs", "--bootstrap", type=int, default=1000, help="Number of bootstrap (default: 1000).")

    parser_4 = subparsers.add_parser("marker", help="Generating list of potential marker resion from tree file.")
    parser_4.add_argument("-sp", "--species_tree", type=str, default="species_tree.tree", help="reference species tree (default: species_tree.tree).")

    parser_5 = subparsers.add_parser("primer", help="Generate multi fasta for primer design.")
    parser_5.add_argument("-g", "--target_gene", required=True, metavar="STR", help="Name of target gene.")
    parser_5.add_argument("-i", "--input_species", type=str, nargs='+', help="Name of species use for desgin primaer.")

    parser_6 = subparsers.add_parser("pipeline", help="Run pipleline prepare, curate, phylo and marker.")
    parser_6.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")
    parser_6.add_argument("-r", "--ref_cds", required=True, metavar="STR", help="reference CDS fasta.")
    parser_6.add_argument("-p", action="store_true", help="Get amino acid seq of single copy genes (default: nucleotide seq).")
    parser_6.add_argument("-bs", "--bootstrap", type=int, default=1000, help="Number of bootstrap (default: 1000).")
    parser_6.add_argument("-sp", "--species_tree", type=str, default="species_tree.tree", help="reference species tree (default: species_tree.tree).")

    parser_7 = subparsers.add_parser("edid", help="Edit fasta id.")
    parser_7.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(), help="Number of threads to use (default: all available cores).")
    parser_7.add_argument("-insert", action="store_true", help="Insert sample_name to fasta_id.")
    parser_7.add_argument("-eject", action="store_true", help="Eject gene_id form fasta_id.")
    parser_7.add_argument("-restore", action="store_true", help="Restore sample_name_gene_id using tsv file.")

    if len(sys.argv) < 2 or sys.argv[1] in ["-h", "--help"]:
        print("cusco  Version=0.1.0\n\n\
Usage:\n\
    {} [command] option\n\n\
Command list:\n\
    prepare       Generate input files.\n\
    curate        Generate multi fasta and table of single copy genes.\n\
    phylo         Generate phylogenetic tree using single copy genes.\n\
    marker        Generate list of potential marker resion.\n\
    primer        Generate multi fasta for primer design.\
    pipeline      All processing is done automatically.\n\
    edid          Edit fasta id.\n".format(sys.argv[0]))
        sys.exit(0)

    args = parser.parse_args()

    if args.command == "prepare":
        prepare(args.threads, args.ref_cds, args.p)
    elif args.command == "curate":
        curate(args.threads, pip=False)
    elif args.command == "phylo":
        phylo(args.threads, args.bootstrap, species_tree="species_tree.tree")
    elif args.command == "marker":
        marker(pip=False, bootstrap=1000, species_tree=args.species_tree)
    elif args.command =="primer":
        primer(args.target_gene, args.input_species)
    elif args.command == "pipeline":
        pipeline(args.threads, args.ref_cds, args.p, args.bootstrap, species_tree=args.species_tree)
    elif args.command == "edid":
        edid(args.threads, args.insert, args.eject, args.restore)
    else:
        print("Error: {} [command] \n\n    [prepare][curate][phylo][marker][primer][pipeline][edid]\n".format(sys.argv[0]))
        sys.exit(1)

if __name__ == "__main__":
    main()
