#!/usr/bin/env python3

import sys
import argparse
import os
import subprocess
import random as rd
import matplotlib.pyplot as plt


RUNID = rd.randint(1000,9999)


def main():
    args = parseArgs()
    labels = []
    for x in range(0,len(args.input)):
        bam = mapReads(args.input[x], args.pair[x], args.fasta, args.cores, args.out)
        bed = bamToBed(bam, args.out)
        total_reads = countReads(bed)
        print("total reads:",total_reads)
        read_counts, feat_counts = getSaturationData(bed, total_reads, args.gff, args.stranded, args.out)
        name = name.split(args.input[x],"/")[-1]
        labels.append(name)
        plt.plot(read_counts, gene_counts)

    plt.legend(labels)
    plt.xlabel("Mapped Reads")
    plt.ylabel("features with >= 10 read(s) mapped")
    plt.savefig(outs+"/saturation_plot.png", dpi=600)



def parseArgs():
    parser = argparse.ArgumentParser(description="Script to produce an RNAseq saturation plot")

    ## required ##
    parser.add_argument("-i", "--input", nargs="+", type=str, help="fastq(.gz) file containing RNAseq reads", required=True)
    parser.add_argument("-f", "--fasta", type=str, help="reference fasta file to use for mapping", required=True)
    parser.add_argument("-g", "--gff", type=str, help="reference gff/bed file that contains features to test for saturation (ex. genes)", required=True)

    ## optional ##
    parser.add_argument("-p", "--pair", nargs="+", type=str, help="second fastq(.gz) if sample is paired-end", required=False)
    parser.add_argument("-o", "--out", type=str, help="directory to output data (default = '.') ", required=False)
    parser.add_argument("-c", "--cores", type=int, help="max cpu cores to use (default = '1') ", required=False)
    parser.add_argument("-s", "--stranded", action="store_true", help="use stranded information when counting reads", required=False)

    args = parser.parse_args()

    # check if required input files exist
    for fq in args.input:
        if not fileExists(fq):
            sys.exit("ERROR: "+fq+" file not found....exiting...\n")

    if not fileExists(args.fasta):
        sys.exit("ERROR: "+args.fasta+" file not found....exiting...\n")
    if not fileExists(args.gff):
        sys.exit("ERROR: "+args.gff+" file not found....exiting...\n")

    # checks if fastq 2 exists
    if args.pair is not None:
        if not fileExists(args.pair):
            sys.exit("ERROR: "+args.pair+" file not found....exiting...\n")
    
    else:
        args.pair = [None] * len(args.input)

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    #setup cores setting
    if args.cores is None:
        args.cores = 1
    
    if args.stranded is None:
        args.stranded = False

    return args


def fileExists(inFile):
    try:
        test = open(inFile, 'r')
    except Exception as e:
        return False
    
    return True


def mapReads(fq, fq2, fasta, cores, out):
    out_sam_file = out+"/"+str(RUNID)+"tmp.sam"

    command = ["bwa", "mem", "-t", str(cores), "-o", out_sam_file, fasta, fq]
    if fq2 is not None:
        command.append(fq2)

    subprocess.call(command)

    out_bam_mapped = out+"/"+str(RUNID)+"tmp.sam.mapped"
    subprocess.call(["samtools","view","-F","4","-b","-o",out_bam_mapped,"-@",str(cores),out_sam_file])

    os.remove(out_sam_file)

    return out_bam_mapped

def bamToBed(bam, out):
    out_bed = out+"/"+str(RUNID)+"tmp.allreads.bed"

    bed = open(out_bed,"w")
    subprocess.call(["bamToBed","-i",bam], stdout=bed)
    bed.close()
    os.remove(bam)

    return out_bed


def countReads(bed):
    read_count = 0
    with open(bed,"r") as bd:
        for line in bd:
            read_count+=1
    
    return read_count


def getSaturationData(bed, total_reads, gff, stranded, out):
    interval = total_reads//30
    threshold = 1

    read_counts = [0]
    gene_counts = [0]

    for x in range(interval,total_reads,interval):
        read_subset = rd.sample(range(0,total_reads),x)
        read_subset.sort()

        sub_bed = out+"/"+str(RUNID)+"tmp.subset.bed"

        with open(bed,"r") as bd:
            with open(sub_bed,"w") as sub:
                lines = bd.readlines()
                for num in read_subset:
                    sub.write(lines[num])
                # for num,line in enumerate(bd):
                #     if num in read_subset:
                #         sub.write(line)
                #         continue
        
        gene_count = countDiscoveredGenes(sub_bed, gff, threshold, stranded, out)
        print(x,"-",gene_count)
        read_counts.append(x)
        gene_counts.append(gene_count)
    
    os.remove(sub_bed)
    os.remove(bed)

    return read_counts, gene_counts


def countDiscoveredGenes(sub_bed, gff, threshold, stranded, out):

    command = ["bedtools","intersect","-c"]

    if stranded is not None:
        command.append("-s")

    command += ["-a",gff,"-b",sub_bed]

    count_bed = out+"/"+str(RUNID)+"tmp.count.bed"
    out_bed = open(count_bed,"w")
    subprocess.call(command,stdout=out_bed)
    out_bed.close()

    gene_count = 0
    with open(count_bed,"r") as bed:
        bed.seek(0)
        for line in bed:
            split_line = line.split("\t")

            count = int(split_line[9])

            if count >= threshold:
                gene_count += 1
    

    os.remove(count_bed)

    return gene_count



if __name__ == "__main__":
    main()