#!/usr/bin/env python3

import mysql.connector
import pysam
import pybedtools
import subprocess

__author__ = 'Clemens Spielvogel'

class Assignment1:
    def __init__(self):
        # EDIT TO MY GENE
        self.gene = "DHCR7"

        # Konvertierung von Bam-Datei in Sam-Datei
        self.sam_file = pysam.AlignmentFile('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam', "rb")


    def fetch_gene_coordinates(self, genome_reference, file_name):
        print("Connecting to UCSC to fetch data")

        # Verbindung oeffnen
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password',
                                      db=genome_reference)

        # Cursor erstellen
        cursor = cnx.cursor()

        # Anfragefelder erstellen
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        # Datenbankabfrage erstellen und ausfuehren
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        cursor.execute(query)

        self.geneinfo = []
        with open(file_name, "w") as filehandle:
            for row in cursor:
                if row[0] == "DHCR7":
                    for attribtute in row:
                        filehandle.write(str(attribtute) + "\n")
                        self.geneinfo.append(attribtute)

        # Cursor und Verbindung schliessen
        cursor.close()
        cnx.close()

        print("Done fetching data")


    def get_sam_header(self):
        # Gibt die Liste mit den Headern des Files zurueck
        header = []
        for line in self.sam_file.header:
            print(line)
            header.append(line)
        return header


    def get_properly_paired_reads_of_gene(self):
        # Aufruf der samtools via bash und Zwischenspeicherung in einem extra File
        cmd = ["samtools flagstat HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam > outputfile.txt"]
        subprocess.call(cmd, shell=True)
        file = open("outputfile.txt", "r")
        lines = []
        for line in file:
            lines.append(line)
        target_line = lines[6]
        properly_paired_reads = target_line.split(" ")[0]

        # Textfile wieder leoschen
        cmd = ["rm outputfile.txt"]
        subprocess.call(cmd, shell=True)
        file.close()

        print(properly_paired_reads)

        return properly_paired_reads


    def get_gene_reads_with_indels(self):
        reads = []
        for read in self.sam_file:
            columns = str(read).split("\t")
            if "I" in str(columns[5]) or "D" in str(columns[5]):
                print(read)
                reads.append(read)

        return reads


    def calculate_total_average_coverage(self):
        a = pybedtools.BedTool('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam')
        coverage = a.genome_coverage(bg=True)
        print(coverage)
        return coverage


    def calculate_gene_average_coverage(self):
        pre_coverage = pybedtools.BedTool(self.geneinfo)
        coverage = pre_coverage.coverage(bg=True)
        print(coverage.head())
        return coverage.head()


    def get_number_mapped_reads(self):
        cmd = ["samtools flagstat HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam > outputfile.txt"]
        subprocess.call(cmd, shell=True)
        file = open("outputfile.txt", "r")
        lines = []
        for line in file:
            lines.append(line)
        target_line = lines[2]
        number_mapped_reads = target_line.split(" ")[0]

        cmd = ["rm outputfile.txt"]
        subprocess.call(cmd, shell=True)
        file.close()

        print(number_mapped_reads)

        return number_mapped_reads


    def get_gene_symbol(self):
        print(self.geneinfo[0])
        return self.geneinfo[0]


    def get_region_of_gene(self):
        print("Chromosome: {}\nStart: {}\nEnd: {}".format(self.geneinfo[2], self.geneinfo[3], self.geneinfo[4]))
        return [self.geneinfo[2], self.geneinfo[3], self.geneinfo[4]]


    def get_number_of_exons(self):
        print(self.geneinfo[6])
        return self.geneinfo[6]


    def print_summary(self):
        self.fetch_gene_coordinates("hg19", "gene.txt")
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        #self.calculate_total_average_coverage() # Pybed Modul funktioniert nicht?
        #self.calculate_gene_average_coverage() # Pybed Modul funktioniert nicht?
        self.get_number_mapped_reads()
        self.get_gene_symbol()
        self.get_region_of_gene()
        self.get_number_of_exons()



if __name__ == '__main__':
    print("Assignment 1")
    assignment1 = Assignment1()
    assignment1.print_summary()
