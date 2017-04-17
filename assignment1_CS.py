#!/usr/bin/env python3

'''
Assignment 1 - Medizinische Genomanalyse
von Clemens Spielvogel

Zusaetzlich zu den importierten Modulen muessen bedtools fuer Linux auf dem System installiert sein.
(sudo apt-get install bedtools)

VOR DEM AUSFUEHREN muss der absolute Pfad der Bam-Datei fuer die Variable hg_bam (Zeile 19) angegeben werden!
'''

import mysql.connector
import pysam
import pybedtools
import subprocess

__author__ = 'Clemens Spielvogel'
hg_bam = '/home/vortex/Bioinformatik/med_genomanalyse/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam'

class Assignment1:
    def __init__(self, hg_bam):
        self.gene = "DHCR7"

        # Absoluten Dateipfad der Bam-Datei angeben
        self.hg_bam = hg_bam

        # Konvertierung von Bam-Datei in Sam-Datei
        self.sam_file = pysam.AlignmentFile(self.hg_bam, "rb")


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
        print("\n---------------\nSam Header:")

        # Gibt die Liste mit den Headern des Files zurueck
        header = []
        for line in self.sam_file.header:
            print(line)
            header.append(line)

        return header


    def get_properly_paired_reads_of_gene(self):
        print("\n---------------\nProperly paired reads of gene:")

        # Aufruf der samtools via bash und Zwischenspeicherung in einem extra File
        cmd = ["samtools flagstat {} > outputfile.txt".format(self.hg_bam)]
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


    def get_gene_reads_with_indels(self):   # Gibt aus uebersichtlichkeitsgruenden die Anzahl aus und returned die reads
        print("\n---------------\nGene reads with indels:")

        read_count = 0
        reads = []
        for read in self.sam_file:
            columns = str(read).split("\t")
            if "I" in str(columns[5]) or "D" in str(columns[5]):
                read_count += 1
                reads.append(read)

        print(read_count)
        return reads


    def calculate_total_average_coverage(self):
        print("\n---------------\nTotal average coverage:")

        a = pybedtools.BedTool(self.hg_bam)
        b = a.genome_coverage(bg=True)

        line_count = 0  # Anzahl der Regionen
        total_coverage = 0  # Summierte Coverage ueber alle Regionen
        for line in b:
            line_count += 1
            total_coverage += int(line[3])

        tot_avg_cov = total_coverage / line_count
        print(round(tot_avg_cov, 2))
        return tot_avg_cov


    def calculate_gene_average_coverage(self):
        print("\n---------------\nGene average coverage:")

        a = pybedtools.BedTool(self.hg_bam)
        b = a.genome_coverage(bg=True)

        line_count = 0  # Anzahl der Regionen
        total_coverage = 0  # Summierte Coverage ueber alle Regionen
        for line in b:
            if int(line[1]) >= int(self.geneinfo[3]) and int(line[2]) <= int(self.geneinfo[4]):
                line_count += 1
                total_coverage += int(line[3])

        gene_avg_cov = total_coverage / line_count
        print(round(gene_avg_cov, 2))
        return gene_avg_cov


    def get_number_mapped_reads(self):
        print("\n---------------\nNumber of mapped reads:")

        cmd = ["samtools flagstat {} > outputfile.txt".format(self.hg_bam)]
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
        print("\n---------------\nGene symbol:")

        print(self.geneinfo[0])
        return self.geneinfo[0]


    def get_region_of_gene(self):
        print("\n---------------\nRegion of gene:")

        print("Chromosome: {}\nStart: {}\nEnd: {}".format(self.geneinfo[2], self.geneinfo[3], self.geneinfo[4]))
        return [self.geneinfo[2], self.geneinfo[3], self.geneinfo[4]]


    def get_number_of_exons(self):
        print("\n---------------\nNumber of exons:")

        print(self.geneinfo[6])
        return self.geneinfo[6]


    def print_summary(self):
        self.fetch_gene_coordinates("hg19", "gene.txt")
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.calculate_total_average_coverage()
        self.calculate_gene_average_coverage()
        self.get_number_mapped_reads()
        self.get_gene_symbol()
        self.get_region_of_gene()
        self.get_number_of_exons()



if __name__ == '__main__':
    print("Assignment 1")
    assignment1 = Assignment1(hg_bam)
    assignment1.print_summary()
