import subprocess
import re
import os
import configparser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF
import textwrap

class augustusPredictor:
    def __init__(self, file):
        self.file = file
        self.config = snakemake.config
        try: 
            dic = SeqIO.to_dict(SeqIO.parse(open(self.file), "fasta"))
        except:
            self.make_unique_headers()
        self.augustus_runnen()
        self.signalPpath = ""
        self.sequentie()
        if self.config["Parameters"]["signalp_version"] == "4":
            self.RunSignalP4()
        else:
            self.RunSignalP5()
        #input()
        self.ExtractOrfToFasta()
        self.write_out()
        self.remove()

    def make_unique_headers(self):
        count = 0
        new_file = open(self.file + "2", "w")
        with open(self.file) as f:
           for line in f:
               if line.startswith(">"):
                   count = count + 1
                   new_header = line.split(" ")
                   extra_info = new_header[-1].split("|")[1]
                   new_file.write(new_header[0] + str(count) + "|" + extra_info)
               else:
                   new_file.write(line)
               if line is None:
                  break
        new_file.close()
        self.file = self.file + "2"


    def augustus_runnen(self):
        """
        Runs Augustus to predict genes

        Parameters:
             - "self.file: The name of the file with to mimps
             -  Self.file[:-6]: Name of the file without .fasta
        """
        self.directory=self.file.split("/")[1].split(".")[0]
        command = ("augustus --species=fusarium --singlestrand=true " +
                   self.file + " > " + self.directory + "_augustus.gff")
        print("running augustus......")
        os.system(command)


    def sub_feature_extract(self):
        """Loops trought the augustus file and checks if a sentace is a CDS
        if that is the case the exonlocation are saved. Returns the CDS sequence
        exon locations and the strand (normal or reverse)

        Parameters:
            - sub_feature: What augustus found, can be CDS, transcript or stop_codon
            - self.rec: the sequence of the founded gene
            - sub_feature.location.start: the start location of the sequence
            - ub_feature.location.end: The end location of the sequence
        """
        self.CDS = ""
        self.exonlocations = []
        self.feature_stand = ""
        for sub_feature in self.feature.sub_features:
            if (sub_feature.type) == "CDS" and sub_feature.strand == 1:
                self.CDS += sub_feature.extract(self.rec).seq
            elif (sub_feature.type) == "CDS" and sub_feature.strand == -1:
                self.CDS += self.feature.extract(self.rec).seq.reverse_complement()
            if (sub_feature.type == "CDS"):
                self.exonlocations.append(
                    [int(sub_feature.location.start), int(self.feature.location.end)])
            self.feature_stand = sub_feature.strand


    def introns_exons_extract(self):
        """
        Finds the intron locations by looping through exonlocations
        intron location is the location between 2 exons. return all_location and intronlocation

        Parameters:
            - self.exonlocations: a list with the exonlocations
        """
        all_locations = []
        intronlocations = []
        for x in range(len(self.exonlocations)):
            all_locations.append(self.exonlocations[x][0])
            all_locations.append(self.exonlocations[x][1])
        for y in range(len(self.exonlocations)-1):
            intronlocations.append(
                [all_locations[(y*2)+1], all_locations[(y+1)*2]])
        self.all_locations= all_locations
        self.intronlocations = intronlocations[1:-1]


    def ORF_finder(self):
        """
        Finds Open reading frame by looping through te exons and introns.
        if there are no intron locations the ORF is between the start exon location
        and end. Returns the ORF

        Parameters:
             - self.intronlocations: List of the intronlocations
             - self.rec.seq: the sequence of a found gene by augustus
             - self.exonlocations: list with the positions of exons
             - self.intronlocations: list with the positions of intros
        """
        ORF = ""
        for x in range(len(self.intronlocations)):
            ORF += self.rec.seq[self.exonlocations[x][0]:self.exonlocations[x][1]].upper()
            ORF += self.rec.seq[self.intronlocations[x][0]:self.intronlocations[x][1]].lower()
        if self.exonlocations != []:
            ORF += self.rec.seq[self.exonlocations[len(self.exonlocations)-1]
                           [0]:self.exonlocations[len(self.exonlocations)-1][1]].upper()
        self.ORF = ORF


    def header_info(self):
        """
        Splits the header, so the useful information is kept

        Parameters:
                - self.rec.id: The ID of a predicted gene
                - self.feature.location.start: The start of the gene
                - self.feature.location.end: The end of the gene
        """
        self.contig = self.rec.id.split("|")[0]
        self.region_start = self.rec.id.split(":")[1].split("-")[0]
        self.region_end = self.rec.id.split(":")[1].split("-")[1].split("_")[0]
        if "rev.transc" in self.rec.id:
            self.strand = "-"
        else:
            self.strand = "+"
        self.mimp_IR_seq =  self.rec.id.split("-")[1].split("^")[0]
        self.mimp_IR_loc = self.rec.id.split("^")[1].split("-")[1]
        self.gene_start_relative = self.feature.location.start
        self.gene_end_relative = self.feature.location.end
        


    def real_orientation_gen(self):
        """checks if the CDS and ORFS are good orientated. If not so it is changed

        Parameters:
                 - self.feature_stand: The side the sequence is read, can be + or -
                 - self.strand: The side the sequence is read, can be + or -
                 - self.region_start: the start of the sequence
                 - self.gene_start_relative: the start of the gene
                 - self.region_end: end of the sequence
                 - self.region_end_relative: end of the gene"""
        self.real_orientation = ""
        self.gene_start_absolute = ""
        self.gene_end_absolute = ""
        if self.feature_stand == 1 and self.strand == "+":
            self.real_orientation = 1  # gene is found on coding strand downstream of MIMP IR
            self.gene_start_absolute = int(self.region_start)+int(self.gene_start_relative)
            self.gene_end_absolute = int(self.region_start)+int(self.gene_end_relative)-1
        elif self.feature_stand == 1 and self.strand == "-":
            self.real_orientation = -1  # gene is found on other strand upstream of MIMP IR
            self.gene_start_absolute = int(self.region_end)-int(self.gene_start_relative)
            self.gene_end_absolute = int(self.region_end)-int(self.gene_end_relative)+1
        elif self.feature_stand == -1:
            self.ORF = self.ORF.reverse_complement()
            self.CDS = self.CDS.reverse_complement()
            if self.strand == "+":
                self.real_orientation = -1
                self.gene_start_absolute = int(self.region_start)+int(self.gene_start_relative)
                self.gene_end_absolute = int(self.region_start)+int(self.gene_end_relative)-1
            elif self.strand == "-":
                self.real_orientation = 1
                self.gene_start_absolute = int(self.region_end)-int(self.gene_start_relative)
                self.gene_end_absolute = int(self.region_end)-int(self.gene_end_relative)+1


    def protein_sequence(self):
        """Checks if CDS is filled (found) if not the case a translated
        CDS is returned

        Parameters:
              - self.CDS: the Coding Sequence of the gene.
        """
        if self.CDS == "":
            self.proteinseq = "incomplete"
        else:
            self.proteinseq = self.CDS.translate()[:-1]


    def met_check(self):
        """Check if the translated CDS starts with a methiodine, and if
        the protein is longer then 25 AA and shorter then 600 AA. Return
        ORFS, CDS, PROT that is in the specific range. The length of protein
        are defined in the config file, and can be changed.

        Parameters:
              - self.proteinseq: CDS sequence translated to a protein sequence
              - Self.gene_start_relative: the start of the gene
              - self.ORF: the Open reading frame of a gene
              - self. CDS: the coding sequence of a gene
              - self.new_fastaheader_id: a new header for the founded CDS
              - self.new_fastaheader_description: a new description for the header
              - min_len: minimum length of protein, defined in the config file
              - max_len: maximum length of protein, defined in the config file
              - d2m: distance that the protein maximum can have from the mimp, defined in the config file
        """
        min_len = self.config["Augustus"]["min_protein_length"]
        max_len = self.config["Augustus"]["max_protein_length"]
        d2m = self.config["Augustus"]["max_distance_mimp"]
        self.new_ORF = ""
        self.new_CDS = ""
        self.new_PROT = ""
        if str(self.proteinseq[:1]) == 'M' and len(self.proteinseq) > int(min_len) and len(self.proteinseq) < int(max_len):
            if int(self.gene_start_relative) < int(d2m):
                self.new_ORF = (SeqRecord(seq=self.ORF, id=self.new_fastaheader_id,
                                     description=self.new_fastaheader_description))
                self.new_CDS = (SeqRecord(seq=self.CDS, id=self.new_fastaheader_id,
                                     description=self.new_fastaheader_description))
                self.new_PROT = (SeqRecord(seq=self.proteinseq, id=self.new_fastaheader_id,
                                      description=self.new_fastaheader_description))


    def feature_extract(self):
        """Loops through the features of each gen. And connects the different functions.
        Returns a list with all the proteins sequence, ORF and CDS

        Parameters:
                - self.rec.features: Seq info about the founded gene
                - self.feature: More info about the founded gene example the location
                - self.contig: contig
                - self.gene_start_absolute: the start of the sequence
                - self.gene_end_absolute: the start of the sequence
                - self.self.gene_start_relative: Distance to mimp
                - self.real_orientation: the strand of the sequence
                - self.proteinseq: protein sequence
                - self.mimp_IR_seq: signalpeptide
                - self.mimp_IR_loc: location
        """
        for self.feature in self.rec.features:
            a = self.feature.sub_features
            self.sub_feature_extract()
            self.introns_exons_extract()
            self.ORF_finder()
            self.header_info()
            self.real_orientation_gen()
            self.protein_sequence()
            self.new_fastaheader_id = self.contig+'|'+str(self.gene_start_absolute)+'-'+str(self.gene_end_absolute)+'|'+'d2m:'+str(
                self.gene_start_relative)+'bp|'+str(self.real_orientation)+'|f?'+'|len:'+str(len(self.proteinseq))
            self.new_fastaheader_description = str(self.mimp_IR_seq)+'^'+self.mimp_IR_loc
            self.met_check()
            if type(self.new_ORF) != str:
                self.CDS_list.append(self.new_CDS)
                self.ORF_list.append(self.new_ORF)
                self.PROT_list.append(self.new_PROT)
        


    def sequentie(self):
        """Opens the augustus output file, gives every gen info to the function
        feature_extract by using a loop. From that function is get a compleet list of CDS, PROT and ORF
        finaly it writes is to files.

        Parameters:
               - self.file: the file with the mimps
               - self.ORF_list: A list with al the ORFS and the headers
               - self.CDS_list: A list with al the CDS and the headers
               - self.PROT_list: A list with al the proteins and the headers
        """
        self.ORF_list = []
        self.CDS_list = []
        self.PROT_list = []
        in_seq_handle = open(self.file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        gff_file = open(self.directory + "_augustus.gff")
        for self.rec in GFF.parse(gff_file, base_dict=seq_dict):
            self.feature_extract()
        #SeqIO.write(self.ORF_list, snakemake.output.augustus_ORFS, 'fasta')
        SeqIO.write(self.CDS_list, "CDS", 'fasta') #Augustus CDS file
        SeqIO.write(self.PROT_list, "PROTEIN", 'fasta') #Augustus PROTEIN FILE
        self.CDS_list

    def RunSignalP4(self):
        SignalPpath = self.config["Parameters"]["signalp_path"]
        probability = self.config["Augustus"]["probability"]
        naamlijst = []
        self.signalpep_list = []
        self.sig = []
        print('// Running SignalP 4.1...')
        print("\n....\n........\n........")
        cline = SignalPpath +' -t euk -f summary -u %s %s > %s' % (probability, "PROTEIN", "result.txt")
        os.system(cline)
        print(cline)
        
        file = open("result.txt", "r")
        id_line = file.readlines()
        self.cleavage = []
        self.id_line = []
        for x in range(len(id_line)):
            try:
                self.cleavage.append(id_line[x].split("between pos. ")[1].split(":")[0].replace(" and ", "-"))
                self.id_line.append(id_line[x].split("Name=")[1].split("\t")[0])
            except:
                pass
        self.CDS = SeqIO.parse("PROTEIN", 'fasta')
        for y in self.CDS:
            for x in range(len(self.id_line)):
                if y.id == self.id_line[x]:
                    self.seq = str(y.seq)
                    self.signalpep_list.append(self.seq[:int(self.cleavage[x].split(
                        "-")[0])].lower() + self.seq[int(self.cleavage[x].split("-")[0]):])
                    self.sig.append(self.seq[:int(self.cleavage[x].split("-")[0])])


    def RunSignalP5(self):
        """Runs signalP5 and parses it so that the signalpeptide sequence is found
        Returns the cleavage positions, signalpeptide sequence, signalpeptide + CDS, and header

        Parameters:
                - self.signalPpath: the location where signalP is installed
        """
        probability = self.config["Augustus"]["probability"]
        naamlijst = []
        self.signalpep_list = []
        self.sig = []
        print ('// Running SignalP 5...')
        os.system(self.signalPpath +
                  "signalp -fasta PROTEIN -format short -mature -org 'euk' -prefix test")
        os.system(
            "cat test_summary.signalp5 | egrep 'CS pos' | awk '{if($NF > " + probability +" ){print($0)}}' > result.txt")
        file = open("result.txt", "r")
        self.id_line = file.readlines()
        self.cleavage = []
        for x in range(len(self.id_line)):
            self.cleavage.append(self.id_line[x].split(" ")[2])
            self.id_line[x] = self.id_line[x].split("\t")[0]
        self.CDS = SeqIO.parse("PROTEIN", 'fasta')
        print(self.cleavage, self.id_line)
        for y in self.CDS:
            for x in range(len(self.id_line)):
                if y.id == self.id_line[x]:
                    self.seq = str(y.seq)
                    self.signalpep_list.append(self.seq[:int(self.cleavage[x].split(
                        "-")[0])].lower() + self.seq[int(self.cleavage[x].split("-")[0]):])
                    self.sig.append(self.seq[:int(self.cleavage[x].split("-")[0])])


    def ExtractOrfToFasta(self):
        """Uses the signalpeptide sequence, CDS and cleavage site to find the putative
        effectors. Return effector header and effector sequence

        Parameters:
             self.id_line: the id of the lines from the file with the mimps
             self.sig: signalpeptide sequence
             self.cleavage: cleavage posistion
        """
        count = 10000
        self.CDS = SeqIO.parse("CDS", 'fasta')
        self.effector_header = []
        self.effector_seq = []
        for y in self.CDS:
            for x in range(len(self.id_line)):
                if y.id == self.id_line[x]:
                    count = count + 1
                    header = self.id_line[x].split("|")
                    if header[3] == "1":
                        nuc_seq = str(y.seq)
                        self.effector_header.append(str(count)[
                                               1:] + "|" + self.sig[x] + "|" + header[2] + "|" + header[-1]+"|" + header[3])
                        self.effector_seq.append(nuc_seq[:int(self.cleavage[x].split(
                            "-")[0])].lower() + nuc_seq[int(self.cleavage[x].split("-")[0]):])
                    else:
                        nuc_seq = str(y.seq.reverse_complement())
                        self.effector_header.append(str(count)[
                                               1:] + "|" + self.sig[x] + "|" + header[2] + "|" + header[-1] + "|" + header[3])
                        self.effector_seq.append(nuc_seq[:int(self.cleavage[x].split(
                            "-")[0])].lower() + nuc_seq[int(self.cleavage[x].split("-")[0]):])


    def write_out(self):
        """Writes Signalpeptide sequence to file and writes the putative effectors to file

        Parameters:
               - self.file: File with the mimps
               - self.id_line: the id of the lines from the file with the mimps
               - self.signalpep_list: list with al the signalpeptides
               - self.effector_seq: list with sequence of al the effectors"""
        SP_file = open(snakemake.output.augustus_effectors_prot, "w")
        EF_file = open(snakemake.output.augustus_effectors_nuc, "w")
        for x in range(len(self.id_line)):
            SP_file.write(">" + self.id_line[x] + "\n")
            SP_file.write(re.sub("(.{60})", "\\1\n", str(
                self.signalpep_list[x]), 0, re.DOTALL) + "\n\n")
            EF_file.write(">" + self.effector_header[x] + "|" + self.id_line[x].split("|")[0] + "\n")
            EF_file.write(re.sub("(.{60})", "\\1\n", str(
                self.effector_seq[x]), 0, re.DOTALL) + "\n\n")

    def remove(self):
        """Removes temporary files"""
        os.remove("CDS")
        os.system("rm *.gff")
        os.remove("PROTEIN")
        if self.config["Parameters"]["signalp_version"] != "4":
            os.remove("result.txt")
            os.remove("test_mature.fasta")
            os.system("mv test_summary.signalp5 " + snakemake.output.augustus_signalp)
        else:
            os.system("mv result.txt " + snakemake.output.augustus_signalp)

if __name__ == '__main__':
    sample = augustusPredictor(snakemake.input.founded_mimps)

