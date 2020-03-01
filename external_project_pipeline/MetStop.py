import os, re, sys
import configparser
import textwrap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


class MetStop:
    """
    The class MetStop uses the content in the attribute 'file', a file containing DNA sequence(s) from the MIMP Finder module.

    Firstly, all three frames are generated.
    Secondly, the open reading frames (ORF) are found with sequences larger than 23 aminoacids and shorter than 300 aminoacids.
    Thirdly, SignalP-5.0 is used to find possible signalpeptide sequences in the ORFs, this results in putative effectors.
    At last, useful headers are generated for each putative effector and the sequences are written down in a FASTA-file.
    """
    def __init__(self, file):
        self.__file = file
        self.datahandler_list = []
        for seq_record in SeqIO.parse(self.get_file(), "fasta", IUPAC.ambiguous_dna):
            self.datahandler_list.append(SeqRecord(
                                         seq=seq_record.seq,
                                         id=seq_record.id,
                                         description=""))
        self.parser()
        self.three_frame_translation()
        self.orf_finder()
        self.orf_writer()
        if self.signalversion == "4":
            self.signalp4()
            os.system("mv results.txt "  + snakemake.output.b)
            os.system("rm ORFS.fasta")
        else:
            self.signalp5()
            self.signalp_writer()
            print("Warning: deleting files: MetStop_mature.fasta, ORFS.fasta and results.txt")
            os.system("rm MetStop_mature.fasta ORFS.fasta results.txt")
            os.system("mv MetStop_summary.signalp5 " + snakemake.output.b)
        

    def get_file(self):
        """
        Returns self.__file
        :return: a string, returning the file name with DNA sequences from the MIMP Finder module
        """
        return self.__file

    def set_file(self, n_file):
        """
        It changes self.__file to n_file.
        :param n_file: a string, containing the new file name
        :return: none
        """
        self.__file = n_file
        return
    
    def parser(self):
        """
        The function parses the config file and returns the contents.
        :return: a string containing the file name
        """
        config = configparser.ConfigParser()
        config.read("config.ini")
        self.directory = config['Parameters']['name_dir']
        self.probability =  config['MetStop']['probability']
        self.min_protein_length =  int(config['MetStop']['min_protein_length'])
        self.max_protein_length =  int(config['MetStop']['max_protein_length'])
        self.signalversion = config["Parameters"]["signalp_version"]
        self.signalp_path = config["Parameters"]["signalp_path"]
        return

    def three_frame_translation(self):
        """
        The function finds every possible frames using three frame translation method.
        The sequences are not translated, yet.
        :return: none
        """
        self.translated_id = []
        self.nuc = []
        for x in range(len(self.datahandler_list)):
            header = self.datahandler_list[x].id
            nucleotide_string = self.datahandler_list[x].seq
            while len(nucleotide_string) % 3 != 0:
                nucleotide_string = nucleotide_string + "N"
            for frame in range(0,4):
                if header.find("strand:+") and frame != 0:
                    self.translated_id.append(">" + header + "|frame:" + "+" + str(frame))
                elif header.find("strand:-") and frame != 0:
                    self.translated_id.append(">" + header + "|frame:" + "-" + str(frame))
                elif frame == 0:
                    self.translated_id.append(">" + header + "|frame:" + str(frame))
                self.nuc.append(nucleotide_string[frame:])
        return


    def orf_finder(self):
        """
        The functions finds ORF sequences larger than 26 and short than 300 amino acids and writes it in a file.
        :return: none
        """
        self.orfs = []
        self.header_list = []
        count = 0
        for x in range(len(self.nuc)):
            for y in range(0, len(self.nuc[x]), 3):
                if self.nuc[x][y:(y+3)] == "ATG":
                    start = y
                    stop = self.nuc[x][start::]
                    for z in range(0, len(stop), 3):
                        if stop[z:(z+3)] == "TAA" or stop[z:(z + 3)] == "TAG" or stop[z:(z + 3)] == "TGA":
                            if len(self.nuc[x][start:start+z+1]) >= (self.min_protein_length * 3)\
                                    and len(self.nuc[x][start:start+z+1]) <= (self.max_protein_length * 3)\
                                    and self.nuc[x][start:start+z]\
                                    not in self.orfs:
                                self.orfs.append(self.nuc[x][start:start+z])
                                count += 1
                                self.header_list.append(">candidate_effector:%d|%s|len:%d\n" % (count,
                                                        self.translated_id[x][1::], len(self.nuc[x][start:start+z])))
                            break
        return

    def orf_writer(self):
        """
        The function writes the ORFs in a file.
        :return: none
        """
        orfs_file = open("ORFS.fasta", "w+")
        count = 0
        for x in range(len(self.orfs)):
            count += 1
            orfs_file.write(self.header_list[x])
            orfs_file.write(textwrap.fill(str(self.orfs[x].translate()), width=60) + "\n")                
        orfs_file.close()
        return

    def signalp4(self):
        """
        The function runs SignalP and finds putative effectors.
        :return: none
        """
        cline = self.signalp_path + ' -t euk -f summary -u %s %s > %s' % (self.probability, "ORFS.fasta", "results.txt")
        print("Running signalp4.......\n\n\n")
        print(cline)
        os.system(cline)

        file = open("results.txt", "r")
        content = file.readlines()
        cleavage = []
        id_list = []
        if len(open("results.txt").readlines()) == 0:
            for file in [snakemake.output.a, snakemake.output.c]:
                proteinseq_file = open(file, "w+")
                proteinseq_file.write("No signalpeptides has been found...")
            return
        for x in range(len(content)):
            if "SP='YES'" in content[x]:
                cleavage.append(content[x].split("Cleavage site between pos. ")[1].split(":")[0].replace(" and ", "-"))
                id_list.append(content[x].split("\t")[0].split("Name=")[1])
        for file in [snakemake.output.a, snakemake.output.c]:
            proteinseq_file = open(file, "w+")
            for x in range(len(self.header_list)):
                for y in range(len(id_list)):
                    if ">"+ id_list[y] == self.header_list[x].strip():
                        self.signalpeptideseq_AA = str(self.orfs[x].translate()[:int(cleavage[y].split("-")[0])])
                        self.signalpeptideseq = str(self.orfs[x][:int(cleavage[y].split("-")[0])])
                        if file == snakemake.output.a:
                            header = self.header_maker(True, self.header_list[x])
                            proteinseq_file.write("%s\n%s\n" % (header, textwrap.fill((self.signalpeptideseq_AA.lower()
                                                                        + str(self.orfs[x].translate())), width=60)))
                        if file == snakemake.output.c:
                            header = self.header_maker(False, self.header_list[x])
                            proteinseq_file.write("%s\n%s\n" % (header, textwrap.fill(self.signalpeptideseq.lower()
                                                                        + str(self.orfs[x]), width=60)))
            proteinseq_file.close()
        return
            
        

    def signalp5(self):
        """
        The function runs SignalP and finds putative effectors.
        :return: none
        """
        os.system("signalp -fasta ORFS.fasta -format short -mature -org 'euk' -prefix 'MetStop'")
        os.system("cat MetStop_summary.signalp5 | egrep 'CS pos' | awk '{if($NF > " + self.probability +
                  "){print($0)}}' > results.txt")
        return

    def signalp_writer(self):
        """
        The function writes the results of SignalP in a file.
        :return: none
        """
        file = open("results.txt", "r")
        content = file.readlines()
        cleavage = []
        print("yolo")
        #return
     
        for x in range(len(content)):
            cleavage.append(content[x].split(" ")[2][:-1]) #the cleavage of the signalpeptide
            content[x] = content[x].split("\t")[0] #id
        
        for file in [snakemake.output.a, snakemake.output.c]:
            proteinseq_file = open(file, "w+")
            for x in range(len(self.header_list)):
                for y in range(len(content)):
                    if ">"+ content[y] == self.header_list[x].strip():
                        self.signalpeptideseq_AA = str(self.orfs[x].translate()[:int(cleavage[y].split("-")[0])])
                        self.signalpeptideseq = str(self.orfs[x][:int(cleavage[y].split("-")[0])])
                        if file == snakemake.output.a:
                            header = self.header_maker(True, self.header_list[x])
                            print(textwrap.fill((self.signalpeptideseq_AA.lower())))
                            proteinseq_file.write("%s\n%s\n" % (header, textwrap.fill((self.signalpeptideseq_AA.lower()
                                                                        + str(self.orfs[x].translate())), width=60)))
                        if file == snakemake.output.c:
                            header = self.header_maker(False, self.header_list[x])
                            print(textwrap.fill(self.signalpeptideseq.lower()))
                            proteinseq_file.write("%s\n%s\n" % (header, textwrap.fill(self.signalpeptideseq.lower()
                                                                        + str(self.orfs[x]), width=60)))
            proteinseq_file.close()
        return

    def header_maker(self, boolean, header):
        """
        The function adds the correct length of the effector to the current header
        :return: returns the new header
        """
        header_list = header.split("|")
        length = int(re.search("len:(.*)", header).group(1))
        if boolean:
            new_header = "%s|len:%d" % ("|".join(header_list[0:4]), ((length / 3) + len(self.signalpeptideseq_AA)))
        elif boolean == False:
            new_header = "%s|len:%d" % ("|".join(header_list[0:4]), (length + len(self.signalpeptideseq)))
        h_list = new_header.split("|")
        if len(h_list[0] + "|" + self.signalpeptideseq_AA + "|" + "|".join(h_list[1:])) > 60:
           sp_header = h_list[0] + "|" + self.signalpeptideseq_AA + "|" + "".join(h_list[1])
        else:
           sp_header = h_list[0] + "|" + self.signalpeptideseq_AA + "|" + "|".join(h_list[1:])
        return sp_header


if __name__ == '__main__':
    sample = MetStop(snakemake.input.a)
