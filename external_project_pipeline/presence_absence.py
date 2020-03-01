import os, re
import subprocess
import os.path
import configparser

class BlastAndCluster:
    def __init__(self):
        self.config = snakemake.config
        self.file = self.config["Parameters"]["name_dir"] + "/" + "all_putative_effectors_concatenated_clustered.fasta"
        self.db_file = self.config["Parameters"]["database_location"]
        self.e_value_thresh = str(self.config["Presence_absence"]["e_value_thresh"])
        genome_folder = self.config['Parameters']['input_folder']
        identitie_thresh = self.config["Presence_absence"]["ident_thresh"] 
        fasta_files = os.popen("ls " + genome_folder +" | egrep 'fasta|fa|fna'").read()
        genome_list = fasta_files.split("\n")[:-1]
        self.scores_for_genome = {}
        for x in range(len(genome_list)):
            self.make_blast_database(genome_list[x], genome_folder)
            self.local_blast(genome_list[x])
            self.parse_blastout2(genome_list[x].split(".")[0], identitie_thresh)
        self.score_to_genome(genome_list)
        self.outwrite()
       

    def make_blast_database(self, genome, genome_folder):
        """This function makes a database per genome

        Parameters:
           - self.db_file: the path where the database are placed
           - genome: 1 genome out of the genome folder """
        database_command = "makeblastdb -in " + genome_folder + "/" + genome +" -out " + self.db_file + genome.split(".")[0] +" -parse_seqids -dbtype nucl"
        os.system(database_command)
        print(database_command)

    def local_blast(self, genome):
        print("local_blast")
        self.outfile = self.file.split(".")[0] + "_VS_" + genome.split(".")[0] + "_blast.out"
        #blastn -db database/Fol4287sample -query Fol4287sample.fasta -evalue 0.001 -out result.out
        os.system("blastn -db " + self.db_file + genome.split(".")[0] + " -query " + self.file +" -evalue " + self.e_value_thresh + " -out " + self.outfile)

    def parse_blastout2(self, genome, identitie_thresh):
        """Checks in if founded effectors are also found in other genomes. If thats the case
        the signalpeptide of the effector is put into a dictionary with the a list of the genome names

        Parameters:
            - genome:  1 genome out of the genome folder
            - self.outfile: The blast file that is created in the local_blast function
            - self.scores_for_genome: dictionary with the signalpeptide of each effector, and value of each genome it is present."""
        outfile = open(self.outfile, "r")
        blastout = outfile.read()
        blast_list = blastout.split("Query=")[1:]
        for x in range(len(blast_list)):
            if "*" not in blast_list[x]:
                scores = blast_list[x].split("\n")
                query = scores[0]
                print(query)
                for y in range(len(scores)):
                    count_hit = 0
                    if "Identities" in scores[y]:
                        identities = int(scores[y].split("%")[0].split("(")[1])
                    if "Expect = " in scores[y]:
                        evalue = (scores[y].split("= ")[2])
                    if "Score =" in blast_list[x]:
                        count_hit = count_hit + 1
                    if "Score =" in blast_list[x] and count_hit == 2:
                        break
                if int(identities) >= int(identitie_thresh) and str(evalue) < str(self.e_value_thresh) or "e" in str(evalue) and int(identities) >= int(identitie_thresh):
                    if query not in self.scores_for_genome.keys():
                           self.scores_for_genome.update({query : [genome]})
                    else:
                           update_value =  self.scores_for_genome.get(query).append(genome)

    def score_to_genome(self, genome_list):
        """this function makes a dictionary, and compares the dic self.scores_for_genome
        if that dicionary has the genome name in it. The key for the self.score_dic is a 1
        if not the key becomes a 0. Eventually there is a key with a list.

        Parameters:
                - genome_list: a list with all the genome names
                - self.scores_for_genome: dictionary with the signalpeptide of each effector, and value of each genome it is present."""
        self.genome_type = [i.split('.', 1)[0] for i in genome_list]
        self.score_dic = {}
        self.peptide_string = ""
        count = 0
        for key, value in self.scores_for_genome.items():
            count = count + 1
            self.peptide_string =  self.peptide_string + "\t" + str(count) + "." + key
            score = []
            for x in range(len(self.genome_type)):
                if self.genome_type[x] in value:
                    score.append("1")
                else:
                    score.append("0")
            self.score_dic.update({key : score})

    def outwrite(self):
        """This funtion makes a tabulair output in the file precence_absence

        Parameters:
                - peptide_string: a string with the signalpeptide and a tab between them.
                - self.score_dic: a dictionary, the key is the signalpeptide and the value is a list with 0 or 1"""
        present_absence = open(snakemake.output.presence_absence, "w")
        present_absence.write(self.peptide_string + "\n")
        count = 0
        total_list = []
        for x in range(len(self.genome_type)):
            total_list.append([])
        for key , value in self.score_dic.items():
            for x in range(len(value)):
                total_list[x].append(value[x])
        for x in range(len(total_list)):
            present_absence.write(self.genome_type[x] + "\t" + "\t".join(total_list[x]))
            present_absence.write("\n")

       
             
               
           







if __name__ == '__main__':
    # moet  BlastAndCluster(genome_list[:-1]) worden zodra de rule ervoor klaar is
    BlastAndCluster()
