import sys

def maak_dictionary():
    """
    De functie maakt een dictionary waarbij de grote van de dictionary
    gelijk is de langste aantal reads. De dictionary wordt gebruikt voor
    latere QC-analyses.

    """
    dic = {}
    file = open(str(snakemake.input), "r")
    lijn = file.readlines()[1].strip()

    for z in range(len(lijn)):
        dic[z] = [0, 0, 0, 0]

    file.close()
    return dic

def QC_analyses(dic):
    """
    De functie doet QC-analyses zoals het berekenen van de totaal aantal reads,
    de kortste read, langste read, etc.

    Parameters:
    - dic: een dictionary waarbij de GC% per positie wordt opgeslagen.
    """
    reads = 0
    nucleotides = 0
    s_read = [1000]
    l_read = [0]
    A = 0
    T = 0
    C = 0
    G = 0
    x = 0
    
    with open(str(snakemake.input), "r") as data:
        for line in data.readlines():
            x += 1
            if x == 2:
                x = -2
                reads += 1
                read = line.strip()
                nucleotides += len(read)
                A += read.count("A")
                T += read.count("T")
                G += read.count("G")
                C += read.count("C")
                if len(read) < s_read[0]:
                    s_read[0] = len(read)
                if len(read) > l_read[0]:
                    l_read[0] = len(read)
                for n in range(len(dic)):
                    dic[n][0] += read.count("A", n, n + 1)
                    dic[n][1] += read.count("T", n, n + 1)
                    dic[n][2] += read.count("G", n, n + 1)
                    dic[n][3] += read.count("C", n, n + 1)
    return reads, nucleotides, s_read, l_read, A, T, C, G, dic

def wegschrijven(reads, nucleotides, s_read, l_read, A, T, C, G, dic):
    """
    De data van de QC-analyse wordt weggeschreven naar een txt-file.

    parameters:
    - reads: aantal reads
    - nucleotides: aantal nucleotides
    - s_read: de aantal van de korste read
    - l_read: de aantal can de langste read
    - A: aantal adenines
    - T: aantal thymines
    - C: aantal cytosines
    - G: aantal guanines
    - dic: een dictionary met GC% per positie
    """
    with open(str(snakemake.output), "w") as output:
        output.write("The number of reads = %d\n" % reads)
        output.write("The average read length = %.0fbp\n" % (nucleotides / reads))
        output.write("The shortest read = %dbp\n" % s_read[0])
        output.write("The longest read = %dbp\n" % l_read[0])
        output.write("The global GC-content = %.1f%%\n" % (100 / (A + T + C + G) * (C + G)))
        for z in range(len(dic)):
            output.write("The GC-content on position %d = %.1f%%\n" % (
                        z + 1, 100 / (dic[z][0] + dic[z][1] + dic[z][2] + dic[z][3]) * (dic[z][2] + dic[z][3])))
    
def main():
    dic = maak_dictionary()
    reads, nucleotides, s_read, l_read, A, T, C, G, dic = QC_analyses(dic)
    wegschrijven(reads, nucleotides, s_read, l_read, A, T, C, G, dic)

def usage_information():
    print("""
    Usage information: Deelopdracht_1.py
    ____________________________________________________________________________________________
    
    i. [Functie van het script]
    Het script doet een QC op de paired-end reads en schrijft de resultaten weg. De analyses die
    uitgevoerd worden zijn:
    - Het aantal reads.
    - De gemiddelde kengte van de reads.
    - De lengte van de kortse read.
    - De lengte van de langste read.
    - Het totale GC-percentage.
    - Het GC-percentage per positie.
   
    ii. [Manier van gebruiken van het script]
    Het script wordt aangeroepen door Snakemake. Dit gaat als volgt:
    snakemake --snakefile Snakefile_s1082644 Deelopdracht_1

    iii. [Voorbeeld van het starten van het script]
    snakemake --snakefile Snakefile_s1082644 Deelopdracht_1

    iv. [Overzicht van de naamgeving van de outputbestanden van het script]
    1. paired_1.QC;
    2. paired_2.QC.
    ____________________________________________________________________________________________
    """)

try:
   if sys.argv[1] == "-h":
      usage_information()
except:
   main() 
