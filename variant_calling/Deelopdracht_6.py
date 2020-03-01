import os
import sys


def fasta():
    """
    Deze functie convert de fastq bestand naar een fasta bestand: "consensus.fa"
    """
    os.system("""cat """ + snakemake.input.consensus + """ | tr -d "\n" | tr "@" "\n" |
               tr "+" "\n" | egrep ^"[0-9]" |
               awk '{print "@"substr($0,0,5), "\t" substr($0,6,$NF)}' | tr "\t" '\n' |
               tr "@" ">" > ./Deelopdracht_6/consensus.fa""")
    return


def variaties():
    """
    De functie telt alle mogelijke puntmutaties uit de VCF-file en returnt een dictionary met de mutaties.
    """
    puntmutaties = {"AT": 0, "AC": 0, "AG": 0,
                    "TA": 0, "TC": 0, "TG": 0,
                    "CA": 0, "CG": 0, "CT": 0,
                    "GA": 0, "GC": 0, "GT": 0}

    os.system("cat " + snakemake.input.vcf + "| awk '{print $4, $5}' | egrep ^[ACTG] | cut -c -3 > "
                                  "./Deelopdracht_6/varianten_enzo.vcf")

    with open("./Deelopdracht_6/varianten_enzo.vcf", "r") as file:
        while True:
            mutatie = file.readline().replace(" ", "")[:2]
            if puntmutaties.get(mutatie) != None:
                value = puntmutaties.get(mutatie) + 1
                puntmutaties[mutatie] = value
            if not mutatie:
                break
    return puntmutaties


def inserties_deleties(inserties, deleties, seq_referentie, seq_assembly):
    """
    De functie bepaalt of er een insertie of deletie is aan de hand van de sequentie lengte in de FASTA- file van de
    referentie en de consenus FASTA-file van de assembly.

    Parameters:
         - inserties: integer, teller voor aantal insterties
         - deleties: integer, teller voor aantal deleties
         - seq_referentie: string, de sequentie van 1 read van de FASTA-file
         - seq_assembly: string, de sequentie van 1 read van de consensus FASTA-file
    """
    if len(seq_referentie) > len(seq_assembly):
        deleties = (len(seq_referentie) - len(seq_assembly)) + deleties
    elif len(seq_referentie) < len(seq_assembly):
        inserties = (len(seq_assembly) - len(seq_referentie)) + inserties
    return inserties, deleties

def wegschrijven(inserties, deleties, dic_mutaties):
    """
    De functie schrijft de resultaten weg naar een file

    Paramters:
         - inserties: integer, teller voor aantal insterties
         - deleties: integer, teller voor aantal deleties
         - dic_mutaties: een dictionary met alle puntmutaties.
    """
    file = open(str(snakemake.output), "w")
    file.write("Aantal insterties: " + str(inserties) + "\n")
    file.write("Aantal deleties: " + str(deleties) + "\n")
    file.write("Ratio deleties/inserties is: " + str(round(deleties / inserties, 2)) + "\n")
    for key, value in dic_mutaties.items():
        file.write("Aantal " + key[0] + "--->" + key[1] + ": " + str(value) + "\n")
    file.close()

def main():
    fasta()
    dic_mutaties = variaties()

    inserties = 0
    deleties = 0
    with open("./Deelopdracht_5/infected_consensus.fasta", "r") as referentie, open("Deelopdracht_6/consensus.fa",
                                                                                    "r") as assembly:
        while True:
            id_referentie = (referentie.readline()[:-1])
            seq_referentie = referentie.readline()[:-1]
            seq_assembly = assembly.readline()[:-1]
            inserties, deleties = inserties_deleties(inserties, deleties, seq_referentie, seq_assembly)
            if not id_referentie:
                break

    wegschrijven(inserties, deleties, dic_mutaties)

def usage_information():
    print("""
    Usage information: Deelopdracht_6.py
    ________________________________________________________________________________________________________________
    i. [Functie van het script]
    Het script creeert een txt-file en berekend:
    - het aantal deleties;
    - het aantal inserties;
    - de ratio deleties en inserties;
    - alle mogelijke puntmutaties.
   
    ii. [Manier van gebruiken van het script]
    Het script wordt aangeroepen met Snakemake. Dit gaat als volgt:
    snakemake --snakefile Snakefile_s1082644 Deelopdracht_6

    iii. [Voorbeeld van het starten van het script]
    snakemake --snakefile Snakefile_s1082644 Deelopdracht_6

    iv. [Overzicht van de naamgeving van de outputbestanden van het script]
    1.  "statistics_variants.txt";
    2.  "consensus.fa";
    3.  "varianten_enzo.vcf"
    ________________________________________________________________________________________________________________
    """)

try:
    if sys.argv[1] == "-h":
        usage_information()
except IndexError:
    main()
