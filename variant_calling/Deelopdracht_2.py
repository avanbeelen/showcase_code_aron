import sys

"""
Zodra de kwaliteitscore, van een read
gedurende 3 opeenvolgende basen, een gemiddelde Q-score (Phred) van 30 heeft dan worden die 3 opeenvolgende basen
getrimd met behulp van de sliding-window methode. De reden van een gemiddelde Q-score van 30 is, 
omdat het 99,99% (30 = -10 log10 (1/1000)) zeker is dat die 3 opeenvolgende basen juiste base calls hebben.
"""

def sliding_window_begin(paired):
    """
    De functie doet een sliding window trim methode aan het begin van de read.
    De read trimt de paren gelijdelijk.

    parameters:
    1. "paired": list, een lijst waar een paired end reads paar in zit, omgeven door een tuple
    """
    score = 0
    read1 = ""
    for x in range(0, len(paired[3][0][:-1])):
        for char in paired[3][0][:-1][x:(x + 3)]:
            score += ord(char) - 33
        if (score / 3) >= 30:
            index_find = paired[3][0][:-1].find(paired[3][0][:-1][x:(x + 3)])
            read1 = paired[1][0][:-1][index_find::]
            read2 = paired[1][1][:-1][index_find::]
            quality1 = paired[3][0][:-1][index_find::]
            quality2 = paired[3][1][:-1][index_find::]
            return read1, read2, quality1, quality2
        score = 0
    return paired[1][0][:-1], paired[1][1][:-1], paired[3][0][:-1], paired[3][1][:-1]

def sliding_window_eind(read1, read2, quality1, quality2):
    """
    De functie doet een sliding window trim methode aan het eind van de read.
    De read trimt de paren gelijdelijk.

    parameters:
    1. "read1": een string met de getrimde read uit paar (1/2)
    2. "read2": een string met de getrimde read uit paar (2/2)
    3. "quality1": de kwaliteitsscore uit paar (1/2)
    4. "quality2": de kwaliteitsscore uit paar (2/2)
    """
    
    n_quality1 = quality1
    score = 0
    
    for x in range(0, len(quality1)):
        for char in quality1[x:(x + 3)]:
            score += ord(char) - 33
        if (score / 3) >= 30:
            index_find = n_quality1.find(quality1[x:(x + 3)])
            read1 = read1[index_find::]
            read2 = read2[index_find::]
            quality1 = quality1[index_find::]
            quality2 = quality2[index_find::]
            return read1[::-1], read2[::-1], quality1[::-1], quality2[::-1]
        score = 0
    return read1[::-1], read2[::-1], quality1[::-1], quality2[::-1]
        

def wegschrijven(id1, id2, read1, read2, quality1, quality2):
    """
    De functie schrijft de paired end reads paren weg naar een file.

    parameters:
    1. "read1": een string met de getrimde read uit paar (1/2)
    2. "read2": een string met de getrimde read uit paar (2/2)
    3. "quality1": de kwaliteitsscore uit paar (1/2)
    4. "quality2": de kwaliteitsscore uit paar (2/2)
    5. "id1": een string met de header uit paar (1/2)
    6. "id2": een string met de header uit paar (2/2)
    """
    
    file1 = open(snakemake.output.trimmed1, "a")
    file2 = open(snakemake.output.trimmed2, "a")
    if len(read1) > 50 and len(read2) > 50:
        file1.write(id1)
        file2.write(id2)
        file1.write(read1 + "\n")
        file2.write(read2 + "\n")
        file1.write("+\n")
        file2.write("+\n")
        file1.write(quality1 + "\n")
        file2.write(quality2 + "\n")
    return

def main():
    paired_end_reads = []
    with open(snakemake.input.pair1) as pair1, open(snakemake.input.pair2) as pair2:
        for x in zip(pair1.readlines(), pair2.readlines()):
            paired_end_reads.append(x)
            if len(paired_end_reads) == 4:
                read1, read2, quality1, quality2 = sliding_window_begin(paired_end_reads)
                read1, read2, quality1, quality2 = sliding_window_eind(read1[::-1], read2[::-1], quality1[::-1], quality2[::-1])
                wegschrijven(paired_end_reads[0][0], paired_end_reads[0][1], read1, read2, quality1, quality2)
                paired_end_reads = []

def usage_information():
    print("""
    Usage information: Deelopdracht_2.py
    ________________________________________________________________________________________________________________
    i. [Functie van het script]
    Het script trimt de sequenties uit de paired end reads bestanden. Zodra de kwaliteitscore, van een read
    gedurende 3 opeenvolgende basen, een gemiddelde Q-score (Phred) van 30 heeft dan worden die 3 opeenvolgende basen
    getrimd met behulp van de sliding-window methode. De reden van een gemiddelde Q-score van 30 is, 
    omdat het gemiddeld 99,99% (30 = -10 log10 (1/1000)) zeker is dat die 3 opeenvolgende basen juiste base calls hebben.
   
    ii. [Manier van gebruiken van het script]
    Het script wordt aangeroepen met Snakemake. Dit gaat als volgt:
    snakemake --snakefile Snakefile_s1082644 Deelopdracht_2

    iii. [Voorbeeld van het starten van het script]
    snakemake --snakefile Snakefile_s1082644 Deelopdracht_2

    iv. [Overzicht van de naamgeving van de outputbestanden van het script]
    1.  trimmed_1.fastq;
    2.  trimmed_2.fastq.
    ________________________________________________________________________________________________________________
    """)

try:
   if sys.argv[1] == "-h":
      usage_information()
except:
   main() 
            
