import re, configparser
from os import path
from Bio import SeqIO
from Bio.Alphabet import IUPAC

def FindSequence(infile):
    '''
    picks the sequence from infile_dc. Than it makes a new prefix (id) for in the fastafile that is made from this sript.
    parameters: infile_dc: The file that is given in the configfile.
    return: sequence and the new id.
    '''
    full_id_sequence_list = []
    for seq_record in SeqIO.parse(infile,'fasta', IUPAC.ambiguous_dna):
        id_seq_list = []
        id_seq_list.append(seq_record.description)
        id_seq_list.append(seq_record.seq)
        full_id_sequence_list.append(id_seq_list)
    return full_id_sequence_list

def MimpFinder(motive, sequence_list, length):
    '''Finds the mimps in the sequence. When it finds a mimp it makes a flank of the length given in  the config file.
    than the flank is put into a list with the header.
    Parameters:
    ‘sequence’: a string, variable that contains only the sequence of the organism, without header
    ‘new_id’: a string, the id that belongs to the sequence
    '''
    header_flank_list = []
    tir_start_stop_list = []
    for id_seq in sequence_list:
        for match in re.finditer(motive, str(id_seq[1])):
            start_stop_mimp = str(match.start() + 1)+ '-' + str(match.end())
            eind_flank = int(match.end() + length)
            if eind_flank > len(id_seq[1]):
                eind_flank = len(id_seq[1])
            start_flank = int(match.end())
            start_stop_flank = str(start_flank + 1) + '-' + str(eind_flank)
            flank = id_seq[1][start_flank:eind_flank]
            header_flank_list.append(str('>' + id_seq[0] + '|mimp_downstreamregion:' + start_stop_flank +
                                '_(strand:+)_' + match.__getitem__(0) + '^' + start_stop_mimp + '\n' + flank + '\n'))
            tir_start_stop_list.append([match.start() + 1, match.end()])
    return header_flank_list, tir_start_stop_list

def MimpFinder_rc(motive, sequence_list, length):
    '''Finds the mimps in the sequence. When it finds a mimp it makes a flank of the length given in  the config file.
    than the flank is put into a list with the header.
    Parameters:
    ‘sequence’: a string, variable that contains only the sequence of the organism, without header
    ‘new_id’: a string, the id that belongs to the sequence
    '''
    header_flank_list = []
    tir_start_stop_list_rc = []
    for id_seq in sequence_list:
        for match in re.finditer(motive, str(id_seq[1])):
            start_stop_mimp = str(match.start() + 1) + '-' + str(match.end())
            eind_flank = int(match.start())
            start_flank = int(match.start() - length)
            if start_flank < 0:
                start_flank = 0
            start_stop_flank = str(start_flank) + '-' + str(eind_flank)
            flank = id_seq[1][start_flank:eind_flank]
            rev_match = reverse_complement(match.__getitem__(0))
            header_flank_list.append(str('>' + id_seq[0] + '|mimp_upstreamregion:' + start_stop_flank +
                                '_(strand:-,rev.transc.)_' + rev_match + '^' + start_stop_mimp + '\n' + flank + '\n'))
            tir_start_stop_list_rc.append([match.start() + 1, match.end()])
    return header_flank_list, tir_start_stop_list_rc

def TIR_finder(sequence_list, tir_start_stop_list, motive_rc):
    #contig in de header en tir compleet maken. inverted repeat vinden
    tir_list = []
    for id_seq in sequence_list:
        for start_stop in tir_start_stop_list:
            tir = id_seq[1][start_stop[0]-400:start_stop[0]]
            for tir_match in re.finditer(motive_rc, str(tir)):
                start = start_stop[0]-tir_match.start()
                stop = start_stop[1]
                tir_list.append(str('>tir_location:' + str(start) + '-' + str(stop) +
                                                 '\n' + id_seq[1][start:stop]+'\n'))
    return tir_list

def TIR_finder_rc(sequence_list, tir_start_stop_list, motive):
    tir_list_rc = []
    for id_seq in sequence_list:
        for start_stop in tir_start_stop_list:
            tir = id_seq[1][start_stop[0]:start_stop[0]+400]
            for tir_match in re.finditer(motive, str(tir)):
                start = start_stop[0]
                stop = start_stop[0]+tir_match.end()
                tir_list_rc.append(str('>tir_location:' + str(start) + '-' + str(stop) +
                                             '\n' + id_seq[1][start:stop] + '\n'))
    return tir_list_rc

def reverse_complement(dna):
    '''makes a reverse complement version sequence
    Parameters:
    ‘dna’: a string, the DNA sequence that needs to be in reverse complement
    '''
    complement = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def WriteToFile(header_flank_lijst, outputfile):
    """imput list of headers with flanks
    function writes each element in the list to a file.
    Parameters:
    ‘header_flank_lijst’: a string, the full header including start/stop flank, mimp match, start/stop mimp
    and the sequence of the flank that has to be written to a file.
    """
    text_output = str(outputfile).split('/')[-1]
    for header_flank in header_flank_lijst:
        with open(outputfile, 'a') as file:
            file.write(header_flank)
        print(header_flank.split('\n')[0] + ' is found and put into ' + text_output)

def main():
    config = configparser.ConfigParser()
    config.read(snakemake.input.a)  # Input a in snakemake file
    # reads the config file and put it into variables
    #file = config['Parameters']['file']
    file = snakemake.input.b
    print(file)
    """
    if path.exists(subfile + ".fasta") == True:
        file = subfile + "fasta"
    elif path.exists(subfile + ".fa") == True:
        file = subfile + "fa"
    elif path.exists(subfile + ".fna") == True:
        file = subfile + "fna"""
    motive = config['Mimpfinder']['motive']
    motive_rc = config['Mimpfinder']['motive_rc']
    length = int(config['Mimpfinder']['length'])
    outputfile = snakemake.output.a  # Output a in snakemake file.
    tir_file = snakemake.output.b

    # makes sequence list
    sequence_list = FindSequence(file)

    #Finds flanks
    header_flank_list, tir_start_stop_list = MimpFinder(motive, sequence_list, length)
    header_flank_list_rc, tir__start_stop_list_rc = MimpFinder_rc(motive_rc, sequence_list, length)

    #Writes flanks to file
    if len(header_flank_list) != 0 or len(header_flank_list_rc) != 0:
        WriteToFile(header_flank_list, outputfile)
        WriteToFile(header_flank_list_rc, outputfile)
    else:
        print("\nNo MIMP's motives have been found.\n\nThe program will stop...\n")
        exit()

    #Finds complete MIMP's. If not, it writes 'No complete MIMP's have been found...'
    tir_list = TIR_finder(sequence_list, tir_start_stop_list, motive_rc)
    if len(tir_list) > 0:
        WriteToFile(tir_list, tir_file)
    else:
        with open(snakemake.output.b, "a") as file:
            file.write("No complete MIMP's have been found in the normal flank...\n")
        file.close()

    # Finds complete MIMP's. If not, it writes 'No complete MIMP's have been found...'
    tir_list_rc = TIR_finder_rc(sequence_list, tir__start_stop_list_rc, motive)
    if len(tir_list_rc) > 0:
        WriteToFile(tir_list_rc, tir_file)
    else:
        with open(snakemake.output.b, "a") as file:
            file.write("No complete MIMP's have been found in the rc flank...\n")
        file.close()

main()

