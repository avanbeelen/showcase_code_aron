import os
from Bio import SeqIO
import configparser
from Bio import AlignIO
from Bio.Align import AlignInfo

def run_blast(file, dbname):
    '''Makes a blast database and then runs a blast. The output of the blast is stored in
        blast_results.out
        Paramaters:
            - file: the inputfile
            - dbname: the location of the blast database
    '''
    outfmt_command = "'6 qseqid sseqid pident qstart qend evalue bitscore qlen'"
    # index in outputfile:  0       1     2      3      4    5       6      7

    #make blastdb
    os_command_makedb = 'makeblastdb -in ' + file + ' -dbtype nucl -out ' + dbname + '/' + file.split('.')[0]
    print('---------Command to make blastdb:---------\n' + os_command_makedb)
    os.system(os_command_makedb)

    #run blastn
    os_command_blastn = "blastn -outfmt " + outfmt_command + " -db " + dbname + '/' + file.split('.')[0] + \
                        ' -query ' + file + ' -out blast_results.out'
    print('\n---------Command to blastn:---------\n' + os_command_blastn + '\n')
    os.system(os_command_blastn)


def find_clusters(PERC_IDENTITY_THRESH, LENGTH_THRESH):
    '''Makes clusters with the blast output. Filters the clusters on pers_thresh and length
        Parameters:
            - PERC_IDENTITY_THRESH: the threshold where a sequence is put into a cluster.
    '''
    dict_effector2homologs = {}
    with open('blast_results.out', 'r') as effector_file:
        lines_in_file = effector_file.readlines()
        for line in lines_in_file:
            tabs = line.strip().split('\t')
            if float(tabs[2]) > PERC_IDENTITY_THRESH and (int(tabs[4]) - (int(tabs[3])-1)) / float(tabs[7]) > LENGTH_THRESH:
                if tabs[0] in dict_effector2homologs.keys():
                    dict_effector2homologs[tabs[0]].add(tabs[1])  # add all BLAST-associated hits to this entry
                else:
                    dict_effector2homologs[tabs[0]] = set([tabs[1]])
                if tabs[1] in dict_effector2homologs.keys():
                    dict_effector2homologs[tabs[1]].add(
                        tabs[0])  # the other way around; also add the BLAST query as an association to each BLAST hit.
                else:
                    dict_effector2homologs[tabs[1]] = set([tabs[0]])
    return single_linkage(dict_effector2homologs)

def single_linkage(node_partners):
    '''Does a singel linkage algorithm on a pair of sequences. If they are in the same cluster they are put into the
    same key in the dictionary.
        Parameter:
            - node_partners: dictionary of the partners found in blast
    '''
    clusters = []
    nodes = list(node_partners.keys())
    while len(nodes) > 0:
        nodes = list(node_partners.keys())
        cluster = update_cluster(node_partners, nodes[0], set([]))
        clusters.append(cluster)
        nodes = node_partners.keys()
    return clusters

def update_cluster(node_partners, todo, cluster = set([])):
    '''Updates the dictionary of the clusters'''
    new_todo = set([])
    for t in node_partners[todo]:
        partners = node_partners[t]
        new_todo = new_todo.union(set(partners))
        cluster.add(t)
        del node_partners[t]
    new_todo = new_todo.difference(cluster)
    if len(new_todo) > 0:
        return update_cluster(node_partners, new_todo, cluster)
    else:
        return cluster

def find_longest_sequence_in_cluster(clusters, infile):
    '''Finds the longet sequence in the cluster en retuns it
        Parameters:
            - clusters: List of the clusters found in singel linkage.
            - infile: the input file
    '''
    #Make a dict of the sequences in the infile. With the headers as keys
    with open(infile, 'r') as handle:
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    list_longest_elements = [] #list is returned with longest sequece of the clusters.
    for c in clusters:
        print('Nodes in this cluster cluster', len(c))
        longest_element = [0, '', ''] #length, id, sequence
        for element in c:
            length_elem = len(record_dict[element].seq)
            if length_elem > longest_element[0]:
                longest_element[0] = length_elem
                longest_element[1] = element
                #longest_element[2] = record_dict[element].seq
        list_longest_elements.append(record_dict[element])
    return list_longest_elements

def consensus_sequence(clusters, infile):
    '''Finds the consensus equence of a cluster.'''
    with open(infile, 'r') as handle:
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    list_consenses_elements = [] #list is returned with consensus sequence of the clusters.
    for clust in clusters:
        print('Nodes in this cluster cluster', len(clust))
        if len(clust) > 1:
            consensus, id = make_consensus(clust, record_dict)
            list_consenses_elements.append([str(consensus), str(id)])
        else:
            one_list = list(clust)
            list_consenses_elements.append([str(record_dict[one_list[0]].seq), str(record_dict[one_list[0]].id)])
    return list_consenses_elements

def make_consensus(clust, record_dict):
    with open('temp.fasta', 'w') as file:
        for elem in clust:
            file.write('>' + str(record_dict[elem].id) + '\n')
            file.write(str(record_dict[elem].seq) + '\n')
            id = elem
    alignment = AlignIO.read('temp.fasta', 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    os.system('rm temp.fasta')
    return consensus, id

def write_long(list_elements, outputfile):
    with open(outputfile, 'w') as file:
        for record in list_elements:
            new_id = record.id.split("|")[1]
            record.id = new_id
            record.description = new_id
            SeqIO.write(record, file, 'fasta')

def write_con(list_elements, outputfile):
    with open(outputfile, 'w') as handle:
        for elem in list_elements:
            seq = elem[0]
            id = elem[1].split('')[1]
            handle.write('>' + id + '\n' + seq + '\n')

def main():
    config = configparser.ConfigParser()
    config.read('config.ini')

    infile = snakemake.input.metstop_eff_nuc
    PERC_IDENTITY_THRESH = int(config['Clustering']['pers_ident_thresh'])
    LENGTH_THRESH = float(config['Clustering']['length_thresh'])
    longest_or_consensus = config['Clustering']['longest_or_consensus']
    dbname = config["Parameters"]["database_location"]
    output_folder = config["Parameters"]["name_dir"]
    outputfile = snakemake.output.eff_clustered
    run_blast(infile, dbname)
    clusters = find_clusters(PERC_IDENTITY_THRESH, LENGTH_THRESH)

    if longest_or_consensus == 'longest':
        list_longest_cluster = find_longest_sequence_in_cluster(clusters, infile)
        write_long(list_longest_cluster, outputfile)
    elif longest_or_consensus == 'consensus':
        list_consensus_cluster = consensus_sequence(clusters, infile)
        write_con(list_consensus_cluster, outputfile)
    else:
        print('Not the right Length_or_Consensus variable in the config file. Must be longenst or consensus')
    os.system("mv blast_results.out " + output_folder + "/clustering_results.blastout")


main()
