from typing import List, Dict, Iterable
from collections import deque
import csv
import sys
sys.setrecursionlimit(20000)
#################
#Read in Data
#################

def fasta_to_dict(file_path):
    fasta_dict = {}
    current_header = ""

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                current_header = line[1:]  # Remove the '>' character
                fasta_dict[current_header] = ""
            else:  # Sequence line
                fasta_dict[current_header] += line

    return fasta_dict


print("begin")
proj_2b_reads = '/Users/AdrianHanson/Downloads/project2b_reads.fasta'
# proj_2b_reads = '/Users/AdrianHanson/Downloads/project2_sample2_reads (1).fasta'


def dict_to_list(dic: dict):
    array = list()
    for i in  dic:
        array.append(dic[i])
    return array
    
proj_2b_reads = fasta_to_dict(proj_2b_reads)
proj_2b_reads_list = dict_to_list(proj_2b_reads)

# print(proj_2b_reads)

#################
#PROJECT
#################
def make_histogram(reads:dict, threshold:int, k:int):
    counts = dict()
    for read in reads: 
        this_read = reads[read]
        for i in range(len(this_read)-k):
            kmer = this_read[i:i+k]



            if kmer not in counts:
                counts[kmer] = 1
            elif kmer in counts:
                counts[kmer] += 1

    frequent_kmers = list()
    for read in counts:
        if counts[read] > threshold:
            frequent_kmers.append(read)

    return frequent_kmers

hist_data = make_histogram(proj_2b_reads, 5, 20)
# print(hist_data)


def de_bruijn_kmers(k_mers: List[str]):
    graph = {}
    for i in range(len(k_mers)):
        try:
            graph[k_mers[i][:-1]].append(k_mers[i][1:])
        except:
            graph[k_mers[i][:-1]] = [k_mers[i][1:]]
    
    # print(graph)
    return graph

def helper(current, G, visited_edges, current_path, all_paths):
    while current in G and G[current]:
        next_node = G[current][-1]
        edge = current + "->" + next_node

        if edge not in visited_edges:
            visited_edges.add(edge)
            G[current].pop()
            helper(next_node, G, visited_edges, current_path, all_paths)
        else:
            G[current].pop()

    current_path.append(current)


def eulerian_paths(G:dict):
    # Create a copy of the graph to modify during the search
    copy = G.copy()
    # Set to keep track of visited edges
    visited_edges = set()
    # List to store all Eulerian paths found
    all_paths = []

    # Iterate through each node in the graph
    for start_node in G:
        # While there are still edges to explore from the current node
        while copy[start_node]:
            # List to store the current path
            current_path = []
            # Call the helper function to explore the path from the current node
            helper(start_node, copy, visited_edges, current_path, all_paths)
            # Reverse the current path to get the correct order (since the helper function adds nodes in reverse order)
            current_path.reverse()
            # Add the current path to the list of all Eulerian paths
            all_paths.append(current_path)

    # Return the list of all Eulerian paths
    return all_paths




def genome_path(paths: List[list]):
    result = ''
    for path in paths:
        this_path = ''.join([e[0] for e in path])+path[-1][1:]
        result = result + this_path
    return result


def recover_read_numbers(text:str, spectrum_dict:dict):
    read_order_list = list()
    for i in range(len(text)):
        for key, value in spectrum_dict.items():
            if text[i:i+len(value)//5] == value[0:len(value)//5]:
                format = f'>{key}'
                format = format.split('/')[0]
                read_order_list.append(format)
    return read_order_list


def find_read_path(patterns: List[str], read_kmer_dict:dict):
    
    db = de_bruijn_kmers(patterns)
    path = eulerian_paths(db)
    text = genome_path(path)
    # print(text)
    # print(len(text))
    read_order = recover_read_numbers(text,read_kmer_dict)
    return read_order


#################
#################



# print("reads from solutions")
# print(sample_1_spectrum['read_210'])
# print(sample_1_spectrum['read_154'])
# print(sample_1_spectrum['read_181'])


read_list = make_histogram(proj_2b_reads, 5, 20)

path = find_read_path(read_list, proj_2b_reads)

print(path)

# Writing to a text file
with open('output.txt', 'w') as txt_file:
    for item in path:
        txt_file.write(item + '\n')

# Alternative approach with a single write operation
with open('output.txt', 'w') as txt_file:
    txt_file.writelines(item + '\n' for item in path)

print("Done")
# # Writing to the CSV file
# with open('predictions.csv', 'w', newline='') as csv_file:
#     # Step 4: Using csv.writer to write the list to the CSV file
#     writer = csv.writer(csv_file)
#     writer.writerow(path) # Use writerow for single list


# with open("its_done.csv", "w", newline='') as csv_file:
#     csv_writer = csv.writer(csv_file)

#     for snp in path:
#      csv_writer.writerow([snp])
