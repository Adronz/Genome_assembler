from typing import List, Dict, Iterable
from collections import deque
import csv
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

annotated_spectrum_path = '/Users/AdrianHanson/Downloads/project2_sample1_spectrum_with_positions.fasta'
annotated_sample_spectrum = fasta_to_dict(annotated_spectrum_path)

sample_1_spectrum = '/Users/AdrianHanson/Downloads/project2_sample1_spectrum.fasta'
sample_1_spectrum = fasta_to_dict(sample_1_spectrum)

proj_2a_spec = '/Users/AdrianHanson/Downloads/project2a_spectrum.fasta'
proj_2a_spec = fasta_to_dict(proj_2a_spec)
# print(sample_1_spectrum)


proj_2b_reads = '/Users/AdrianHanson/Downloads/project2_sample2_reads (1).fasta'

def dict_to_list(dic: dict):
    array = list()
    for i in  dic:
        array.append(dic[i])
    return array
    
proj_2b_reads = fasta_to_dict(proj_2b_reads)
proj_2b_reads_list = dict_to_list(proj_2b_reads)
read_list = dict_to_list(proj_2a_spec)


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
            frequent_kmers.append(counts[read])

    return frequent_kmers

hist_data = make_histogram(proj_2b_reads, 5, 20)
# print(hist_data)


def de_bruijn_kmers(k_mers: List[str]):
    """Forms the de Bruijn graph of a collection of k-mers."""
    k = len(k_mers[0])
    graph = {}

    for i in range(len(k_mers)):
        try:
            graph[k_mers[i][:-1]].append(k_mers[i][1:])
        except:
            graph[k_mers[i][:-1]] = [k_mers[i][1:]]
    return graph



def ensure_all_vertices_present(g: Dict[int, List[int]]):
    # Find all vertices mentioned in the graph
    all_vertices = set(g.keys())
    for neighbors in g.values():
        all_vertices.update(neighbors)
    
    # Ensure each vertex has an entry in the graph dictionary
    for vertex in all_vertices:
        g.setdefault(vertex, [])

def find_start_vertex(g: Dict[int, List[int]]) -> int:
    # Calculate in-degrees
    in_degrees = {u: 0 for u in g}
    for u in g:
        for v in g[u]:
            in_degrees[v] = in_degrees.get(v, 0) + 1

    # Determine the start vertex
    start_vertex = None
    for u in g:
        if len(g[u]) - in_degrees.get(u, 0) == 1:
            return u
    for u in g:
        if len(g[u]) > in_degrees.get(u, 0):
            return u
    return next(iter(g))



def eulerian_path(g: Dict[int, List[int]]) -> Iterable[int]:
    ensure_all_vertices_present(g)
    start_vertex = find_start_vertex(g)
    all_paths = list()
    path = deque([start_vertex])
    cycle = []

    while path:
        vertex = path[-1]
        if g[vertex]:
            next_vertex = g[vertex].pop(0)
            path.append(next_vertex)
        else:
            cycle.append(path.pop())
        
    return cycle[::-1]

def genome_path(path: List[str]) -> str:
    """Forms the genome path formed by a collection of patterns."""
    return ''.join([e[0] for e in path])+path[-1][1:]

# Insert your string_reconstruction function here, along with any subroutines you need
def find_read_path(patterns: List[str], read_kmer_dict:dict):
    
    db = de_bruijn_kmers(patterns)
    path = eulerian_path(db)
    text = genome_path(path)
    read_order = recover_read_numbers(text,read_kmer_dict)
    return read_order

def recover_read_numbers(text:str, spectrum_dict:dict):
    read_order_list = list()
    for i in range(len(text)):
        for key, value in spectrum_dict.items():
            if text[i:i+20] == value:
                format = f'>{key}'
                read_order_list.append(format)
    return read_order_list

#################
#################



# print("reads from solutions")
# print(sample_1_spectrum['read_210'])
# print(sample_1_spectrum['read_154'])
# print(sample_1_spectrum['read_181'])




path = find_read_path(read_list, proj_2a_spec)

# Writing to the CSV file
with open('last_proj_test.csv', 'w', newline='') as csv_file:
    # Step 4: Using csv.writer to write the list to the CSV file
    writer = csv.writer(csv_file)
    writer.writerow(path) # Use writerow for single list


with open("last_proj_test.csv", "w", newline='') as csv_file:
    csv_writer = csv.writer(csv_file)

    for snp in path:
     csv_writer.writerow([snp])
