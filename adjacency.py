import random
import numpy as np
import itertools


def debruijnize(reads):
    # Initialize sets to store nodes, nodes that are not starts, and list to store edges
    nodes = set()
    not_starts = set()
    edges = []
    
    # Iterate over each read in the input list of reads
    for r in reads:
        # Extract the substrings by removing the last and first characters
        r1 = r[:-1]
        r2 = r[1:]
        
        # Add the substrings to the set of nodes
        nodes.add(r1)
        nodes.add(r2)
        
        # Create an edge between the substrings and append it to the list of edges
        edges.append((r1, r2))
        
        # Add the second substring to the set of nodes that are not starting nodes
        not_starts.add(r2)
    
    # Return a tuple containing the set of nodes, list of edges, and list for start node if exists
    return (nodes, edges, list(nodes - not_starts))

def build_k_mer(str, k):
    # Initialize an empty list to store the k-mers
    k_mers = []
    
    # Iterate over the string with a sliding window of size k
    for i in range(0, len(str) - k + 1):
        # Extract a k-mer substring starting from index i
        k_mer = str[i:k+i]
        
        # Append the extracted k-mer to the list
        k_mers.append(k_mer)
    
    # Return the list of k-mers
    return k_mers

def make_node_edge_map(edges):
    # Initialize an empty dictionary to store node-edge mappings
    node_edge_map = {}
    
    # Iterate over each edge in the list of edges
    for e in edges:
        # Extract the source node from the edge
        n = e[0]
        
        # Check if the source node already exists in the mapping
        if n in node_edge_map:
            # If it exists, append the destination node to its list of connected nodes
            node_edge_map[n].append(e[1])
        else:
            # If it doesn't exist, create a new entry with the source node as key
            # and a list containing the destination node as its value
            node_edge_map[n] = [e[1]]
    
    # Return the node-edge mapping dictionary
    return node_edge_map

def eulerian_trail(m, v):
    # Copy the input node-edge mapping to a local variable
    nemap = m
    
    # Initialize an empty list to store the result trail
    result_trail = []
    
    # Set the starting node for the trail
    start = v
    
    # Append the starting node to the result trail
    result_trail.append(start)
    
    # Loop until the trail is completed
    while(True):
        # Initialize an empty list to store the current trail
        trail = []
        
        # Initialize the previous node to the starting node
        previous = start
        
        # Traverse the current trail until it loops back to the starting node
        while(True):
            # Check if the previous node has any outgoing edges
            if(previous not in nemap):
                break
            
            # Pop the next node from the list of outgoing edges of the previous node
            next = nemap[previous].pop()
            
            # Remove the previous node from the mapping if it has no more outgoing edges
            if(len(nemap[previous]) == 0):
                nemap.pop(previous, None)
            
            # Append the next node to the current trail
            trail.append(next)
            
            # Break the loop if the next node is the starting node, indicating completion of the trail
            if(next == start):
                break
            
            # Update the previous node to the next node for the next iteration
            previous = next
        
        # Combine the current trail with the result trail
        index = result_trail.index(start)
        result_trail = result_trail[0:index + 1] + trail + result_trail[index + 1:len(result_trail)]
        
        # Choose a new start node if there are remaining edges
        if(len(nemap) == 0):
            break
        
        # Find a new start node from the result trail
        found_new_start = False
        for n in result_trail:
            if n in nemap:
                start = n
                found_new_start = True
                break
        
        # Break the loop if no new start node is found
        if not found_new_start:
            break
    
    # Return the result trail
    return result_trail


def assemble_trail(trail):
    # Check if the trail is empty
    if len(trail) == 0:
        return ""
    
    # Initialize the result string with the first node's substring excluding its last character
    result = trail[0][:-1]
    
    # Concatenate the last character of each node in the trail to the result string
    for node in trail:
        result += node[-1]
    
    # Return the assembled result string
    return result

def test_assembly_debruijn(str, k):
    # Generate all k-mers from the input string
    reads = build_k_mer(str, k)
    
    # Construct the De Bruijn graph from the k-mers
    G = debruijnize(reads)
    
    # Create a node-edge mapping from the edges of the De Bruijn graph
    nemap = make_node_edge_map(G[1])
    
    # Choose a starting node for the Eulerian trail
    start = next(iter(G[2])) if (len(G[2]) > 0) else next(iter(G[0]))
    
    # Find the Eulerian trail in the De Bruijn graph
    trail = eulerian_trail(nemap, start)
    
    # Assemble the Eulerian trail into a string
    return assemble_trail(trail)

def generate_random_genome(length):
    nucleotides = ['A', 'C', 'G', 'T']
    genome = ''.join(random.choices(nucleotides, k=length))
    return genome

class Adjacency:
    def __init__(self, kmer_length):
        # Initialize the length of k-mers and generate all possible k-mers of that length
        self.kmer_length = kmer_length
        self.all_kmers = self.get_all_kmers()
        self.reads_length = kmer_length + 1
        # Create a mapping from k-mers to their indices in the list of all k-mers
        self.kmer_to_index = {kmer: i for i, kmer in enumerate(self.all_kmers)}
    
    def graph2adjacency(self, k_mers):
        # Determine the size of the adjacency matrix based on the number of all possible k-mers
        matrix_size = len(self.all_kmers)
        
        # Initialize the adjacency matrix with zeros
        adjacency_matrix = np.zeros((matrix_size, matrix_size), dtype=int)
        
        # Populate the adjacency matrix based on the input k-mers
        for kmer, neighbors in k_mers.items():
            row_index = self.kmer_to_index[kmer]
            for neighbor in neighbors:
                col_index = self.kmer_to_index[neighbor]
                adjacency_matrix[row_index, col_index] += 1
        
        # Return the constructed adjacency matrix
        return adjacency_matrix
    
    def genome2adjacency(self, genome):
        reads = build_k_mer(genome, self.reads_length)
        nodes, edges, start_node = debruijnize(reads)
        node_edge_map = make_node_edge_map(edges)
        adjacency_matrix = self.graph2adjacency(node_edge_map)
        return (adjacency_matrix, start_node)
    
    def adjacency2genome(self, adjacency_matrix, start_node):
        node_edge_map_from_adjacency = self.adjacency2graph(adjacency_matrix)

        nodes = set()
        for value in node_edge_map_from_adjacency.values():
            nodes.update(value)

        start = start_node[0] if (len(start_node) > 0) else next(iter(nodes))
        trail = eulerian_trail(node_edge_map_from_adjacency, start)
        rec_genome = assemble_trail(trail)
        return rec_genome

    def get_all_kmers(self):
        # Generate all possible k-mers of the specified length
        nucleotides = ['A', 'C', 'G', 'T']
        kmers = [''.join(p) for p in itertools.product(nucleotides, repeat=self.kmer_length)]
        return kmers
    
    def adjacency2graph(self, adjacency_matrix):
        # Initialize an empty dictionary to represent the graph
        graph = {}
        
        # Iterate over each k-mer and its index in the kmer_to_index mapping
        for kmer, index in self.kmer_to_index.items():
            neighbors = []
            
            # Iterate over each entry in the row of the adjacency matrix corresponding to the current k-mer
            for i, val in enumerate(adjacency_matrix[index]):
                # If the value is greater than 0, add the corresponding k-mer to the neighbors list
                if val > 0:
                    neighbors.extend([self.all_kmers[i]] * val)
            
            # If the neighbors list is not empty, add the k-mer and its neighbors to the graph
            if neighbors:
                graph[kmer] = neighbors
        
        # Return the constructed graph
        return graph