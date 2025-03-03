''' Implementation of the Dijkstra algorithm. Code skeleton taken from www.datacamp.com which
    I extended with the making of the graph (make_graph() funciton) and tracking of reactions.'''
#%%
import numpy as np
from heapq import heapify, heappop, heappush
import sys
import pickle
import copy
import os
from contextlib import redirect_stdout

output_folder = '/scratch/s2555875/output'
pathway_folder = '/scratch/s2555875/pathways'
conv_file = '/scratch/s2555875/converged.txt'
with open(conv_file, 'r') as f:
    conv_text = f.read()
#%%
class Graph:
    def __init__(self, graph: dict = {}, g_reactions: dict = {}, dat: dict = {}):
        self.graph = graph
        self.g_reactions = g_reactions
        if dat:
            self.dat = dat
            self.spec_list = self.dat['variable']['species']
            self.rea_dict = self.create_reaction_dict()
            self.graph = {key: {} for key in self.spec_list} # dictionary to hold the graph in which keys are the nodes (species) and items are edges (lifetime, i.e. n/reactionrate)
            self.g_reactions = {key: {} for key in self.spec_list} # dictionary to hold information on reaction ID, the reaction rate of the chosen path and the total reaction rate of all the paths between two nodes; and the lifetime
        
    
    def make_graph(self, n_layer):
        # given the  multiplicity of the network, we need to track which reactions give the weights
        # and what percantage that reaction contribute to that edge by choosing the most prominent reaction (hence the need for the total reaction rate)
        
        # iterate over reactions and fill dictionary value
        for re_id,rea in self.rea_dict.items():
            reactant_list = rea[0]
            product_list = rea[1]
            # excluding loops
            if set(reactant_list) == set(product_list):
                continue
            # first calculate reaction rate
            rate = self.dat['variable']['k'][re_id][n_layer]
            if rate > 0: # save resource by calculating rate only if rate coefficient is non-zero (i.e. the reaction takes place)
                for reactant in reactant_list:
                    if reactant == 'M': # when third body is present
                        rate *= self.dat['atm']['n_0'][n_layer]
                    else: # mono- and bimolecular reactions
                        rate *= self.dat['variable']['y'][n_layer, self.spec_list.index(reactant)]
            # then fill up the dictionary
            duplicate_reactant = False
            for reactant in reactant_list:
                # avoid repetition and third body
                if duplicate_reactant:
                    duplicate_reactant = False # in case there would be a double team up of species
                    continue
                if reactant == 'M':
                    continue
                # if reactant is in reaction twice, need to double the rate so the weight would be quicker (shorter lifetime)
                rate_sp = copy.copy(rate)
                if reactant_list.count(reactant) == 2:
                    rate_sp *= 2
                    duplicate_reactant = True
                # setting weight, default is inf when unconnected, i.e. rate = 0, or it's a loop
                # otherwise weight is lifetime tau_i = n/(k*n_i*nr) where nr is the number density oh the rest of reactants and rate = k*n*nr
                lifetime_reactant = np.inf
                if rate_sp > 0:
                    lifetime_reactant = self.dat['variable']['y'][n_layer, self.spec_list.index(reactant)] / rate_sp
                duplicate_product = False
                for product in product_list:
                    # avoiding third body and repetation
                    if product == 'M':
                        continue
                    if duplicate_product:
                        duplicate_product = False
                        continue
                    lifetime_product = copy.copy(lifetime_reactant)
                    rate_sp_product = copy.copy(rate_sp)
                    if product_list.count(product) == 2: # if twice as much is produced, this route is favoured, takes half the time for the same amount
                        duplicate_product = True
                        lifetime_product /= 2
                        rate_sp_product *= 2
                    if product not in self.graph[reactant]: # if first entry for this edge, simply assign values
                        self.graph[reactant][product] = copy.copy(lifetime_product)
                        self.g_reactions[reactant][product] = [copy.copy(re_id), copy.copy(rate_sp_product), copy.copy(rate_sp_product), copy.copy(lifetime_product)] # reaction ID, rate, total reaction rate and lifetime
                    else: # if there is already a value for this edge, could override if new reaction is faster (shorter lifetime)
                        self.g_reactions[reactant][product][2] += copy.copy(rate_sp_product) # adding to the total reaction rate
                        if lifetime_product < self.graph[reactant][product]: # override if faster
                            self.graph[reactant][product] = copy.copy(lifetime_product)
                            self.g_reactions[reactant][product][0] = copy.copy(re_id)
                            self.g_reactions[reactant][product][1] = copy.copy(rate_sp_product) 
                            self.g_reactions[reactant][product][3] = copy.copy(lifetime_product) 
        
    def get_species(self, eq_side):
        ''' Returns the species in a given reaction side in a list.'''
        side_split = eq_side.split('+')
        side_split = [r.strip() for r in side_split] # stripping them from white spaces
        return side_split
    
    def create_reaction_dict(self):
        re_dict_sim = {}
        for re_id,rea in self.dat['variable']['Rf'].items():
            reagents_products = rea.split('->') # separating production and destruction
            # assigning reactant and product for re_id and reverse for the reverse reaction (re_id+1)
            re_dict_sim[re_id] = [self.get_species(reagents_products[0]), self.get_species(reagents_products[1])]
            re_dict_sim[re_id+1] = [self.get_species(reagents_products[1]), self.get_species(reagents_products[0])]
        return re_dict_sim
        
    def shortest_distances(self, source: str):
        # Initialize the values of all nodes with infinity
        distances = {node: float("inf") for node in self.graph}
        distances[source] = 0  # Set the source value to 0

        # Initialize a priority queue
        pq = [(0, source)]
        heapify(pq)

        # Create a set to hold visited nodes
        visited = set()
        while pq:  # While the priority queue isn't empty
            current_distance, current_node = heappop(pq)  # Get the node with the min distance
            
            if current_node in visited:
                continue  # Skip already visited nodes
            visited.add(current_node)  # Else, add the node to visited set
            
            for neighbour, weight in self.graph[current_node].items():
                # Calculate the distance from current_node to the neighbor
                tentative_distance = current_distance + weight
                if tentative_distance < distances[neighbour]:
                    distances[neighbour] = tentative_distance
                    heappush(pq, (tentative_distance, neighbour))

        predecessors = {node: None for node in self.graph}

        for node, distance in distances.items():
            for neighbour, weight in self.graph[node].items():
                if distances[neighbour] == distance + weight:
                    predecessors[neighbour] = node
                    
        return distances, predecessors

    def shortest_path(self, source: str, target: str):
        # Generate the predecessors dict
        _, predecessors = self.shortest_distances(source)

        path = []
        reactions = {}
        current_node = target

        # Backtrack from the target node using predecessors
        while current_node:
            path.append(current_node)
            current_node = predecessors[current_node]

        # Reverse the path and reactions and their contributions and rates (as dict -> id: [reaction, contribution, rate])
        path.reverse()
        for i,node in enumerate(path):
            if i == len(path) - 1:
                break
            neighbour = path[i+1]
            # reaction is zeroth element in the graph, 1st element is the rate of that reaction and 2nd is the total reaction rate between the two species
            reactions[self.g_reactions[node][neighbour][0]] = [self.rea_dict[self.g_reactions[node][neighbour][0]], '{:.3f} %'.format(100*self.g_reactions[node][neighbour][1]/self.g_reactions[node][neighbour][2]), 'k = {:.3e} cm-3s-1'.format(self.g_reactions[node][neighbour][1]), 'tau = {:.4e} s'.format(self.g_reactions[node][neighbour][3])]
        
        return path, reactions
    
    def get_nodes(self):
        "Returns the nodes of the graph."
        return self.graph.keys()
    
    def get_outgoing_edges(self, node):
        "Returns the neighbors of a node."
        connections = []
        for out_node in self.get_nodes():
            if self.graph[node].get(out_node, False) != False:
                connections.append(out_node)
        return connections
    
    def value(self, node1, node2):
        "Returns the value of an edge between two nodes."
        return self.graph[node1][node2]

    def dijkstra_algorithm(self, start_node):
        unvisited_nodes = list(self.get_nodes())
    
        # We'll use this dict to save the cost of visiting each node and update it as we move along the graph   
        shortest_path = {}
    
        # We'll use this dict to save the shortest known path to a node found so far
        previous_nodes = {}
    
        # We'll use max_value to initialize the "infinity" value of the unvisited nodes   
        max_value = np.inf
        for node in unvisited_nodes:
            shortest_path[node] = max_value
        # However, we initialize the starting node's value with 0
        shortest_path[start_node] = 0
        
        # The algorithm executes until we visit all nodes
        while unvisited_nodes:
            # The code block below finds the node with the lowest score
            current_min_node = None
            for node in unvisited_nodes: # Iterate over the nodes
                if current_min_node == None:
                    current_min_node = node
                elif shortest_path[node] < shortest_path[current_min_node]:
                    current_min_node = node
                    
            # The code block below retrieves the current node's neighbors and updates their distances
            neighbors = self.get_outgoing_edges(current_min_node)
            for neighbor in neighbors:
                tentative_value = shortest_path[current_min_node] + self.value(current_min_node, neighbor)
                if tentative_value < shortest_path[neighbor]:
                    shortest_path[neighbor] = tentative_value
                    # We also update the best path to the current node
                    previous_nodes[neighbor] = current_min_node
    
            # After visiting its neighbors, we mark the node as "visited"
            unvisited_nodes.remove(current_min_node)
        
        self.previous_nodes = previous_nodes
        self.path = shortest_path
        
    def print_result(self, start_node, target_node):
        path = []
        reactions = {}
        node = target_node
        
        while node != start_node:
            path.append(node)
            node = self.previous_nodes[node]
    
        # Add the start node manually
        path.append(start_node)
        path = path[::-1]
        print("We found the following best path with a value of {:.4e} s.".format(self.path[target_node]))
        print(" -> ".join(path))
        for i,node in enumerate(path):
            if i == len(path) - 1:
                break
            neighbour = path[i+1]
            print('{} -> {}'.format(node, neighbour))
            print('ID: {}, reaction: {}, rate = {:.4e} cm-3s-1, percentage: {:.3f}%, tau = {} s'.format(self.g_reactions[node][neighbour][0], self.rea_dict[self.g_reactions[node][neighbour][0]], self.g_reactions[node][neighbour][1], 100*self.g_reactions[node][neighbour][1]/self.g_reactions[node][neighbour][2], self.g_reactions[node][neighbour][3]))
# %%
def get_layer(dat, pressure):
    ''' Returns the layer that is the closest to the given pressure (in bar) in the given VULCAN data.'''
    return np.argmin(abs(dat['atm']['pco']/1e6 - pressure))

def patway_analysis(pressures, start_sp, end_sp, sim_type, network, nsim):
    ''' Performs the pathway analysis for a given simulation type and pressures.
        The results are saved in a txt the pathways folder.'''
    for p in pressures:
        pathways_file = os.path.join(pathway_folder, '{}->{}_p{:.1e}_{}{}.txt'.format(start_sp, end_sp, p, sim_type, network))
        f = open(pathways_file, 'w')
        f.write('# This file contains the pathway analysis at p = {:.1e} bar, starting from {} and ending in {} throughout the {} simulations\n'.format(p, start_sp, end_sp, sim_type))
        for i in range(nsim):
            sim = 'sim_'
            if i < 10:
                sim += '0{}_{}{}.vul'.format(i, sim_type, network)
            else:
                sim += '{}_{}{}.vul'.format(i, sim_type, network)
            convergence = '\t(NOT converged)'
            if sim in conv_text:
                convergence = '\t(converged)'
            # if rerun was needed, combiune the results
            sim_rerun_file = sim[:-4] + '_rerun.vul'
            if os.path.exists(os.path.join(output_folder, sim_rerun_file)):
                with open(os.path.join(output_folder, sim_rerun_file), 'rb') as handle:
                    data = pickle.load(handle)
                if sim_rerun_file in conv_text:
                    convergence = '\t(converged after rerun)'
            else:
                with open(os.path.join(output_folder, sim), 'rb') as handle:
                    data = pickle.load(handle)

            n_layer = get_layer(data, p)
            G = Graph(dat = data)
            G.make_graph(n_layer)
            G.dijkstra_algorithm(start_node = start_sp)
            f.write('\n{}{}\n'.format(sim, convergence))
            with redirect_stdout(f):
                G.print_result(start_node = start_sp, target_node = end_sp)
        f.close()
#%%
network = ''
#network = '_ncho'
pressures = [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
start_sp = 'CH4'
end_sp = 'HCN'
patway_analysis(pressures = pressures, start_sp = start_sp, end_sp = end_sp, sim_type = 'meteor', network = network, nsim = 15)
patway_analysis(pressures = pressures, start_sp = start_sp, end_sp = end_sp, sim_type = 'CtoO', network = network, nsim = 15)
patway_analysis(pressures = pressures, start_sp = start_sp, end_sp = end_sp, sim_type = 'dist', network = network, nsim = 15)
patway_analysis(pressures = pressures, start_sp = start_sp, end_sp = end_sp, sim_type = 'star', network = network, nsim = 13)
# %%
