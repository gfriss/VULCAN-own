''' Implementation of the Dijkstra algorithm. Code skeleton taken from www.datacamp.com which
    I extended with the making of the graph (in __init__ funciton) and tracking of reactions.'''
#%%
import numpy as np
from heapq import heapify, heappop, heappush

class Graph:
    def __init__(self, dat, n_layer):
        # prepare variable used to construct matrix
        self.spec_list = dat['variable']['species']
        self.rea_dict = self.create_reaction_dict(dat)
            
        # create the dictionary to hold the graph in which keys are the nodes (species)
        # and items are edges (reactions) as dictionaries, such that if A is connected to B and C with
        # weigths of 1 and 2 through reactions b and c and total reaction rate for all reaction contributing to this edge b_t and c_t, 
        # respectively then it is 'A': {'B':[b,b_t,1], 'C':[c,c_t,2]}
        graph = {key: {} for key in self.spec_list}
        # given the  multiplicity of the network, we need to track which reactions give the weights
        # and what percantage that reaction contribute to that edge by choosing the most prominent reaction (hence the need for the total reaction rate)
        
        # iterate over reactions and fill dictionary value
        for i,rea in self.rea_dict.items():
            # first calculate reaction rate
            rate = dat['variable']['k'][i][n_layer]
            for reactant in rea[0]: # first get the reaction rate
                if reactant != 'M':
                    rate *= dat['variable']['y'][n_layer, self.spec_list.index(reactant)]
            if rate == 0: # in case there is a zero value change it so no calculation errors when later taking 1/rate
                rate = 1e-99
            # then fill up the dictionary
            for reactant in rea[0]:
                if reactant != 'M':
                    for product in rea[1]:
                        if product != 'M' and product in graph[reactant]: # if there is already a value for this edge, could override
                            chosen_reaction = np.argmax([graph[reactant][product][1], rate]) # decide whether rate of new reaction is greater than previously found
                            graph[reactant][product][0] = [graph[reactant][product][0], i][chosen_reaction] # use chosen reaction to be saved
                            graph[reactant][product][1] += rate # adding to the total reaction rate, how to avoid duplication???
                            graph[reactant][product][2] = [graph[reactant][product][2], 1/rate][chosen_reaction]
                        elif product != 'M' and product not in graph[reactant]: # if first entry for this edge, simply assign:
                            graph[reactant][product] = [0, 0, 0]
                            graph[reactant][product][0] = i # reaction ID
                            graph[reactant][product][1] = rate # total reaction rate
                            graph[reactant][product][2] = 1/rate # weight
        
        self.graph = graph
        
    def get_species(self, eq_side):
        ''' Returns the species in a given reaction in an array.'''
        side_split = eq_side.split('+')
        if len(side_split) == 1: # stripping them from white spaces
            side_split = np.array([side_split[0].strip()]) # with array length 1 the other method fails so doing it separately
        else:
            side_split = np.array([r.strip() for r in side_split])
        return side_split
    
    def create_reaction_dict(self, dat):
        re_dict_sim = {}
        for re in dat['variable']['k'].keys():
            if re % 2 == 0:
                reagents_products = dat['variable']['Rf'][re-1].split('->') # separating production and destruction
                reagents = list(self.get_species(reagents_products[1])) # reversed because reaction is reversed
                products = list(self.get_species(reagents_products[0]))
            else:
                reagents_products = dat['variable']['Rf'][re].split('->') # separating production and destruction
                reagents = list(self.get_species(reagents_products[0]))
                products = list(self.get_species(reagents_products[1]))
            re_dict_sim[re] = [reagents, products]
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
            current_distance, current_node = heappop(
                pq
            )  # Get the node with the min distance

            if current_node in visited:
                continue  # Skip already visited nodes
            visited.add(current_node)  # Else, add the node to visited set
            
            for neighbour, weight in self.graph[current_node].items():
                # Calculate the distance from current_node to the neighbor
                tentative_distance = current_distance + weight[2]
                if tentative_distance < distances[neighbour]:
                    distances[neighbour] = tentative_distance
                    heappush(pq, (tentative_distance, neighbour))

        predecessors = {node: None for node in self.graph}

        for node, distance in distances.items():
            for neighbour, weight in self.graph[node].items():
                if distances[neighbour] == distance + weight[2]:
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
            # reaction is zeroth element in the graph, contribution is rate/total_rate but I have 1/rate (2nd item) and total rate (first item)
            # so it becomes 1/(total_rate * 1/rate)
            reactions[self.graph[node][neighbour][0]] = [self.rea_dict[self.graph[node][neighbour][0]], 1/(self.graph[node][neighbour][1] * self.graph[node][neighbour][2]), 1/self.graph[node][neighbour][2]]
        
        return path, reactions
# %%
import pickle
vul_data = 'output/a10_conv.vul'
with open(vul_data, 'rb') as handle:
    data = pickle.load(handle)
# %%
G = Graph(data, 0)
#%%
G.shortest_path('CH4', 'HCN')
# %%
