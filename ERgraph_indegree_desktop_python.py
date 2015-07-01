import matplotlib.pyplot as plt
import random
###############################################
"""
Algorithm ER
"""
###############################################
def compute_in_degrees(digraph):
    """
     return a dictionary with the same set of keys (nodes) 
     as digraph whose corresponding values are the number 
     of edges whose head matches a particular node
    """
    degree = dict()
    for values in digraph.values():
        for node in values:
            if degree.has_key(node):
                degree[node] += 1
            else:
                degree[node] = 1
    for node in digraph.keys():
        if not degree.has_key(node):
            degree[node] = 0
    return degree

def in_degree_distribution(digraph):
    """
    computes the unnormalized distribution of the in-degrees of the graph
    """
    in_degree_dis = dict()
    in_gress = compute_in_degrees(digraph)
    for value in in_gress.values():
        if in_degree_dis.has_key(value):
            in_degree_dis[value] += 1
        else:    
            in_degree_dis[value] = 1
    return in_degree_dis

def algorithm_ER(nnodes, p):
    """
    generate ER graph
    """
    graph = dict()
    for index_1 in range(nnodes):
        if graph.has_key(index_1):
            neighbour = graph[index_1]
        else:
            neighbour = set()
        for index_2 in range(nnodes):
            if index_1 != index_2:             
                a = random.random()
                if a < p:
                    if not graph.has_key(index_2):
                        graph[index_2] = set([index_1])
                    else:
                        graph[index_2].add(index_1)
                    neighbour.add(index_2)           
        graph[index_1] = neighbour 
    return graph

def in_degree_distribution_ER(in_degree_dis, noden):
    """
    computes the normalized distribution of the in-degrees of the ER graph
    """
    in_degree_dis_normal = dict()
    total = 0
    for value in in_degree_dis.values():
        total += value
    for node in in_degree_dis.keys():
        in_degree_dis_normal[node] = float(in_degree_dis[node])/total
    for index in range(noden):
        if index not in in_degree_dis_normal:
            in_degree_dis_normal[index] = 0               
    return in_degree_dis_normal

def generate_ingree_ER(node,p):
    """
    generate ER graph based on the number of nodes and p
    Then get the in-degree distribution
    """
    graph = algorithm_ER(node,p)
    indegree_R = in_degree_distribution(graph)
    indegree_dis_ER = in_degree_distribution_ER(indegree_R, node)
    return indegree_dis_ER

graph_in = generate_ingree_ER(10000, 0.5)
plt.xlabel('In-degree number')
plt.ylabel('Normalized quantity')
plt.title('ER graph in-degree distribution with nodeno = 10000 p = 0.5')
plt.plot(graph_in.keys(), graph_in.values(),'bo')
plt.show()
