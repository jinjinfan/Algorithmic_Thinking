"""
Degree distributions for graphs
"""
from alg_load_graph import *
import simpleplot
import math
import random
# Set timeout for CodeSkulptor if necessary
import codeskulptor
codeskulptor.set_timeout(100)
EX_GRAPH0 = {0:set([1,2]),1:set([]),2:set([])}
EX_GRAPH1 = {0:set([1,4,5]),1:set([2,6]),2:set([3]),3:set([0]),4:set([1]),5:set([2]),6:set([])}
EX_GRAPH2 = {0:set([1,4,5]),1:set([2,6]),2:set([3,7]),3:set([7]),4:set([1]),5:set([2]),6:set([]),\
             7:set([3]),8:set([1,2]),9:set([0,3,4,5,6,7])}

def make_complete_graph(num_nodes):
    """
    returns a dictionary corresponding to a complete directed graph
    """
    digraph = dict()
    if num_nodes > 0:
        for index1 in range(num_nodes):
            neighbour = []
            for index2 in range(num_nodes):
                if index2 != index1:
                    neighbour.append(index2)
            digraph[index1] = set(neighbour)
    return digraph
def compute_out_degrees(digraph):
    """
     return a dictionary with the same set of keys (nodes)
     as digraph whose corresponding values are out degree
     of each node
    """
    out_degree = dict()
    for key in digraph.keys():
        out_degree[key] = len(digraph[key])
    return out_degree

def compute_out_degrees_average(digraph):
    """
     return the average out-degree of the graph
    """
    out_degree = compute_out_degrees(digraph)
    totalnode = len(out_degree.keys())
    print out_degree.values()
    total = 0
    for value in out_degree.values():
        total += value
    return int(math.ceil(float(total)/totalnode))


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

def in_degree_distribution_normal(in_degree_dis):
    """
    computes the normalized distribution of the in-degrees of the graph
    """
    in_degree_dis_normal = dict()
    total = 0
    for value in in_degree_dis.values():
        total += value
    for node in in_degree_dis.keys():
        in_degree_dis_normal[node] = float(in_degree_dis[node])/total
    return in_degree_dis_normal

def build_plot(plot_size, digraph):
    """
    Build plot for citation graph
    """
    print digraph
    total = 0
    plot = []
    for key in digraph.keys():
#       using linear scale
       plot.append([key, digraph[key]])
#       using log/log scale
#        plot.append([math.log(key), math.log(digraph[key])])
    return plot

# plottting code
plot_size = 40

#citation_graph = load_graph(CITATION_URL)
#indegre_citation= in_degree_distribution(citation_graph)
#indegree = in_degree_distribution_normal(indegre_citation)
##
#plot1 = build_plot(plot_size, indegree)
#simpleplot.plot_scatter ("citation graph in-degree distribution", 600, 600,
#                      "In-degree number", "Normalized quantity", [plot1])

# calculate average out degree for the use in DPA algorithm
#compute_out_degrees = compute_out_degrees(citation_graph)
#print compute_out_degrees_average(citation_graph)
###############################################
"""
Algorithm ER
"""
###############################################
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

plot1 = build_plot(plot_size, generate_ingree_ER(2000,0.5))
simpleplot.plot_lines("ER graph in-degree distribution", 600, 600,
                     "In-degree number", "Normalized quantity", [plot1], True)

###############################
#Algorithm DPA
##############################
import alg_dpa_trial

def DPA_algorithm(n, m):
    nodes = range(m)
    dpa_trial_instance=alg_dpa_trial.DPATrial(m)
    subset_graph = make_complete_graph(m)
    for index in range(m, n):
        nodes_add = set()
        nodes_add.update(dpa_trial_instance.run_trial(m))
        nodes.append(index)
        subset_graph[index] = nodes_add
    return subset_graph

def DPA_algorithm_indegree(n,m):
    subset_graph = DPA_algorithm(n,m)
    in_degree = in_degree_distribution(subset_graph)
    in_degree_normal = in_degree_distribution_normal(in_degree)
    print in_degree_normal
    return in_degree_normal

#plot1 = build_plot(plot_size, DPA_algorithm_indegree(27770,13))
#simpleplot.plot_scatter("The normalized in-degree distribution for the DPA graph", 600, 600,
#                      "log10 of in-degrees", "log10 of fraction of nodes with in-degrees", [plot1])

#