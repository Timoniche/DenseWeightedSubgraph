# #####################################################
#   Authors  #  Varun Gohil (gohil.varun@iitgn.ac.in) #
#######################################################

import maxflow


def WFind_Density(answer, filepath):
    ''' Finds the density of the returned subgraph.'''
    degree = 0
    file = open(filepath, "r")
    while True:
        edge = file.readline()
        # print "Reading an edge"
        if not edge:
            break
        from_node, to_node, w = edge.split()
        if (int(from_node) in answer and int(to_node) in answer):
            degree += (2 * float(w))
    file.close()
    return degree / (2.0 * (len(answer)))



def WFind_Densest_Subgraph(number_of_nodes, number_of_edges, filepath):
    ''' This function performs the binary search of the density of subgraph and finds the densest subgraph.'''
    min_degree = 0

    file = open(filepath, "r")
    sum_w = 0.0
    while True:
        edge = file.readline()
        if not edge:
            break
        from_node, to_node, w = edge.split()
        w = float(w)
        sum_w += w
    file.close()

    max_degree = sum_w

    subgraph = []
    difference = 1.0 / (number_of_nodes * (number_of_nodes + 1))
    # print difference, " diff"
    while (max_degree - min_degree >= difference):
        print("...")
        # print "max - min = ", max_degree - min_degree
        least_density = (max_degree + min_degree) / 2.0
        # print "ld--->", least_density
        source_segment = Wmake_graph(number_of_nodes, number_of_edges, least_density, filepath)
        if (source_segment == []):
            max_degree = least_density
        else:
            min_degree = least_density
            subgraph = source_segment
        # print subgraph
    return subgraph


def Wmake_graph(number_of_nodes, number_of_edges, least_density, filepath):
    ''' Constructs the network as per the specifications given by Goldberg'''
    graph = maxflow.Graph[float](number_of_nodes, number_of_edges)
    nodes = graph.add_nodes(number_of_nodes)
    # print nodes
    degrees = {}
    # print degrees
    file = open(filepath, "r")
    sum_w = 0.0
    while True:
        edge = file.readline()
        if not edge:
            break
        from_node, to_node, w = edge.split()
        w = float(w)
        sum_w += w
        # print edge.split()
        graph.add_edge(nodes[int(from_node)], nodes[int(to_node)], w, w)
        if from_node in degrees:
            degrees[from_node] += w
        else:
            degrees[from_node] = w
        if to_node in degrees:
            degrees[to_node] += w
        else:
            degrees[to_node] = w
    file.close()
    for i in range(number_of_nodes):
        if str(i) not in degrees:
            degrees[str(i)] = 0
        graph.add_tedge(nodes[i], sum_w, sum_w + 2 * least_density - degrees[str(i)])
        # print "s -- ",number_of_edges,"-->", nodes[i], "--",number_of_edges + 2*least_density - degrees[str(i)], "-->t\n"
    source_segment = []
    '''Computes the max-flow in the graph'''
    print(f'max flow is {graph.maxflow()}')
    print(f'least_density is {least_density}')
    '''The following section of code finds which node belongs to which cutset.'''
    for i in nodes:
        # print nodes[i] ,"--->", graph.get_segment(nodes[i])
        if (graph.get_segment(nodes[i]) == 0):
            source_segment.append(nodes[i])
    # print degrees
    print(source_segment)
    return source_segment


def main():
    number_of_nodes = int(input("Enter number of nodes: "))
    number_of_edges = int(input("Enter number of edges: "))

    answer = WFind_Densest_Subgraph(number_of_nodes, number_of_edges, "edges_w.txt")
    print(answer)
    print(WFind_Density(answer, "edges_w.txt"))


if __name__ == '__main__':
    main()
