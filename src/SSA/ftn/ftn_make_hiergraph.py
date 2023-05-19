import numpy as np
import networkx as nx
def make_hiermat(limitset_list):
    """return a matrix representing the hierarchy of limitset_list
    Each column and row represents a limitset.
    
    Parameters
    ----------
    limitset_list : list
        list of limitsets(buffering structures)
    
    Returns
    -------
    hiermat : numpy.ndarray
        matrix representing the hierarchy of limitset_list
    """
    # return a matrix representing the hierarchy of limitset_list

    limitset_list_name = []
    for lset in limitset_list:
        lset_name = lset[0]+lset[1]  # reaction_idを追加
        limitset_list_name.append(lset_name)

    l = len(limitset_list_name)

    relmat = np.zeros((l, l), dtype=int)  # matrix representing inclusion relation of limitsets
    # relmat[i,j]=1 if limitset[j] is in limitset[i].
    for i, s in enumerate(limitset_list_name):
        for j, t in enumerate(limitset_list_name):
            if set(t) < set(s):
                relmat[i, j] = 1

    relmat2 = []  # inclusion relation passing more than two edges in hierarchy graph
    for i in range(l):
        row = np.dot(relmat[i], relmat)
        row2 = np.array([1 if p > 0 else 0 for p in row])
        relmat2.append(row2)
    relmat2 = np.array(relmat2)

    hiermat = relmat-relmat2
    if np.min(hiermat) < 0:
        raise Exception('inclusive relation is not correct')
        
    return hiermat

def make_hieredge(limitset_list):
    """return a hierarchy graph as a list of edges
    
    Parameters
    ----------
    limitset_list : list
        list of limitsets(buffering structures)
    
    Returns
    -------
    node_list : list
        list of nodes of hierarchy graph
    edge_list : list
        list of edges of hierarchy graph
    """
    # return a hierarchy graph as a list of edges

    limitset_list_all = []  # reactionとcompoundを一緒のリストにする
    for lset in limitset_list:
        lset_all = lset[0]+lset[1]  # reaction_idを追加
        limitset_list_all.append(lset_all)

    l = len(limitset_list_all)
    hiermat = make_hiermat(limitset_list)

    # make nodes of hierarchy graph
    node_list = []
    for i in range(l):
        elim = set()
        # if j is in i, eliminate j from i
        for j in range(l):
            if hiermat[i, j] == 1:
                elim |= set(limitset_list_all[j])
        node = sorted(list(set(limitset_list_all[i])-elim))

        node_list.append(node)

    # ヒエラルキーのエッジを作成
    edge_list = []
    for i, s in enumerate(node_list):  # source
        for j, t in enumerate(node_list):  # target
            if hiermat[i, j] == 1:
                edge_list.append((s, t))

    return node_list, edge_list

def make_hiergraph(limitset_list):
    """return a hierarchy graph as a Agraph
    
    Parameters
    ----------
    limitset_list : list
        list of limitsets(buffering structures)
        
    Returns
    -------
    hier_agraph : Agraph
        hierarchy graph
    """
    #return an Agraph of the hierarchy
    hiernode_list, hieredge_list=make_hieredge(limitset_list)
    # add nodes with short names
    nodes = [short_name(' '.join(node)) for node in hiernode_list]
    hier = []
    for edge in hieredge_list:
        s = ' '.join(edge[0])
        t = ' '.join(edge[1])
        hier.append([short_name(s), short_name(t)])
    hier_graph = nx.DiGraph(hier)
    hier_graph.add_nodes_from(nodes)
    hier_agraph = nx.nx_agraph.to_agraph(hier_graph)

    # when node degree is large, node hight is fixed
    for node in hier_agraph.nodes():
        if hier_agraph.degree(node) > 30:
            hier_agraph.get_node(node).attr['height'] = 5
            hier_agraph.get_node(node).attr['width'] = 5

    hier_agraph.layout(prog='dot', args="-Nshape=box")
    return hier_agraph

def short_name(name):
    l = len(name)
    N = 40  # N文字に一回改行
    if l<N:
        return name
    else: 
        short_name = name [:N]
        for i in range(l//N):
            short_name += '\n'+name[(i+1)*N:(i+2)*N]
            # too long name is truncated
            if i == 5:
                short_name += '...'
                break
        return short_name