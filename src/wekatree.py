#!/usr/bin/python3
# -*- coding: utf-8 -*-

import re

class CompleteNodeException(Exception):
    def __init__(self, message, node_id):
        super(CompleteNodeException, self).__init__(message)
        self.node_id = node_id

class UndefinedEdgeException(Exception):
    def __init__(self, message, node_id):
        super(UndefinedEdgeException, self).__init__(message)
        self.node_id = node_id

class UnparsableTreeException(Exception):
    def __init__(self, message):
        super(UnparsableTreeException, self).__init__(message)


class DecisionTree(object):
    """ A decision tree is a (binary) tree whose internal nodes are labelled with a
    feature and the edges splits the feature on a threshold.
    Leaves are labelled with the predicted class.
    Features are genes and feature values represent gene expression data (credo)"""

    internal_regex = re.compile(r"^(\w+) (\S+) (\d+)$")
    leaf_regex = re.compile(r"^(\w+) (\S+) (\d+): (\w+) \((\S+)\)$")

    def __init__(self):
        """ The tree is a recursive data structure. The root is the entire tree.
        __nodes attribute is a dictionary that uses gene names as keys and
        a list of nodes as the corresponding value. """
        self.__root = Node()
        self.__nodes = dict()

        self.__num_leaves = 0

    def __get_candidates(self, node_id, threshold_value):
        #ultimo nodo incompleto con tale id e tale threshold
        candidates = filter(lambda node: node.threshold == threshold_value and not node.is_complete(), self.__nodes[node_id])
        #[node for node in self.__nodes[node_id] if node.threshold == threshold_value]
        return list(candidates) #candidate.pop() if candidate else False

    @staticmethod
    def parse(filename, verbose=False):
        """ Parsing method for a single decision tree. Return the parsed tree"""

        current_tree = DecisionTree()

        parsing_failed = False
        prev_node, prev_threshold = None, None

        with open(filename) as f:
            for line in f:
                try:
                    stripped = line.replace("|", "").strip()
                    is_leaf = ":" in stripped
                    result = DecisionTree.leaf_regex.match(stripped) if is_leaf \
                        else DecisionTree.internal_regex.match(stripped)

                    #extract gene id
                    gene = result.group(1)

                    #initialize the edge: extract edge label
                    relation, threshold = result.group(2), result.group(3)
                    edge = Edge(relation, threshold)

                    if is_leaf:
                        #extract predicted class
                        edge.set_leaf(label = result.group(4), ratio = result.group(5))
                        current_tree.__num_leaves += 1

                except AttributeError:
                    parsing_failed = True
                    break

                #initialize a new node
                new_node = Node(gene).set_edge(edge)

                if len(current_tree.__nodes) > 0:
                    if new_node.node_id in current_tree.__nodes:
                        #current gene is already present in the tree:
                        #obtaining incomplete nodes
                        candidates = current_tree.__get_candidates(new_node.node_id, edge.threshold)

                        if len(candidates) > 0:
                            #link edge to the most recent incomplete node
                            node = candidates[-1].set_edge(edge)

                            if verbose:
                                print("linking {} to father {}".format(gene, node.father), end="")
                        else:
                            #no available nodes: adding new node
                            father = current_tree.__get_candidates(prev_node, prev_threshold)[-1].add_child(new_node)
                            current_tree.__nodes[new_node.node_id].append(new_node)

                            if verbose:
                                print("duplicating node {} (father {})".format(gene, father), end="")
                    else:
                        #adding a new node
                        father = current_tree.__get_candidates(prev_node, prev_threshold)[-1].add_child(new_node)
                        current_tree.__nodes[new_node.node_id] = [new_node]

                        if verbose:
                            print("new node {} w/ father {}".format(gene, father), end="")
                else:
                    #tree initialization: setting the root
                    current_tree.__root = new_node
                    current_tree.__nodes[new_node.node_id] = [new_node]

                    if verbose:
                        print("Root: {}".format(new_node.node_id), end="")

                if verbose:
                    print("  ({})".format(edge.target.label) if is_leaf else "")

                prev_node = new_node.node_id
                prev_threshold = threshold

        if parsing_failed:
            raise UnparsableTreeException("Cannot parse {} file".format(filename))

        return current_tree


    def bfs(self, max_depth=None):
        """Breadth-first search (BFS) explores the tree in a level-wise fashion, starting
        from the root to the leaves."""

        open_nodes = [(self.__root, 0)]
        visited = list()

        while len(open_nodes) > 0:
            node, depth = open_nodes.pop(0)

            if not node.is_leaf():
                if max_depth is None or depth <= max_depth:
                    open_nodes.extend([(node.left_child, depth + 1), (node.right_child, depth + 1)])
                    visited.append((node, depth))

        return visited

    def __len__(self):
        """Returns the total number of nodes (internals and leaves) of the tree"""
        return sum([len(x) for x in self.__nodes.values()]) + self.__num_leaves

    def num_leaves(self):
        return self.__num_leaves



class Node(object):
    """An internal node is identified by a gene id """

    def __init__(self, node_id = None):
        self.__node_id = node_id
        self.__father = None
        self.__left = None
        self.__right = None

    def is_leaf(self):
        return False

    def is_complete(self):
        return bool(self.__left and self.__left.target and self.__right and self.__right.target)

    def set_edge(self, edge):
        """Add a edge to a node and returns the node"""
        if self.__left is None:
            self.__left = edge
        elif self.__right is None:
            self.__right = edge
        else:
            raise CompleteNodeException("aiutooo", self.node_id)

        return self

    def add_child(self, node):
        """Add a child node to another node and returns the father node """
        edge = self.__right if self.__right is not None else self.__left

        if edge is None:
            raise UndefinedEdgeException("OLBIA GRAN TURISMO", node.node_id)

        edge.target = node
        node.father = self

        return self

    @property
    def father(self):
        return self.__father

    @property
    def threshold(self):
        return self.__left.threshold

    @father.setter
    def father(self, father_node):
        self.__father = father_node

    @property
    def left_child(self):
        return self.__left.target

    @property
    def right_child(self):
        return self.__right.target

    @property
    def node_id(self):
        return self.__node_id

    def __str__(self):
        return self.node_id

    def __repr__(self):
        return str(self)

class Edge(object):
    """ Edges link the father node with its subtrees. Blabla """

    #arco: label + nodo target o label
    def __init__(self, relation, threshold):    # = None, target = None):
        self.__relation = relation
        self.__label = threshold
        self.__target = None #Node o LeafNode

    @property
    def threshold(self):
        return self.__label

    @property
    def target(self):
        return self.__target

    @target.setter
    def target(self, node):
        self.__target = node


    def set_leaf(self, label, ratio=None):
        self.__target = LeafNode(label, ratio)

class LeafNode(Node):
    """ A leaf node is a terminal node. A leaf has no subtrees and it is labelled
    with a predicted class """

    def __init__(self, label, ratio=None):
        super().__init__()
        self.__label = label
        self.__ratio = ratio if ratio else "?"

    def is_leaf(self):
        return True

    @property
    def label(self):
        return self.__label

    @property
    def ratio(self):
        return self.__ratio

    def __str__(self):
        return self.__label
