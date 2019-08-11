#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from math import log2
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


class CoveredExamples(object):
    def __init__(self, n_pos, n_neg):
        self.__pos = n_pos
        self.__neg = n_neg
        self.__tot = n_pos + n_neg

    @property
    def num_pos(self):
        return self.__pos

    @property
    def num_neg(self):
        return self.__neg

    @property
    def num_covered(self):
        return self.__tot

    @property
    def entropy(self):
        p = self.num_pos / self.num_covered
        if p == 0 or p == 1:
            return 0

        return -p*log2(p) - (1-p)*log2(1-p) if p not in (0,1) else 0

    def __str__(self):
        return "[{}+, {}-]".format(self.num_pos, self.num_neg)



class DecisionTree(object):
    """ A decision tree is a (binary) tree whose internal nodes are labelled with a
    feature and the edges splits the feature on a threshold.
    Leaves are labelled with the predicted class.
    Features are genes and feature values represent gene expression data (credo)"""

    split_regex = re.compile(r"[:\ ]")


    def __init__(self):
        """ The tree is a recursive data structure. The root is the entire tree.
        __nodes attribute is a dictionary that uses gene names as keys and
        a list of nodes as the corresponding value. """
        self.__root = Node()
        self.__nodes = dict()

        self.__num_leaves = 0
        self.__filename = None

    def __get_candidates(self, node_id, threshold_value):
        #ultimo nodo incompleto con tale id e tale threshold
        candidates = filter(lambda node: node.threshold == threshold_value and not node.is_complete(), self.__nodes[node_id])
        #[node for node in self.__nodes[node_id] if node.threshold == threshold_value]
        return list(candidates) #candidate.pop() if candidate else False

    @property
    def root(self):
        return self.__root

    @property
    def filename(self):
        return self.__filename

    @filename.setter
    def filename(self, value):
        self.__filename = value.split("/")[-1].split(".")[0]

    @staticmethod
    def parse(filename, verbose=False):
        """ Parsing method for a single decision tree. Return the parsed tree"""

        current_tree = DecisionTree()
        current_tree.filename = filename

        parsing_failed = False
        prev_node, prev_threshold = None, None

        with open(filename) as f:
            for index, line in enumerate(f):
                try:
                    tokens = [token.strip() for token in DecisionTree.split_regex.split(line.replace("|", "")) if token != ""]

                    #extract gene id, relation and threshold 
                    gene, relation, threshold = tokens[:3]
                    #initialize edge 
                    edge = Edge(relation, threshold)
     
                    if ":" in line: #is leaf?
                        #extract predicted class 
                        label, ratio = tokens[3], tokens[4][1:-1] #remove '(' and ')'
                        edge.set_leaf(label, ratio) 
                        current_tree.__num_leaves += 1

                except IndexError:
                    print("Parsing failed at line {}: {}".format(index, line)) #print random bestemmie
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

        current_tree.__root.set_coverage()

        return current_tree


    @staticmethod
    def get_entropies(exploration):
        weighted_entropies = dict()
        tot_elements = exploration[0][0].covered.num_covered

        for node, depth in exploration:
            if depth not in weighted_entropies:
                weighted_entropies[depth] = list()
            weighted_entropies[depth].append((node.entropy*node.covered.num_covered, node.covered.num_covered))

        max_depth = max(weighted_entropies.keys())
        result = [(
            sum([ws for ws, c in weighted_entropies[depth]])/ tot_elements,
            sum([ws for ws, _ in weighted_entropies[depth]]) / sum([c for _, c in weighted_entropies[depth]]))
            for depth in range(max_depth+1)]
        return result

        # result = [for depth in range(len(weighted_entropies))]
        #
        # d = {n:list() for n in range(max_depth+1)}
        # for node, depth in explorations[0]:
        #     d[depth].append((node.entropy*node.covered.num_covered, node.covered.num_covered))
        # else:
        #     asd = (explorations[0][0][0].covered.num_covered)
        #     for depth in d.keys():
        #         weighted, elements = zip(*d[depth])
        #         d[depth] = sum(weighted) / asd#sum(elements)
        #         print("{}: {}".format(depth, d[depth]))

    def bfs(self, max_depth=None):
        """Breadth-first search (BFS) explores the tree in a level-wise fashion,
        starting from the root to the leaves."""

        open_nodes = [(self.__root, 0)]
        visited = list()

        while len(open_nodes) > 0:
            node, depth = open_nodes.pop(0)

            if not node.is_leaf():
                if max_depth is None or depth <= max_depth:
                    open_nodes.extend([(node.left_child, depth + 1), (node.right_child, depth + 1)])

            visited.append((node, depth))

        return visited

    def get_node(self, target):
        #dirty code?
        def get_info(node, depth, which = "sx"):
            child, edge = node.left_child, node.left_edge
            if which == "dx":
                child, edge = node.right_child, node.right_edge
            
            return [
                depth, 
                edge.relation, 
                edge.threshold, 
                child.covered.num_covered,
                child.covered.num_pos,
                child.covered.num_neg
            ]

        data = list()

        for nodes in self.__nodes.values():
            for node in nodes:
                if node.node_id == target:
                    #obtain node depth 
                    depth, tmp = 0, node.father
                    while tmp is not None: 
                        tmp = tmp.father 
                        depth += 1

                    data.append(get_info(node, depth))
                    data.append(get_info(node, depth, "dx"))
    
        data.sort(key=lambda l: l[1])
        return data


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
        #details about covered examples in training set
        self.covered = None

    @property
    def entropy(self):
    #    return self.covered.entropy

        if self.is_leaf():
            return self.covered.entropy
        else:
            #calculate weighted entropy for mutually esclusive nodes
            left = self.left_child.covered
            right = self.right_child.covered

            return (left.num_covered * left.entropy + right.num_covered * right.entropy) / self.covered.num_covered


    # def visit(self, depth=0):
    #     if not self.is_leaf():
    #         self.left_child.visit(depth+1)
    #         print("{}> {}, {}".format(depth, self.entropy, self.covered))
    #         self.right_child.visit(depth+1)
    #     else:
    #         print("{}) {}, {}".format(depth, self.entropy, self.covered))
    #

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

    def set_coverage(self):
        if not self.is_leaf():
            left_coverage = self.left_child.set_coverage()
            right_coverage = self.right_child.set_coverage()

            pos = left_coverage.num_pos + right_coverage.num_pos
            neg = left_coverage.num_neg + right_coverage.num_neg

            self.covered = CoveredExamples(pos, neg)
        return self.covered

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
    def left_edge(self):
        return self.__left
    
    @property
    def right_edge(self):
        return self.__right

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
    def relation(self):
        return self.__relation

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

        if ratio is not None:
            bigger, smaller = [float(x) for x in ratio.split("/")] if "/" in ratio \
                else (float(ratio), 0)
            self.covered = CoveredExamples(smaller, bigger) if label == "CX" \
                else CoveredExamples(bigger, smaller)


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
