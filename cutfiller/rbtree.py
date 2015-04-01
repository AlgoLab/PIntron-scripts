#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import getrefcount
import sys

import copy

class Leaf(object):
    """Leaf of an red-black tree."""
    def __init__(self):
        self.color = "BLACK"
        self.parent = None

    def insert(self, interval):
        node = RBIntervalNode(interval, "RED", self.parent)
        if self.parent.right == self:
            self.parent.right = node
        else:
            self.parent.left = node
        return [node]
    
    def s(self, depth):
        """DEBUG ONLY: Print informations about this leaf."""
        if not self.parent:
            dire = "ROOT"
        elif self is self.parent.left:
            dire = "L"
        else:
            dire = "R"
        s = "{0}{3} -{4}- : {1} - {2}\n".format("   "*depth, "LEAF", self.color, depth, dire)
        return s        

    def search(self, interval):
        return []

    def read_hit(self, interval):
        pass

class RBIntervalNode(object):
    """Base element of a red-black tree. Describe a genomic interval."""
    def __init__(self, root=None, color="BLACK", parent=None, left=None, right=None):
        self.root = root
        self.color = color
        self.parent = parent
        if not left:
            self.left = Leaf()
        else:
            self.left = left
        self.left.parent = self
        if not right:
            self.right = Leaf()
        else:
            self.right = right
        self.right.parent = self
        self.nprev = None # predecessor
        self.previnit = None # end position of prev
        self.nnext = None # successor
        self.nextinit = None # start position of next
        self.cutnext = False # Multiple next node? 
        self.cutprev = False # Multiple prev node?
        self.int_cover = [] # intervals covered by at least one read
        self.num_hit = 0 # read mapped to this block
        self.int_len = 0 # length of this node (without intron)

    def read_hit(self, interval):
        """Test the parameter falls within this interval and update values accordingly."""
        #print >> sys.stderr, "INTERVALLO:", interval
        if interval[0] < self.root[0]:
            self.left.read_hit(interval)
        if interval[1] > self.root[1]:
            self.right.read_hit(interval)
        if interval[1] > self.root[0] and interval[0] < self.root[1]:
            cover_interval = [interval[0], interval[1]]
            if cover_interval[0] < self.root[0]:
                cover_interval[0] = self.root[0]
            if cover_interval[1] > self.root[1]:
                cover_interval[1] = self.root[1]
            self.int_cover.append(cover_interval)
            self.int_cover.sort(key=lambda tup : tup[0])
            y = [self.int_cover[0]]
            for x in self.int_cover[1:]:
                if y[-1][1] < x[0]:
                    y.append(x)
                elif y[-1][1] == x[0]:
                    y[-1][1] = x[1]
                elif y[-1][1] < x[1]:
                    y[-1][1] = x[1]
            self.int_cover = y
            self.num_hit += 1

    def search(self, interval):
        """Test if this interval falls within the parameter and return the correct
nodes describing it."""
        ret_int = []
        if interval[1] < self.root[0]:
            ret_int += self.left.search(interval)
        elif self.root[1] < interval[0]:
            ret_int += self.right.search(interval)
        else:
            if interval[0] < self.root[0]:
                new_interval = [interval[0], self.root[0]-1]
                interval[0] = self.root[0]
                ret_int += self.left.search(new_interval)
                ret_int += self.search(interval)
            elif interval[1] > self.root[1]:
                new_interval = [self.root[1]+1, interval[1]]
                interval[1] = self.root[1]
                ret_int += self.search(interval)
                ret_int += self.right.search(new_interval)
            elif interval == self.root:
                ret_int.append(self)
        return ret_int

    def insert(self, interval):
        """Insert interval in this subtree."""
        ret_int = []
        if self.root is None:
            # First node
            self.root = interval
            return [self.root]
        elif self.root == interval:
            return [self]
        elif self.root[0] > interval[1]:
            # interval < self
            if self.left.__class__ is Leaf().__class__:
                self.left = RBIntervalNode(interval, "RED", self)
                return [self.left]
            else:
                return self.left.insert(interval)
        elif self.root[1] < interval[0]:
            # interval > self
            if self.right.__class__ is Leaf().__class__:
                self.right = RBIntervalNode(interval, "RED", self)
                return [self.right]
            else:
                return self.right.insert(interval)
        else:
            # interval e self is sovrappongono in qualche modo
            if interval[0] < self.root[0]:
                new_interval = [self.root[0], interval[1]]
                interval[1] = self.root[0]-1
                
                # Questo inserimento ricade nella prima casistica.
                # self.prev viene modificato.
                ret_int += self.insert(interval)
                
                # Questo inserimento ricade nell'ultima casistica.
                # se new_interval != self.root self.prev e self.succ 
                # vengono modificati
                ret_int += self.insert(new_interval)

                return ret_int

            elif interval[1] > self.root[1]:
                new_interval = [interval[0], self.root[1]]
                interval[0] = self.root[1]+1

                # Questo inserimento non dovrebbe inserire nulla
                ret_int += self.insert(new_interval)
                
                ret_int += self.insert(interval)

                return ret_int

            else:
                if interval[0] == self.root[0]:                
                    # interval inizia dove inizia self.root e finisce prima di 
                    # self.root
                    self.root[0] = interval[1]
                    interval[1] = interval[1] -1
                    if interval[0] > interval[1]:
                        return [self]
                    if self.left.__class__ is Leaf().__class__:
                        self.left = RBIntervalNode(interval, "RED", self)
                        if self.prev:
                            self.prev.succ = self.left
                            self.left.prev = self.prev
                            self.left.cutprev = self.cutprev
                        self.left.succ = self
                        self.prev = None
                        self.cutprev = False
                        ret_int = [self.left]
                    else:
                        self.prev = None
                        ret_int = self.left.insert(interval)
                        ret_int[0].cutprev = self.cutprev
                        ret_int[-1].succ = self
                        self.cutprev = False
                    
                if interval[1] == self.root[1]:
                    # interval finisce dove finisce self.root ma inizia dopo di
                    # self.root
                    self.root[1] = interval[0]
                    interval[0] = interval[0] +1
                    if interval[1] < interval[0]:
                        return [self]
                    if self.right.__class__ is Leaf().__class__:
                        self.right = RBIntervalNode(interval, "RED", self)
                        if self.succ:
                            self.succ.prev = self.right
                            self.right.succ = self.succ
                            self.right.cutnext = self.cutnext
                        self.right.prev = self
                        self.succ = None
                        self.cutnext = False
                        ret_int.append(self.right)
                    else:
                        self.succ = None
                        ret_int = self.right.insert(interval)
                        ret_int[-1].cutnext = self.cutnext
                        ret_int[0].prev = self
                        self.cutnext = False
                return ret_int

            if self.root[0] < interval[0] < interval[1] < self.root[1]:
                # interval cade in mezzo a self.root
                left_interval = [self.root[0], interval[0]-1]
                right_interval = [interval[1]+1, self.root[1]]
                
                # non è necessario modificare self.root in quanto
                # le seguenti chiamate ricorsive vengono effettuate sullo
                # stesso nodo e ne modificano il valore
                ret_int += self.insert(left_interval)
                ret_int += self.insert(right_interval)
                return ret_int
        return ret_int

    def __str__(self):
        return "ELEMENT: {0}".format(self.root)

    def s(self, depth):
        """DEBUG ONLY: Print informations about this leaf."""
        if not self.parent:
            dire = "ROOT"
        elif self is self.parent.left:
            dire = "L"
        else:
            dire = "R"
        s = "{0}{3} -{4}- : {1} - {2}\n".format("   "*depth, self, self.color, depth, dire)
        s += self.left.s(depth+1)
        s += self.right.s(depth+1)
        return s

    def __eq__(self, other):
        return self.root == other.root
    
    def __ne__(self, other):
        return self.root != other.root

    @property
    def prevsucc(self):
        return [self.prev, self.succ]

    @property
    def succ(self):
        """Return the successor node."""
        return self.nnext
    
    @succ.setter
    def succ(self, value):
        """Set the successor node. """
        if value is not None and self.nnext is not None:
            if self.nextinit != value.root[0]:
                self.cutnext = True
        self.nnext = value
        if self.nnext:
            self.nextinit = self.nnext.root[0]
        else:
            self.nextinit = None

    @property
    def prev(self):
        """Return the predecessor node."""
        return self.nprev
    
    @prev.setter
    def prev(self, value):
        """Set the predecessor node."""
        if value is not None and self.nprev is not None:
            if self.previnit != value.root[1]:
                self.cutprev = True
        self.nprev = value
        if self.nprev:
            self.previnit = self.nprev.root[1]
        else:
            self.previnit = None

    @property
    def grandparent(self):
        """Return the grandparent, if present."""
        if self.parent:
            return self.parent.parent
        else:
            return None

    def get_max( self ):
        if self.right.__class__ is Leaf( ).__class__:
            return self.root[ 1 ]
        else:
            return self.right.get_max( )

class RBIntervalTree(object):
    """Red-black interval tree."""
    def __init__(self):
        self.root = None

    def read_hit(self, interval):
        """Test if the parameter is present in this tree and update values accordingly.

        """
        if not self.root:
            return
        else:
            self.root.read_hit(interval)
            pass
        
    def rbinsert(self, interval):
        """Insert interval in this tree."""
        added_nodes = []
        if not self.root:
            self.root = RBIntervalNode(interval, "BLACK")
            added_nodes = [self.root]
        else:
            added_nodes = self.root.insert(interval)
            for node in added_nodes:
                while node is not self.root and node and node.parent.color is "RED":
                    if node.parent is node.grandparent.left:
                        y = node.grandparent.right
                        if y.color is "RED":
                            # Recolor
                            node.parent.color = "BLACK"
                            y.color = "BLACK"
                            node.grandparent.color = "RED"
                            node = node.grandparent
                        else:
                            if node is node.parent.right:
                                # Rotate left
                                node.parent.right = node.left
                                node.parent.right.parent = node.parent
                                node.left = node.parent
                                node.parent = node.grandparent
                                node.parent.left = node
                                node.left.parent = node
                            else:
                                node = node.parent
                            
                            p = node.parent
                            if p.parent is None:
                                self.root = node

                            # Rotate right
                            p.left = node.right
                            p.left.parent = p
                            node.right = p
                            node.parent = p.parent
                            p.parent = node
                            if node.parent is not None:
                                if p is node.parent.right:
                                    node.parent.right = node
                                else:
                                    node.parent.left = node
                            node.color = "BLACK"
                            p.color = "RED"
                            node = node.parent
                    else:
                        y = node.grandparent.left
                        if y.color is "RED":
                            # Recolor
                            node.parent.color = "BLACK"
                            y.color = "BLACK"
                            node.grandparent.color = "RED"
                            node = node.grandparent
                        else:
                            if node is node.parent.left:
                                # Rotate right
                                node.parent.left = node.right
                                node.parent.left.parent = node.parent
                                node.right = node.parent
                                node.parent = node.grandparent
                                node.parent.right = node
                                node.right.parent = node
                            else:
                                node = node.parent

                            p = node.parent
                            if p.parent is None:
                                self.root = node

                            # Rotate left
                            p.right = node.left
                            p.right.parent = p
                            node.left = p
                            node.parent = p.parent
                            p.parent = node
                            if node.parent is not None:
                                if p is node.parent.right:
                                    node.parent.right = node
                                else:
                                    node.parent.left = node
                            node.color = "BLACK"
                            p.color = "RED"
                            node = node.parent
                self.root.color = "BLACK"
        return added_nodes

    def search(self, interval):
        """Search interval in this tree and return the nodes that describe it."""
        return self.root.search(interval)

    def __str__(self):
        if not self.root:
            return ""
        return self.root.s(0)

    def get_max( self ):
        return self.root.get_max( )
