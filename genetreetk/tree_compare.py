###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2017'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import logging

import dendropy
from dendropy.calculate import treecompare

from biolib.newick import parse_label, create_label

class TreeCompare(object):
    """Compare pairs of trees."""
    
    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
    def _check_fractional_bootstraps(self, tree):
        """Check if bootstrap values are between [0, 1] and change to [0, 100]."""
        
        fractional_bootstrap = True
        for n in tree.preorder_node_iter():
            support, label, aux_info = parse_label(n.label) 
            if support is not None and support > 1.0:
                fractional_bootstrap = False
                break
        
        if fractional_bootstrap:
            for n in tree.preorder_node_iter():
                support, label, aux_info = parse_label(n.label) 
                if support is not None:
                    n.label = create_label(int(support*100 + 0.5), label, aux_info)
        
    def _prune(self, tree1, tree2):
        """Prune trees to common set of taxa."""
        
        # prune both trees to the set of common taxa
        taxa1 = set([t.taxon for t in tree1.leaf_node_iter()])  
        taxa2 = set([t.taxon for t in tree2.leaf_node_iter()])

        if len(taxa1) != len(taxa2) or taxa1 != taxa2:
            taxa_in_common = taxa1.intersection(taxa2)
            self.logger.info('Tree 1 contains %d taxa.' % len(taxa1))
            self.logger.info('Tree 2 contains %d taxa.' % len(taxa2))
            self.logger.info('Pruning trees to the %d taxa in common.' % len(taxa_in_common))
            
            if not taxa_in_common:
                self.logger.error('No taxa in common.')
                sys.exit(-1)
            
            tree1.retain_taxa(taxa_in_common)
            tree2.retain_taxa(taxa_in_common)
        
    def _read_trees(self, tree1_file, tree2_file, taxa_list, prune=True):
        """Read trees from file."""
        
        tns = dendropy.TaxonNamespace()
        tree1 = dendropy.Tree.get_from_path(tree1_file, 
                                            schema='newick',
                                            rooting='force-unrooted', 
                                            preserve_underscores=True,
                                            taxon_namespace=tns)
                      
        tree2 = dendropy.Tree.get_from_path(tree2_file, 
                                            schema='newick',
                                            rooting='force-unrooted', 
                                            preserve_underscores=True,
                                            taxon_namespace=tns)
                                            
        # check if bootstrap values are fractional
        self._check_fractional_bootstraps(tree1)
        self._check_fractional_bootstraps(tree2)
                                            
        # prune trees to specified taxa
        if taxa_list:
            taxa_to_keep = set()
            for line in open(taxa_list):
                taxa_to_keep.add(line.strip().split('\t')[0])
            tree1.retain_taxa_with_labels(taxa_to_keep)
            tree2.retain_taxa_with_labels(taxa_to_keep)
                                                            
        # verify trees were defined over the same set of taxa
        if prune:
            self._prune(tree1, tree2)
 
        tree1.encode_bipartitions()
        tree2.encode_bipartitions()
        
        assert tree1.taxon_namespace is tree2.taxon_namespace

        return tree1, tree2

    def euclidean(self, tree1, tree2, taxa_list):
        """Calculate Euclidean distance between two trees."""
        
        tree1, tree2 = self._read_trees(tree1, tree2, taxa_list) 
            
        return treecompare.euclidean_distance(tree1, tree2)
        
    def robinson_foulds(self, tree1, tree2, taxa_list):
        """Calculate Robinson-Foulds (i.e., symmetric_difference) distance between two trees."""
        
        tree1, tree2 = self._read_trees(tree1, tree2, taxa_list)
            
        rf = treecompare.symmetric_difference(tree1, tree2)
        
        num_taxa = len([t for t in tree1.leaf_node_iter()])
        normalized_rf = float(rf) / (2*(num_taxa - 3))
            
        return rf, normalized_rf
        
    def weighted_robinson_foulds(self, tree1, tree2, taxa_list):
        """Calculate weighted Robinson-Foulds distance between two trees."""

        tree1, tree2 = self._read_trees(tree1, tree2, taxa_list) 
            
        return treecompare.weighted_robinson_foulds_distance(tree1, tree2)
        
    def report_all(self, tree1_file, tree2_file):
        """Report results for all tree comparison statistics."""
        
        tree1, tree2 = self._read_trees(tree1_file, tree2_file)

        print('Euclidean: %f' % self.euclidean(tree1, tree2))
        print('Robinson-Foulds: %f' % self.robinson_foulds(tree1, tree2))
        print('Weighted Robinson-Foulds: %f' % self.weighted_robinson_foulds(tree1, tree2))
        
    def report_missing_splits(self, ref_tree, compare_tree, min_support, taxa_list):
        """Report supported bipartitions in reference tree not in comparison tree."""
        
        ref_tree, compare_tree = self._read_trees(ref_tree, compare_tree, taxa_list)
        
        incongruent = 0
        print('Missing splits with support >= %f:' % min_support)
        for n in ref_tree.preorder_node_iter(lambda n: not n.is_leaf()):
            support, label, aux_info = parse_label(n.label)
            
            if support >= min_support:
                if n.bipartition not in compare_tree.bipartition_encoding:
                    incongruent += 1
                    if label:
                        print(label, n.edge.length)
                    else:
                        print(','.join([t.taxon.label for t in n.leaf_iter()]))
                    
        print('Missing splits: %d' % incongruent)
        
    def _avg_descendant_rate(self, tree):
        """Calculate average rate of divergence for each nodes in a tree.

        The average rate is the arithmetic mean of the
        branch length to all descendant taxa.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
        """

        # calculate the mean branch length to extant taxa
        for node in tree.postorder_node_iter():
            avg_div = 0
            if node.is_leaf():
                node.mean_dist = 0.0
                node.num_taxa = 1
            else:
                node.num_taxa = sum([1 for _ in node.leaf_iter()])
                for c in node.child_node_iter():
                    num_tips = c.num_taxa
                    avg_div += (float(c.num_taxa) / node.num_taxa) * (c.mean_dist + c.edge_length)

            node.mean_dist = avg_div

    def _rel_dist(self, tree):
        """Calculate relative distance to each internal node.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
          rel_dists: relative distance of node between root and extant organisms
        """
        
        self._avg_descendant_rate(tree)

        for node in tree.preorder_node_iter():
            if node == tree.seed_node:
                node.rel_dist = 0.0
            elif node.is_leaf():
                node.rel_dist = 1.0
            else:
                a = node.edge_length
                b = node.mean_dist
                x = node.parent_node.rel_dist

                if (a + b) != 0:
                    rel_dist = x + (a / (a + b)) * (1.0 - x)
                else:
                    # internal node has zero length to parent,
                    # so should have the same relative distance
                    # as the parent node
                    rel_dist = x

                node.rel_dist = rel_dist
        
    def _supported(self, ref_tree, compare_tree, min_support, max_depth):
        """Determine supported bipartitions in reference tree not in comparison tree."""
        
        congruent = 0
        congruent_w = 0
        incongruent = 0
        incongruent_w = 0
        nontrivial_splits = 0
        nontrivial_splits_w = 0
        congruent_splits = {}
        incongruent_splits = {}
        for n in ref_tree.preorder_node_iter(lambda n: not n.is_leaf()):
            if not n.parent_node:
                continue

            nontrivial_splits += 1
            nontrivial_splits_w += n.edge.length
            
            support, label, aux_info = parse_label(n.label)
            if support is None or (support >= min_support and n.rel_dist <= max_depth):
                split_lca = n.child_nodes()[0].leaf_nodes()[0].taxon.label
                split_lca += '|'
                split_lca += n.child_nodes()[1].leaf_nodes()[0].taxon.label
                if n.bipartition not in compare_tree.bipartition_encoding:
                    incongruent += 1
                    incongruent_w += n.edge.length
                    incongruent_splits[split_lca] = (n.edge.length, support)
                else:
                    congruent += 1
                    congruent_w += n.edge.length
                    congruent_splits[split_lca] = (n.edge.length, support)

        return congruent, congruent_w, incongruent, incongruent_w, nontrivial_splits, nontrivial_splits_w, congruent_splits, incongruent_splits
        
    def _shared_support(self, tree1, tree2, min_support, max_depth):
        """Determine supported bipartitions common to a pair of trees."""
        
        assert tree1.taxon_namespace is tree2.taxon_namespace
        
        common_supported_splits = 0
        common_supported_splits_w = 0
        null_support = False
        for n in tree1.preorder_node_iter(lambda n: not n.is_leaf()):
            if not n.parent_node:
                continue

            support, label, aux_info = parse_label(n.label)
            if support is None:
                null_support = True
                
            if support is None or (support >= min_support and n.rel_dist <= max_depth):
                if n.bipartition in tree2.bipartition_encoding:
                    edge2 = tree2.bipartition_edge_map[n.bipartition]
                    support2, label2, aux_info2 = parse_label(edge2.head_node.label)

                    if support2 is None or (support2 >= min_support and edge2.head_node.rel_dist <= max_depth):
                        common_supported_splits += 1
                        common_supported_splits_w += n.edge.length + edge2.length
  
        if null_support:
            self.logger.warning('Some internal nodes lack support values and were treated as supported.')
        
        return common_supported_splits, common_supported_splits_w
        
    def supported_splits(self, tree1_file, tree2_file, split_file, min_support, max_depth, taxa_list):
        """Supported bipartitions of common taxa shared between two trees."""
        
        tree1, tree2 = self._read_trees(tree1_file, tree2_file, taxa_list)
        
        self._rel_dist(tree1)
        self._rel_dist(tree2)

        r = self._supported(tree1, tree2, min_support, max_depth)
        (congruent12, congruent12_w, 
            incongruent12, incongruent12_w, 
            nontrivial_splits1, nontrivial_splits1_w,
            congruent_splits12, incongruent_splits12) = r
        
        r = self._supported(tree2, tree1, min_support, max_depth)
        (congruent21, congruent21_w, 
            incongruent21, incongruent21_w, 
            nontrivial_splits2, nontrivial_splits2_w,
            congruent_splits21, incongruent_splits21) = r
        
        common_supported_splits, common_supported_splits_w = self._shared_support(tree1, tree2, min_support, max_depth)

        print('')
        print('T1 -> T2')
        print('Non-trivial splits in T1: %d' % nontrivial_splits1)
        print('Weight of non-trivial splits in T1: %g' % nontrivial_splits1_w)
        print('Well-supported splits in T1: %d' % (congruent12 + incongruent12))
        print('Weight of well-supported splits in T1: %g' % (congruent12_w + incongruent12_w))
        print('Well-supported splits in T1 present in T2: %d' % congruent12)
        print('Weight of well-support splits in T1 present in T2: %g' % congruent12_w)

        print('')
        print('T2 -> T1')
        print('Non-trivial splits in T2: %d' % nontrivial_splits2)
        print('Weight of non-trivial splits in T2: %g' % nontrivial_splits2_w)
        print('Well-supported splits in T2: %d' % (congruent21 + incongruent21))
        print('Weight of well-supported splits in T2: %g' % (congruent21_w + incongruent21_w))
        print('Well-supported splits in T2 present in T1: %d' % congruent21)
        print('Weight of well-support splits in T2 present in T1: %g' % congruent21_w)
        
        print('')
        print('T1 <-> T2')
        print('Non-trivial splits in T1 and T2: %d' % (nontrivial_splits1 + nontrivial_splits2))
        print('Weight of non-trivial splits in T1 and T2: %g' % (nontrivial_splits1_w + nontrivial_splits2_w))
        print('Well-supported splits in T1 and T2: %d' % (congruent12 + incongruent12 + congruent21 + incongruent21))
        supported_splits_t12 = congruent12_w + incongruent12_w + congruent21_w + incongruent21_w
        print('Weight of well-supported splits in T1 and T2: %g' % supported_splits_t12)
        print('Well-supported splits in common between T1 and T2: %d' % common_supported_splits)
        print('Weight of well-support splits in common between T1 and T2: %g' % common_supported_splits_w)
        print('Fraction of common well-supported splits: %g' % (common_supported_splits_w / supported_splits_t12))
        
        if split_file:
            fout = open(split_file, 'w')
            fout.write('Split\tStatus\tTree Comparison\tBranch Length\tSupport\n')
            for split_lca, (length, support) in congruent_splits12.items():
                if support is not None:
                    fout.write('%s\t%s\t%s\t%f\t%.2f\n' % (split_lca, 'Congruent', 'T1 -> T2', length, support))
            for split_lca, (length, support) in congruent_splits21.items():
                if support is not None:
                    fout.write('%s\t%s\t%s\t%f\t%.2f\n' % (split_lca, 'Congruent', 'T2 -> T1', length, support))
            for split_lca, (length, support) in incongruent_splits12.items():
                if support is not None:
                    fout.write('%s\t%s\t%s\t%f\t%.2f\n' % (split_lca, 'Incongruent', 'T1 -> T2', length, support))
            for split_lca, (length, support) in incongruent_splits21.items():
                if support is not None:
                    fout.write('%s\t%s\t%s\t%f\t%.2f\n' % (split_lca, 'Incongruent', 'T2 -> T1', length, support))
            fout.close()
            
           
