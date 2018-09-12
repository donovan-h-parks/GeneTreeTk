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

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2014"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import os
import sys
import ntpath
import logging

import biolib.seq_io as seq_io
from biolib.external.execute import check_dependencies
from biolib.external.fasttree import FastTree
from biolib.external.raxml import RAxML

import dendropy


class TreeWorkflow():
    """Workflow for inferring tree."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use during homology search.
        """

        check_dependencies(['FastTreeMP', 
                            'raxmlHPC-PTHREADS-SSE3'])

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def run(self, msa_file, 
                tree_program, 
                prot_model,
                output_dir):
        """Infer tree.

        Parameters
        ----------
        msa_file : str
          Multiple sequence alignment in fasta format.
        tree_program : str
          Program to use for tree inference ['fasttree', 'raxml'].
        prot_model : str
          Protein substitution model for tree inference ['WAG', 'LG', 'AUTO'].
        output_dir : str
          Directory to store results.
        """
        
        num_seqs = sum([1 for _, _ in seq_io.read_seq(msa_file)])
        if num_seqs <= 2:
            self.logger.error('Insufficient number of sequences in MSA to infer tree.')
            raise SystemExit('Tree inference failed.')
        
        output_file = ntpath.basename(msa_file)
        prefix = output_file[0:output_file.rfind('.')]
        suffix = output_file[output_file.rfind('.') + 1:]

        if tree_program == 'fasttree':
            self.logger.info('Inferring gene tree with FastTree using %s+GAMMA.' % prot_model)
            fasttree = FastTree(multithreaded=(self.cpus > 1))

            tree_unrooted_output = os.path.join(output_dir, prefix + '.unrooted.tree')
            tree_log = os.path.join(output_dir, prefix + '.tree.log')
            tree_output_log = os.path.join(output_dir, 'fasttree.log')
            fasttree.run(msa_file, 'prot', prot_model, tree_unrooted_output, tree_log, tree_output_log)
        elif tree_program == 'raxml':
            self.logger.info('Inferring gene tree with RAxML using PROTGAMMA%s.' % prot_model)
            
            # create phylip MSA file
            phylip_msa_file = msa_file.replace('.faa', '.phyx')
            cmd = 'seqmagick convert %s %s' % (msa_file, phylip_msa_file)
            os.system(cmd)
            
            # run RAxML
            raxml_dir = os.path.abspath(os.path.join(output_dir, 'raxml'))
            tree_output_log = os.path.join(output_dir, 'raxml.log')
            
            raxml = RAxML(self.cpus)
            tree_unrooted_output = raxml.run(phylip_msa_file, prot_model, raxml_dir)

        # root tree at midpoint
        seqs = seq_io.read(msa_file)
        if len(seqs) > 2:
            self.logger.info('Rooting tree at midpoint.')
            tree = dendropy.Tree.get_from_path(tree_unrooted_output, schema='newick', rooting="force-rooted", preserve_underscores=True)
            tree.reroot_at_midpoint(update_bipartitions=False)
        
        tree_output = os.path.join(output_dir, prefix + '.rooted.tree')
        tree.write_to_path(tree_output, schema='newick', suppress_rooting=True, unquoted_underscores=True)
        
        return tree_output