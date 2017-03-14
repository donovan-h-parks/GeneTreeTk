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
import logging

from genetreetk.msa_workflow import MsaWorkflow
from genetreetk.tree_workflow import TreeWorkflow

import biolib.seq_io as seq_io

import dendropy

class Reduce():
    """Workflow for inferring a tree over a reduced set of genes."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use during homology search.
        """

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def read_ids(self, gene_id_file):
        """Read gene id file.

        Read file with gene ids to retain. Each id
        should be put on a separate line.

        Parameters
        ----------
        gene_id_file : str
            File with gene ids.

        Returns
        -------
        set
          Gene or genome ids to retain.
        """

        genes_to_retain = set()
        for line in open(gene_id_file):
            genes_to_retain.add(line.split()[0].strip())

        return genes_to_retain

    def run(self, 
                homolog_file, 
                gene_id_file, 
                taxonomy_file, 
                min_per_taxa,
                consensus,
                min_per_bp,
                use_trimAl,
                msa_program,
                tree_program, 
                prot_model, 
                output_dir):
        """Infer a tree over a reduced set of genes.

        Filter a set of homolgs to a specified set of gene ids, 
        and infer tree over this reduced set of proteins.

        Parameters
        ----------
        homolog_file : str
            Fasta file containing homologs.
        gene_ids : str
            File with gene ids to retain in tree.
        taxonomy_file : str
            Taxonomic assignment of each reference genomes.
        min_per_taxa : float
            Minimum percentage of taxa required to retain a column.
        consensus : float
            Minimum percentage of the same amino acid required to retain column.
        min_per_bp : float
            Minimum percentage of base pairs required to keep trimmed sequence.
        use_trimAl : boolean
            Filter columns using trimAl.
        msa_program : str
            Program to use for multiple sequence alignment ['mafft', 'muscle'].
        tree_program : str
            Program to use for tree inference ['fasttree', 'raxml'].
        prot_model : str
            Protein substitution model for tree inference ['WAG', 'LG', 'AUTO'].
        output_dir: str
            Output directory.
        """

        # generate msa with reduced sequences
        self.logger.info('Extracting sequences to retain.')
        genes_to_retain = self.read_ids(gene_id_file)

        seqs = seq_io.read_fasta(homolog_file)
        reduced_seqs = {}
        for seq_id, seq in seqs.iteritems():
            if seq_id in genes_to_retain:
                reduced_seqs[seq_id] = seq

        reduced_homolog_file = homolog_file[0:homolog_file.rfind('.')]
        reduced_homolog_file += '.reduced.' + homolog_file[homolog_file.rfind('.') + 1:]
        seq_io.write_fasta(reduced_seqs, reduced_homolog_file)

        self.logger.info('Retained %d sequences.' % len(reduced_seqs))

        # infer multiple sequence alignment
        msa = MsaWorkflow(self.cpus)
        trimmed_msa_output = msa.run(reduced_homolog_file,
                                        min_per_taxa, 
                                        consensus, 
                                        min_per_bp, 
                                        use_trimAl, 
                                        msa_program,
                                        output_dir)

        # infer tree
        tw = TreeWorkflow(self.cpus)
        tree_output = tw.run(trimmed_msa_output,
                                tree_program,
                                prot_model,
                                output_dir)
                                
        # create tax2tree consensus map and decorate tree
        self.logger.info('Decorating internal tree nodes with tax2tree.')
        t2t_tree = tree_output.replace('.tree', '.tax2tree.tree')
        os.system('t2t decorate -m %s -t %s -o %s' % (taxonomy_file, tree_output, t2t_tree))
