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
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import logging

from biolib.external.fasttree import FastTree
from biolib.external.raxml import RAxML


class Bootstrap():
    """Calculate bootstrap support."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def run(self, 
                input_tree, 
                msa_file, 
                tree_program,
                prot_model,
                num_replicates,
                output_dir):
        """Calculate bootstraps.

        Calculate support for tree using  the non-parametric
        bootstrap methods.

        Parameters
        ----------
        input_tree : str
            Tree requiring bootstrap support values.
        msa_file : str
            Multiple sequence alignment used to infer input tree (fasta format).
        tree_program : str
            Program to use for tree inference ['fasttree', 'raxml'].
        prot_model : str
            Protein substitution model for tree inference ['WAG', 'LG'].
        num_replicates : str
            Number of bootstrap replicates to perform.
        output_tree : float
            Output tree with bootstrap values.
        """

        if tree_program == 'fasttree':
            self.logger.info('Calculating bootstraps with FastTree under %s+GAMMA.' % prot_model)
            ft = FastTree(multithreaded=False)
            ft.bootstrap(input_tree,
                         msa_file,
                         'prot', prot_model,
                         num_replicates,
                         output_dir,
                         self.cpus)
        elif tree_program == 'raxml':
            self.logger.info('Calculating bootstraps with RAxML under PROTGAMMA%s.' % prot_model)
            raxml = RAxML(cpus=1)
            raxml.bootstrap(input_tree,
                            msa_file,
                            prot_model,
                            num_replicates,
                            output_dir,
                            self.cpus)
