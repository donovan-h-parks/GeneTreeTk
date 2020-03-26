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

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2020"
__credits__ = ["Donovan Parks","Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft@gmail.com"
__status__ = "Development"

import os
import sys
import logging

import biolib.seq_io as seq_io
import biolib.seq_tk as seq_tk
from biolib.common import concatenate_files
from biolib.taxonomy import Taxonomy
from biolib.external.blast import Blast
from biolib.external.diamond import Diamond
from biolib.external.execute import check_dependencies

from genetreetk.common import validate_seq_ids, extract_homologs_and_context, create_arb_metadata, remove_stop_codons, gene_distribution, extract_sequences
from genetreetk.msa_workflow import MsaWorkflow
from genetreetk.tree_workflow import TreeWorkflow

from numpy import (percentile as np_percentile,
                    mean as np_mean)

import dendropy
import extern

class ArbDbCreator():
    """TODO"""

    def __init__(self):
        """Initialization.

        Parameters
        ----------
        """

        check_dependencies(['mfqe'])

        self.logger = logging.getLogger('timestamp')

    def create_from_protein_alignment(self, **kwargs):
        # """Infer a gene tree for each query gene identified by blast.

        # Workflow for inferring a gene tree from sequences identified as being
        # homologs to a set of query proteins. Homologs are identified using BLASTP
        # and a set of user-defined parameters.

        # Parameters
        # ----------
        # query_proteins : str
        #     Fasta file containing query proteins.
        # db_file : str
        #     BLAST database of reference proteins.
        # custom_db_file : str
        #     Custom database of proteins.
        # taxonomy_file : str
        #     Taxonomic assignment of each reference genomes.
        # # custom_taxonomy_file : str
        # #     Taxonomic assignment of genomes in custom database.
        # evalue : float
        #     E-value threshold used to define homolog.
        # per_identity : float
        #     Percent identity threshold used to define a homolog.
        # per_aln_len : float
        #     Alignment length threshold used to define a homolog.
        # max_matches : int
        #     Maximum matches per query protein.
        # metadata : dict[genome_id] -> metadata dictionary
        #     Metadata for genomes.
        # homology_search : str
        #     Type of homology search to perform.
        # min_per_taxa : float
        #     Minimum percentage of taxa required to retain a column.
        # consensus : float
        #     Minimum percentage of the same amino acid required to retain column.
        # min_per_bp : float
        #     Minimum percentage of base pairs required to keep trimmed sequence.
        # use_trimAl : boolean
        #     Filter columns using trimAl.
        # restrict_taxon : str
        #     Restrict alignment to specific taxonomic group (e.g., k__Archaea).
        # msa_program : str
        #     Program to use for multiple sequence alignment ['mafft', 'muscle'].
        # tree_program : str
        #     Program to use for tree inference ['fasttree', 'raxml'].
        # prot_model : str
        #     Protein substitution model for tree inference ['WAG', 'LG', 'AUTO'].
        # skip_rooting : boolean
        #     Skip midpoint rooting if True.
        # output_dir : str
        #     Directory to store results.
        # """

    # req_args.add_argument('-i', '--protein_ids', help='IDs of proteins to be included in the database (newline separated)', required=True)
    # req_args.add_argument('-d', '--db_file', help='BLAST database of reference proteins', required=True)
    # req_args.add_argument('-t', '--taxonomy_file', help='taxonomic assignment of each reference genomes', required=True)
    # req_args.add_argument('-o', '--output_file', help='output file', required=True)
        alignment_file = kwargs.pop('alignment_file')
        db_file = kwargs.pop('db_file')
        taxonomy_file = kwargs.pop('taxonomy_file')
        output_file = kwargs.pop('output_file')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        # Read in protein IDs. Take only the bit before the space.
        # self.logger.info("Reading protein IDs ..")
        # ids_set = set()
        # with open(protein_ids_file) as f:
        #     for line in f:
        #         ids_set.add(line.strip())
        # ids = list(ids_set)
        # self.logger.info("Read in {} IDs to create arb DB for".format(len(ids)))

        # read taxonomy file
        self.logger.info('Reading taxonomy file.')
        taxonomy = Taxonomy().read(taxonomy_file)

        # setup metadata for ARB file
        src_dir = os.path.dirname(os.path.realpath(__file__))
        version_file = open(os.path.join(src_dir, 'VERSION'))

        metadata = {}
        metadata['genetreetk_version'] = version_file.read().strip()
        # metadata['genetreetk_query_proteins'] = query_proteins
        # metadata['genetreetk_db_file'] = db_file
        # metadata['genetreetk_taxonomy_file'] = taxonomy_file
        # metadata['genetreetk_blast_evalue'] = str(evalue)
        # metadata['genetreetk_blast_per_identity'] = str(per_identity)
        # metadata['genetreetk_blast_per_aln_len'] = str(per_aln_len)
        # metadata['genetreetk_blast_max_matches'] = str(max_matches)
        # metadata['genetreetk_homology_search'] = homology_search

        # metadata['genetreetk_msa_min_per_taxa'] = str(min_per_taxa)
        # metadata['genetreetk_msa_consensus'] = str(consensus)
        # metadata['genetreetk_msa_min_per_bp'] = str(min_per_bp)
        # metadata['genetreetk_msa_program'] = msa_program

        # metadata['genetreetk_tree_program'] = tree_program
        # metadata['genetreetk_tree_prot_model'] = prot_model

        # create ARB metadata file
        self.logger.info('Creating ARB metadata file.')
        create_arb_metadata({}, alignment_file, taxonomy,
                                 metadata,
                                 {}, {},
                                 output_file)

    def _blast_result_to_dictionary(self, blast_object, blast_output):
        # Assumes format is custom, as used above. Return a dictionary of query
        # to list of BlastHitCustom. Only take the first hit from each query/hit
        # pair.

        query_to_hits = {}
        for res in blast_object.read_hit(blast_output, 'custom'):
            query = res.query_id
            if query not in query_to_hits:
                query_to_hits[query] = []
            
            if len(query_to_hits[query]) == 0 or query_to_hits[query][-1].subject_id != res.subject_id:
                query_to_hits[query].append(res)
        return query_to_hits

    def _extract_sequences(self, sequence_ids, db_file, output_file):
        extern.run("mfqe --fasta-read-name-lists \/dev/stdin --input-fasta '{}' --output-fasta-files '{}' --output-uncompressed".format(
            db_file, output_file),
            stdin="\n".join(sequence_ids))