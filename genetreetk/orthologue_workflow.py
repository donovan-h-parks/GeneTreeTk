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

from genetreetk.common import validate_seq_ids, extract_homologs_and_context, create_arb_metadata, remove_stop_codons, gene_distribution
from genetreetk.msa_workflow import MsaWorkflow
from genetreetk.tree_workflow import TreeWorkflow

from numpy import (percentile as np_percentile,
                    mean as np_mean)

import dendropy


class OrthologueWorkflow():
    """Blast-based workflow for building gene trees after finding orthologues."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use during homology search.
        """

        check_dependencies(['blastp',
                            'mafft',
                            'muscle',
                            'FastTreeMP',
                            'raxmlHPC-PTHREADS-SSE3',
                            't2t',
                            'seqmagick',
                            'trimal'])

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def run(self, **kwargs):
        """Infer a gene tree for each query gene identified by blast.

        Workflow for inferring a gene tree from sequences identified as being
        homologs to a set of query proteins. Homologs are identified using BLASTP
        and a set of user-defined parameters.

        Parameters
        ----------
        query_proteins : str
            Fasta file containing query proteins.
        db_file : str
            BLAST database of reference proteins.
        custom_db_file : str
            Custom database of proteins.
        taxonomy_file : str
            Taxonomic assignment of each reference genomes.
        # custom_taxonomy_file : str
        #     Taxonomic assignment of genomes in custom database.
        evalue : float
            E-value threshold used to define homolog.
        per_identity : float
            Percent identity threshold used to define a homolog.
        per_aln_len : float
            Alignment length threshold used to define a homolog.
        max_matches : int
            Maximum matches per query protein.
        metadata : dict[genome_id] -> metadata dictionary
            Metadata for genomes.
        homology_search : str
            Type of homology search to perform.
        min_per_taxa : float
            Minimum percentage of taxa required to retain a column.
        consensus : float
            Minimum percentage of the same amino acid required to retain column.
        min_per_bp : float
            Minimum percentage of base pairs required to keep trimmed sequence.
        use_trimAl : boolean
            Filter columns using trimAl.
        restrict_taxon : str
            Restrict alignment to specific taxonomic group (e.g., k__Archaea).
        msa_program : str
            Program to use for multiple sequence alignment ['mafft', 'muscle'].
        tree_program : str
            Program to use for tree inference ['fasttree', 'raxml'].
        prot_model : str
            Protein substitution model for tree inference ['WAG', 'LG', 'AUTO'].
        skip_rooting : boolean
            Skip midpoint rooting if True.
        output_dir : str
            Directory to store results.
        """
        query_proteins = kwargs.pop('query_proteins')
        db_file = kwargs.pop('db_file')
        # custom_db_file = kwargs.pop('custom_db_file')
        taxonomy_file = kwargs.pop('taxonomy_file')
        # custom_taxonomy_file = kwargs.pop('custom_taxonomy_file')
        evalue = kwargs.pop('evalue')
        per_identity = kwargs.pop('per_identity')
        per_aln_len = kwargs.pop('per_aln_len')
        max_matches = kwargs.pop('max_matches')
        homology_search = kwargs.pop('homology_search')
        min_per_taxa = kwargs.pop('min_per_taxa')
        consensus = kwargs.pop('consensus')
        min_per_bp = kwargs.pop('min_per_bp')
        use_trimAl = kwargs.pop('use_trimAl')
        restrict_taxon = kwargs.pop('restrict_taxon')
        msa_program = kwargs.pop('msa_program')
        tree_program = kwargs.pop('tree_program')
        prot_model = kwargs.pop('prot_model')
        skip_rooting = kwargs.pop('skip_rooting')
        output_dir = kwargs.pop('output_dir')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        # validate query sequence names for use with GeneTreeTk
        validate_seq_ids(query_proteins)

        # read taxonomy file
        self.logger.info('Reading taxonomy file.')
        taxonomy = Taxonomy().read(taxonomy_file)

        # report distribution of query genes
        mean_len, max_len, min_len, p10, p50, p90 = gene_distribution(query_proteins)
        self.logger.info('Query gene lengths: min, mean, max = %d, %.1f, %d | p10, p50, p90 = %.1f, %.1f, %.1f' % (
                                                                                        min_len, mean_len, max_len,
                                                                                        p10, p50, p90))

        # identify homologs using BLASTP
        self.logger.info('Identifying homologs using %s.' % homology_search)
        blast = Blast(self.cpus)
        blast_output = os.path.join(output_dir, 'reference_hits.tsv')
        # Commented out below for debug
        # if homology_search == 'diamond':
        #     raise Exception("DIAMOND is not implemented all the way through")
        #     diamond = Diamond(self.cpus)
        #     diamond.blastp(query_proteins, db_file, evalue, per_identity, per_aln_len, max_matches, blast_output, output_fmt='custom')
        # else:
        #     blast.blastp(query_proteins, db_file, blast_output, evalue, max_matches, output_fmt='custom', task=homology_search)
        homologs = blast.identify_homologs(blast_output, evalue, per_identity, per_aln_len)
        self.logger.info('Identified %d homologs in reference database.' % len(homologs))

        if len(homologs) == 0:
            self.logger.error('No homologs were identified at all. Gene tree cannot be inferred.')
            sys.exit()

        # extract homologs
        self.logger.info('Extracting homologs and determining local gene context.')
        db_homologs_tmp = os.path.join(output_dir, 'homologs_db.tmp')

        gene_precontext, gene_postcontext = extract_homologs_and_context(set(homologs.keys()), db_file, db_homologs_tmp)

        # report gene length distribution of homologs
        mean_len, max_len, min_len, p10, p50, p90 = gene_distribution(db_homologs_tmp)
        self.logger.info('Homolog gene lengths: min, mean, max = %d, %.1f, %d | p10, p50, p90 = %.1f, %.1f, %.1f' % (
                                                                                        min_len, mean_len, max_len,
                                                                                        p10, p50, p90))

        # Split up blast hits into [{query_id->[hits]}] structure
        query_to_hits = {}
        for hit_name, blast_hit in homologs.items():
            query = blast_hit.query_id
            if query not in query_to_hits:
                query_to_hits[query] = []
            query_to_hits[query].append(blast_hit)
        self.logger.info("Found {} query sequences with 1 or more homologs".format(len(query_to_hits)))

        # For each query, collect the top 5 hits' IDs
        # Remove duplicates from top hit ID list
        top5_hit_ids = set()
        for query, hits in query_to_hits.items():
            if len(hits) < 5:
                for h in hits:
                    top5_hit_ids.add(h.subject_id)
            else:
                for i in range(5):
                    top5_hit_ids.add(hits[i].subject_id)
        self.logger.info("Extracting homologous sequences for {} top hits ..".format(len(top5_hit_ids)))

        # BLAST all the top query sequences against the DB itself
        db_top_homologs_tmp = os.path.join(output_dir, 'top_homologs_db.tmp')
        extract_homologs_and_context(top5_hit_ids, db_homologs_tmp, db_top_homologs_tmp)
        self.logger.info("BLASTing top homologues against DB ..")
        tophit_blast_output = os.path.join(output_dir, 'tophit_hits.tsv')
        if homology_search == 'diamond':
            diamond = Diamond(self.cpus)
            diamond.blastp(db_top_homologs_tmp, db_file, evalue, per_identity, per_aln_len, max_matches, tophit_blast_output, output_fmt='custom')
        else:
            blast.blastp(db_top_homologs_tmp, db_file, tophit_blast_output, evalue, max_matches, output_fmt='custom', task=homology_search)
        top_homologs = blast.identify_homologs(tophit_blast_output, evalue, per_identity, per_aln_len)

        # Identify bitscores of blast against the query sequences
        back_blast_output = os.path.join(output_dir, 'back_blast.tsv')
        if homology_search == 'diamond':
            raise Exception("Not implemented")
        else:
            query_db = os.path.join(output_dir, 'query_blast_db')
            logging.info("Creating BLAST DB for query sequences ..")
            # Currently a bug in biolib: https://github.com/dparks1134/biolib/pull/3
            os.system(blast.create_blastp_db_cmd(query_proteins, query_db))
            logging.info("Querying top homologues against query sequences")
            blast.blastp(db_top_homologs_tmp, query_db back_blast_output, evalue, max_matches, output_fmt='custom', task=homology_search)

        # For each query:
        # Find the sequences that are in all lists
        # Output a FASTA file for that group
        # If there are >= 4 sequences, then add it to the list for tree building

        # TODO: In parallel MSA and tree making?
        # Run MSA workflow
        # Run tree making software
        # Make Arb DB



        # concatenate homologs with initial query genes
        homolog_ouput_tmp = os.path.join(output_dir, 'homologs.faa.tmp')
        if False:#custom_homologs:
            # custom_db_homologs_tmp = os.path.join(output_dir, 'custom_homologs_db.tmp')
            # custom_gene_precontext, custom_gene_postcontext = self.extract_homologs_and_context(custom_homologs.keys(), custom_db_file, custom_db_homologs_tmp)
            # gene_precontext.update(custom_gene_precontext)
            # gene_postcontext.update(custom_gene_postcontext)
            # homologs.update(custom_homologs)
            # concatenate_files([query_proteins, db_homologs_tmp, custom_db_homologs_tmp], homolog_ouput_tmp)
            # os.remove(custom_db_homologs_tmp)
            pass
        else:
            concatenate_files([query_proteins, db_homologs_tmp], homolog_ouput_tmp)

        os.remove(db_homologs_tmp)

        # remove stop codons
        homolog_ouput = os.path.join(output_dir, 'homologs.faa')
        remove_stop_codons(homolog_ouput_tmp, homolog_ouput)
        os.remove(homolog_ouput_tmp)

        # infer multiple sequence alignment
        msa = MsaWorkflow(self.cpus)
        trimmed_msa_output = msa.run(homolog_ouput,
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
                                skip_rooting,
                                output_dir)

        # create tax2tree consensus map and decorate tree
        self.logger.info('Decorating internal tree nodes with tax2tree.')
        output_taxonomy_file = os.path.join(output_dir, 'taxonomy.tsv')
        fout = open(output_taxonomy_file, 'w')
        for homolog_id in homologs.keys():
            genome_id = homolog_id.split('~')[0]
            t = taxonomy.get(genome_id, None)
            if t:
                fout.write(homolog_id + '\t' + ';'.join(t) + '\n')
        fout.close()

        t2t_tree = os.path.join(output_dir, 'homologs.tax2tree.tree')
        cmd = 't2t decorate -m %s -t %s -o %s' % (output_taxonomy_file, tree_output, t2t_tree)
        os.system(cmd)

        # create tree with leaf nodes given as genome accessions
        tree = dendropy.Tree.get_from_path(t2t_tree,
                                            schema='newick',
                                            rooting='force-rooted',
                                            preserve_underscores=True)

        for leaf in tree.leaf_node_iter():
            leaf.taxon.label = leaf.taxon.label.split('~')[0]

        genome_tree = os.path.join(output_dir, 'homologs.tax2tree.genome_accessions.tree')
        tree.write_to_path(genome_tree,
                            schema='newick',
                            suppress_rooting=True,
                            unquoted_underscores=True)

        # setup metadata for ARB file
        src_dir = os.path.dirname(os.path.realpath(__file__))
        version_file = open(os.path.join(src_dir, 'VERSION'))

        metadata = {}
        metadata['genetreetk_version'] = version_file.read().strip()
        metadata['genetreetk_query_proteins'] = query_proteins
        metadata['genetreetk_db_file'] = db_file
        metadata['genetreetk_taxonomy_file'] = taxonomy_file
        metadata['genetreetk_blast_evalue'] = str(evalue)
        metadata['genetreetk_blast_per_identity'] = str(per_identity)
        metadata['genetreetk_blast_per_aln_len'] = str(per_aln_len)
        metadata['genetreetk_blast_max_matches'] = str(max_matches)
        metadata['genetreetk_homology_search'] = homology_search

        metadata['genetreetk_msa_min_per_taxa'] = str(min_per_taxa)
        metadata['genetreetk_msa_consensus'] = str(consensus)
        metadata['genetreetk_msa_min_per_bp'] = str(min_per_bp)
        metadata['genetreetk_msa_program'] = msa_program

        metadata['genetreetk_tree_program'] = tree_program
        metadata['genetreetk_tree_prot_model'] = prot_model

        # create ARB metadata file
        self.logger.info('Creating ARB metadata file.')
        arb_metadata_file = os.path.join(output_dir, 'arb.metadata.txt')
        create_arb_metadata(homologs, trimmed_msa_output, taxonomy,
                                 metadata,
                                 gene_precontext, gene_postcontext,
                                 arb_metadata_file)
