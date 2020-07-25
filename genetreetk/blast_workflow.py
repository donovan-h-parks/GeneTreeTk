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


class BlastWorkflow():
    """Blast-based workflow for building a gene tree."""

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

    def run(self, query_proteins,
            db_file, custom_db_file,
            taxonomy_file, custom_taxonomy_file,
            evalue, per_identity, per_aln_len, max_matches, homology_search,
            min_per_taxa, consensus, min_per_bp, use_trimAl, restrict_taxon,
            msa_program, tree_program, prot_model, skip_rooting,
            output_dir):
        """Infer a gene tree for homologs genes identified by blast.

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
        custom_taxonomy_file : str
            Taxonomic assignment of genomes in custom database.
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

        # validate query sequence names for use with GeneTreeTk
        validate_seq_ids(query_proteins)

        # read taxonomy file
        self.logger.info('Reading taxonomy file.')
        taxonomy = Taxonomy().read(taxonomy_file)

        if custom_taxonomy_file:
            custom_taxonomy = Taxonomy().read(custom_taxonomy_file)
            taxonomy.update(custom_taxonomy)
            
        # report distribution of query genes
        mean_len, max_len, min_len, p10, p50, p90 = self._gene_distribution(query_proteins)
        self.logger.info('Query gene lengths: min, mean, max = %d, %.1f, %d | p10, p50, p90 = %.1f, %.1f, %.1f' % (
                                                                                        min_len, mean_len, max_len, 
                                                                                        p10, p50, p90))

        # identify homologs using BLASTP
        self.logger.info('Identifying homologs using %s.' % homology_search)
        blast = Blast(self.cpus)
        blast_output = os.path.join(output_dir, 'reference_hits.tsv')
        if homology_search == 'diamond':
            diamond = Diamond(self.cpus)
            diamond.blastp(query_proteins, db_file, evalue, per_identity, per_aln_len, max_matches, blast_output, output_fmt='custom')
        else:
            blast.blastp(query_proteins, db_file, blast_output, evalue, max_matches, output_fmt='custom', task=homology_search)
        homologs = blast.identify_homologs(blast_output, evalue, per_identity, per_aln_len)
        self.logger.info('Identified %d homologs in reference database.' % len(homologs))

        custom_homologs = None
        if custom_db_file:
            custom_blast_output = os.path.join(output_dir, 'custom_hits.tsv')
            if homology_search == 'diamond':
                diamond = Diamond(self.cpus)
                diamond.blastp(query_proteins, custom_db_file, evalue, per_identity, per_aln_len, max_matches, custom_blast_output, output_fmt='custom')
            else:
                blast.blastp(query_proteins, custom_db_file, custom_blast_output, evalue, max_matches, output_fmt='custom', task=homology_search)
            custom_homologs = blast.identify_homologs(custom_blast_output, evalue, per_identity, per_aln_len)
            self.logger.info('Identified %d homologs in custom database.' % len(custom_homologs))
            
        # restrict homologs to specific taxonomic group
        if restrict_taxon:
            self.logger.info('Restricting homologs to %s.' % restrict_taxon)
            restricted_homologs = {}
            for query_id, hit in homologs.items():
                genome_id = hit.subject_id.split('~')[0]
                if restrict_taxon in taxonomy[genome_id]:
                    restricted_homologs[query_id] = hit

            self.logger.info('%d of %d homologs in reference database are from the specified group.' % (len(restricted_homologs), len(homologs)))
            homologs = restricted_homologs

        if len(homologs) == 0:
            self.logger.error('Too few homologs were identified. Gene tree cannot be inferred.')
            sys.exit()

        # extract homologs
        self.logger.info('Extracting homologs and determining local gene context.')
        db_homologs_tmp = os.path.join(output_dir, 'homologs_db.tmp')
        gene_precontext, gene_postcontext = extract_homologs_and_context(homologs.keys(), db_file, db_homologs_tmp)

        # report gene length distribution of homologs
        mean_len, max_len, min_len, p10, p50, p90 = gene_distribution(db_homologs_tmp)
        self.logger.info('Homolog gene lengths: min, mean, max = %d, %.1f, %d | p10, p50, p90 = %.1f, %.1f, %.1f' % (
                                                                                        min_len, mean_len, max_len, 
                                                                                        p10, p50, p90))
        
        # concatenate homologs with initial query genes
        homolog_ouput_tmp = os.path.join(output_dir, 'homologs.faa.tmp')
        if custom_homologs:
            custom_db_homologs_tmp = os.path.join(output_dir, 'custom_homologs_db.tmp')
            custom_gene_precontext, custom_gene_postcontext = self.extract_homologs_and_context(custom_homologs.keys(), custom_db_file, custom_db_homologs_tmp)
            gene_precontext.update(custom_gene_precontext)
            gene_postcontext.update(custom_gene_postcontext)
            homologs.update(custom_homologs)
            concatenate_files([query_proteins, db_homologs_tmp, custom_db_homologs_tmp], homolog_ouput_tmp)
            os.remove(custom_db_homologs_tmp)
        else:
            concatenate_files([query_proteins, db_homologs_tmp], homolog_ouput_tmp)

        os.remove(db_homologs_tmp)
        
        # remove stop codons
        homolog_ouput = os.path.join(output_dir, 'homologs.faa')
        self._remove_stop_codons(homolog_ouput_tmp, homolog_ouput)        
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
        self.create_arb_metadata(homologs, trimmed_msa_output, taxonomy,
                                 metadata,
                                 gene_precontext, gene_postcontext,
                                 arb_metadata_file)
