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

from genetreetk.arb_parser import ArbParser
from genetreetk.common import validate_seq_ids
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

    def extract_homologs_and_context(self, homologs, db_file, output_file):
        """Extract homologs sequences from database file, and local gene context.

        This function extract sequences information for each
        homolog and writes this to file for downstream processing.
        In addition, it determines the local gene context for each
        gene. Specifically, it saves the annotations for the
        3 genes prior to and after a given gene.

        This function assumes the database is sorted according
        to the order genes are identified on each contig.

        Parameters
        ----------
        homologs : iterable
            Unique identifiers of sequences to extract
        db_file : str
            Fasta file with sequences.
        output_file : str
            File to write homologs.

        Returns
        -------
        dict
            d[seq_id] -> list of annotations for pre-context genes
        dict
            d[seq_id] -> list of annotations for post-context genes
        """

        gene_precontext = {}
        gene_postcontext = {}

        if len(homologs) == 0:
            return gene_precontext, gene_postcontext

        if type(homologs) is not set:
            homologs = set(homologs)

        fout = open(output_file, 'w')
        local_context = [('unknown~unknown_x', None)] * 3
        post_context_counter = {}
        for seq_id, seq, annotation in seq_io.read_fasta_seq(db_file, keep_annotation=True):
            if seq_id in homologs:
                fout.write('>' + seq_id + ' ' + annotation + '\n')
                fout.write(seq + '\n')

                gene_precontext[seq_id] = list(local_context)
                post_context_counter[seq_id] = 3

            # record 3 precontext genes
            local_context[0] = local_context[1]
            local_context[1] = local_context[2]
            local_context[2] = (seq_id, annotation)

            # record 3 postcontext genes
            if len(post_context_counter):
                key_to_remove = None
                for seq_id, count in post_context_counter.iteritems():
                    count -= 1
                    if count == -1:
                        gene_postcontext[seq_id] = list(local_context)
                        key_to_remove = seq_id
                    else:
                        post_context_counter[seq_id] = count

                if key_to_remove:
                    post_context_counter.pop(key_to_remove)

        fout.close()

        # filter gene context to contain only genes on the same scaffold
        gene_precontext = self._filter_gene_context(gene_precontext)
        gene_postcontext = self._filter_gene_context(gene_postcontext)

        return gene_precontext, gene_postcontext

    def _filter_gene_context(self, gene_context):
        """Filter gene context to contain only genes on the same scaffold.

        This function assumes sequence identifies have the following format:
            <genome_id>~<scaffold_id>_<gene_#> [gtdb_taxonomy] [NCBI organism name] [annotation]

        Parameters
        ----------
        gene_context : d[seq_id] -> [(seq_id, annotation), ..., (seq_id, annotation)]
            Gene context.

        Returns
        -------
        dict: d[seq_id] -> [annotation, ..., annotation]
            Filtered to contain only annotations from the same scaffold.
        """

        filtered_gene_context = {}
        for seq_id, context in gene_context.iteritems():
            _genome_id, gene_id = seq_id.split('~')
            scaffold_id = gene_id[0:gene_id.rfind('_')]

            filtered_context = []
            for local_seq_id, annotation in context:
                _local_genome_id, local_gene_id = local_seq_id.split('~')
                local_scaffold_id = local_gene_id[0:local_gene_id.rfind('_')]

                # strip organism name and IMG gene id
                annotation = annotation[0:annotation.rfind('[')]
                annotation = annotation[0:annotation.rfind('[')].strip()

                if scaffold_id == local_scaffold_id:
                    filtered_context.append(annotation)

            filtered_gene_context[seq_id] = filtered_context

        return filtered_gene_context

    def create_arb_metadata(self,
                            homologs, 
                            msa_output, 
                            taxonomy,
                            metadata,
                            gene_precontext, 
                            gene_postcontext,
                            output_file):
        """Create metadata file suitable for import into ARB.

        Parameters
        ----------
        homologs : d[seq_id] -> namedtuple of BlastHit information
            BLAST results for identified homologs.
        msa_output : str
            Fasta file with aligned homologs.
        taxonomy : d[genome_id] -> list of taxa
            Taxonomic information for genomes.
        metadata : d[key] - string
            Additional metadata to write to ARB file.
        gene_precontext : d[seq_id] -> list of annotations for pre-context genes
            Annotation for genes preceding a gene.
        gene_postcontext: d[seq_id] -> list of annotations for post-context genes
            Annotation for genes following a gene.
        output_file : str
            File to write metadata information.
        """

        arb_metadata_list = []
        for seq_id, seq, annotation in seq_io.read_seq(msa_output, keep_annotation=True):
            if '~' in seq_id:
                genome_id, scaffold_gene_id = seq_id.split('~')
            else:
                scaffold_gene_id = seq_id
                genome_id = ''

            arb_metadata = {}
            arb_metadata['db_name'] = seq_id
            arb_metadata['genome_id'] = genome_id
            arb_metadata['scaffold_id'] = scaffold_gene_id[0:scaffold_gene_id.rfind('_')]
            arb_metadata['scaffold_gene_id'] = scaffold_gene_id
            arb_metadata['gtdb_tax_string'] = ';'.join(taxonomy.get(genome_id, ''))
            arb_metadata['aligned_seq'] = seq

            for k, v in metadata.iteritems():
                arb_metadata[k] = v

            arb_metadata['gene_precontext'] = ' -> '.join(gene_precontext.get(seq_id, []))
            arb_metadata['gene_postcontext'] = ' <- '.join(gene_postcontext.get(seq_id, []))

            hit_info = homologs.get(seq_id, None)
            if hit_info:
                arb_metadata['blast_evalue'] = '%.1g' % hit_info.evalue
                arb_metadata['blast_bitscore'] = '%.1f' % hit_info.bitscore
                arb_metadata['blast_perc_identity'] = '%.1f' % hit_info.perc_identity
                arb_metadata['blast_subject_perc_alignment_len'] = '%.1f' % hit_info.subject_perc_aln_len
                arb_metadata['blast_query_perc_alignment_len'] = '%.1f' % hit_info.query_perc_aln_len
                arb_metadata['blast_query_id'] = hit_info.query_id

            if annotation:
                annotation_split = annotation.split('] [')
                if len(annotation_split) == 3:
                    # assume format is [gtdb_taxonomy] [NCBI organism name] [annotation]
                    gtdb_taxonomy, organism_name, gene_annotation = annotation_split
                    gtdb_taxonomy = gtdb_taxonomy.replace('[', '')
                    gene_annotation = gene_annotation.replace(']', '')
                else:
                    # no idea what the format is, so just save the annotation
                    gene_annotation = annotation
                    organism_name = ''
                    gtdb_taxonomy = ''

                arb_metadata['gene_annotation'] = gene_annotation
                arb_metadata['organism'] = organism_name
                arb_metadata['full_name'] = organism_name

            arb_metadata_list.append(arb_metadata)

        fout = open(output_file, 'w')
        arb_parser = ArbParser()
        arb_parser.write(arb_metadata_list, fout)
        fout.close()
        
    def _gene_distribution(self, seq_file):
        """Calculate length distribution of sequences."""
        
        gene_lens = []
        for seq_id, seq in seq_io.read_seq(seq_file):
            gene_lens.append(len(seq))
            
        p10, p50, p90 = np_percentile(gene_lens, [10, 50, 90])
        
        return np_mean(gene_lens), max(gene_lens), min(gene_lens), p10, p50, p90
        
    def _remove_stop_codons(self, input_file, output_file):
        """Remove stop codons at end of sequences."""
        
        fout = open(output_file, 'w')
        for seq_id, seq, annotation in seq_io.read_seq(input_file, keep_annotation=True):
            fout.write('>%s %s\n' % (seq_id, annotation))
            
            if seq[-1] == '*':
                seq = seq[0:-1]
            fout.write('%s\n' % seq)
        fout.close()

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
            for query_id, hit in homologs.iteritems():
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
        gene_precontext, gene_postcontext = self.extract_homologs_and_context(homologs.keys(), db_file, db_homologs_tmp)

        # report gene length distribution of homologs
        mean_len, max_len, min_len, p10, p50, p90 = self._gene_distribution(db_homologs_tmp)
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
