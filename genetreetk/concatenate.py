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
import math
import logging
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.taxonomy import Taxonomy
from biolib.external.fasttree import FastTree
from biolib.external.raxml import RAxML
from biolib.external.execute import check_dependencies

from genetreetk.arb_parser import ArbParser

import dendropy


class Concatenate():
    """Concate MSA for genes and infer concatenated gene tree."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use during homology search.
        """

        check_dependencies(['FastTreeMP', 
                            'raxmlHPC-PTHREADS-SSE3', 
                            't2t', 
                            'seqmagick'])

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def create_arb_metadata(self, 
                            msa_output, 
                            taxonomy,
                            metadata,
                            output_file):
        """Create metadata file suitable for import into ARB.

        Parameters
        ----------
        msa_output : str
            Fasta file with aligned homologs.
        taxonomy : d[genome_id] -> list of taxa
            Taxonomic information for genomes.
        metadata : d[key] - string
            Additional metadata to write to ARB file.
        output_file : str
            File to write metadata information.
        """

        arb_metadata_list = []
        for seq_id, seq in seq_io.read_seq(msa_output):
            arb_metadata = {}
            arb_metadata['db_name'] = seq_id
            arb_metadata['genome_id'] = seq_id
            arb_metadata['gtdb_tax_string'] = ';'.join(taxonomy.get(seq_id, ''))
            arb_metadata['aligned_seq'] = seq

            for k, v in metadata.iteritems():
                arb_metadata[k] = v

            arb_metadata_list.append(arb_metadata)

        fout = open(output_file, 'w')
        arb_parser = ArbParser()
        arb_parser.write(arb_metadata_list, fout)
        fout.close()
        
    def _split_ids(self, label, split_chars):
        """Split leaf label into taxon and gene identifiers."""
        
        taxon_id = None
        gene_id = None
        for ch in split_chars:
            if ch in label:
                taxon_id, gene_id = label.split(ch, 1)
                break
                    
        return taxon_id, gene_id
                
    def run(self, 
            gene_dirs,
            min_per_gene,
            min_per_bps,
            tree_program,
            prot_model,
            split_chars,
            output_dir):
        """Infer concatenated gene tree.

        Parameters
        ----------
        gene_dirs : list
            GeneTreeTk output directories with information for individual genes.
        min_per_gene : float
            Minimum percentage of genes required to retain taxa.
        min_per_bps : float
            Minimum percentage of base pairs required to retain taxa.
        tree_program : str
            Program to use for tree inference ['fasttree', 'raxml'].
        prot_model : str
            Protein substitution model for tree inference ['WAG', 'LG', 'AUTO'].
        output_dir : str
            Directory to store results.
        """

        # read MSA files
        concat = defaultdict(lambda: defaultdict(list))
        msa_length = 0
        gene_lengths = {}
        for gene_dir in gene_dirs:
            homologs = os.path.join(gene_dir, 'homologs.trimmed.aligned.faa')

            for seq_id, seq in seq_io.read_seq(homologs):
                taxon_id, gene_id = self._split_ids(seq_id, split_chars)
                if not taxon_id:
                    self.logger.error('Failed to split identifier: %s' % seq_id)
                    sys.exit(-1)
  
                concat[taxon_id][gene_dir].append(seq)
                
            msa_length += len(seq)
            gene_lengths[gene_dir] = len(seq)
                
        # filter taxon
        mc_filter = set()
        min_per_gene_filter = set()
        min_per_bps_filter = set()
        for taxon_id in concat:
            # check if multiple copy
            missing = 0
            taxon_msa_len = 0
            for gene_id in gene_dirs:
                if gene_id not in concat[taxon_id]:
                    missing += 1
                    continue
                    
                if len(concat[taxon_id][gene_id]) > 1:
                    mc_filter.add(taxon_id)
                    break
 
                taxon_msa_len += len(concat[taxon_id][gene_id][0])
            
            if taxon_id not in mc_filter:
                if missing > len(gene_dirs)*(1.0-float(min_per_gene)/100.0):
                    min_per_gene_filter.add(taxon_id)
                elif taxon_msa_len < msa_length*float(min_per_bps)/100.0:
                    min_per_bps_filter.add(taxon_id)
                    
        min_req_genes = math.ceil(len(gene_dirs) * float(min_per_gene)/100.0)
            
        filtered_taxa = mc_filter.union(min_per_gene_filter).union(min_per_bps_filter)
        remaining_taxa = set(concat) - filtered_taxa
        self.logger.info('No. genes: %d' % len(gene_dirs))
        self.logger.info('No. taxa across all genes: %d' % len(concat))
        self.logger.info('Total filtered taxa: %d' % len(filtered_taxa))
        self.logger.info('  Due to multi-copy genes: %d' % len(mc_filter))
        self.logger.info('  Due to having <%d of the genes: %d' % (min_req_genes, len(min_per_gene_filter)))
        self.logger.info('  Due to an insufficient number of base pairs: %d' % len(min_per_bps_filter))
        self.logger.info('Remaining taxa: %d' % len(remaining_taxa))
        self.logger.info('Length of concatenated MSA: %d' % msa_length)

        # create the multiple sequences alignment
        msa_file = os.path.join(output_dir, 'concatenated.faa')
        fout = open(msa_file, 'w')
        for taxon_id in remaining_taxa:
            msa = ''
            for gene_id in gene_dirs:
                if gene_id not in concat[taxon_id]:
                    msa += '-'*gene_lengths[gene_id]
                else:
                    msa += concat[taxon_id][gene_id][0]
        
            fout.write('>%s\n' % taxon_id)
            fout.write('%s\n' % msa)
        fout.close()
        
        # read all taxonomy files 
        # (assumes taxonomy is the same for taxa across all genes)
        taxonomy = {}
        for gene_id in gene_dirs:
            taxonomy_file = os.path.join(gene_id, 'taxonomy.tsv')
            t = Taxonomy().read(taxonomy_file)
            for label, taxa_str in t.iteritems():
                taxon_id, gene_id = self._split_ids(label, split_chars)
                taxonomy[taxon_id] = taxa_str
        
        # create taxonomy file for retained taxa
        self.logger.info('Creating taxonomy file for retained taxa.') 
        output_taxonomy_file = os.path.join(output_dir, 'taxonomy.tsv')
        fout = open(output_taxonomy_file, 'w')
        for taxon_id in remaining_taxa:
            if taxon_id in taxonomy: # query genomes will generally be missing
                fout.write('%s\t%s\n' % (taxon_id, ';'.join(taxonomy[taxon_id])))
        fout.close()
        
        # infer tree
        if tree_program == 'fasttree':
            self.logger.info('Inferring gene tree with FastTree using %s+GAMMA.' % prot_model)
            fasttree = FastTree(multithreaded=(self.cpus > 1))

            tree_unrooted_output = os.path.join(output_dir, 'concatenated.unrooted.tree')
            tree_log = os.path.join(output_dir, 'concatenated.tree.log')
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
        self.logger.info('Rooting tree at midpoint.')
        tree = dendropy.Tree.get_from_path(tree_unrooted_output, schema='newick', rooting="force-rooted", preserve_underscores=True)
        if len(remaining_taxa) > 2:
            tree.reroot_at_midpoint(update_bipartitions=False)
        tree_output = os.path.join(output_dir, 'concatenated.rooted.tree')
        tree.write_to_path(tree_output, schema='newick', suppress_rooting=True, unquoted_underscores=True)

        # create tax2tree consensus map and decorate tree
        t2t_tree = os.path.join(output_dir, 'concatenated.tax2tree.tree')
        cmd = 't2t decorate -m %s -t %s -o %s' % (output_taxonomy_file, tree_output, t2t_tree)
        os.system(cmd)

        # setup metadata for ARB file
        src_dir = os.path.dirname(os.path.realpath(__file__))
        version_file = open(os.path.join(src_dir, 'VERSION'))

        metadata = {}
        metadata['genetreetk_version'] = version_file.read().strip()
        
        metadata['genetreetk_tree_program'] = tree_program
        metadata['genetreetk_tree_prot_model'] = prot_model

        # create ARB metadata file
        self.logger.info('Creating ARB metadata file.')
        arb_metadata_file = os.path.join(output_dir, 'arb.metadata.txt')
        self.create_arb_metadata(msa_file, 
                                    taxonomy,
                                    metadata,
                                    arb_metadata_file)
