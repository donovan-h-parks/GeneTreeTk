#!/usr/bin/env python

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

__prog_name__ = 'make_database'
__prog_desc__ = 'create dereplicated database of genes from reference genomes'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import string
import random
import tempfile
import shutil
import logging
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.common import make_sure_path_exists, remove_extension
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies
from biolib.misc.time_keeper import TimeKeeper

from genetreetk.common import validate_seq_ids


class CreateDatabase(object):
    """Make a dereplicated database of genes.

    Dereplication is done between genes within a named taxonomic
    group (e.g., genomes in the same genus) and is based on the
    average amino acid identity (AAI) between genes. Groups with large
    numbers of taxa can take an excessive amount of time to
    process so are subsampled to a specific number of taxa.
    Subsampling is done in a manor which aims to retain
    phylogenetic diversity and thus helps ensures a good
    distribution of genes within the group. Care is taken
    to ensure type strains are retained during dereplication.
    """

    def __init__(self, cpus):
        """Initialize.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        
        self.logger = logging.getLogger('timestamp')

        check_dependencies(['comparem', 'diamond', 'makeblastdb'])

        self.underclassified = 'underclassified'

        self.rank_prefixes = Taxonomy.rank_prefixes
        self.rank_index = Taxonomy.rank_index
        self.rank_labels = Taxonomy.rank_labels

        self.cpus = cpus

    def read_type_strain(self, type_strain_file):
        """Read type strain file.

        The type strain file should have the following format:
            <genome_id>\t<genome_name>

        Parameters
        ----------
        type_strain_file : str
            File specifying type strains.

        Returns
        -------
        set
            Set of all genome ids specified as type strains.
        """

        type_strains = set()
        for line in open(type_strain_file):
            line_split = line.split('\t')
            type_strains.add(line_split[0])

        return type_strains

    def select_taxa(self, genome_list, taxonomy, type_strains, max_taxa):
        """Select subset of genomes with a good distribution across named groups.

        Groups genomes into named groups and subsamples evenly across
        these groups. Ideally, genomes would be grouped into species, but
        some genomes may not have a species identifier. Such genomes are
        assigned to the most specific named group possible. Any genome
        marked as a type strain will be retained.

        Parameters
        ----------
        genome_list : iterable of genome ids
            Genomes to subsample.
        taxonomy : d[genome_id] -> [domain, ..., species]
            Taxonomy of each genome.
        type_strains : iterable
            Genome identifiers of type strains.
        max_taxa : int
            Number of genomes to retain.

        Returns
        -------
        iterable
            Subsampled list of genomes.
        """

        if len(genome_list) <= max_taxa:
            return genome_list

        reduced_genome_list = []

        # group genomes into the most specific named groups possible
        groups = defaultdict(set)
        for genome_id in genome_list:
            # add in type strains regardless of taxonomy
            if genome_id in type_strains:
                reduced_genome_list.append(genome_id)
                continue

            # get first classified rank
            for rank_index in xrange(self.rank_index['s__'], -1, -1):
                taxa = taxonomy[genome_id][rank_index]
                if taxa != self.rank_prefixes[rank_index]:
                    break

            groups[taxa].add(genome_id)

        # sample genomes from each named group
        while len(reduced_genome_list) < max_taxa:
            genomes_to_select = max_taxa - len(reduced_genome_list)
            genomes_per_group = max(genomes_to_select / len(groups), 1)
            for taxa, genome_ids in groups.iteritems():
                selected_genomes = random.sample(genome_ids, min(len(genome_ids), genomes_per_group))
                groups[taxa] = genome_ids.difference(selected_genomes)

                reduced_genome_list.extend(selected_genomes)

                if len(reduced_genome_list) == max_taxa:
                    break  # special case where we are adding single genomes from each group

        return reduced_genome_list

    def write_gene_file(self, gene_out, gene_dir, genome_list, taxonomy, genes_to_ignore):
        """Write genes to output stream.

        Parameters
        ----------
        gene_out : stream
            Output stream.
        gene_dir : str
            Directory containing called genes in amino acid space.
        genome_list : iterable
            Genomes to process.
        genes_to_ignore : set
            Genes which should not be written to file.
        """

        genes_kept = 0
        for genome_id in genome_list:
            genome_gene_file = os.path.join(gene_dir, genome_id + '.faa')
            if not os.path.exists(genome_gene_file):
                print '[WARNING] Missing gene file for genome %s.' % genome_gene_file
                continue

            if os.stat(genome_gene_file).st_size == 0:
                print '[WARNING] Gene file is empty for genome %s.' % genome_gene_file
                continue

            for gene_id, seq, annotation in seq_io.read_fasta_seq(genome_gene_file, keep_annotation=True):
                if gene_id in genes_to_ignore:
                    continue

                gene_out.write('>' + gene_id + ' ' + annotation + '\n')
                gene_out.write(seq + '\n')
                genes_kept += 1

        return genes_kept

    def reformat_gene_id_to_scaffold_id(self, gene_file, gff_file, taxonomy, output_file):
        """Reformat gene ids to format which explicitly gives scaffold names.

        <genome_id>~<scaffold_id>_<gene_#> [gtdb_taxonomy] [NCBI organism name] [annotation]

        Parameters
        ----------
        gene_file : str
            Gene file for genome.
        gff_file : str
            General feature file (GFF) for genome.
        output_file : float
            File to contain modified gene fasta file.
        """

        # determine source scaffold for each gene
        gene_id_to_scaffold_id = {}
        gene_number = defaultdict(int)
        for line in open(gff_file):
            if line.startswith('##FASTA'):
                # start of FASTA section with individual sequences
                break

            if line[0] == '#':
                continue

            line_split = line.split('\t')
            scaffold_id = line_split[0]
            info = line_split[8]
            if info != '':  # this will be empty for non-protein coding genes
                gene_id = info.split(';')[0].replace('ID=', '')

                gene_number[scaffold_id] += 1
                gene_id_to_scaffold_id[gene_id] = scaffold_id + '_' + str(gene_number[scaffold_id])

        # write out gene file with modified identifiers
        fout = open(output_file, 'w')
        for gene_id, seq, annotation in seq_io.read_fasta_seq(gene_file, keep_annotation=True):
            genome_id = remove_extension(gene_file)
            fout.write('>%s [%s] [%s] [%s]\n' % (gene_id_to_scaffold_id[gene_id],
                                                    ';'.join(taxonomy.get(genome_id, ['none'])),
                                                    'none',
                                                    annotation))
            fout.write(seq + '\n')
        fout.close()

    def amend_gene_identifies(self, gene_dir, output_dir):
        """Modify gene ids to include source genome id.

        The following format is used:
          <genome_id>~<gene_id>

        Parameters
        ----------
        gene_dir : str
            Directory with fasta files containing protein sequences.
        output_dir : float
            Directory to contain modified fasta files.
        """

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for f in os.listdir(gene_dir):
            gf = os.path.join(gene_dir, f)
            genome_id = remove_extension(gf)

            aa_file = os.path.join(output_dir, genome_id + '.faa')
            fout = open(aa_file, 'w')
            for seq_id, seq, annotation in seq_io.read_fasta_seq(gf, keep_annotation=True):
                fout.write('>%s~%s %s\n' % (genome_id, seq_id, annotation))
                if seq[-1] == '*':
                    seq = seq[0:-1]
                fout.write(seq + '\n')
            fout.close()

    def filter_aai(self, tmp_dir, gene_dir, ammended_gene_dir, per_identity, per_aln_len, cpus):
        """Filter genes with similar amino acid identity.

        Parameters
        ----------
        tmp_dir : str
            Temporary directory for storing results.
        gene_dir : str
            Directory with fasta files containing protein sequences.
        ammended_gene_dir : str
            Directory to store protein sequences with ammended gene ids.
        per_identity : float
            Percent identity for subsampling similar genes.
        per_aln_len : float
            Percent alignment length for subsampling similar genes.
        cpus : int
            Number of cpus to use.

        Returns
        -------
        genes_to_remove : set
            Unique identifiers of genes to filter.
        """

        rblast_dir = os.path.join(tmp_dir, 'rblast')
        os.system('comparem rblast -e 1e-10 -p %d -c %d %s %s' % (per_identity, cpus, gene_dir, rblast_dir))
        aai_dir = os.path.join(tmp_dir, 'aai')
        os.system('comparem aai -p %d -a %d -c %d %s %s' % (per_identity, per_aln_len, cpus, rblast_dir, aai_dir))

        # identify homologs to be filtered
        print ''
        print '  Identifying homologs to be filtered.'
        shared_genes_dir = os.path.join(aai_dir, 'shared_genes')
        files = os.listdir(shared_genes_dir)

        homologs = defaultdict(set)
        for f in files:
            with open(os.path.join(shared_genes_dir, f)) as fin:
                fin.readline()

                for line in fin:
                    line_split = line.split('\t')

                    gene_idA = line_split[0]
                    gene_idB = line_split[1]

                    homologs[gene_idA].add(gene_idB)
                    homologs[gene_idB].add(gene_idA)

        genes_to_remove = set()
        genes_to_keep = set()
        sorted_keys = sorted(homologs, key=lambda k: len(homologs[k]), reverse=True)
        for gene_id in sorted_keys:
            gene_set = homologs[gene_id]

            if len(gene_set.intersection(genes_to_keep)) > 0:
                genes_to_remove.update(gene_set - genes_to_keep)
                genes_to_remove.add(gene_id)
            else:
                genes_to_keep.add(gene_id)
                genes_to_remove.update(gene_set - genes_to_keep)

        # The CompareM call to rblast creates fasta files where gene ids are modified to
        # also contain genome ids. This is just a hack so to point to the directory with
        # these amended fasta files.
        os.system('ln -s %s %s' % (os.path.join(rblast_dir, 'genes'), ammended_gene_dir))

        return genes_to_remove

    def run(self,
                taxonomy_file, type_strains_file,
                genome_prot_dir, extension,
                max_taxa, rank,
                per_identity, per_aln_len,
                genomes_to_process, keep_all_genes,
                no_reformat_gene_ids,
                output_dir):
        """ Create dereplicate set of genes.

        Taxonomy file should have the following format:
            <genome_id>\t<taxonomy_str>

            where taxonomy_str is in GreenGenes format:
                d__Bacteria;p__Proteobacteria;...;s__Escherichia coli

        Type strain file should have the following format:
            <genome_id>\t<genome name>

        Parameters
        ----------
        taxonomy_file : str
            File indicating taxonomy string for all genomes of interest
        type_strains_file : str
            File indicating type strains.
        genome_prot_dir : str
            Directory containing amino acid genes for each genome.
        extension : str
            Extension of files with called genes.
        max_taxa : int
            Maximum taxa to retain in a named group.
        rank : int
            Taxonomic rank to perform dereplication (0 = domain, ..., 6 = species).
        per_identity : float
            Percent identity for subsampling similar genes.
        per_aln_len : float
            Percent alignment length for subsampling similar genes.
        genomes_to_process : str
            File with list of genomes to retain instead of performing taxon subsampling.
        keep_all_genes : boolean
            Flag indicating that no gene subsampling should be performed.
        no_reformat_gene_ids : boolean
            Flag indicating if gene ids should be reformatted to include scaffold names given by the GFF file.
        output_dir : str
            Desired output directory for storing results.
        """

        make_sure_path_exists(output_dir)
        self.logger.info('Dereplicating at the rank of %s.' % self.rank_labels[rank])

        # get taxonomy string for each genome
        taxonomy = {}
        if taxonomy_file:
            self.logger.info('Reading taxonomy file.')
            taxonomy = Taxonomy().read(taxonomy_file)
            self.logger.info('There are %d genomes with taxonomy strings.' % len(taxonomy))

        # get type strains; genomes which should never be dereplicated
        type_strains = set()
        if type_strains_file:
            self.logger.info('Reading type strain file.')
            type_strains = self.read_type_strain(type_strains_file)
            self.logger.info('There are %d type strains.' % len(type_strains))

        # get specific list of genomes to process
        genomes_to_retain = set()
        if genomes_to_process:
            self.logger.info('Reading genomes to retain.')
            for line in open(genomes_to_process):
                line_split = line.split()
                genomes_to_retain.add(line_split[0])
            self.logger.info('Retaining %d genomes.' % len(genomes_to_retain))
            
        # make sure extension filter starts with a '.'
        if not extension.startswith('.'):
            extension = '.' + extension

        # identify unique genes in each named group
        fout = open(os.path.join(output_dir, 'genomes_without_called_genes.tsv'), 'w')
        rank_genomes = defaultdict(list)
        genome_files = os.listdir(genome_prot_dir)
        underclassified_genomes = 0
        genomes_with_missing_data = 0
        for genome_file in genome_files:
            genome_id = remove_extension(genome_file, extension)

            if not genome_file.endswith(extension):
                continue

            if genomes_to_process and genome_id not in genomes_to_retain:
                continue

            genome_file = os.path.join(genome_prot_dir, genome_file)
            if not os.path.exists(genome_file):
                genomes_with_missing_data += 1
                fout.write(genome_id + '\t' + ';'.join(taxonomy[genome_id]) + '\n')
                continue

            t = taxonomy.get(genome_id, self.rank_prefixes)
            taxa = t[rank]
            if taxa[3:] == '':
                underclassified_genomes += 1
                rank_genomes[self.underclassified].append(genome_id)
            else:
                rank_genomes[taxa].append(genome_id)

            validate_seq_ids(genome_file)

        fout.close()

        total_genomes_to_process = sum([len(genome_list) for genome_list in rank_genomes.values()])
        if total_genomes_to_process == 0:
            self.logger.error('No genomes found in directory: %s. Check the --extension flag used to identify genomes.' % genome_prot_dir)
            sys.exit(-1)

        self.logger.info('Under-classified genomes automatically placed into the database: %d' % underclassified_genomes)
        self.logger.info('Genomes with missing sequence data: %d' % genomes_with_missing_data)
        self.logger.info('Total named groups: %d' % len(rank_genomes))
        self.logger.info('Total genomes to process: %d' % total_genomes_to_process)

        # process each named group
        gene_file = os.path.join(output_dir, 'custom_db.faa')
        gene_out = open(gene_file, 'w')

        taxonomy_out = open(os.path.join(output_dir, 'custom_taxonomy.tsv'), 'w')

        tmp_dir = tempfile.mkdtemp()
        total_genes_removed = 0
        total_genes_kept = 0
        total_genomes_kept = 0
        processed_genomes = 0
        for taxa, genome_list in rank_genomes.iteritems():
            processed_genomes += len(genome_list)

            print '-------------------------------------------------------------------------------'
            self.logger.info('Processing %s | Finished %d of %d (%.2f%%) genomes.' % (taxa, processed_genomes, total_genomes_to_process, processed_genomes * 100.0 / total_genomes_to_process))

            # create directory with selected genomes
            taxon_dir = os.path.join(tmp_dir, 'taxon')
            os.mkdir(taxon_dir)

            reduced_genome_list = genome_list
            if not genomes_to_process and taxa != self.underclassified:  # perform taxon subsampling
                reduced_genome_list = self.select_taxa(genome_list, taxonomy, type_strains, max_taxa)
            total_genomes_kept += len(reduced_genome_list)

            gene_dir = os.path.join(taxon_dir, 'genes')
            os.mkdir(gene_dir)
            for genome_id in reduced_genome_list:
                taxonomy_out.write(genome_id + '\t' + ';'.join(taxonomy.get(genome_id, self.rank_prefixes)) + '\n')

                genome_gene_file = os.path.join(genome_prot_dir, genome_id + extension)
                gff_file = os.path.join(genome_prot_dir, genome_id + '.gff')
                output_gene_file = os.path.join(gene_dir, genome_id + '.faa')
                if not no_reformat_gene_ids:
                    self.reformat_gene_id_to_scaffold_id(genome_gene_file, gff_file, taxonomy, output_gene_file)
                else:
                    os.system('cp %s %s' % (genome_gene_file, output_gene_file))

            # filter genes based on amino acid identity
            genes_to_remove = []
            amended_gene_dir = os.path.join(taxon_dir, 'amended_genes')
            if keep_all_genes or taxa == self.underclassified:
                # modify gene identifiers to include genome ids
                self.amend_gene_identifies(gene_dir, amended_gene_dir)
            else:
                # filter genes on AAI
                genes_to_remove = self.filter_aai(taxon_dir, gene_dir, amended_gene_dir, per_identity, per_aln_len, self.cpus)

            self.logger.info('Writing unique genes from genomes in %s.' % taxa)
            genes_kept = self.write_gene_file(gene_out, amended_gene_dir, reduced_genome_list, taxonomy, genes_to_remove)

            self.logger.info('Retain %d of %d taxa.' % (len(reduced_genome_list), len(genome_list)))
            self.logger.info('Genes to keep: %d' % genes_kept)
            self.logger.info('Genes removed: %d' % len(genes_to_remove))

            total_genes_kept += genes_kept
            total_genes_removed += len(genes_to_remove)

            shutil.rmtree(taxon_dir)

        taxonomy_out.close()
        gene_out.close()

        self.logger.info('Retain %d of %d (%.1f%%) genomes' % (total_genomes_kept, total_genomes_to_process, total_genomes_kept * 100.0 / (total_genomes_to_process)))
        self.logger.info('Total genes kept: %d' % total_genes_kept)
        self.logger.info('Total genes removed: %d (%.1f%%)' % (total_genes_removed, total_genes_removed * 100.0 / (total_genes_kept + total_genes_removed)))

        self.logger.info('Creating BLAST database.')
        os.system('makeblastdb -dbtype prot -in %s' % gene_file)

        shutil.rmtree(tmp_dir)
