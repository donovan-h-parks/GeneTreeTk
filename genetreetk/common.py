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


import sys
import logging

import biolib.seq_io as seq_io

from numpy import (percentile as np_percentile,
                    mean as np_mean)
import extern

from genetreetk.arb_parser import ArbParser


def validate_seq_ids(query_proteins):
    """Ensure all sequence identifiers contain only acceptable characters.

    Parameters
    ----------
    query_proteins : str
        Fasta file containing query proteins.
    """

    invalid_chars = set('()[],;=')
    for seq_id, _seq in seq_io.read_seq(query_proteins):
        if any((c in invalid_chars) for c in seq_id):
            logging.getLogger('no_timestamp').error('Invalid sequence header in file %s' % query_proteins)
            logging.getLogger('no_timestamp').error('Sequence contains an invalid character: %s' % seq_id)
            logging.getLogger('no_timestamp').error('Sequence identifiers must not contain the following characters: ' + ''.join(invalid_chars))
            sys.exit()

def extract_homologs_and_context(homologs, db_file, output_file):
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
            for seq_id, count in post_context_counter.items():
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
    gene_precontext = _filter_gene_context(gene_precontext)
    gene_postcontext = _filter_gene_context(gene_postcontext)

    return gene_precontext, gene_postcontext

def _filter_gene_context(gene_context):
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
    for seq_id, context in gene_context.items():
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

def create_arb_metadata(homologs, 
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

        for k, v in metadata.items():
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
    
def gene_distribution(seq_file):
    """Calculate length distribution of sequences."""
    
    gene_lens = []
    for _seq_id, seq in seq_io.read_seq(seq_file):
        gene_lens.append(len(seq))
        
    p10, p50, p90 = np_percentile(gene_lens, [10, 50, 90])
    
    return np_mean(gene_lens), max(gene_lens), min(gene_lens), p10, p50, p90
    
def remove_stop_codons(input_file, output_file):
    """Remove stop codons at end of sequences."""
    
    fout = open(output_file, 'w')
    for seq_id, seq, annotation in seq_io.read_seq(input_file, keep_annotation=True):
        fout.write('>%s %s\n' % (seq_id, annotation))
        
        if seq[-1] == '*':
            seq = seq[0:-1]
        fout.write('%s\n' % seq)
    fout.close()

def extract_sequences(self, sequence_ids, db_file, output_file):
    extern.run("mfqe --fasta-read-name-lists \/dev/stdin --input-fasta '{}' --output-fasta-files '{}' --output-uncompressed".format(
        db_file, output_file),
        stdin="\n".join(sequence_ids))