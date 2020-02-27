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
import biolib.seq_tk as seq_tk
from biolib.external.muscle import Muscle
from biolib.external.mafft import Mafft
from biolib.external.execute import check_dependencies

class MsaWorkflow():
    """Workflow for creating multiple sequence alignment."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use during homology search.
        """

        check_dependencies(['mafft',
                            'muscle', 
                            'seqmagick',
                            'trimal'])

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def run(self, homolog_file,
                min_per_taxa, 
                consensus, 
                min_per_bp, 
                use_trimAl, 
                msa_program,
                output_dir):
        """Create multiple sequence alignment.

        Parameters
        ----------
        homolog_file : str
            File containing sequences to align
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
        output_dir : str
            Directory to store results.
        """

        # infer multiple sequence alignment
        self.logger.info('Inferring multiple sequence alignment with %s.' % msa_program)
        
        output_file = ntpath.basename(homolog_file)
        prefix = output_file[0:output_file.rfind('.')]
        suffix = output_file[output_file.rfind('.') + 1:]

        msa_output = os.path.join(output_dir,
                                    prefix + '.aligned.' + suffix)
        if msa_program == 'mafft':
            mafft = Mafft(self.cpus)
            msa_log = os.path.join(output_dir, 'mafft.log')
            mafft.run(homolog_file, msa_output, msa_log)
        elif msa_program == 'muscle':
            muscle = Muscle()
            msa_log = os.path.join(output_dir, 'muscle.log')
            muscle.run(homolog_file, msa_output, msa_log)

        # trim multiple sequence alignment
        trimmed_msa_output = os.path.join(output_dir,
                                            prefix + '.trimmed.aligned.' + suffix)
        if use_trimAl:
            self.logger.info('Using trimAl to filter poorly represented columns from alignment.')
            
            # convert MSA to relaxed phylip format
            phylip_msa_output = msa_output.replace('.faa', '.phyx')
            cmd = 'seqmagick convert %s %s' % (msa_output, phylip_msa_output)
            os.system(cmd)
            
            tmp_output = os.path.join(output_dir, 'tmp.faa')
            cmd = 'trimal -in %s -out %s -automated1 -fasta' % (phylip_msa_output, tmp_output)
            os.system(cmd)
            
            cmd = 'trimal -in %s -out %s -resoverlap 0.75 -seqoverlap %f' % (tmp_output, trimmed_msa_output, min_per_bp)
            os.system(cmd)
            
            seqs = seq_io.read_fasta(msa_output)
            tmp_seqs = seq_io.read_fasta(tmp_output)
            trimmed_seqs = seq_io.read_fasta(trimmed_msa_output)
            self.logger.info('Trimmed alignment from %d to %d AA.' % (len(list(seqs.values())[0]), len(list(trimmed_seqs.values())[0])))
            self.logger.info('%d of %d taxa were deemed to be too short and removed.' % (len(tmp_seqs)-len(trimmed_seqs), len(seqs)))
            os.remove(tmp_output)
        else:
            self.logger.info('Trimming poorly represented columns from alignment.')
            seqs = seq_io.read_fasta(msa_output, keep_annotation=True)
            trimmed_seqs, pruned_seqs, min_taxa_filtered, consensus_filtered = seq_tk.trim_seqs(seqs, 
                                                                                                min_per_taxa / 100.0, 
                                                                                                consensus / 100.0, 
                                                                                                min_per_bp / 100.0)
            
            self.logger.info('Trimmed alignment from %d to %d AA (%d by minimum taxa percent, %d by consensus).' % (len(list(seqs.values())[0]), 
                                                                                                                    len(list(trimmed_seqs.values())[0]),
                                                                                                                    min_taxa_filtered, 
                                                                                                                    consensus_filtered))
            self.logger.info('%d of %d taxa were deemed to be too short and removed.' % (len(pruned_seqs), len(seqs)))

            if len(pruned_seqs) > 0:
                prune_seqs_out = os.path.join(output_dir, 'filtered_seqs.too_short.txt')
                self.logger.info('Pruned sequences written to %s.' % prune_seqs_out)
                seq_io.write_fasta(pruned_seqs, prune_seqs_out)

            if len(pruned_seqs) == len(seqs):
                self.logger.error('Too many sequences were pruned. Gene tree cannot be inferred.')
                sys.exit()

            seq_io.write_fasta(trimmed_seqs, trimmed_msa_output)
            
        return trimmed_msa_output
