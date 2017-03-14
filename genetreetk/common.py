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
