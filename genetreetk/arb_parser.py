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
__copyright__ = "Copyright 2015"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = ""
__status__ = "Development"

import logging


class ArbParser:
    """Provides methods for reading and writing files for use with ARB.

    This class currently supports the GreenGenes format file:

    BEGIN
    db_name=A00000001
    organism=Korebacter versatilis Ellin345
    genome_tree_description=None
    prokMSA_id=A00000001
    owner=root
    genome_tree_tax_string=k__Bacteria; p__Acidobacteria; c__Acidobacteriia; o__Acidobacteriales; f__; g__; s__;
    greengenes_tax_string=k__Bacteria;p__Acidobacteria;c__;o__;f__Koribacteraceae;g__CandidatusKoribacter;s__CandidatusKoribacterversatilis
    blast_hits_16s=1144/1144 (100%) -- 157743 k__Bacteria; p__Acidobacteria; c__Acidobacteriia; o__Acidobacteriales; f__Koribacteraceae; g__Candidatus Koribacter; s__versatilis
    img_tax_string=k__Bacteria; p__Acidobacteria; c__unclassified; o__unclassified; f__unclassified; g__Candidatus Koribacter; s__versatilis;
    img_tax_tax2tree_corrected=
    checkm_completeness=1.0
    checkm_contamination=0.0
    core_list_status=public
    warning=
    aligned_seq=-KRTHKCGELRAADANKNVVLMGWVNRRRDLGGLIFIDLRDRTGITQIVFDNSSELQAKAGDLRSEYCIAVIGTVAKREANTVNKNLPTGEIEVVAKEMRLFNDSKVLPFSIANSNVNEEVRLKYRYLDLRRPEMQANVQMRHDVTFAIRNYLASQNFLEVETPIMTRSTPEGARDYLVPSRVHPGEFYALPQSPQIFKQ
    END

    BEGIN
    ...

    The particular keys are not necessarily as above, just the = BEGIN/END and general layout.

    Information about ARB import filters can be found at:
      http://help.arb-home.de/importift.html
    """

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')

    def read(self, greengenes_file, public):
        """Read records from a GreenGenes file.

        Parameters
        ----------
        greengenes_file : str
            Name of GreenGenes file to read.
        public : boolean
            Flag indicating if only public records should be read.

        Returns
        -------
        dict : dict[genome_id] -> metadata dictionary
            Meatadata for genomes.
        """

        logging.info('Reading taxonomy information from ARB GreenGenes file.')
        genome_metadata = {}
        for entry in self.each(open(greengenes_file)):
            db_name = entry['db_name']

            try:
                if not public:  # read every record
                    genome_metadata[db_name] = entry
                elif entry['core_list_status'] == 'public':  # only retrieve public records
                    genome_metadata[db_name] = entry
            except KeyError:
                logging.warn("Metadata record not found for ID: %s, skipping" % db_name)

        return genome_metadata

    def each(self, file_handle):
        """Generator function for reading GreenGenes formatted files.

        Parameters
        ----------
        file_handle : input stream
            Stream to GreenGenes formatted data.

        Returns
        -------
        dict
            Metadata described in a single BEGIN/END block.
        """

        # state machine - are we inside or outside a BEGIN/END block
        state = 'outside'
        current_hash = {}

        for line_number, line in enumerate(file_handle):
            line = line.strip()

            if len(line) == 0:
                continue

            if line == 'BEGIN':
                if state == 'outside':
                    state = 'inside'
                    current_hash = {}
                else:
                    raise Exception("Badly formatted file type 1 on line %s, cannot continue, error on line: %s" % (line_number + 1, line))

            elif line == 'END':
                if state == 'inside':
                    yield current_hash
                    state = 'outside'
                else:
                    raise Exception("Badly formatted file type 2 on line %s, cannot continue, error on line: %s" % (line_number + 1, line))

            elif state == 'inside':
                splits = line.split('=', 1)
                if len(splits) == 2:
                    current_hash[splits[0]] = splits[1]
                else:
                    raise Exception("Badly formatted file type 3 on line %s, cannot continue, error on line: %s" % (line_number + 1, line))

            else:
                raise Exception("Badly formatted file type 4 on line %s, cannot continue, error on line: %s" % (line_number + 1, line))

        if state != 'outside':
            raise Exception("Badly formatted file type 5 on line %s, cannot continue, error at the end of the file (ended while expecting an END)" % (line_number + 1))

    def write(self, hashes, io):
        """Write data to a GreenGenes formatted files.

        Parameters
        ----------
        hashes : list of dict[feature] -> value
            GreenGenes style metadata.
        io : output stream
            Handle to output stream.
        """
        is_first = True
        for dahash in hashes:
            if is_first:
                is_first = False
            else:
                io.write('\n')

            io.write('BEGIN\n')
            for key in sorted(dahash.keys()):
                if key != 'warning' and key != 'aligned_seq':
                    io.write('='.join([key, dahash[key]]) + '\n')

            # the warning field must be the second to last as it is used
            # to indicates that the aligned sequence is to follow
            io.write('='.join(['warning', dahash.get('warning', '')]) + '\n')

            # the aligned sequence must be the last field
            io.write('='.join(['aligned_seq', dahash.get('aligned_seq', '')]) + '\n')

            io.write('END\n')
