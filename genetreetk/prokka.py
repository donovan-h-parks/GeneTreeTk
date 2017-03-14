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

__prog_name__ = 'runProkka'
__prog_desc__ = 'run Prokka on a set of genome bins'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os


class Prokka(object):
    def __init__(self, cpus):
        """Initialization."""
        self.cpus = cpus

    def run(self, genome_dir, kingdom, extension, output_dir):
        """Run Prokka on each genome bin."""

        files = os.listdir(genome_dir)
        for f in files:
            genome_file = os.path.join(genome_dir, f)

            genome_id = os.path.basename(genome_file)
            genome_id = os.path.splitext(genome_id)[0]

            os.system('prokka --force --cpus %d --prefix %s --outdir %s --kingdom %s %s' % (self.cpus, genome_id, output_dir, kingdom, genome_file))
