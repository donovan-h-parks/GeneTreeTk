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

import os
import sys
import logging

import dendropy

from biolib.common import check_file_exists, make_sure_path_exists
from biolib.external.execute import check_dependencies

from genetreetk.blast_workflow import BlastWorkflow
from genetreetk.concatenate import Concatenate
from genetreetk.reduce_workflow import Reduce
from genetreetk.bootstrap import Bootstrap
from genetreetk.prune import Prune
from genetreetk.prokka import Prokka
from genetreetk.create_database import CreateDatabase
from genetreetk.tree_compare import TreeCompare
from genetreetk.orthologue_workflow import OrthologueWorkflow


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger('timestamp')

    def blast(self, options):
        """Infer gene tree using BLAST."""
        
        check_file_exists(options.query_proteins)
        check_file_exists(options.db_file)
        check_file_exists(options.taxonomy_file)
        
        # sanity check arguments
        if options.prot_model == 'AUTO' and options.tree_program != 'raxml':
            self.logger.error("The 'AUTO' protein model can only be used with RAxML.")
            sys.exit(-1)

        blast_workflow = BlastWorkflow(options.cpus)
        blast_workflow.run(options.query_proteins,
                           options.db_file,
                           options.custom_db_file,
                           options.taxonomy_file,
                           options.custom_taxonomy_file,
                           options.evalue,
                           options.per_identity,
                           options.per_aln_len,
                           options.max_matches,
                           options.homology_search,
                           options.min_per_taxa,
                           options.consensus,
                           options.min_per_bp,
                           options.use_trimAl,
                           options.restrict_taxon,
                           options.msa_program,
                           options.tree_program,
                           options.prot_model,
                           options.skip_rooting,
                           options.output_dir)
                           
    def concat(self, options):
        """Infer concatenated gene tree."""
        
        make_sure_path_exists(options.output_dir)
        
        c = Concatenate(options.cpus)
        c.run(options.gene_dirs,
                options.min_per_gene,
                options.min_per_bps,
                options.tree_program,
                options.prot_model,
                options.split_chars,
                options.output_dir)
             
    def reduce(self, options):
        """Infer tree for reduced set of genes."""
        
        check_file_exists(options.homolog_file)
        check_file_exists(options.gene_ids)
        check_file_exists(options.taxonomy_file)
        
        make_sure_path_exists(options.output_dir)
        
        r = Reduce(options.cpus)
        r.run(options.homolog_file, 
                options.gene_ids, 
                options.taxonomy_file,
                options.min_per_taxa,
                options.consensus,
                options.min_per_bp,
                options.use_trimAl,
                options.msa_program,
                options.tree_program,
                options.prot_model,
                options.output_dir)
                
    def bootstrap(self, options):
        """Calculate bootstrap support for tree."""
        
        check_file_exists(options.tree)
        
        bootstrap = Bootstrap(options.cpus)
        bootstrap.run(options.tree, 
                    options.msa_file,
                    options.tree_program,
                    options.prot_model,
                    options.num_replicates,
                    options.output_dir)
                    
    def prune(self, options):
        """Prune tree."""
        
        check_file_exists(options.tree)
        check_file_exists(options.taxa_to_retain)
        
        prune = Prune()
        prune.run(options.tree,
                    options.taxa_to_retain,
                    options.output_tree)
    
    def prokka(self, options):
        """Run Prokka across multiple genome bins."""
        
        prokka = Prokka(options.cpus)
        prokka.run(options.genome_dir, 
                    options.kingdom, 
                    options.extension, 
                    options.output_dir)
                    
    def create_db(self, options):      
        """Create dereplicated GeneTreeTk-compatible database."""
        
        create_db = CreateDatabase(options.cpus)
        create_db.run(options.taxonomy,
                         options.type_strains,
                         options.genome_prot_dir,
                         options.extension,
                         options.max_taxa,
                         options.rank,
                         options.per_identity,
                         options.per_aln_len,
                         options.genomes_to_process,
                         options.keep_all_genes,
                         options.no_reformat_gene_ids,
                         options.output_dir)
                         
    def robinson_foulds(self, options):
        """Compare unrooted trees using common statistics."""
        
        check_file_exists(options.tree1)
        check_file_exists(options.tree2)
        
        tc = TreeCompare()
        if options.weighted:
            wrf = tc.weighted_robinson_foulds(options.tree1, 
                                                options.tree2,
                                                options.taxa_list)
            print(('Weighted Robinson-Foulds: %.3f' % wrf))
        else:
            rf, normalized_rf = tc.robinson_foulds(options.tree1, 
                                                    options.tree2,
                                                    options.taxa_list)
            print(('Robinson-Foulds: %d' % rf))
            print(('Normalized Robinson-Foulds: %.3f' % normalized_rf))
                         
    def supported_splits(self, options):
        """Supported bipartitions of common taxa shared between two trees."""
        
        check_file_exists(options.tree1)
        check_file_exists(options.tree2)
        
        tc = TreeCompare()
        tc.supported_splits(options.tree1, 
                            options.tree2,
                            options.split_file,
                            options.min_support,
                            options.max_depth,
                            options.taxa_list)
        
    def missing_splits(self, options):
        """Report supported bipartitions in reference tree not in comparison tree."""
        
        check_file_exists(options.ref_tree)
        check_file_exists(options.compare_tree)
        
        tc = TreeCompare()
        tc.report_missing_splits(options.ref_tree, 
                                    options.compare_tree,
                                    options.min_support,
                                    options.taxa_list)
                                    
    def midpoint(self, options):
        """"Midpoint root tree."""
        
        check_file_exists(options.in_tree)
        
        tree = dendropy.Tree.get_from_path(options.in_tree, 
                                            schema='newick', rooting='force-rooted', 
                                            preserve_underscores=True)
        tree.reroot_at_midpoint()
        
        tree.write_to_path(options.out_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)

    def orthologue(self, options):
        """Infer gene tree using BLAST after Orthologue clustering."""
        
        check_file_exists(options.query_proteins)
        check_file_exists(options.db_file)
        check_file_exists(options.taxonomy_file)

        # sanity check arguments
        if options.prot_model == 'AUTO' and options.tree_program != 'raxml':
            self.logger.error("The 'AUTO' protein model can only be used with RAxML.")
            sys.exit(-1)

        workflow = OrthologueWorkflow(options.cpus)
        workflow.run(
            query_proteins=options.query_proteins,
            db_file=options.db_file,
            #custom_db_file=options.custom_db_file,
            taxonomy_file=options.taxonomy_file,
            custom_taxonomy_file=options.custom_taxonomy_file,
            evalue=options.evalue,
            per_identity=options.per_identity,
            per_aln_len=options.per_aln_len,
            max_matches=options.max_matches,
            homology_search=options.homology_search,
            min_per_taxa=options.min_per_taxa,
            consensus=options.consensus,
            min_per_bp=options.min_per_bp,
            use_trimAl=options.use_trimAl,
            restrict_taxon=options.restrict_taxon,
            msa_program=options.msa_program,
            tree_program=options.tree_program,
            prot_model=options.prot_model,
            skip_rooting=options.skip_rooting,
            output_dir=options.output_dir)



    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if options.subparser_name == 'blast':
            self.blast(options)
        elif options.subparser_name == 'concat':
            self.concat(options)
        elif options.subparser_name == 'orthologue':
            self.orthologue(options)
        elif options.subparser_name == 'reduce':
            self.reduce(options)
        elif options.subparser_name == 'bootstrap':
            self.bootstrap(options)
        elif options.subparser_name == 'prune':
            self.prune(options)
        elif options.subparser_name == 'prokka':
            self.prokka(options)
        elif options.subparser_name == 'create_db':
            self.create_db(options)
        elif options.subparser_name == 'robinson_foulds':
            self.robinson_foulds(options)   
        elif options.subparser_name == 'supported_splits':
            self.supported_splits(options)
        elif options.subparser_name == 'missing_splits':
            self.missing_splits(options)
        elif options.subparser_name == 'midpoint':
            self.midpoint(options)   
        else:
            self.logger.error('  [Error] Unknown GeneTreeTk command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
