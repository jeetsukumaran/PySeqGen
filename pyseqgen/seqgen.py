#! /usr/bin/env python

############################################################################
##  seqgen.py
##
##  Part of the PySeqGen wrapper library.
##
##  Copyright 2007 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

import subprocess
import StringIO
import tempfile
import random
import os
import sys

from optparse import OptionGroup
from optparse import OptionParser

import dendropy

LONG_MAX = 2**32

class SeqGen(object):
    """
    This class wraps all attributes and input needed to make a call to SeqGen.
    """

    class SubstitutionModel(object):
        def __init__(self, idstr):
            self.idstr = idstr
        def __str__(self):
            return self.idstr

    F84 = SubstitutionModel("F84")
    HKY = SubstitutionModel("HKY")
    GTR = SubstitutionModel("GTR")
    JTT = SubstitutionModel("JTT")
    WAG = SubstitutionModel("WAG")
    PAM = SubstitutionModel("PAM")
    BLOSUM = SubstitutionModel("BLOSUM")
    MTREV = SubstitutionModel("MTREV")
    CPREV = SubstitutionModel("CPREV")
    GENERAL = SubstitutionModel("GENERAL")
    MODELS = [F84, HKY, GTR, JTT, WAG, PAM, BLOSUM, MTREV, CPREV, GENERAL]
    MODEL_IDS = [str(m) for m in MODELS]


    def get_model(idstr):
        for model in SeqGen.MODELS:
            if idstr.upper() == model.idstr.upper():
                return model
        return None
    get_model = staticmethod(get_model)

    def __init__(self):
        """
        Sets up all properties, which (generally) map directly to command
        parameters of Seq-Gen.
        """

        # python object specific attributes
        self.seqgen_path = '/usr/local/bin/seq-gen'
        self.num_replicates = 1
        self.rng_seed = None
        self.__rng = None
        self.trees = None
        self.num_sites_per_tree = None
        self.rel_rates_per_tree = None
        self.overwrite = False
        self.output_dir = '.'
        self.output_filename_prefix = 'data'
        self.output_filename_suffix = None
        self.output_nexus = False
        self.simple_nexus = True
        self.dump_tree = False

        # following are passed to seq-gen in one form or another
        self.char_model = 'HKY'
        self.seq_len = None
        self.num_partitions = None
        self.scale_branch_lens = None
        self.scale_tree_len = None
        self.codon_pos_rates = None
        self.gamma_shape = None
        self.gamma_cats = None
        self.prop_invar = None
        self.state_freqs = None
        self.ti_tv = None
        self.general_rates = None
        self.ancestral_seq = None
        self.output_text_append = None
        self.write_ancestral_seqs = False
        self.write_site_rates = False
        self.quiet = False
        self.dry_run = False
        self.dump_command = True

    def compose_output_filepath(self, filename_index=1):
        output_dir = os.path.expanduser(os.path.expandvars(self.output_dir))
        filename_prefix = os.path.expanduser(os.path.expandvars(self.output_filename_prefix))
        if self.output_filename_suffix is not None:
            suffix = self.output_filename_suffix
        else:
            suffix = ""
        if self.output_nexus:
            filename_ext = "nex"
        else:
            filename_ext = "xml"
        output_filepath = os.path.join(output_dir, "%s%03d%s.%s" % (filename_prefix, filename_index, suffix, filename_ext))
        return output_filepath

    def _get_rng(self):
        if self.__rng is None:
            self.__rng = random.Random(self.rng_seed)
        return self.__rng

    def _set_rng(self, rng):
        self.__rng = rng

    rng = property(_get_rng, _set_rng)

    def compose_input_string(self):
        """
        Given a list of trees (as PyNexml/DendroPy tree objects),
        this will compose the input `file` for the Seq-Gen call. If more than
        one tree is given, `self.num_sites_per_tree` should be list of integers specifying the
        number of characters/base-pairs to be generated on each tree.
        `self.rel_rate_per_tree`, if specified, should be a list of reals specifying the relative
        rate by which to scale each tree.
        """
        if isinstance(self.trees, list):
            input = []
            for idx, tree in enumerate(self.trees):
                if self.num_sites_per_tree or self.rel_rates_per_tree:
                    parts = []
                    if self.num_sites_per_tree:
                        parts.append(self.num_sites_per_tree[idx])
                    if self.rel_rates_per_tree:
                        parts.append(self.rel_rates_per_tree[idx])
                    parts = '[%s]' % (','.join(parts))
                else:
                    parts = ""
                input.append("%s%s" % (parts, tree.as_string(schema="newick", write_rooting=False, internal_labels=False)))
            return '\n'.join(input)
        else:
            return "".join([t.as_string(schema="newick", write_rooting=False, internal_labels=False) for t in self.trees])
#
#     def compose_input_file(self, trees, num_sites=None, rel_rates=None):
#         return StringIO.StringIO(self.compose_input_string(trees, num_sites, rel_rates))

    def compose_arguments(self):
        """
        Composes and returns a list of strings that make up the arguments to a Seq-Gen
        call, based on the attribute values of the object.
        """

        args = []
        args.append(self.seqgen_path)

        if self.char_model:
            args.append("-m%s" % str(self.char_model))

        if self.seq_len:
            args.append("-l%s" % self.seq_len)

        if self.num_partitions:
            args.append("-p%s" % self.num_partitions)

        if self.scale_branch_lens:
            args.append("-s%s" % self.scale_branch_lens)

        if self.scale_tree_len:
            args.append("-d%s" % self.scale_tree_len)

        if self.codon_pos_rates:
            args.append("-c%s" % (",".join(self.codon_pos_rates)))

        if self.gamma_shape:
            args.append("-a%s" % self.gamma_shape)

        if self.gamma_cats:
            args.append("-g%s" % self.gamma_cats)

        if self.prop_invar:
            args.append("-i%s" % self.prop_invar)

        if self.state_freqs:
            if isinstance(self.state_freqs, str):
                args.append("-f%s" % self.state_freqs)
            else:
                args.append("-f%s" % (",".join([str(s) for s in self.state_freqs])))

        if self.ti_tv:
            args.append("-t%s" % self.ti_tv)

        if self.general_rates:
            if isinstance(self.general_rates, str):
                args.append("-r%s" % self.general_rates)
            else:
                args.append("-r%s" % (",".join(self.general_rates)))

        if self.ancestral_seq:
            args.append("-k%s" % self.ancestral_seq)

        if self.output_text_append:
            args.append("-x'%s'" % self.output_text_append)

        if self.write_ancestral_seqs:
            args.append("-wa")

        if self.write_site_rates:
            args.append("-wr")

        if self.quiet:
            args.append("-q")

        # following are controlled directly by the wrapper

        # we explicitly pass a random number seed on each call
        args.append("-z%s" % self.rng.randint(0, LONG_MAX))

        # force nexus
        args.append("-on")

        # force one dataset at a time
        args.append("-n1")

        return args

#     def get_trees_from_nexus_file(self, filepath):
#         filepath = os.path.expanduser(os.path.expandvars(filepath))
#         dataset = dataio.from_nexus(filepath)
#         if dataset.tree_lists:
#             trees = dataset.tree_lists[0]
#
#     def get_trees_from_newick_file(self, filepath):
#         filepath = os.path.expanduser(os.path.expandvars(filepath))
#         dataset = dataio.from_newick(filepath)
#         if dataset.tree_lists:
#             trees = dataset.tree_lists[0]
#
#     def get_trees_from_nexml_file(self, filepath):
#         filepath = os.path.expanduser(os.path.expandvars(filepath))
#         dataset = dataio.from_nexml(filepath)
#         if dataset.tree_lists:
#             trees = dataset.tree_lists[0]

    def generate_and_save_nexus(self):
        index_offset = 1
        if not os.path.exists(self.seqgen_path):
            raise Exception('seq-gen not found at "%s"' % self.seqgen_path)
        input_string = self.compose_input_string()
        if not self.quiet and self.dump_tree:
            print >>sys.stderr, input_string
        for rep in range(self.num_replicates):
            output_filepath = self.compose_output_filepath(rep+index_offset)
            args = self.compose_arguments()
            if not self.quiet and self.dump_command:
                print >>sys.stderr, "--- INVOKING Seq-Gen: %s" % (' '.join(args))
            if not self.dry_run:
                input = tempfile.NamedTemporaryFile()
                input.write(input_string)
                input.flush()
                output = open(output_filepath, "w")
                args.append(input.name)
                run = subprocess.Popen(args, shell=False, stdin=None, stdout=output, stderr=subprocess.PIPE)
                stdout, stderr = run.communicate()
                if stderr.upper().count("ERROR"):
                    sys.stderr.write("--- Seq-Gen ERROR DETECTED:\n")
                    sys.stderr.write(stderr)
                    sys.stderr.write("--- Seq-Gen EXITED WITH ERROR\n\n")
                    sys.stderr.write("Failed input tree was:\n%s\n\n" % input_string)
                elif not self.quiet:
                    print >>sys.stderr, "--- Seq-Gen EXITED WITH STATUS CODE %d\n" % run.returncode

    def generate_and_save_datasets(self):
        index_offset = 1
        input_string = self.compose_input_string()
        if not self.quiet and self.dump_tree:
            print >>sys.stderr, input_string
        for rep in range(self.num_replicates):
            dataset = self.generate_dataset(input_string=input_string)
            if dataset:
                output_filepath = self.compose_output_filepath(rep+index_offset)
                if not self.overwrite:
                    while os.path.exists(output_filepath):
                        index_offset = index_offset + 1
                        output_filepath = self.compose_output_filepath(rep+index_offset)
                if self.output_nexus:
                    schema = "nexus"
                else:
                    schema = "nexml"
                dataset.write(open(output_filepath, "w"), schema)

    def generate_dataset(self, input_string=None, dataset=None):
        #if not os.path.exists(self.seqgen_path):
        #    raise Exception('seq-gen not found at "%s"' % self.seqgen_path)
        if input_string is None:
            input_string = self.compose_input_string()
        args=self.compose_arguments()
        if not self.quiet and self.dump_command:
            print >>sys.stderr, ' '.join(args)
            #print >>sys.stderr, input_string
        if not self.dry_run:
            #inputf = open('/tmp/input', 'w')
            inputf = tempfile.NamedTemporaryFile()
            inputf.write(input_string)
            inputf.flush()
            #outputf = open('/tmp/output', 'w')
            outputf = tempfile.NamedTemporaryFile()
            args.append(inputf.name)

            run = subprocess.Popen(args,
                                   stdin=None,
                                   stdout=outputf)
            run.communicate()

            if not self.quiet:
                print >>sys.stderr, "\n--(SEQ-GEN PROCESS EXITED WITH STATUS %d)--" % run.returncode
            if not self.quiet:
                print >>sys.stderr, "\nReading generated file ..."

            if dataset is None:
                dataset = dendropy.DataSet()
            dataset.read(open(outputf.name, "rU"), "nexus")
            return dataset


def parse_list_option(list_str, option_name, number_type=float):
    if list_str.count(','):
        list_str = list_str.replace(' ', '').split(',')
    else:
        list_str = list_str.split(' ')
    values = []
    for idx, item in enumerate(list_str):
        try:
            value = item #number_type(item)
        except ValueError:
            print >>sys.stderr, "%s: item %d in '%s' is not a valid %s value" % (option_name, idx+1, list_str, str(number_type))
            sys.exit(1)
        values.append(value)
    return values

if __name__ == "__main__":
    usage = '%prog [options]'
    description = 'Seq-Gen Wrapper'
    version = 'PySeqGen Version 1.0'

    seqgen = SeqGen()

    parser = OptionParser(usage=usage, add_help_option=True, version = version, description=description)

    sim_optgroup = OptionGroup(parser, "Simulation Options")
    sim_optgroup.add_option('--seq-gen', action='store', dest='seqgen_path', default='/usr/bin/seq-gen', metavar='SEQ-GEN PATH', help='path to seq-gen binary [\'/usr/bin/seq-gen\']')
    sim_optgroup.add_option('--random-seed', action='store', dest='rng_seed', default=None, metavar='RANDOM SEED', help='integer seed for random number generator')
    sim_optgroup.add_option('--seq-len', action='store', dest='seq_len', default=None, metavar='SEQ-LEN', help='length of sequences')
    sim_optgroup.add_option('--num-partitions', action='store', dest='num_partitions', default=None, metavar='NUM-PARTS', help='number of data partitions')
    sim_optgroup.add_option('--ancestral-seq', action='store', dest='ancestral_seq', default=None, metavar='ANCESTRAL SEQUENCE', help='set ancestral character sequence')
    sim_optgroup.add_option('--num-replicates', action='store', dest='num_replicates', default=1, metavar='NUM REPLICATES', help='number of independent datasets to generate')
    sim_optgroup.add_option('--quiet', action='store_true', dest='quiet', default=False, help='suppress progress messages')
    sim_optgroup.add_option('--dry-run', action='store_true', dest='dry_run', default=False, help='do not actually do anything')
    parser.add_option_group(sim_optgroup)

    input_tree_optgroup = OptionGroup(parser, "Model Tree Options")
    input_tree_optgroup.add_option('--from-nexus', default=None, dest='from_nexus', metavar='TREE-FILE', help="[DEFAULT] tree file in NEXUS schema")
    input_tree_optgroup.add_option('--from-newick', default=None, dest='from_newick', metavar='TREE-FILE', help="tree file in Newick schema")
    input_tree_optgroup.add_option('--from-nexml', default=None, dest='from_nexml', metavar='TREE-FILE', help="tree file in nexml schema")
    input_tree_optgroup.add_option('--scale-branches', default=None, dest='scale_branch_lens', metavar='SCALE-FACTOR', help="rescale branch lengths by SCALE-FACTOR")
    input_tree_optgroup.add_option('--scale-tree', default=None, dest='scale_tree_len', metavar='SCALE-FACTOR', help="rescale total tree length by SCALE-FACTOR")
    parser.add_option_group(input_tree_optgroup)

    char_model_optgroup = OptionGroup(parser, "Character Substitution Model Options")
    char_model_optgroup.add_option("-m", "--model", default=None, dest="char_model", metavar='MODEL', help="character substitution model (%s)" % (', '.join(SeqGen.MODEL_IDS)))
    char_model_optgroup.add_option("--ti-tv", default=None, dest="ti_tv", metavar="RATIO", help="transition/transversion ratio for HKY model")
    char_model_optgroup.add_option("--rates", default=None, dest="general_rates", metavar="'#.##, #.##, #.## ...'", help="substitution rates for general model")
    char_model_optgroup.add_option("--state-freqs", default=None, dest="state_freqs", metavar="'#.##, #.##, #.## ...' or 'e'", help="relative frequencies of character states, or 'e' to force equal")
    char_model_optgroup.add_option("--codon-pos-rates", default=None, dest="codon_pos_rates", metavar="'#.##, #.##, #.##'", help="relative rates of first, second and third codon positions")
    char_model_optgroup.add_option("--gamma-shape", default=None, dest="gamma_shape", metavar="ALPHA", help="shape parameter for gamma-distributed rates")
    char_model_optgroup.add_option("--gamma-cats", default=None, dest="gamma_cats", metavar="NUM-CATS", help="number of rate categories if using discretized gama")
    char_model_optgroup.add_option("--prop-invar", default=None, dest="prop_invar", metavar="PROP", help="proportion of invariant sites")
    parser.add_option_group(char_model_optgroup)

    output_optgroup = OptionGroup(parser, "Output Options")
    default_to_schema = 'SIMPLE-NEXUS'
    output_optgroup.add_option('--output-directory', default='.', dest='output_dir', help='directory to store simulated dataset files (default=current)')
    output_optgroup.add_option('--output-prefix', default='data', dest='output_filename_prefix', help='common prefix for data files generated')
    output_optgroup.add_option('--to-nexus', action='store_const', default=default_to_schema, const='NEXUS', dest='to_schema', help="[DEFAULT] saves the data in standard NEXUS schema (will remove comments and advanced script blocks)")
    output_optgroup.add_option('--to-simple-nexus', action='store_const', default=default_to_schema, const='SIMPLE-NEXUS', dest='to_schema', help="saves the data in the older NEXUS schema, using a DATA block instead of TAXA/CHARACTER blocks")
    output_optgroup.add_option('--to-nexml', action='store_const', default=default_to_schema, const='NEXML', dest='to_schema', help="saves the data in NEXML schema")
    output_optgroup.add_option('--write-ancestral-seqs', action='store_true', dest='write_ancestral_seqs', default=False, help='write ancestral sequences')
    output_optgroup.add_option('--write-site-rates', action='store_true', dest='write_site_rates', default=False, help='write site rates')
    output_optgroup.add_option('--overwrite', action='store_true', default=False, dest='overwrite', help="overwrite files in destination directoy if already existing if false [Default = True]")
    parser.add_option_group(output_optgroup)


    (opts, args) = parser.parse_args()

    seqgen.seqgen_path = os.path.expanduser(os.path.expandvars(opts.seqgen_path))
    if not os.path.exists(seqgen.seqgen_path):
        print >>sys.stderr, 'seq-gen binary "%s" not found' % seqgen.seqgen_path
        sys.exit(1)

    if opts.from_nexus:
        tree_filepath = opts.from_nexus
        schema_desc = "NEXUS"
    elif opts.from_newick:
        tree_filepath = opts.from_newick
        schema_desc = "NEWICK"
    elif opts.from_nexml:
        tree_filepath = opts.from_nexml
        schema_desc = "NEXML"
    else:
        print >>sys.stderr, "source of model tree(s) not specified"
        sys.exit(1)

    tree_filepath = os.path.expanduser(os.path.expandvars(tree_filepath))
    if not os.path.exists(tree_filepath):
        print >>sys.stderr, 'tree file "%s" not found' % tree_filepath
        sys.exit(1)

#     source_dataset = dataio.dataset_from_file(tree_filepath, schema=schema_desc)
    try:
        source_dataset = dendropy.DataSet()
        source_dataset.read(open(tree_filepath, "rU"), schema=schema_desc)
    except Exception, e:
        print >>sys.stderr, e
        print >>sys.stderr, '"%s" is not a valid %s file' % (tree_filepath, schema_desc)
        raise e

    if source_dataset and len(source_dataset.tree_lists) > 0 and len(source_dataset.tree_lists[0]) > 0:
        seqgen.trees = source_dataset.tree_lists[0]
    else:
        print >>sys.stderr, 'no trees found in tree file "%s"' % tree_filepath
        sys.exit(1)

    if not opts.char_model:
        print >>sys.stderr, "character substitution model not specified"
        sys.exit(1)

    seqgen.char_model = SeqGen.get_model(opts.char_model)
    if seqgen.char_model is None:
        print >>sys.stderr, "'%s' no a recognized character model" % opts.char_model
        sys.exit(1)

    seqgen.ti_tv = opts.ti_tv
    if opts.general_rates:
        seqgen.general_rates = parse_list_option(opts.general_rates, "rates", float)
    if opts.state_freqs:
        if opts.state_freqs.lower() != 'e':
            seqgen.state_freqs = parse_list_option(opts.state_freqs, "state-freqs", float)
        else:
            seqgen.state_freqs = 'e'
    if opts.codon_pos_rates:
        seqgen.codon_pos_rates = parse_list_option(opts.codon_pos_rates, "codon-pos-rates", float)
    seqgen.gamma_shape = opts.gamma_shape
    seqgen.gamma_cats = opts.gamma_cats
    seqgen.prop_invar = opts.prop_invar

    seqgen.scale_branch_lens = opts.scale_branch_lens
    seqgen.scale_tree_len = opts.scale_tree_len

    seqgen.dry_run = opts.dry_run
    if opts.rng_seed:
        seqgen.rng_seed = long(opts.rng_seed)
    if opts.num_replicates:
        seqgen.num_replicates = long(opts.num_replicates)
    seqgen.seq_len = opts.seq_len
    seqgen.num_partitions = opts.num_partitions
    seqgen.ancestral_seq = opts.ancestral_seq
    seqgen.write_ancestral_seqs = opts.write_ancestral_seqs
    seqgen.write_site_rates = opts.write_site_rates
    seqgen.quiet = opts.quiet

    seqgen.output_dir = os.path.expanduser(os.path.expandvars(opts.output_dir))
    seqgen.output_filename_prefix = os.path.expandvars(opts.output_filename_prefix)
    seqgen.output_nexus = opts.to_schema == 'NEXUS' or opts.to_schema == 'SIMPLE-NEXUS'
    seqgen.simple_nexus = opts.to_schema == 'SIMPLE-NEXUS'
    seqgen.overwrite = opts.overwrite

    seqgen.generate_and_save_datasets()



