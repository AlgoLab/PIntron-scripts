#!/usr/bin/env python3
####
#
#
#                              PIntron - scripts
#
# Scripts accompanying PIntron, a novel pipeline for computational
# gene-structure prediction based on spliced alignment of expressed sequences
# (ESTs and mRNAs).
#
# Copyright (C) 2011  Yuri Pirola
#
# Distributed under the terms of the GNU Affero General Public License (AGPL)
#
#
# This file is part of PIntron-scripts.
#
# The scripts contained in PIntron-scripts are free software: you can
# redistribute them and/or modify them under the terms of the GNU Affero
# General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# The scripts contained in PIntron-scripts are distributed in the hope that
# they will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with PIntron-scripts.  If not, see <http://www.gnu.org/licenses/>.
#
####


import json
import logging
import optparse
import os
import os.path
import sys

class AnalysisResults:
    pass

def extract_exons(isoforms):
    exons= []
    for isoform in isoforms:
        for exon in isoform["exons"]:
            assert "relative start" in exon and "relative end" in exon
            exons.append((exon["relative start"], exon["relative end"]))
    return frozenset(exons)

def extract_isoforms(isoforms):
    exons_list= []
    for isoform in isoforms:
        exons_list.append(tuple(((exon["relative start"], exon["relative end"])
                                 for exon in isoform["exons"])))
    return frozenset(exons_list)

def analyze_generic(gold_isoforms, new_isoforms, extract_fun):
    results= AnalysisResults()
    gold_objs= extract_fun(gold_isoforms)
    results.no_of_true= len(gold_objs)
    new_objs= extract_fun(new_isoforms)
    results.no_of_predicted= len(new_objs)

    results.FN_objs= gold_objs - new_objs
    results.FP_objs= new_objs - gold_objs
    results.TP= len(gold_objs & new_objs)
    results.FN= len(results.FN_objs)
    results.FP= len(results.FP_objs)
    if results.TP+results.FP > 0:
        results.precision= results.TP/(results.TP+results.FP)
    else:
        results.precision= float('nan')
    if results.TP+results.FN > 0:
        results.recall= results.TP/(results.TP+results.FN)
    else:
        results.recall= float('nan')
    return results



def analyze_and_print_generic(name, desc, gold_isoforms, new_isoforms, extract_fun):
    res= analyze_generic(gold_isoforms, new_isoforms, extract_fun)

    name_cap= name.capitalize()
    logging.info("TRUE %ss      (%s): %5d",
                 name, desc, res.no_of_true)
    logging.info("PREDICTED %ss (%s): %5d",
                 name, desc, res.no_of_predicted)
    logging.info("%s TP         (%s): %5d", name_cap, desc, res.TP)
    logging.info("%s FN         (%s): %5d", name_cap, desc, res.FN)
    logging.info("%s FP         (%s): %5d", name_cap, desc, res.FP)
    logging.info("%s Precision  (%s): %.3f",
                 name_cap, desc, res.precision)
    logging.info("%s Recall     (%s): %.3f",
                 name_cap, desc, res.recall)
    logging.debug("FN %ss        (%s): [%s]",
                  name, desc, ";   ".join([str(x) for x in res.FN_objs]))
    logging.debug("FP %ss        (%s): [%s]",
                  name, desc, ";   ".join([str(x) for x in res.FP_objs]))
    return res

def analyze_exons(desc, gold_isoforms, new_isoforms):
    return analyze_and_print_generic("exon", desc, gold_isoforms, new_isoforms, extract_exons)

def analyze_isoforms(desc, gold_isoforms, new_isoforms):
    return analyze_and_print_generic("isoform", desc, gold_isoforms, new_isoforms, extract_isoforms)




parser= optparse.OptionParser(usage="usage: "
                              "%prog -g <GOLD RESULTS> -n <NEW RESULTS> [-v]")
parser.add_option("-g", "--gold",
                  action="store", dest="gold",
                  type="string", default=None,
                  help="the file containing the gold-standard results.",
                  metavar="FILE")
parser.add_option("-n", "--new",
                  action="store", dest="new",
                  type="string", default=None,
                  help="the file containing the new results.",
                  metavar="FILE")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose",
                  default=False,
                  help="print additional log messages")

(options, args)= parser.parse_args()



log_level= logging.DEBUG if options.verbose else logging.INFO

logging.basicConfig(level=log_level,
                    format='%(levelname)-8s [%(asctime)s]  %(message)s')


logging.info("PIntron -- Compare PIntron predictions")

if (  options.gold is None
      or not os.path.isfile(options.gold)
      or not os.access(options.gold, os.R_OK)):
    logging.fatal("Input file for gold-standard results not specified or invalid. Given: '%s'.",
                  options.gold)
    sys.exit("File '{}' not found.".format(options.gold))

if (  options.new is None
      or not os.path.isfile(options.new)
      or not os.access(options.new, os.R_OK)):
    logging.fatal("Input file for new results not specified or invalid. Given: '%s'.",
                  options.new)
    sys.exit("File '{}' not found.".format(options.new))

rgold= {}
rnew= {}
logging.info("Gold results: %s", options.gold)
with open(options.gold, "r") as fgold:
    rgold= json.load(fgold)
logging.info("New results:  %s", options.new)
with open(options.new, "r") as fnew:
    rnew= json.load(fnew)

if "program_version" in rgold and "program_version" in rnew:
    logging.info("Comparing results computed by %s with %s",
                 rgold["program_version"].strip(), rnew["program_version"].strip())

if "version" in rgold and "version" in rnew:
    assert rgold["version"]==rnew["version"]
    logging.info("Prediction result format: %s", rgold["version"])

assert "version" not in rgold or int(rgold["version"])>=3
assert "version" not in rnew or int(rnew["version"])>=3

logging.info("Checking consistency...")
logging.debug("  ...of genomic sequence...")
assert rgold["genome"]["sequence_id"]==rnew["genome"]["sequence_id"]

logging.info("The two files are consistent!")


rgold_iso_all= rgold["isoforms"].values() 
rnew_iso_all= rnew["isoforms"].values()
rgold_iso_rs= [ isoform for isoform in rgold["isoforms"].values() if isoform["from RefSeq?"] ]
rnew_iso_rs= [ isoform for isoform in rnew["isoforms"].values() if isoform["from RefSeq?"] ]

logging.info("###  Analyzing the exons of all the isoforms...")
res_exons_all= analyze_exons("all", rgold_iso_all, rnew_iso_all)
logging.info("###  Analyzing the exons of the isoforms from RefSeqs...")
res_exons_rs= analyze_exons("only RefSeq", rgold_iso_rs, rnew_iso_rs)
logging.info("###  Analyzing all the isoforms...")
res_isoforms_all= analyze_isoforms("all", rgold_iso_all, rnew_iso_all)
logging.info("###  Analyzing the isoforms from RefSeqs...")
res_isoforms_rs= analyze_isoforms("only RefSeq", rgold_iso_rs, rnew_iso_rs)

for x in [ res_exons_all, res_exons_rs,
           res_isoforms_all, res_isoforms_rs ]:
    x.precision= "{:.3}".format(x.precision)
    x.recall= "{:.3}".format(x.recall)

print("\t".join([str(x.__dict__[y])
                 for x in [ res_exons_all, res_exons_rs,
                            res_isoforms_all, res_isoforms_rs ]
                 for y in ['no_of_true', 'no_of_predicted',
                           'TP', 'FN', 'FP',
                           'precision', 'recall']]))


