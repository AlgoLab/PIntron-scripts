#!/usr/bin/env python

from __future__ import print_function

import argparse
import logging
import glob
import os.path

##


VALID_SEQNAMES = [str(x) for x in range(1, 22)] + ["X", "Y", "MT"]


class GTFRecord:
    def __init__(self, line, attribute_sep):
        pline = line.split("\t")
        self.seqname = pline[0]
        if self.seqname.startswith("chr"):
            self.seqname = self.seqname[3:]
        self.source = pline[1]
        self.feature = pline[2]
        self.start = int(pline[3])
        self.end = int(pline[4])
        self.score = pline[5]
        self.strand = pline[6]
        self.frame = pline[7]
        aline = pline[8].split("; ")
        self.attributes = dict([(key, value.strip("\""))
                                for key, value in
                                [elem.strip().split(attribute_sep, 1)
                                 for elem in aline]])

    def __str__(self):
        return ("{self.seqname}\t{self.source}\t{self.feature}\t" +
                "{self.start}\t{self.end}\t{self.score}\t" +
                "{self.strand}\t{self.frame}\t").format(self=self)


def gtf_parser(stream, attribute_sep=" "):
    for line in stream:
        line = line.strip()
        if line.startswith("#"):
            continue
        record = GTFRecord(line, attribute_sep)
        yield record


def filter_all(filter_list):
    return lambda record: all([f(record) for f in filter_list])


def filter_sources(valid_sources=("havana", "ensembl_havana")):
    return lambda record: record.source in valid_sources


def filter_seqnames(valid_seqnames=VALID_SEQNAMES):
    return lambda record: record.seqname in valid_seqnames


def filter_features(valid_features=("transcript", "gene", "exon")):
    return lambda record: record.feature in valid_features


def filter_attribute(attr_name, attr_value):
    return lambda record: record.attributes[attr_name] == attr_value


def filter_attributes(attr_name, attr_values):
    return lambda record: record.attributes[attr_name] in attr_values


def exons2introns(exons):
    return [(e1.end+1, e2.start-1) for e1, e2 in zip(exons[:-1], exons[1:])]


class Transcript:
    def __init__(self, gene, transcript_id=None):
        self.gene = gene
        self.transcript_id = transcript_id
        self.exons = []

    def __str__(self):
        return "{transcript_id}:\t [\n  {exons}]".format(
            transcript_id=self.transcript_id,
            exons="\n  ".join(["{x[0]}, {x[1]}".format(x=x)
                               for x in self.exons])
        )

    def introns(self):
        return exons2introns(self.exons)

    def is_predicted(self, exon_list):
        logging.debug("Testing if %s is predicted.", self.transcript_id)
        cond_a = True
        cond_b = all([(exon[0], exon[1]) in exon_list
                      for exon in self.exons[1:-1]])
        first_exon_idx = 0 if self.gene.strand == "+" else -1
        last_exon_idx = -1 if self.gene.strand == "+" else 0
        cond_c = ((self.exons[first_exon_idx] in exon_list) or
                  (self.exons[first_exon_idx][1] in [e[1] for e in exon_list]))
        cond_d = ((self.exons[last_exon_idx] in exon_list) or
                  (self.exons[last_exon_idx][0] in [e[0] for e in exon_list]))
        logging.debug("is_predicted cond_a? %s", cond_a)
        logging.debug("is_predicted cond_b? %s", cond_b)
        logging.debug("is_predicted cond_c? %s", cond_c)
        logging.debug("is_predicted cond_d? %s", cond_d)
        return cond_a and cond_b and cond_c and cond_d


class Gene:
    def __init__(self, gene_id, gene_name):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.seqname = None
        self.strand = None
        self.transcripts = {}

    def add_transcript(self, transcript_id):
        self.transcripts[transcript_id] = Transcript(self, transcript_id)

    def __str__(self):
        return "{gene_name} [{gene_id}]:\n{transcripts}".format(
            gene_name=self.gene_name,
            gene_id=self.gene_id,
            transcripts="\n".join([str(x) for x in self.transcripts.values()])
        )


import itertools
import gzip


def main():
    pred_folder = "predictions"
    pred_file = "predicted-exons-filtered.gtf"
    annot_gtf = "Homo_sapiens.GRCh38.78.gtf.gz"

    parser = argparse.ArgumentParser()
    parser.add_argument('gene',
                        nargs='+')
    parser.add_argument('-p', '--prediction-root-dir',
                        default=pred_folder)
    parser.add_argument('-a', '--annotations',
                        default=annot_gtf)
    parser.add_argument('--prediction-filename',
                        default=pred_file)
    parser.add_argument('-v', '--verbose',
                        help='increase output verbosity',
                        action='count', default=0)

    args = parser.parse_args()
    args = vars(args)
    if args['verbose'] == 0:
        log_level = logging.INFO
    elif args['verbose'] == 1:
        log_level = logging.DEBUG
    else:
        log_level = logging.DEBUG

    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt="%y%m%d %H%M%S")

    pred_folder = args["prediction_root_dir"]
    pred_file = args["prediction_filename"]
    annot_gtf = args["annotations"]

    gene_names = args["gene"]

    filter_fn = filter_all((filter_attributes("gene_name", gene_names),
                            filter_seqnames(),
                            filter_sources(),
                            filter_features(("exon",))))

    logging.info("Reading annotations from '%s'...", annot_gtf)
    genes = {}
    with gzip.open(annot_gtf) as gtfin:
        for record in itertools.ifilter(filter_fn,
                                        gtf_parser(gtfin)):
            if record.attributes["gene_id"] not in genes:
                gene = Gene(record.attributes["gene_id"],
                            record.attributes["gene_name"])
                gene.seqname = record.seqname
                gene.strand = record.strand
                genes[record.attributes["gene_id"]] = gene
            else:
                gene = genes[record.attributes["gene_id"]]

            if record.attributes["transcript_name"] not in gene.transcripts:
                gene.add_transcript(record.attributes["transcript_name"])
            transcript = gene.transcripts[record.attributes["transcript_name"]]

            transcript.exons.append((record.start, record.end,
                                     record.attributes["exon_id"]))

    logging.info("Performing comparisons...")
    for gene in genes.values():
        logging.info("Gene: %s", gene.gene_name)
        gene_name = gene.gene_name
        with open("report-"+gene_name+".csv", "w") as repfout:
            for pred_name in glob.glob(os.path.join(pred_folder,
                                                    gene_name,
                                                    "*",
                                                    pred_file)):
                pred_exons = []
                dataset_name = os.path.basename(os.path.dirname(pred_name))
                dataset_name = dataset_name.replace("_", ",")
                logging.info("  Prediction: %s", pred_name)
                with open(pred_name) as predfin:
                    for record in gtf_parser(predfin, "="):
                        assert record.seqname == gene.seqname, \
                            "Gene {} is on sequence {} but predictions are on sequence {}." + \
                            "Read: {}".format(gene_name,
                                              gene.seqname,
                                              record.seqname,
                                              str(record))
                        pred_exons.append((record.start, record.end))

                for transcript in gene.transcripts.values():
                    print(",".join([gene.gene_id, gene.gene_name,
                                    transcript.transcript_id,
                                    dataset_name,
                                    str(transcript.is_predicted(pred_exons))]),
                          file=repfout)


if __name__ == "__main__":
    main()
