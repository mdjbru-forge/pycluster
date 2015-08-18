### * Description

# Entry points for the command line scripts

### * Wishlist

# pycluster ABC blastp.out -o blastp.out.ABC
# pycluster mcl blastp.out.ABC -I 1.4 -o blastp.mcl-1.4
# pycluster extract gene.table blastp.mcl-1.4

### * Setup

### ** Import

import sys
import os
import argparse
import hashlib
import pygenes as pygenes
import pycluster as pycluster

### * Parser

def makeParser() :
    """Prepare the parser

    Returns:
        ArgumentParser: An argument parser

    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help = "One of the available actions")
    # Convert blastp output to ABC format
    sp_ABC = subparsers.add_parser("ABC",
                                   help = "Convert blastp output to ABC format")
    sp_ABC.add_argument("blastpOut", metavar = "BLASTP_TAB",
                          type = str,
                          help = "Blastp output file (tabular format)")
    sp_ABC.add_argument("-o", "--output", metavar = "BLASTP_ABC",
                        help = "Converted output file (ABC format) "
                        "(default: input file name with ABC extension)")
    sp_ABC.add_argument("-q", "--query", metavar = "FIELD_NUM", type = int,
                        default = 1,
                        help = "Field number for query (default: 1)")
    sp_ABC.add_argument("-s", "--subject", metavar = "FIELD_NUM", type = int,
                        default = 2,
                        help = "Field number for subject (default: 2)")
    sp_ABC.add_argument("-e", "--evalue", metavar = "FIELD_NUM", type = int,
                        default = 11,
                        help = "Field number for e-value (default: 11)")
    sp_ABC.set_defaults(action = "ABC")
    # Perform mcl
    sp_mcl = subparsers.add_parser("mcl",
                                    help = "Run mcl (Markov cluster algorithm)")
    sp_mcl.add_argument("input", metavar = "ABC_FILE", type = str,
                        help = "Blastp output in ABC format")
    sp_mcl.add_argument("-I", "--inflation", metavar = "FLOAT", type = float,
                        default = 1.4,
                        help = "Inflation value used by MCL (default: 1.4)")
    sp_mcl.add_argument("-n", "--nthreads", type = int, default = 1,
                        help = "Number of threads (default: 1)")
    sp_mcl.add_argument("-o", "--output", metavar = "FILE", type = str,
                        help = "Output file name (default: input file name "
                        "with clusters extension)")
    sp_mcl.set_defaults(action = "mcl")
    # Make mapping table
    sp_table = subparsers.add_parser("table",
                                     help = "Convert mcl output to a tabular file")
    sp_table.add_argument("clusters", metavar = "MCL_OUTPUT", type = str,
                          help = "Mcl output file")
    sp_table.add_argument("-o", "--output", metavar = "TAB_OUTPUT", type = str,
                          help = "Cluster table (default: input file with "
                          "table extension)")
    sp_table.set_defaults(action = "table")
    # Extract sequences
    sp_extract = subparsers.add_parser("extract",
                                       help = "Extract peptide sequences "
                                       "corresponding to clusters")
    sp_extract.add_argument("geneTable", metavar = "GENE_TABLE", type = str,
                          help = "Gene table file (produced by pygenes)")
    sp_extract.add_argument("clusters", metavar = "CLUSTER_TABLE", type = str,
                            help = "Cluster table file (columns: gene "
                            "identifier, cluster identifier)")
    sp_extract.add_argument("outDir", metavar = "OUTPUT_DIR", type = str,
                            help = "Output folder")
    sp_extract.set_defaults(action = "extract")
    # Return
    return parser

### * Mains

### ** Main entry point

def main(args = None, stdout = None, stderr = None) :
    """Main entry point

    Args:
        args (namespace): Namespace with script arguments, parse the command 
          line arguments if None
        stdout (file): Writable stdout stream (if None, use `sys.stdout`)
        stderr (file): Writable stderr stream (if None, use `sys.stderr`)

    """
    if args is None :
        parser = makeParser()
        args = parser.parse_args()
    if stdout is None :
        stdout = sys.stdout
    if stderr is None :
        stderr = sys.stderr
    dispatch = dict()
    dispatch["ABC"] = main_ABC
    dispatch["mcl"] = main_mcl
    dispatch["table"] = main_table
    dispatch["extract"] = main_extract
    dispatch[args.action](args, stdout, stderr)

### ** Main ABC

def main_ABC(args, stdout, stderr) :
    pycluster.prepareBlastpABC(args.blastpOut, args.output, args.query,
                               args.subject, args.evalue)

### ** Main mcl

def main_mcl(args, stdout, stderr) :
    pycluster.runMcxload(args.input, args.input + ".mci",
               args.input + ".tab")
    if args.output is None :
        args.output = args.input + ".cluster"
    pycluster.runMcl(args.input + ".mci", args.inflation,
           args.input + ".tab",
           args.output, args.nthreads)
    os.remove(args.input + ".mci")
    os.remove(args.input + ".tab")

### ** Main table

def main_table(args, stdout, stderr) :
    if args.output is None :
        args.output = args.clusters + ".table"
    n = 0
    with open(args.clusters, "r") as fi :
        with open(args.output, "w") as fo :
            for l in fi :
                l = l.strip()
                if l != "" :
                    items = l.split("\t")
                    for i in items :
                        fo.write(i + "\t" + "cl" + str(n) + "\n")
                    n += 1
    
### ** Main extract

def main_extract(args, stdout, stderr) :
    # Load gene table
    geneTable = pygenes.GeneTable()
    geneTable.loadTable(args.geneTable)
    # Load cluster mapping
    mergedPep2cl = dict()
    with open(args.clusters, "r") as fi :
        for l in fi :
            l = l.strip()
            if l != "" :
                l = l.split("\t")
                geneId = l[0]
                clusterId = l[1]
                assert not mergedPep2cl.get(geneId, False)
                mergedPep2cl[geneId] = clusterId
    # Check output directory            
    if not os.path.isdir(args.outDir) :
        if os.path.isfile(args.outDir) :
            raise Exception("Output directory name already used for a file")
        os.makedirs(args.outDir)
    # Build the extracted clusters
    extClusters = dict()
    for gene in geneTable :
        clId = mergedPep2cl.get(gene.mergedPeptideHash, "NA")
        extClusters[clId] = extClusters.get(clId, []) + [(gene.geneId,
                                                          gene.peptideSeq)]
    # Write the clusters
    for (clId, seqs) in extClusters.items() :
        if clId != "NA" :
            path = os.path.join(args.outDir, clId + ".fa")
            with open(path, "w") as fo :
                for seq in seqs :
                    fo.write(">" + seq[0] + "\n")
                    fo.write(seq[1] + "\n")

        

