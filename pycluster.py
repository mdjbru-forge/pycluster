### * Description

# Run mcl and build unaligned clusters from blastp output and gene table

### * Set up

### ** Import

import os
import subprocess

### * Functions

### ** prepareBlastpABC(results, out = None, q = 1, s = 2, e = 11)

def prepareBlastpABC(results, out = None, q = 1, s = 2, e = 11) :
    """Convert a result file to ABC format

    Args:
        results (str): Name of the blastp result file
        out (str): Name of the output file. If None, append the '.ABC' suffix 
          to the result file
        q (int): Query field
        s (int): Subject field
        e (int): E-value field

    Returns:
        Name of the output file

    """
    command = "cut -f " + str(q) + "," + str(s) + "," + str(e) + " "
    command += results
    if out is None :
        out = results + ".ABC"
    command += " > " + out
    os.system(command)
    return out

### ** runMcxload(ABC_file, mci_file, tab_file)

def runMcxload(ABC_file, mci_file, tab_file) :
    """Run Mcxload on an ABC file

    Args:
        ABC_file (str): Name of the input ABC file
        mci_file (str): Name of the output mci file
        tab_file (str): Name of the output tab file

    Returns:
        A tuple (mci_file, tab_file)

    """
    command = ["mcxload"]
    command += ["-abc", ABC_file]
    command += ["--stream-mirror"]
    command += ["--stream-neg-log10"]
    command += ["-stream-tf",  "ceil(200)"]
    command += ["-o", mci_file]
    command += ["-write-tab", tab_file]
    p = subprocess.Popen(command, stdout = subprocess.PIPE,
                         stderr = subprocess.PIPE)
    p.wait()
    return (mci_file, tab_file)

### ** runMcl(mci_file, inflation, mcl_tab, out = None, n_threads = 1)

def runMcl(mci_file, inflation, mcl_tab, out = None, n_threads = 1) :
    """Run MCL on a mci file (network file)

    Args:
        mci_file (str): Name of the input mci file
        inflation (float): Inflation value
        mcl_tab (str): Name of the tab file
        out (str): Output file name. If None, append '.MCL' to the input mci 
          file
        n_threads (int): Number of threads

    Returns:
        Output file name

    """
    command = ["mcl"]
    command += [mci_file]
    command += ["-I", str(inflation)]
    command += ["-use-tab", mcl_tab]
    if out is None :
        out = mci_file + ".cluster"
    command += ["-o", out]
    command += ["-te", str(n_threads)]
    p = subprocess.Popen(command)
    p.wait()
    return out    
