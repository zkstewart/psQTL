#!/usr/bin/env python3
# psQTL_view.py
# Provides an overview of a psQTL working directory at any stage of the pipeline, including
# the inputs prepared by psQTL_prep.py and each analysis run produced by psQTL_proc.py
# (and post-processed by psQTL_post.py). Nothing in the working directory is modified.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.locations import Locations
from modules.parameters import ParameterCache
from modules.parsing import parse_metadata
from modules.runs import get_run, find_matching_run, list_runs, describe_runs, parse_run_id
from modules.viewing import view_inputs, view_runs, view_run_details, get_input_samples
from _version import __version__

def main():
    usage = """%(prog)s presents an overview of a psQTL working directory. By default, it
    shows the input files prepared by psQTL_prep.py alongside a summary of each analysis
    run, where every run holds the psQTL_proc.py statistics and psQTL_post.py outputs
    generated for one particular metadata grouping. Specify --run or --meta to instead
    see the samples, statistics, and outputs of a single run.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="workingDirectory",
                   required=True,
                   help="Specify the location where the analysis is being performed")
    p.add_argument("--run", dest="runID", metavar="RUN",
                   required=False,
                   help="""Optionally, specify an analysis run to view in detail by its
                   number (e.g., '2' or 'run2')""",
                   default=None)
    p.add_argument("--meta", dest="runMetadataFile", metavar="METADATAFILE",
                   required=False,
                   help="""Optionally, specify a metadata file to view the analysis run
                   whose groups match it in detail""",
                   default=None)
    p.add_argument("-v", "--version",
                   action="version",
                   version="psQTL_view.py {version}".format(version=__version__))

    args = p.parse_args()
    locations = Locations(args.workingDirectory) # validates directory existence
    if args.runID != None and args.runMetadataFile != None:
        raise ValueError("--run and --meta cannot be specified together; provide one or the other")

    paramsCache = ParameterCache(locations.workingDirectory)
    paramsCache.merge(args) # raises FileNotFoundError if cache does not exist

    # View a single run in detail
    if args.runID != None or args.runMetadataFile != None:
        print("## psQTL_view.py - View Analysis Run ##")
        if args.runID != None:
            run = get_run(locations.workingDirectory, parse_run_id(args.runID))
        else:
            if not os.path.isfile(args.runMetadataFile):
                raise FileNotFoundError(f"--meta file '{args.runMetadataFile}' does not exist or is not a file!")
            run = find_matching_run(locations.workingDirectory, parse_metadata(args.runMetadataFile))
            if run == None:
                raise FileNotFoundError(f"No analysis run matches the groups in metadata file '{args.runMetadataFile}'.\n" +
                                        describe_runs(list_runs(locations.workingDirectory)))
        view_run_details(locations.workingDirectory, run, get_input_samples(locations.workingDirectory))
    # View the whole working directory
    else:
        print("## psQTL_view.py - View Working Directory ##")
        view_inputs(locations, args)
        view_runs(locations.workingDirectory, args.metadataFile, "psQTL_view.py")

    print("Program completed successfully!")

if __name__ == "__main__":
    main()
