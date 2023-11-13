import pandas as pd
import os
import argparse


"""
checkv_filter Module

This module contains functions to filter a CheckV report based on completeness percent and write selected contig nodes to a file.

Functions:
- filter_checkv_completeness(checkv_report, selected_nodes_filename, completeness_cutoff=50): Filter CheckV report based on completeness.
- main(): Entry point for the script to filter CheckV report.
"""

def filter_checkv_completeness(checkv_report, selected_nodes_filename, completeness_cutoff = 50):
    """
    Filter a CheckV report based on completeness percent and write selected contig nodes to a file.

    Parameters:
    checkv_report (str): Path to the CheckV report in TSV format.
    selected_nodes_filename (str): Path to the file to save selected contig nodes.
    completeness_cutoff (int, optional): The threshold of the CheckV completeness score (default is 50).

    Returns:
    None: This function writes selected contig nodes to the specified file.
    """

    checkvReportDir = checkv_report

    checkvReport = pd.read_csv(checkvReportDir, sep="\t")


    completeness = checkvReport.query("completeness > @completeness_cutoff")

    with open(f"{selected_nodes_filename}", "w") as fn:
        for node in completeness["contig_id"]:
            fn.write(f"{node}\n")
    return None


def main():

    #parse paths
    parser = argparse.ArgumentParser()
    parser.add_argument("checkv_report_path", help="path to the checkv report")
    parser.add_argument("completeness_cutoff", help="threshold of the checkv completeness score", type=int)
    parser.add_argument("selected_nodes_filename", help="path to a file with contig headers (nodes) selected based on the specified completeness cutoff")
    args = parser.parse_args()


    checkv_report_path = args.checkv_report_path
    completeness_cutoff = args.completeness_cutoff
    selected_nodes_filename = args.selected_nodes_filename

    #run the
    filter_checkv_completeness(checkv_report_path, selected_nodes_filename, completeness_cutoff)



if __name__ == "__main__":
    main()

