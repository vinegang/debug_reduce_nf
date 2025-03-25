#!/bin/env python3

import csv
from collections import defaultdict
import sys
import os

# This script takes two inputs
# 1. input samplesheet
# 2. path to output directory

samplesheet = sys.argv[1]
outdir = sys.argv[2]
# Define the expected column names
expected_columns = [
    "sample",
    "library",
    "read1",
    "read2",
    "sample_captures",
    "Diagnosis",
    "Matched_RNA",
    "Matched_normal",
    "casename",
    "type",
]

# Read the sample sheet file
with open(samplesheet, "r") as file:
    reader = csv.DictReader(file)
    columns = reader.fieldnames

    # Check if all expected columns are present
    if not all(col in columns for col in expected_columns):
        print("Error: Missing columns in the sample sheet.")
        exit(0)

    # Store the contents of the sample sheet in a list
    samplesheet_data = list(reader)

for row in samplesheet_data:
    if row["type"] == "xeno_DNA":
        row["type"] = "tumor_DNA"

for row in samplesheet_data:
    if row["type"] == "xeno_RNA":
        row["type"] = "tumor_RNA"

# Remove spaces from columns that don't have any values
for row in samplesheet_data:
    for col in columns:
        if not row[col].strip():  # Check if the column has no values
            row[col] = row[col].strip()

# Create a list to store the matching rows
matching_rows = []
tumor_rnaseq_rows = []
tumor_rnaseq_normal_rows = []
tumor_normal_rows = []

# Create a dictionary to track the occurrences of each sample and casename combination
sample_casename_counts = defaultdict(int)

# Create a set to track processed samples to avoid duplicate entries
processed_samples = set()

for row in samplesheet_data:
    sample = row["sample"]
    casename = row["casename"]
    # Get the values of Matched_RNA and Matched_normal
    matched_rna = row["Matched_RNA"]
    matched_normal = row["Matched_normal"]
    matched_casename = row["casename"]
    row_type = row["type"]
    

    # Increment the count for the current sample and casename combination
    sample_casename_counts[(sample, casename)] += 1

    # Check if the row is a Tumor type
    if row["type"] in ["tumor_DNA", "cell_line_DNA"]:
        # Check additional conditions
        if matched_rna and not matched_normal:
            tumor_rnaseq_rows.append(row)
            # Find other rows with the same library as matched RNA
            for other_row in samplesheet_data:
                if (
                    other_row["library"] == matched_rna
                    and other_row["casename"] == matched_casename
                ):
                    tumor_rnaseq_rows.append(other_row)
                    processed_samples.add(other_row["sample"])

        elif matched_rna and matched_normal:
            tumor_rnaseq_normal_rows.append(row)
            # Find other rows with the same library as matched normal and matched RNA
            for other_row in samplesheet_data:
                if (matched_rna and other_row["library"] == matched_rna) or (
                    matched_normal and other_row["library"] == matched_normal
                ):
                    if other_row["casename"] == matched_casename:
                        tumor_rnaseq_normal_rows.append(other_row)
                        processed_samples.add(other_row["sample"])

        elif matched_normal and not matched_rna:
            tumor_normal_rows.append(row)
            # Find other rows with the same library as matched normal
            for other_row in samplesheet_data:
                if (
                    other_row["library"] == matched_normal
                    and other_row["casename"] == matched_casename
                ):
                    tumor_normal_rows.append(other_row)
                    processed_samples.add(other_row["sample"])
    elif row_type == "normal_DNA":
        patient_id = row["sample"]
        has_tumor_rna = False
        has_tumor_dna = False

        # Check all entries for this patient and casename
        for other_row in samplesheet_data:
            if other_row["sample"] == patient_id and other_row["casename"] == casename:
                if other_row["type"] in ["tumor_RNA", "cell_line_RNA"]:
                    has_tumor_rna = True
                if other_row["type"] == "tumor_DNA":
                    has_tumor_dna = True

        if has_tumor_rna and not has_tumor_dna:
            print(f"Appending Normal DNA {sample} (only Tumor RNA found, no Tumor DNA)")
            tumor_rnaseq_rows.append(row)
            for other_row in samplesheet_data:
                if (
                    other_row["sample"] == patient_id
                    and other_row["casename"] == casename
                    and other_row["type"] in ["tumor_RNA", "cell_line_RNA"]
                ):
                    print(f"Appending RNA sample {other_row['sample']} for {casename}")
                    tumor_rnaseq_rows.append(other_row)
                    processed_samples.add(other_row["sample"])

# Write the matching rows to the respective output files
if tumor_rnaseq_rows:
    tumor_rnaseq_rows = list(
        {tuple(row.items()): row for row in tumor_rnaseq_rows}.values()
    )
    output_file_rnaseq = os.path.join(outdir, "Tumor_RNAseq.csv")
    with open(output_file_rnaseq, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(tumor_rnaseq_rows)

if tumor_rnaseq_normal_rows:
    output_file_rnaseq_normal = os.path.join(outdir, "Tumor_RNAseq_Normal.csv")
    with open(output_file_rnaseq_normal, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(tumor_rnaseq_normal_rows)

if tumor_normal_rows:
    output_file_normal = os.path.join(outdir, "Tumor_Normal.csv")
    with open(output_file_normal, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(tumor_normal_rows)


# Create separate file writers for RNAseq, Normal, and Tumor types
rna_file = None
exome_file = None
normal_file = None
tumor_file = None
rna_writer = None
normal_writer = None
tumor_writer = None
tumor_lib_rows = []
RNAseq_lib_rows = []
rna_lib_rows = []

# Check if there are unique values for RNAseq type
if any(
    sample_casename_counts[(row["sample"], row["casename"])] == 1
    for row in samplesheet_data
    if row["type"] in ["tumor_RNA", "cell_line_RNA", "xeno_RNA"]
):
    output_file = os.path.join(outdir, "RNAseq.csv")
    rna_file = open(output_file, "w", newline="")
    rna_writer = csv.DictWriter(rna_file, fieldnames=columns)
    rna_writer.writeheader()

# Check if there are unique values for Normal type
# if any(
#     sample_casename_counts[(row["sample"], row["casename"])] == 1
#     for row in samplesheet_data
#     if row["type"] == "normal_DNA"
# ):
#     output_file = os.path.join(outdir, "Normal.csv")
#     normal_file = open(output_file, "w", newline="")
#     normal_writer = csv.DictWriter(normal_file, fieldnames=columns)
#     normal_writer.writeheader()

# # Check if there are unique values for Tumor type
# if any(
#     sample_casename_counts[(row["sample"], row["casename"])] == 1
#     for row in samplesheet_data
#     if row["type"] == "tumor_DNA"
# ):
#     output_file = os.path.join(outdir, "Tumor.csv")
#     tumor_file = open(output_file, "w", newline="")
#     tumor_writer = csv.DictWriter(tumor_file, fieldnames=columns)
#     tumor_writer.writeheader()

has_unique_normal = any(
    sample_casename_counts[(row["sample"], row["casename"])] == 1
    for row in samplesheet_data
    if row["type"] == "normal_DNA"
)

# Check if there are unique values for Tumor type
has_unique_tumor = any(
    sample_casename_counts[(row["sample"], row["casename"])] == 1
    for row in samplesheet_data
    if row["type"] in ["tumor_DNA", "cell_line_DNA"]
)

if has_unique_normal or has_unique_tumor:
    output_file = os.path.join(outdir, "Exome.csv")
    exome_file = open(output_file, "w", newline="")
    exome_writer = csv.DictWriter(exome_file, fieldnames=columns)
    exome_writer.writeheader()
    # exome_writer.writerows(row for row in samplesheet_data if row["type"] in ["normal_DNA", "tumor_DNA"])

    # Iterate over the sample sheet data again to write to respective files
for row in samplesheet_data:
    sample = row["sample"]
    casename = row["casename"]
    row_type = row["type"]
    library = row["library"]

    # Check if the sample and casename combination is not unique, type is Tumor, and Matched_RNA and Matched_normal are empty
    if (
        sample_casename_counts[(sample, casename)] > 1
        and row_type
        in [
            "tumor_DNA",
            "cell_line_DNA",
        ]
        and not row["Matched_RNA"]
        and not row["Matched_normal"]
    ):
        tumor_lib_rows.append(row)

    if sample_casename_counts[(sample, casename)] > 1 and row_type in [
        "tumor_RNA",
        "cell_line_RNA",
    ]:
        match_found = False
        has_normal_dna = False
        # Iterate over the rows again to find matches
        for other_row in samplesheet_data:
            if other_row["Matched_RNA"] == library:
                match_found = True
                break
            if (other_row["casename"] == casename and other_row["type"] == "normal_DNA"):
                has_normal_dna = True

        # If no match is found, add the current row to rna_lib_rows
        if not match_found and not has_normal_dna:
            rna_lib_rows.append(row)


if rna_lib_rows:
    output_file = os.path.join(outdir, "RNA_lib.csv")
    with open(output_file, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rna_lib_rows)


# Write the qualifying rows to tumor_lib.csv
if tumor_lib_rows:
    output_file = os.path.join(outdir, "Tumor_lib.csv")
    with open(output_file, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(tumor_lib_rows)


# Iterate over the sample sheet data again to write to respective files
for row in samplesheet_data:
    sample = row["sample"]
    casename = row["casename"]
    row_type = row["type"]

    # Check if the sample and casename combination is unique
    if sample_casename_counts[(sample, casename)] == 1:
        if row_type in ["tumor_RNA", "cell_line_RNA", "xeno_RNA"] and rna_writer:
            rna_writer.writerow(row)
        elif row_type in ["normal_DNA", "tumor_DNA", "cell_line_DNA"] and exome_writer:
            exome_writer.writerow(row)

# Close the file writers
if rna_file:
    rna_file.close()
# if normal_file:
#     normal_file.close()
# if tumor_file:
#     tumor_file.close()
if exome_file:
    exome_file.close()
