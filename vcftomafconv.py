import pandas as pd

# Input and output files
vcf_file = "Annotated.vcf"
maf_file = "Annotated_with_Classificationtest1.maf"

# Read VCF, skipping metadata lines (starting with #)
vcf = pd.read_csv(vcf_file, sep="\t", comment="#", header=None)

# Print the number of columns in the VCF to debug
print(f"VCF file has {vcf.shape[1]} columns")

# Adjust the vcf_columns list to match the actual column count
vcf_columns = [
    "Chromosome", "Start_Position", "ID", "Reference_Allele",
    "Tumor_Seq_Allele2", "Quality", "Filter", "Info"
]

# If the VCF has more columns (e.g., sample-specific fields), add them
if vcf.shape[1] > len(vcf_columns):
    vcf_columns.extend([f"Extra_Column_{i}" for i in range(len(vcf_columns), vcf.shape[1])])

# Assign column names to the VCF
vcf.columns = vcf_columns[:vcf.shape[1]]

def classify_variant(ref, alt, info):
    """Classify the variant type based on reference and alternate alleles."""
   
 #   print(f"Checking variant: REF={ref}, ALT={alt}")  # Debugging print to show REF and ALT values
#    print(f"REF length: {len(ref)}, ALT length: {len(alt)}")  # Debugging print for lengths of ref and alt

    if len(ref) < len(alt):  # Insertion
        if "frameshift" in info.lower():
            return "Frame_Shift_Ins"
        return "In_Frame_Ins"
    elif len(ref) > len(alt):  # Deletion
        if "frameshift" in info.lower():
            return "Frame_Shift_Del"
        return "In_Frame_Del"
    elif len(ref) == len(alt):  # SNP or substitution
        if "ANN=" in info.lower():
            print("Identified SNP, checking ANN field...")  # Debugging print to confirm entering SNP condition
            print(f"ANN Field: {info}")  # Debugging print to show the full ANN field
            ann_entries = info.split("ANN=")[1].split(",")
            for ann_entry in ann_entries:
                ann_fields = ann_entry.split("|")
                if len(ann_fields) > 1:
                    variant_type = ann_fields[1].lower()  # Mutation type in the 2nd field
                    print(f"Variant type from ANN: {variant_type}")  # Print variant type for debugging
                    # Check for missense, nonsense, synonymous, etc.
                    if "missense_variant" in variant_type:
                        return "Missense_Mutation"
                    elif "nonsense_variant" in variant_type:
                        return "Nonsense_Mutation"
                    elif "synonymous_variant" in variant_type:
                        return "Silent"
                    elif "splice_site_variant" in variant_type:
                        return "Splice_Site"
                    elif "intron_variant" in variant_type:
                        return "Intron_Variant"
            return "missense_variant"  # If no known mutation type is found in the annotations
        else:
            return "Missense_Mutation"  # If no ANN field is present



# Function to extract Hugo_Symbol from the ANN field
def extract_hugo_symbol(info_field):
    """Extract the Hugo_Symbol from the ANN field in VCF."""
    if "ANN=" in info_field:
        ann_entries = info_field.split("ANN=")[1].split(",")
        for ann_entry in ann_entries:
            ann_fields = ann_entry.split("|")
            if len(ann_fields) > 3:
                return ann_fields[3]  # Gene name is typically in the 4th field
    return "Unknown"  # Default if Hugo_Symbol is not found

# Extract relevant MAF columns and apply Hugo_Symbol extraction and classification
maf_data = pd.DataFrame({
    "Hugo_Symbol": vcf["Info"].apply(extract_hugo_symbol),  # Hugo_Symbol from ANN field
    "Chromosome": vcf["Chromosome"],
    "Start_Position": vcf["Start_Position"],
    "End_Position": vcf["Start_Position"],  # Default for SNPs; will adjust for indels below
    "Reference_Allele": vcf["Reference_Allele"],
    "Tumor_Seq_Allele1": vcf["Reference_Allele"],  # For indels, this may remain empty or match Ref
    "Tumor_Seq_Allele2": vcf["Tumor_Seq_Allele2"]  # Alt allele; for insertions, this contains the inserted sequence
})

# Adjust Start_Position and End_Position for indels
def adjust_positions(ref, alt, start):
    if len(ref) < len(alt):  # Insertion
        return start, start
    elif len(ref) > len(alt):  # Deletion
        return start, start + len(ref) - 1
    else:  # SNP or substitution
        return start, start

maf_data[["Start_Position", "End_Position"]] = maf_data.apply(
    lambda row: adjust_positions(row["Reference_Allele"], row["Tumor_Seq_Allele2"], row["Start_Position"]),
    axis=1, result_type="expand"
)

# Add Variant_Classification based on the reference, alternate alleles, and INFO field
maf_data["Variant_Classification"] = maf_data.apply(
    lambda row: classify_variant(row["Reference_Allele"], row["Tumor_Seq_Allele2"], row.get("info", "")),
    axis=1
)

# Add additional MAF columns
maf_data["Variant_Type"] = maf_data.apply(
    lambda row: "INS" if len(row["Reference_Allele"]) < len(row["Tumor_Seq_Allele2"]) else
                "DEL" if len(row["Reference_Allele"]) > len(row["Tumor_Seq_Allele2"]) else
                "SNP",
    axis=1
)
maf_data["Tumor_Sample_Barcode"] = "Cancer"  # Replace with actual sample barcode
maf_data["Matched_Norm_Sample_Barcode"] = "Normal"  # Replace if applicable

# Save to MAF format
maf_data.to_csv(maf_file, sep="\t", index=False)

print(f"MAF file saved to {maf_file}")

