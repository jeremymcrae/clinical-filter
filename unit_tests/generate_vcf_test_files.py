""" Generates VCF files for testing variant filtering outputs.
"""

import os
import sys
import sys
import copy
import itertools

class writeVcf(object):
    """ creates a test vcf file
    """
    
    def __init__(self, sample_IDs, filename):
        """ initiate the vcf object, by making a VCF header and opening a file
        
        Args:
            sample_IDs: list of sample IDs for the individual/s
            filename: path to write VCF file to
        """
        
        self.sample_IDs = sample_IDs
        self.filename = filename
        
        # make a complete (but slim) VCF header
        header = self.make_vcf_header(self.sample_IDs)
        
        # write the header to a file
        self.output = open(self.filename, 'w')
        self.output.writelines(header)
    
    def make_vcf_header(self, IDs):
        """ makes a VCF header
        
        Args:
            IDs: the IDs for the indivdual/s in the VCF file
        
        Returns:
            list of VCF header lines
        """
        
        header = []
        
        header.append("##fileformat=VCFv4.1\n")
        header.append("##source=DDD clinical filtering unit testing\n")
        header.append("##ALT=<ID=DEL,Description=\"Deletion\">\n")
        header.append("##FILTER=<ID=temp_filtername,Description=\"Filler field for FILTER\">\n")
        header.append("##INFO=<ID=AMR_AF,Number=.,Type=Float,Description=\"Allele Frequency for AMR samples\">\n")
        header.append("##INFO=<ID=ASN_AF,Number=.,Type=Float,Description=\"Allele Frequency for ASN samples\">\n")
        header.append("##INFO=<ID=AFR_AF,Number=.,Type=Float,Description=\"Allele Frequency for AFR samples\">\n")
        header.append("##INFO=<ID=EUR_AF,Number=.,Type=Float,Description=\"Allele Frequency for EUR samples\">\n")
        header.append("##INFO=<ID=MAX_AF,Number=.,Type=Float,Description=\"Maximum Allele Frequency across 1000G, UK10K, ESP and DDD MAFs\">\n")
        header.append("##INFO=<ID=UK10K_cohort_AF,Number=.,Type=Float,Description=\"Maximum Allele Frequency from UK10K twins\">\n")
        header.append("##INFO=<ID=ESP_AF,Number=.,Type=Float,Description=\"Maximum Allele Frequency from the Exome Sequencing Project (Washington edu)\">\n")
        header.append("##INFO=<ID=DDD_AF,Number=.,Type=Float,Description=\"DDD internal allele frequency.\">\n")
        header.append("##INFO=<ID=CQ,Number=1,Type=String,Description=\"Highest consequences (from ensembl VEP)\">\n")
        header.append("##INFO=<ID=DENOVO-INDEL,Number=0,Type=Flag,Description=\"De-novo indel\">\n")
        header.append("##INFO=<ID=DENOVO-SNP,Number=0,Type=Flag,Description=\"De-novo SNP\">\n")
        header.append("##INFO=<ID=ENST,Number=1,Type=String,Description=\"Transcript id (from ensembl VEP)\">\n")
        header.append("##INFO=<ID=HGNC,Number=1,Type=String,Description=\"HGNC gene identifer (from ensembl VEP)\">\n")
        header.append("##INFO=<ID=PolyPhen,Number=1,Type=String,Description=\"PolyPhen prediction (from ensembl VEP)\">\n")
        header.append("##INFO=<ID=SIFT,Number=1,Type=String,Description=\"SIFT prediction (from ensembl VEP)\">\n")
        header.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        header.append("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality\">\n")
        header.append("##FORMAT=<ID=TEAM29_FILTER,Number=1,Type=String,Description=\"either PASS or the name of the filter that failed: AF_MAX (Population MAF > 1%); inVCF (Variant not present in child VCF, or is present in parent VCF); maxAltInParentFlag (Maximum alternate frequency in parent > 10%); segmentaldup (Variant overlaps segmental duplication); TRF (Variant overlaps tandem repeat)\">\n")
        header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(IDs) + "\n")
        
        return header
    
    def write(self, vcf_entry):
        """ writes a vcf variant entry
        
        Args:
            vcf_entry: dict for VCF line
        """
        
        # write a single vcf variant
        vcf_line = [vcf_entry["CHROM"], vcf_entry["POS"], vcf_entry["ID"], \
                    vcf_entry["REF"], vcf_entry["ALT"], vcf_entry["QUAL"], \
                    vcf_entry["FILTER"]]
        
        # add in the INFO entries
        info = []
        for key in vcf_entry["INFO"]:
            if vcf_entry["INFO"][key] == True:
                info.append(key)
            else:
                info.append(key + "=" + str(vcf_entry["INFO"][key]))
        vcf_line.append(";".join(info))
        
        # now add the FORMAT entries
        format = list(vcf_entry["FORMAT"].keys())
        samples = None
        for key in format:
            
            for position in range(len(vcf_entry["FORMAT"][key])):
                if samples == None:
                    samples = [[]] * len(vcf_entry["FORMAT"][key])
                value = vcf_entry["FORMAT"][key][position]
                
                if value == None:
                    value == ""
                samples[position].append(value)
            
        vcf_line.append(":".join(format))
        
        samples = list(map(":".join, samples))
        vcf_line.append("\t".join(samples))
        
        self.output.write("\t".join(vcf_line) + "\n")


def define_vep_consequences():
    """ define variant effect predictor (VEP) consequences for variants.
    
    From: http://www.ensembl.org/info/genome/variation/predicted_data.html
    All we really need is the list of VEP terms (eg stop_gained), the
    definitions and sequence ontology IDs for the terms are given at the
    above website.
    
    Returns:
        list of sequence ontology terms
    """
    
    vep_consequences = ["transcript_ablation", "splice_donor_variant", \
        "splice_acceptor_variant", "stop_gained", "frameshift_variant", \
        "stop_lost", "initiator_codon_variant", "inframe_insertion", \
        "inframe_deletion", "missense_variant", "transcript_amplification", \
        "splice_region_variant", "incomplete_terminal_codon_variant", \
        "synonymous_variant", "stop_retained_variant",\
        "coding_sequence_variant", "mature_miRNA_variant", 
        "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_exon_variant", \
        "nc_transcript_variant", "intron_variant", "NMD_transcript_variant", \
        "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", \
        "TFBS_amplification", "TF_binding_site_variant", \
        "regulatory_region_variant", "regulatory_region_ablation", \
        "regulatory_region_amplification", "feature_elongation", \
        "feature_truncation", "intergenic_variant"]
    
    return vep_consequences

def create_default_vcf_dict(position, ID):
    """ creates a dict for a default VCF entry, which passes filters
    
    Args:
        position: nucleotide position for the variant
    
    Returns:
        a dictionary for the VCF variant
    """
    
    position = str(position)
    
    default = {}
    
    default["CHROM"] = "0"
    default["POS"] = position
    default["ID"] = ID
    default["REF"] = "C"
    default["ALT"] = "T"
    default["QUAL"] = "50"
    default["FILTER"] = "PASS"
    default["INFO"] = {}
    default["FORMAT"] = {}
    
    # now set up the default INFO fields
    pops = ["AMR_AF", "ASN_AF", "AFR_AF", "EUR_AF", "MAX_AF", \
            "UK10K_cohort_AF", "ESP_AF", "DDD_AF"]
    for pop in pops:
        default["INFO"][pop] = 0.001
    
    # define polyphen and sift values, even though we currently don't use them
    polyphen_types = ["benign", "unknown", "possibly_damaging", \
                      "probably_damaging"]
    sift_types = ["tolerated", "deleterious"]
    
    default["INFO"]["CQ"] = "missense_variant"
    default["INFO"]["ENST"] = ID
    default["INFO"]["HGNC"] = "TEMPGENENAME_" + position + "_" + ID
    default["INFO"]["PolyPhen"] = "probably_damaging(0.998)"
    default["INFO"]["SIFT"] = "deleterious(0)"
    
    default["FORMAT"]["GT"] = ["0/0"]
    default["FORMAT"]["GQ"] = ["50"]
    
    return default

class createTestVcfs(object):
    """ creates VCF files for all filters and genotypes combinations
    
    We need to test a variety of scenarios, from checking that the filtering
    for functional variants, or low MAF variants, or ones passing the FILTER
    field. We also need to test different genotype combinations across the
    family trio, as well as the compound het possibilities.
    
    Each of these types of test is kept separate from the others. I've grouped
    the types of test by nucleotide position, so we can easily find why
    variants do not pass. Note that most positions in these ranges are unused,
    and that we start at 10,000,001 in order to avoid the pseudoautosomal
    regions on the X chrom:
        functional annotation filter:  10000001 - 10009999
        FILTER field values:           10010000 - 10019999
        minor allele frequency filter: 10020000 - 10029999
        single nucleotide genotypes:   10030000 - 10039999
        compound het genotypes:        10040000 - 10049999
    
    The functional annotation, FILTER field and MAF filtering tests use trio 
    genotypes from de novo variants, which should be included in reporting 
    output. The only reason why they are not included is because they do not
    pass the filtering criteria. In contrast, the tests that examine different
    genotype combinations have default filter values that should pass the 
    variant for examination, and only the inheritance checks determine whether
    the variant is included in the reporting output.
    
    Each test is partially identified by the position string, eg 
    20005_MAF:EUR_AF (variant 20005, examining the MAF, specifically the EUR 
    MAF), or 30015_snv:121 (variant 30015, examining a single nucleotide
    variant, with genotype of 1, 2, 1 for child, mother, father.)
    Each of the filtering and snv tests are run on independent (fake) genes, to
    prevent compound het checks. The compound hets are set up as pairs of 
    variants on separate genes, some of which will be reported under single 
    variant checks.
    """
    
    def __init__(self):
        """ initiate the class by starting VCF files for family members
        """
        
        self.maf_pops = ["AMR_AF", "ASN_AF", "AFR_AF", "EUR_AF", "MAX_AF", \
                         "UK10K_cohort_AF", "ESP_AF", "DDD_AF"]
        self.cur_pos = 10000000
        
        # start VCF files for family members
        self.child_vcf = writeVcf(["child"], "child.vcf")
        self.mother_vcf = writeVcf(["mother"], "mother.vcf")
        self.father_vcf = writeVcf(["father"], "father.vcf")
        
        # create variants across different scenarios
        self.vary_consequences()
        self.vary_filter()
        self.vary_maf_populations()
        self.vary_snvs()
        self.vary_compound_snvs()
    
    def round_position_up(self):
        """ group variants with similar filtering purposes by position
        """
        
        increment = 10000
        
        if self.cur_pos % increment == 0:
            pass
        else:
            self.cur_pos = self.cur_pos + increment - self.cur_pos % increment
    
    def write_filtering_varint(self, vcf_entry):
        """ make a variant for filtering purposes
        """
        
        self.mother_vcf.write(vcf_entry)
        self.father_vcf.write(vcf_entry)
        
        # now make the child a de novo, to pass all inheritance checks
        vcf_entry["FORMAT"]["GT"] = ["0/1"]
        vcf_entry["INFO"]["DENOVO-SNP"] = True
        vcf_entry["INFO"]["TEAM29_FILTER"] = "PASS"
        self.child_vcf.write(vcf_entry)
        
        self.cur_pos += 1
    
    def vary_consequences(self):
        """ create variants that cover all VEP consequences
        """
        
        vep_consequences = define_vep_consequences()
        
        for cq in sorted(vep_consequences):
            entry = create_default_vcf_dict(self.cur_pos,  "consequence:" + cq)
            entry["INFO"]["CQ"] = cq
            self.write_filtering_varint(entry)
        
        self.round_position_up()
    
    def vary_filter(self):
        """ create variants that fail for different FILTER values
        """
        
        # list of filter values that a) pass ("PASS"), b) explicitly fail 
        # ("FAIL"), c) implicitly fail ("gtak_pl"), d) fail due to null value 
        # ("."), e) fail due to no value (""), and f) fail even though "PASS" 
        # is part of the value ("NOT_PASSED")
        filter_types = ["PASS", "FAIL", "gtak_pl", ".", "", "NOT_PASSED"]
        
        for filt in sorted(filter_types):
            entry = create_default_vcf_dict(self.cur_pos, "filter:" + filt)
            entry["FILTER"] = filt
            self.write_filtering_varint(entry)
        
        self.round_position_up()
    
    def vary_maf_populations(self):
        """ create variants that fail for each MAF population
        """
        
        # check that we catch each population failing the MAF
        for pop in sorted(self.maf_pops):
            entry = create_default_vcf_dict(self.cur_pos, "MAF:" + pop)
            entry["INFO"][pop] = 0.1
            self.write_filtering_varint(entry)
        
        # and make an entry that lacks any MAF info for the populations
        entry = create_default_vcf_dict(self.cur_pos, "MAF:none")
        for pop in sorted(self.maf_pops):
            del entry["INFO"][pop]
        
        self.write_filtering_varint(entry)
        
        # and make an entry that fails the de novo filter
        entry = create_default_vcf_dict(self.cur_pos, "fail_de_novo")
        self.mother_vcf.write(entry)
        self.father_vcf.write(entry)
        entry["FORMAT"]["GT"] = ["0/1"]
        self.child_vcf.write(entry)
        
        self.round_position_up()
    
    def write_trio_snvs(self, vcf_entry, genotypes):
        """ write vcf lines for each family member, given trio genotypes
        
        Args:
            vcf_entry: vcf dictionary for a variant
            genotypes: list of genotypes for a trio
        """
        
        child = copy.deepcopy(vcf_entry)
        mother = copy.deepcopy(vcf_entry)
        father = copy.deepcopy(vcf_entry)
        
        child["FORMAT"]["GT"] = [genotypes[0]]
        mother["FORMAT"]["GT"] = [genotypes[1]]
        father["FORMAT"]["GT"] = [genotypes[2]]
        
        # annotate de novos as passing filters
        if genotypes == ("0/1", "0/0", "0/0"):
            child["INFO"]["DENOVO-SNP"] = True
            child["INFO"]["TEAM29_FILTER"] = "PASS"
        # and allow for male X-chrom de novos
        if child["CHROM"] == "X" and genotypes == ("1/1", "0/0", "0/0"):
            child["INFO"]["DENOVO-SNP"] = True
            child["INFO"]["TEAM29_FILTER"] = "PASS"
        
        self.child_vcf.write(child)
        self.mother_vcf.write(mother)
        self.father_vcf.write(father)
        
        self.cur_pos += 1
    
    def convert_geno_list_to_string(self, genotypes):
        """ converts genotype list to string with single char values eg "210"
        
        Args:
            genotypes: tuple of vcf encoded genotypes eg ("0/0", "0/1", "1/1")
        
        Returns:
            string with each genotype as single characters of alt allele count
        """
        
        geno_dict = {"0/0": "0", "0/1": "1", "1/1": "2"}
        
        string = ""
        
        for geno in genotypes:
            string += geno_dict[geno]
        
        return string
    
    def vary_snvs(self):
        """ create all possible genotypes amongst the family
        """
        
        genotypes = ["0/0", "0/1", "1/1"]
        combos = itertools.product(genotypes, repeat=3)
        for combo in combos:
            # don't use ref genotype "0/0" for the child
            if combo[0] == "0/0":
                continue
            
            geno = self.convert_geno_list_to_string(combo)
            
            entry = create_default_vcf_dict(self.cur_pos, "snv:" + geno)
            self.write_trio_snvs(entry, combo)
            
            # and do a set for the x-chromosome
            entry = create_default_vcf_dict(self.cur_pos, "snv:" + geno)
            entry["CHROM"] = "X"
            self.write_trio_snvs(entry, combo)
        
        self.round_position_up()
    
    def vary_compound_snvs(self):
        """ create all possible compound het scenarios amongs the family
        """
        
        genotypes = ["0/0", "0/1", "1/1"]
        first_combos = itertools.product(genotypes, repeat=3)
        for first in first_combos:
            second_combos = itertools.product(genotypes, repeat=3)
            for second in second_combos:
                # don't use ref genotype "0/0" for the child
                if first[0] == "0/0" or second[0] == "0/0":
                    continue
                
                geno1 = self.convert_geno_list_to_string(first)
                geno2 = self.convert_geno_list_to_string(second)
                
                first_entry = create_default_vcf_dict(self.cur_pos, "compound:" + geno1 + "_" + geno2)
                self.write_trio_snvs(first_entry, first)
                
                # both variants need to be on the same gene
                second_entry = create_default_vcf_dict(self.cur_pos, "compound:" + geno1 + "_" + geno2)
                second_entry["INFO"]["HGNC"] = "TEMPGENENAME_" + str(self.cur_pos - 1) + "_compound:"+ geno1 + "_" + geno2
                self.write_trio_snvs(second_entry, second)
                
                # and do a set for the x-chromosome
                first_entry = create_default_vcf_dict(self.cur_pos, "compound:" + geno1 + "_" + geno2)
                first_entry["CHROM"] = "X"
                self.write_trio_snvs(first_entry, first)
                
                # both variants need to be on the same gene
                second_entry = create_default_vcf_dict(self.cur_pos, "compound:" + geno1 + "_" + geno2)
                second_entry["CHROM"] = "X"
                second_entry["INFO"]["HGNC"] = "TEMPGENENAME_" + str(self.cur_pos - 1)  + "_compound:"+ geno1 + "_" + geno2
                self.write_trio_snvs(second_entry, second)
        
        self.round_position_up()


def main():
    createTestVcfs()

if __name__ == '__main__':
    main()


