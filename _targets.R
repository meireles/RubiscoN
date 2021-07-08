################################################################################
# Setup
################################################################################

########################################
# Load libraries
########################################

library("targets")     # manages the whole thing
library("visNetwork")  # Explore target pipeline status
library("conflicted")  # complain if function names are duplicated

library("httr")        # talk to the web
library("vroom")       # read tables. fast
library("readr")       # table IO
library("zip")         # cross platform zipping
library("stringr")     # deal with strings

library("taxize")      # taxonomic cleanup
library("rgnparser")   # parse scientific names
library("ape")         # phylogeny stuff
library("Biostrings")  # deal with alignment file
library("rentrez")     # interact with ncbi
library("CHNOSZ")      # parse chemical formulas
library("kmer")        # kmer aa counting


########################################
# Load functions
########################################

invisible(sapply(dir("R/data_processing/", full.names = TRUE), source))
invisible(sapply(dir("R/analyses/", full.names = TRUE), source))

########################################
# Set options
########################################

# Not including packages here. It will be a pain to have to change things just
# because of a minor package update.
tar_option_set(garbage_collection = TRUE)

################################################################################
# Plan
################################################################################

list(
    ############################################################################
    # Find Data
    ############################################################################

    ## rubisco
    tar_target(rubisco_uniprot_raw_data_file,
               download_uniprot_raw_data(path_out   = "raw_data/UniProt/rubisco_uniprot_raw.tsv",
                                         gene       = "rbcl",
                                         min_length = 400,
                                         max_length = 505),
               format = "file"),

    ## TRY
    tar_target(try_raw_data_file,
               "raw_data/TRY/13403.zip",
               format = "file"),


    ## phylogeny
    tar_target(phy_raw_data_file,
               "raw_data/phylo/v0.1/ALLMB.tre",
               format = "file"),

    ############################################################################
    # Read data and update taxonomy
    ############################################################################

    ## rubisco
    tar_target(rubisco_raw_data,
               read.delim(file = rubisco_uniprot_raw_data_file,
                          check.names = FALSE)),

    tar_target(rubisco_species_lookup,
               get_accepted_plant_name(x = rubisco_raw_data$Organism)),

    tar_target(rubisco_newname,
               update_names(x        = rubisco_raw_data,
                            sp_field = "Organism",
                            lookup   = rubisco_species_lookup)),

    tar_target(rubisco_raw_aa_fasta,
               parse_uniprot_fasta(x        = rubisco_newname,
                                   path_out = "raw_data/UniProt/rubisco_raw_aa.fasta"),
               format = "file"),

    tar_target(rubisco_raw_aa_features,
               parse_uniprot_features(x        = rubisco_newname,
                                      path_out = "raw_data/UniProt/rubisco_raw_aa_features.tsv"),
               format = "file"),


    ## TRY
    tar_target(try_raw_data,
               read_and_trim_raw_try_data(path_in = try_raw_data_file)),

    tar_target(try_species_lookup,
               get_accepted_plant_name(x = try_raw_data$AccSpeciesName)),

    tar_target(try_newname,
               update_names(x        = try_raw_data,
                            sp_field = "AccSpeciesName",
                            lookup   = try_species_lookup)),

    ## phylogeny
    tar_target(phy_raw_data,
               ape::read.tree(phy_raw_data_file)),

    tar_target(phy_species_lookup,
               get_accepted_plant_name(x = gsub(pattern     = "_",
                                                replacement = " ",
                                                x           = phy_raw_data$tip.label),
                                       chunk = 5000,
                                       sleep = 2)),

    tar_target(phy_newname,
               update_names(x        = phy_raw_data,
                            sp_field = "tip.label",
                            lookup   = phy_species_lookup)),

    ############################################################################
    # Process data
    ############################################################################

    ## rubisco
    tar_target(rubisco_scrubbed_aa_fasta,
               initial_seq_cleanup(path_in       = rubisco_raw_aa_fasta,
                                   path_out      = "data/rubisco/rubisco_scrubbed_aa.fasta",
                                   uncertain_tol = 0.05,
                                   sigma         = 7,
                                   min_length    = 400),
               format = "file"),

    tar_target(rubisco_aligned_aa_fasta,
               mafft(rubisco_scrubbed_aa_fasta,
                     "data/rubisco/rubisco_aligned_aa.fasta"),
               format = "file"),

    tar_target(rubisco_trimmed_alignment_aa,
               trim_alignment_ends(path_in  = rubisco_aligned_aa_fasta,
                                   path_out = "data/rubisco/rubisco_trimmed_alignment_aa.fasta",
                                   prop_gap = 0.1),
               format = "file"),


    tar_target(rubisco_trimmed_alignment_aa_by_sp,
               clean_tax_and_pick_long_seq(path_in = rubisco_trimmed_alignment_aa,
                                           path_out = "data/rubisco/rubisco_trimmed_alignment_aa_by_sp.fasta"),
               format = "file"),

    ## TRY
    tar_target(try_combined_traits_by_sp,
               aggregate_traits_by_sp(try_newname)),

    ## phylogeny
    tar_target(phy_sp_only,
               remove_genera_and_duplicated_tips(phy_newname)),

    ######################################
    # extract data from rubisco
    ######################################

    tar_target(rubisco_aa_composition,
               get_aa_composition(path_in = rubisco_trimmed_alignment_aa_by_sp)),


    tar_target(rubisco_element_composition,
               get_element_composition(path_in = rubisco_trimmed_alignment_aa_by_sp)),


    ######################################
    # Intersect datasets
    ######################################

    tar_target(sp_intersect,
               Reduce(intersect, list(phy_sp_only$tip.label,
                                      try_combined_traits_by_sp$AccSpeciesName,
                                      rubisco_aa_composition$organism,
                                      rubisco_element_composition$organism))),

    tar_target(phy_sp_only_intersect,
               subset_phylo(tree = phy_sp_only, sp_keep = sp_intersect, brlens_tol = 1e-9)),

    tar_target(try_combined_traits_by_sp_intersect,
               {try_combined_traits_by_sp[try_combined_traits_by_sp$AccSpeciesName %in% sp_intersect, ]}),

    tar_target(rubisco_aa_composition_intersect,
               {rubisco_aa_composition[rubisco_aa_composition$organism %in% sp_intersect, ]}),

    tar_target(rubisco_element_composition_intersect,
               {rubisco_element_composition[rubisco_element_composition$organism %in% sp_intersect, ]}),


    ######################################
    # Repackage datasets
    ######################################

    tar_target(final_data,
               split_data_by_clade_and_repackage(phy   = phy_sp_only_intersect,
                                                 try   = try_combined_traits_by_sp_intersect,
                                                 aa    = rubisco_aa_composition_intersect,
                                                 ele   = rubisco_element_composition_intersect))

    ############################################################################
    # Analyses
    ############################################################################

    # TODO
    # I'm keeping some the exploratory code in the kitchensink diretory

)
