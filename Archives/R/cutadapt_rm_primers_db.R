#' Remove pairs of primers and flanking region from a fasta file using [cutadapt](https://github.com/marcelm/cutadapt/)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
# Function to cutadapt a fasta files using pair of primer.
# Make usage of https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters.

#' @param ref_fasta (required) A link to a database in fasta format
#'   (can be gzip).
#' @param output A path to the output fasta file.
#' @param primer_fw (Required, String) The forward primer DNA sequence.
#' @param primer_rev (Required String)  The reverse primer DNA sequence.
#' @param discard_untrimmed (logical, default TRUE) Do we add the param
#'   '--discard-untrimmed' to cutadapt function.
#' @param nproc (default 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param verbose (logical). If TRUE, print additional information.
#' @param cmd_is_run (logical, default TRUE) Do the cutadapt command is run.
#'   If set to FALSE, the only effect of the function is to return a list of
#'   command to manually run in a terminal.
#' @param return_file_path (logical, default FALSE) If TRUE, the function return
#'   a path to the file instead of returning the command. Useful in {targets}
#'   pipeline.
#' @param start_with_fw (logical, default FALSE) If TRUE, the forward sequences
#'   must be at the left hand without any nucleotides before.
#' @param output_json (logical, default FALSE) If TRUE, a json file summarizing
#'   cutadapt process is write in the working directory.
#' @param args_before_cutadapt (String) A one line bash command to run before
#' to run cutadapt. For examples, "source ~/miniconda3/etc/profile.d/conda.sh && conda activate cutadaptenv &&" allow to bypass the conda init which asks to restart the shell
#'
#' @returns
#' @export
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' \dontrun{
#' cutadapt_rm_primers_db(system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
#'                                output= "unite_cutadapted.fasta",
#'                            primer_fw="GCATCGATGAAGAACGCAGC",
#'                            primer_rev= "TCCTCCGCTTATTGATATGC")
#'
#'  unlink("unite_cutadapted.fasta")
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite cutadapt (\doi{doi:10.14806/ej.17.1.200}).
cutadapt_rm_primers_db <-
  function(ref_fasta,
           output = NULL,
           primer_fw = NULL,
           primer_rev = NULL,
           discard_untrimmed = TRUE,
           nproc = 1,
           verbose = TRUE,
           cmd_is_run = TRUE,
           return_file_path = FALSE,
           start_with_fw = FALSE,
           output_json = FALSE,
           error_tolerance = 0.1,
           args_before_cutadapt =
             "source ~/miniforge3/etc/profile.d/conda.sh && conda activate cutadaptenv && ") {
    if (is.null(output)) {
      output = paste0(basename(file.path(ref_fasta)), "_cutadapted.fasta")
    }

    # primer_fw_RC <- dada2::rc(primer_fw)
    # primer_rev_RC <- dada2::rc(primer_rev)

    cmd <-
      paste0(
        args_before_cutadapt,
        "cutadapt --cores=",
        nproc,
        " -e ",
        error_tolerance,
        " -a '",
        ifelse(start_with_fw, "^", ""),
        primer_fw,
        "...",
        primer_rev,
        "' -o ",
        output
      )

    if (output_json) {
      cmd <- paste0(cmd,
                    " --json=",
                    basename(ref_fasta),
                    "_cutadapt.json")
    }

    if (discard_untrimmed) {
      cmd <- paste0(cmd, " --discard-untrimmed")
    }

    cmd <-   paste0(cmd,
                    " ",
                    normalizePath(ref_fasta))

    if (cmd_is_run) {
      writeLines(cmd, paste0(tempdir(), "/script_cutadapt.sh"))
      system2("bash", paste0(tempdir(), "/script_cutadapt.sh"))
      if (verbose) {
        message(paste0("Output file is available: ", normalizePath(output)))
      }
      unlink(paste0(tempdir(), "/script_cutadapt.sh"))
    } else {
      return(cmd)
    }

    if (verbose) {
      nseq_initial <- count_seq(ref_fasta)
      nseq_final <- count_seq(output)

      message(
        "The cutadapt process trimmed ",
        nseq_initial - nseq_final,
        " (",
        round((nseq_initial - nseq_final) / nseq_initial * 100,2),
        "%)",
        " references sequences, for a final number of ",
        nseq_final,
        " references sequences."
      )

      n_nuc_initial <- sum(Biostrings::width(Biostrings::readDNAStringSet(ref_fasta)))
      n_nuc_final <- sum(Biostrings::width(Biostrings::readDNAStringSet(output)))

      message(
        "The cutadapt process trimmed ",
        n_nuc_initial - n_nuc_final,
        " (",
        round((n_nuc_initial - n_nuc_final) / n_nuc_initial*100,2),
        "%)",
        " nucleotides, for a final number of ",
        n_nuc_final,
        " nucleotides.\n The mean width of references sequences is now ",
        round(mean(Biostrings::width(Biostrings::readDNAStringSet(output))),2),
        " vs ",
        round(mean(Biostrings::width(Biostrings::readDNAStringSet(ref_fasta))),2),
        " in the original database."
      )
    }

    if (return_file_path) {
      return(normalizePath(output))
    } else {
      return(cmd)
    }
  }





#' remove_file_extension("data/data_raw/refseq/dada2_format/Unite.fasta")
#' remove_file_extension("data/data_raw/refseq/dada2_format/Unite.fasta", full.names = TRUE)
remove_file_extension <- function(file_path, full.names=FALSE) {
  ext <- get_file_extension(file_path)
  if(full.names) {
    new_path <- gsub(paste0(".", ext, "$"),"", file_path)
  } else {
    new_path <- gsub(paste0(".", ext, "$"),"", base::basename(file_path))
  }
  return(new_path)
}


#' add_suffix_before_ext("data/data_raw/refseq/dada2_format/Unite.fasta", "_cut")
#' add_suffix_before_ext("data/data_raw/refseq/dada2_format/Unite.fasta", "_cut", full.names=TRUE)
add_suffix_before_ext <- function(file_path, suffix, full.names=FALSE){
  ext <- get_file_extension(file_path)
  new_path <- paste0(remove_file_extension(file_path, full.names=full.names),
                     suffix,
                     ".",
                     ext)
  return(new_path)
}



