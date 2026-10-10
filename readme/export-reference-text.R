args <- commandArgs(trailingOnly = TRUE)
destination <- if (length(args)) args[[1]] else file.path("output", ".reference-text")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)

rd_files <- sort(list.files("man", pattern = "[.]Rd$", full.names = TRUE))
for (path in rd_files) {
  rd <- tools::parse_Rd(path)
  rendered <- capture.output(tools::Rd2txt(rd, fragment = FALSE,
                                            options = list(underline_titles = FALSE)))
  writeLines(enc2utf8(rendered), file.path(destination,
    paste0(tools::file_path_sans_ext(basename(path)), ".txt")), useBytes = TRUE)
}

desc <- read.dcf("DESCRIPTION")[1, ]
metadata <- c(Package = desc[["Package"]], Version = desc[["Version"]],
              Date = format(Sys.Date()), Documents = length(rd_files))
writeLines(paste(names(metadata), metadata, sep = "\t"),
           file.path(destination, "metadata.tsv"), useBytes = TRUE)
message("Exported ", length(rd_files), " help topics to ", destination)
