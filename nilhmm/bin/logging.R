# Shared logger setup for the nilhmm R scripts (lab convention; mirrors zealhmm/scripts/logging.R).
# Source once near the top of each script (finds logging.R next to the running script under Rscript):
#   .bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
#   source(file.path(.bin, "logging.R"))
#
# Emits timestamped, sprintf-style logs to stderr (-> Nextflow .command.err, tail -f'able):
#   [HH:MM:SS] INFO: message
# Use log_info()/log_warn()/log_error() with "%s"/"%d"/"%.1f" formats (NOT paste/sprintf inside --
# the formatter is sprintf). Tag each script's messages, e.g. log_info("[build_hd] ...").
# For any loop, log a running ETA after each iteration:
#   t0 <- Sys.time()                                              # before the loop
#   el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
#   log_info(">>> %d/%d done | elapsed %.1f min | avg %.1f min | ETA ~%.1f min remaining",
#            i, N, el, el/i, (el/i)*(N-i))
suppressMessages(library(logger))
log_layout(layout_glue_generator(format = '[{format(time, "%H:%M:%S")}] {level}: {msg}'))
log_formatter(formatter_sprintf)
log_threshold(INFO)
