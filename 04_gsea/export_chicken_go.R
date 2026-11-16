#!/usr/bin/env Rscript
# Chicken GO biological-process gene sets as GMT, keyed by gene symbol (input to gsea.py).
# GO2ALLEGS gives each term the genes annotated to it and to all of its descendants.
# Set names are "GO:id|term". Sizes 5-2000 here; gsea.py applies the 15-1000 filter.
#   Rscript export_chicken_go.R [GOBP_chicken.gmt]
suppressPackageStartupMessages({ library(org.Gg.eg.db); library(GO.db) })

out   <- if (length(a <- commandArgs(TRUE))) a[1] else "GOBP_chicken.gmt"
go2eg <- AnnotationDbi::as.list(org.Gg.egGO2ALLEGS)
info  <- suppressMessages(AnnotationDbi::select(GO.db, keys = names(go2eg), keytype = "GOID",
                                                columns = c("ONTOLOGY", "TERM")))
go2eg <- go2eg[names(go2eg) %in% info$GOID[info$ONTOLOGY == "BP"]]

eg2sym <- unlist(AnnotationDbi::as.list(org.Gg.egSYMBOL))
sets   <- lapply(go2eg, function(eg) unique(stats::na.omit(eg2sym[unique(eg)])))
sets   <- sets[lengths(sets) >= 5 & lengths(sets) <= 2000]

term <- info$TERM[match(names(sets), info$GOID)]
writeLines(vapply(seq_along(sets), function(i)
  paste(c(paste0(names(sets)[i], "|", term[i]), names(sets)[i], sets[[i]]), collapse = "\t"),
  character(1)), out)
