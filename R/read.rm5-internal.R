##
##
## Definition of auxiliary functions to import RevMan 5 analysis datasets
##
##


extract_outcomes <- function(txt, outcome.type, res,
                             comp.no, complab,
                             debug = 0) {
  ##
  pattern1 <- paste0("<", outcome.type, "_OUTCOME")
  pattern2 <- paste0("</", outcome.type, "_OUTCOME")
  sel1 <- grep(pattern1, txt)
  sel2 <- grep(pattern2, txt) 
  ##
  if (length(sel1) > 0 | length(sel2) > 0) {
    ##
    if (debug)
      cat(paste0("\n*** ", outcome.type, "_OUTCOME ***\n"))
    ##
    if (length(sel1) != length(sel2))
      stop("Malformed XML file (tag: ", outcome.type, "_OUTCOME)")
    ##
    for (j in seq(along = sel1)) {
      txt.j <- txt[sel1[j]:sel2[j]]
      sel.data.j <- c(grep(paste0("<", outcome.type, "_DATA"), txt.j) - 1,
                      grep(paste0("<", outcome.type, "_SUBGROUP"), txt.j) - 1)
      if (outcome.type %in% c("IV", "IPD"))
        sel.data.j <- c(sel.data.j, grep("</EFFECT_MEASURE>", txt.j))
      else
        sel.data.j <- c(sel.data.j, grep("</GRAPH_LABEL_2>", txt.j))
      if (length(sel.data.j) > 0)
        sel.j <- c(1:(max(c(6, min(sel.data.j)))), length(txt.j))
      else
        sel.j <- c(1:6, length(txt.j))
      ##
      xml.j <- as_xml_document(paste(txt.j[sel.j], collapse = " "))
      ##
      outcome.no <- as.numeric(xml_attr(xml.j, "NO"))
      outclab <- xml_text(xml_find_all(xml.j, "//NAME"))
      ##
      if (debug) {
        cat(paste0("Outcome.no: ", outcome.no, "\n"))
        cat(paste0("outclab: ", outclab, "\n"))
        ##
        if (debug == 3)
          print(txt.j[sel.j])
      }
      ##
      overall <- xml_attr(xml.j, "TOTALS") == "YES"
      test.subgroup <- xml_attr(xml.j, "SUBGROUP_TEST") == "YES"
      ##
      n.studies <- as.numeric(xml_attr(xml.j, "DF")) + 1
      ##
      type <- switch(outcome.type,
                     DICH = "D",
                     CONT = "C",
                     IV = "I",
                     IPD = "I")
      ##
      method <- xml_attr(xml.j, "METHOD")
      if (is.na(method) & outcome.type == "IPD")
        method <- "Peto"
      else if (is.na(method))
        method <- "Inverse"
      else if (method == "IV")
        method <- "Inverse"
      ##
      sm <- xml_attr(xml.j, "EFFECT_MEASURE")
      if (is.na(sm)) {
        sm <- xml_text(xml_find_all(xml.j, "//EFFECT_MEASURE"))
        sm[tolower(sm) == "odds ratio"] <- "OR"
        sm[tolower(sm) == "peto odds ratio"] <- "OR"
        sm[tolower(sm) == "peto_or"] <- "OR"
        sm[tolower(sm) == "odds ratio (non-event)"] <- "OR"
        sm[tolower(sm) == "risk ratio"] <- "RR"
        sm[tolower(sm) == "risk difference"] <- "RD"
        sm[tolower(sm) == "mean difference"] <- "MD"
        sm[tolower(sm) == "standardized mean difference"] <- "SMD"
        sm[tolower(sm) == "std. mean difference"] <- "SMD"
        sm[tolower(sm) == "hazard ratio"] <- "HR"
      }
      ##
      logscale <- xml_attr(xml.j, "LOG_DATA") == "NO"
      ##
      random <- xml_attr(xml.j, "RANDOM")
      if (is.na(random) & outcome.type == "IPD")
        random <- "NO"
      model <- ifelse(random == "NO", "Fixed", "Random")
      common <- model == "Fixed"
      random <- model == "Random"
      ##
      TE.pooled <- as.numeric(xml_attr(xml.j, "EFFECT_SIZE"))
      lower.pooled <- as.numeric(xml_attr(xml.j, "CI_START"))
      upper.pooled <- as.numeric(xml_attr(xml.j, "CI_END"))
      weight.pooled <- as.numeric(xml_attr(xml.j, "WEIGHT"))
      level <- as.numeric(xml_attr(xml.j, "CI_STUDY")) / 100
      level.ma <- as.numeric(xml_attr(xml.j, "CI_TOTAL")) / 100
      ##
      Z.pooled <- as.numeric(xml_attr(xml.j, "Z"))
      pval.TE.pooled <- as.numeric(xml_attr(xml.j, "P_Z"))
      Q <- as.numeric(xml_attr(xml.j, "Q"))
      pval.Q <- as.numeric(xml_attr(xml.j, "P_Q"))
      I2 <- as.numeric(xml_attr(xml.j, "I2"))
      tau2 <- as.numeric(xml_attr(xml.j, "TAU2"))
      ##
      event.e.pooled <- as.numeric(xml_attr(xml.j, "EVENTS_1"))
      n.e.pooled <- as.numeric(xml_attr(xml.j, "TOTAL_1"))
      event.c.pooled <- as.numeric(xml_attr(xml.j, "EVENTS_2"))
      n.c.pooled <- as.numeric(xml_attr(xml.j, "TOTAL_2"))
      ##
      swap.events <- xml_attr(xml.j, "SWAP_EVENTS") != "NO"
      ##
      label.e <- xml_text(xml_find_all(xml.j, "//GROUP_LABEL_1"))
      label.c <- xml_text(xml_find_all(xml.j, "//GROUP_LABEL_2"))
      label.left <- xml_text(xml_find_all(xml.j, "//GRAPH_LABEL_1"))
      label.right <- xml_text(xml_find_all(xml.j, "//GRAPH_LABEL_2"))
      ##
      ## Subgroups
      ##
      sel.subgroup1 <- grep(paste0("<", outcome.type, "_SUBGROUP"), txt.j)
      sel.subgroup2 <- grep(paste0("</", outcome.type, "_SUBGROUP"), txt.j)
      ##
      if (length(sel.subgroup1 > 0) | length(sel.subgroup2) > 0) {
        if (debug)
          cat(paste0("\n** ", outcome.type, "_SUBGROUP **\n"))
        ##
        if (length(sel.subgroup1) != length(sel.subgroup2))
          stop("Malformed XML file (tag: ", outcome.type, "_SUBGROUP)")
        ##
        for (k in seq(along = sel.subgroup1)) {
          txt.jk <- txt.j[sel.subgroup1[k]:sel.subgroup2[k]]
          sel.data.jk <- c(grep(paste0("<", outcome.type, "_DATA"), txt.jk) - 1,
                           grep("</NAME>", txt.jk))
          if (length(sel.data.jk) > 0)
            sel.jk <- unique(c(1:(max(c(2, min(sel.data.jk)))), length(txt.jk)))
          else
            sel.jk <- unique(c(1:2, length(txt.jk)))
          ##
          xml.jk <- as_xml_document(paste(txt.jk[sel.jk], collapse = " "))
          ##
          group.no <- as.numeric(xml_attr(xml.jk, "NO"))
          grplab <- xml_text(xml_find_all(xml.jk, "//NAME"))
          ##
          if (debug) {
            cat(paste0("Group.no: ", group.no, "\n"))
            cat(paste0("grplab: ", grplab, "\n"))
            ##
            if (debug == 3)
              print(txt.jk[sel.jk])
          }
          ##
          Q.w <- as.numeric(xml_attr(xml.jk, "NO"))
          pval.Q.w <- as.numeric(xml_attr(xml.jk, "P_CHI2"))
          I2.w <- as.numeric(xml_attr(xml.jk, "I2"))
          ##
          ## Data
          ##
          sel.data1 <- grep(paste0("<", outcome.type, "_DATA"), txt.jk)
          ##
          sel1.data2 <- grepl("/>", txt.jk)
          sel2.data2 <- grepl(paste0("</", outcome.type, "_DATA"), txt.jk)
          sel3.data2 <-
            !(grepl("<GROUP_LABEL", txt.jk) |
              grepl("<GRAPH_LABEL", txt.jk) |
              grepl("<EFFECT_MEASURE", txt.jk))
          sel.data2 <-
            seq_along(txt.jk)[(sel1.data2 | sel2.data2) & sel3.data2]
          ##
          if (length(sel.data1) != length(sel.data2)) {
            if (debug == 3) {
              print(txt.jk)
              cat(paste0("sel.data1: ", sel.data1, "\n"))
              cat(paste0("sel.data2: ", sel.data2, "\n"))
              ##
              print(data.frame(sel1.data2, sel2.data2, sel3.data2))
            }
            stop("Malformed XML file (tag: ", outcome.type, "_DATA)")
          }
          ##
          event.e <- n.e <- event.c <- n.c <-
            mean.e <- sd.e <- mean.c <- sd.c <-
              O.E <- V <-
                TE <- seTE <- lower.TE <- upper.TE <-
                  weight <- order <- 
                    rep(NA, length(sel.data1))
          ##
          id <- rep("", length(sel.data1))
          ##
          if (length(sel.data1) > 0) {
            if (debug)
              cat(paste0("\n* ", outcome.type, "_DATA *\n"))
            for (l in seq(along = sel.data1)) {
              txt.jkl <- paste(txt.jk[sel.data1[l]:sel.data2[l]],
                               collapse = "")
              ##
              if (debug == 3)
                print(txt.jk[sel.data1[l]:sel.data2[l]])
              ##
              xml.jkl <- as_xml_document(txt.jkl)
              ##
              id[l] <- xml_attr(xml.jkl, "STUDY_ID")
              ##
              event.e[l] <- as.numeric(xml_attr(xml.jkl, "EVENTS_1"))
              n.e[l] <- as.numeric(xml_attr(xml.jkl, "TOTAL_1"))
              mean.e[l] <- as.numeric(xml_attr(xml.jkl, "MEAN_1"))
              sd.e[l] <- as.numeric(xml_attr(xml.jkl, "SD_1"))
              ##
              event.c[l] <- as.numeric(xml_attr(xml.jkl, "EVENTS_2"))
              n.c[l] <- as.numeric(xml_attr(xml.jkl, "TOTAL_2"))
              mean.c[l] <- as.numeric(xml_attr(xml.jkl, "MEAN_2"))
              sd.c[l] <- as.numeric(xml_attr(xml.jkl, "SD_2"))
              ##
              O.E[l] <- as.numeric(xml_attr(xml.jkl, "O_E"))
              V[l] <- as.numeric(xml_attr(xml.jkl, "VAR"))
              ##
              TE[l] <- as.numeric(xml_attr(xml.jkl, "EFFECT_SIZE"))
              seTE[l] <- as.numeric(xml_attr(xml.jkl, "SE"))
              lower.TE[l] <- as.numeric(xml_attr(xml.jkl, "CI_START"))
              upper.TE[l] <- as.numeric(xml_attr(xml.jkl, "CI_END"))
              ##
              weight[l] <- as.numeric(xml_attr(xml.jkl, "WEIGHT"))
              order[l] <- as.numeric(xml_attr(xml.jkl, "ORDER"))
            }
            ##
            res.new <-
              data.frame(comp.no, outcome.no, group.no,
                         id = id,
                         event.e = event.e, n.e = n.e,
                         event.c = event.c, n.c = n.c,
                         mean.e = mean.e, sd.e = sd.e,
                         mean.c = mean.c, sd.c = sd.c,
                         O.E = O.E, V = V,
                         TE = TE, seTE = seTE,
                         lower.TE = lower.TE, upper.TE = upper.TE,
                         weight = weight, order = order,
                         level = level,
                         grplab = grplab,
                         overall = overall,
                         test.subgroup = test.subgroup,
                         type = type, method = method, sm = sm, model = model,
                         common = common, random = random,
                         outclab = outclab,
                         k = n.studies,
                         event.e.pooled = event.e.pooled,
                         n.e.pooled = n.e.pooled,
                         event.c.pooled = event.c.pooled,
                         n.c.pooled = n.c.pooled,
                         TE.pooled = TE.pooled,
                         lower.pooled = lower.pooled,
                         upper.pooled = upper.pooled,
                         weight.pooled = weight.pooled,
                         level.ma = level.ma,
                         Z.pooled = Z.pooled, pval.TE.pooled = pval.TE.pooled,
                         Q = Q, pval.Q = pval.Q,
                         I2 = I2, tau2 = tau2,
                         Q.w = Q.w, pval.Q.w = pval.Q.w, I2.w = I2.w,
                         swap.events = swap.events, logscale = logscale,
                         label.e = label.e, label.c = label.c,
                         label.left = label.left,
                         label.right = label.right,
                         complab = complab)
            ##
            if (debug >= 2) {
              cat("\n* New data: *\n")
              print(res.new)
            }
            res <- rbind(res, res.new)
          }
        }
      }
      else {
        group.no <- NA
        grplab <- ""
        Q.w <- NA
        pval.Q.w <- NA
        I2.w <- NA
        ##
        ## Data
        ##
        sel.data1 <- grep(paste0("<", outcome.type, "_DATA"), txt.j)
        ##
        sel1.data2 <- grepl("/>", txt.j)
        sel2.data2 <- grepl(paste0("</", outcome.type, "_DATA"), txt.j)
        sel3.data2 <-
          !(grepl("<GROUP_LABEL", txt.j) |
            grepl("<GRAPH_LABEL", txt.j) |
            grepl("<EFFECT_MEASURE", txt.j))
        sel.data2 <-
          seq_along(txt.j)[(sel1.data2 | sel2.data2) & sel3.data2]
        ##
        if (length(sel.data1) != length(sel.data2)) {
          if (debug == 3) {
            print(txt.j)
            cat(paste0("sel.data1: ", sel.data1, "\n"))
            cat(paste0("sel.data2: ", sel.data2, "\n"))
            ##
            print(data.frame(sel1.data2, sel2.data2, sel3.data2))
          }
          stop("Malformed XML file (tag: ", outcome.type, "_DATA)")
        }
        ##
        event.e <- n.e <- event.c <- n.c <-
          mean.e <- sd.e <- mean.c <- sd.c <-
            O.E <- V <-
              TE <- seTE <- lower.TE <- upper.TE <-
                weight <- order <- 
                  rep(NA, length(sel.data1))
        ##
        id <- rep("", length(sel.data1))
        ##
        if (length(sel.data1) > 0) {
          if (debug) {
            cat(paste0("\n* ", outcome.type, "_DATA *\n"))
            print(sel.data1)
          }
          ##
          for (k in seq(along = sel.data1)) {
            txt.jk <- paste(txt.j[sel.data1[k]:sel.data2[k]], collapse = "")
            ##
            if (debug == 3)
              print(txt.j[sel.data1[k]:sel.data2[k]])
            ##
            xml.jk <- as_xml_document(txt.jk)
            ##
            id[k] <- xml_attr(xml.jk, "STUDY_ID")
            ##
            event.e[k] <- as.numeric(xml_attr(xml.jk, "EVENTS_1"))
            n.e[k] <- as.numeric(xml_attr(xml.jk, "TOTAL_1"))
            mean.e[k] <- as.numeric(xml_attr(xml.jk, "MEAN_1"))
            sd.e[k] <- as.numeric(xml_attr(xml.jk, "SD_1"))
            ##
            event.c[k] <- as.numeric(xml_attr(xml.jk, "EVENTS_2"))
            n.c[k] <- as.numeric(xml_attr(xml.jk, "TOTAL_2"))
            mean.c[k] <- as.numeric(xml_attr(xml.jk, "MEAN_2"))
            sd.c[k] <- as.numeric(xml_attr(xml.jk, "SD_2"))
            ##
            O.E[k] <- as.numeric(xml_attr(xml.jk, "O_E"))
            V[k] <- as.numeric(xml_attr(xml.jk, "VAR"))
            ##
            TE[k] <- as.numeric(xml_attr(xml.jk, "EFFECT_SIZE"))
            seTE[k] <- as.numeric(xml_attr(xml.jk, "SE"))
            lower.TE[k] <- as.numeric(xml_attr(xml.jk, "CI_START"))
            upper.TE[k] <- as.numeric(xml_attr(xml.jk, "CI_END"))
            ##
            weight[k] <- as.numeric(xml_attr(xml.jk, "WEIGHT"))
            order[k] <- as.numeric(xml_attr(xml.jk, "ORDER"))
          }
          ##
          res.new <-
            data.frame(comp.no, outcome.no, group.no,
                       id = id,
                       event.e = event.e, n.e = n.e,
                       event.c = event.c, n.c = n.c,
                       mean.e = mean.e, sd.e = sd.e,
                       mean.c = mean.c, sd.c = sd.c,
                       O.E = O.E, V = V,
                       TE = TE, seTE = seTE,
                       lower.TE = lower.TE, upper.TE = upper.TE,
                       weight = weight, order = order,
                       level = level,
                       grplab = grplab,
                       overall = overall,
                       test.subgroup = test.subgroup,
                       type = type, method = method, sm = sm, model = model,
                       common = common, random = random,
                       outclab = outclab,
                       k = n.studies,
                       event.e.pooled = event.e.pooled,
                       n.e.pooled = n.e.pooled,
                       event.c.pooled = event.c.pooled,
                       n.c.pooled = n.c.pooled,
                       TE.pooled = TE.pooled,
                       lower.pooled = lower.pooled,
                       upper.pooled = upper.pooled,
                       weight.pooled = weight.pooled,
                       level.ma = level.ma,
                       Z.pooled = Z.pooled, pval.TE.pooled = pval.TE.pooled,
                       Q = Q, pval.Q = pval.Q,
                       I2 = I2, tau2 = tau2,
                       Q.w = Q.w, pval.Q.w = pval.Q.w, I2.w = I2.w,
                       swap.events = swap.events, logscale = logscale,
                       label.e = label.e, label.c = label.c,
                       label.left = label.left,
                       label.right = label.right,
                       complab = complab)
          ##
          if (debug >= 2) {
            cat("\n* New data: *\n")
            print(res.new)
          }
          ##
          res <- rbind(res, res.new)        
        }
      }
    }
  }
  ##
  res
}


oct2txt <- function(txt) {
  txt <- gsub("\200", "EUR", txt, useBytes = TRUE)
  txt <- gsub("\202", "'", txt, useBytes = TRUE)
  txt <- gsub("\203", "f", txt, useBytes = TRUE)
  txt <- gsub("\204", '"', txt, useBytes = TRUE)
  txt <- gsub("\205", "...", txt, useBytes = TRUE)
  txt <- gsub("\206", "+", txt, useBytes = TRUE)
  txt <- gsub("\207", "+", txt, useBytes = TRUE)
  txt <- gsub("\210", "^", txt, useBytes = TRUE)
  txt <- gsub("\211", "%%", txt, useBytes = TRUE)
  txt <- gsub("\212", "S", txt, useBytes = TRUE)
  txt <- gsub("\213", "'", txt, useBytes = TRUE)
  txt <- gsub("\214", "OE", txt, useBytes = TRUE)
  txt <- gsub("\216", "Z", txt, useBytes = TRUE)
  txt <- gsub("\221", "'", txt, useBytes = TRUE)
  txt <- gsub("\222", "'", txt, useBytes = TRUE)
  txt <- gsub("\223", '"', txt, useBytes = TRUE)
  txt <- gsub("\224", '"', txt, useBytes = TRUE)
  txt <- gsub("\225", "-", txt, useBytes = TRUE)
  txt <- gsub("\226", "-", txt, useBytes = TRUE)
  txt <- gsub("\227", "-", txt, useBytes = TRUE)
  txt <- gsub("\230", "~", txt, useBytes = TRUE)
  txt <- gsub("\231", "TM", txt, useBytes = TRUE)
  txt <- gsub("\232", "S", txt, useBytes = TRUE)
  txt <- gsub("\233", "'", txt, useBytes = TRUE)
  txt <- gsub("\234", "oe", txt, useBytes = TRUE)
  txt <- gsub("\237", "Y", txt, useBytes = TRUE)
  txt <- gsub("\240", " ", txt, useBytes = TRUE)
  txt <- gsub("\241", "!", txt, useBytes = TRUE)
  txt <- gsub("\242", "cent", txt, useBytes = TRUE)
  txt <- gsub("\243", "Pound", txt, useBytes = TRUE)
  txt <- gsub("\244", "C", txt, useBytes = TRUE)
  txt <- gsub("\245", "Yen", txt, useBytes = TRUE)
  txt <- gsub("\246", "|", txt, useBytes = TRUE)
  txt <- gsub("\247", "S", txt, useBytes = TRUE)
  txt <- gsub("\250", "", txt, useBytes = TRUE)
  txt <- gsub("\251", "(C)", txt, useBytes = TRUE)
  txt <- gsub("\252", "", txt, useBytes = TRUE)
  txt <- gsub("\253", "'", txt, useBytes = TRUE)
  txt <- gsub("\254", "-", txt, useBytes = TRUE)
  txt <- gsub("\255", "-", txt, useBytes = TRUE)
  txt <- gsub("\256", "(R)", txt, useBytes = TRUE)
  txt <- gsub("\257", "-", txt, useBytes = TRUE)
  txt <- gsub("\260", "o", txt, useBytes = TRUE)
  txt <- gsub("\261", "+/-", txt, useBytes = TRUE)
  txt <- gsub("\262", "2", txt, useBytes = TRUE)
  txt <- gsub("\263", "3", txt, useBytes = TRUE)
  txt <- gsub("\264", "", txt, useBytes = TRUE)
  txt <- gsub("\265", "mu", txt, useBytes = TRUE)
  txt <- gsub("\266", "", txt, useBytes = TRUE)
  txt <- gsub("\267", ".", txt, useBytes = TRUE)
  txt <- gsub("\270", " ", txt, useBytes = TRUE)
  txt <- gsub("\271", "1", txt, useBytes = TRUE)
  txt <- gsub("\272", "o", txt, useBytes = TRUE)
  txt <- gsub("\273", '"', txt, useBytes = TRUE)
  txt <- gsub("\274", "1/4", txt, useBytes = TRUE)
  txt <- gsub("\275", "1/2", txt, useBytes = TRUE)
  txt <- gsub("\276", "3/4", txt, useBytes = TRUE)
  txt <- gsub("\277", "?", txt, useBytes = TRUE)
  txt <- gsub("\300", "A", txt, useBytes = TRUE)
  txt <- gsub("\301", "A", txt, useBytes = TRUE)
  txt <- gsub("\302", "A", txt, useBytes = TRUE)
  txt <- gsub("\303", "A", txt, useBytes = TRUE)
  txt <- gsub("\304", "A", txt, useBytes = TRUE)
  txt <- gsub("\305", "A", txt, useBytes = TRUE)
  txt <- gsub("\306", "AE", txt, useBytes = TRUE)
  txt <- gsub("\307", "C", txt, useBytes = TRUE)
  txt <- gsub("\310", "E", txt, useBytes = TRUE)
  txt <- gsub("\311", "E", txt, useBytes = TRUE)
  txt <- gsub("\312", "E", txt, useBytes = TRUE)
  txt <- gsub("\313", "E", txt, useBytes = TRUE)
  txt <- gsub("\314", "I", txt, useBytes = TRUE)
  txt <- gsub("\315", "I", txt, useBytes = TRUE)
  txt <- gsub("\316", "I", txt, useBytes = TRUE)
  txt <- gsub("\317", "I", txt, useBytes = TRUE)
  txt <- gsub("\320", "E", txt, useBytes = TRUE)
  txt <- gsub("\321", "N", txt, useBytes = TRUE)
  txt <- gsub("\322", "O", txt, useBytes = TRUE)
  txt <- gsub("\323", "O", txt, useBytes = TRUE)
  txt <- gsub("\324", "O", txt, useBytes = TRUE)
  txt <- gsub("\325", "O", txt, useBytes = TRUE)
  txt <- gsub("\326", "O", txt, useBytes = TRUE)
  txt <- gsub("\327", "x", txt, useBytes = TRUE)
  txt <- gsub("\330", "O", txt, useBytes = TRUE)
  txt <- gsub("\331", "U", txt, useBytes = TRUE)
  txt <- gsub("\332", "U", txt, useBytes = TRUE)
  txt <- gsub("\333", "U", txt, useBytes = TRUE)
  txt <- gsub("\334", "U", txt, useBytes = TRUE)
  txt <- gsub("\335", "Y", txt, useBytes = TRUE)
  txt <- gsub("\336", "T", txt, useBytes = TRUE)
  txt <- gsub("\337", "ss", txt, useBytes = TRUE)
  txt <- gsub("\340", "a", txt, useBytes = TRUE)
  txt <- gsub("\341", "a", txt, useBytes = TRUE)
  txt <- gsub("\342", "a", txt, useBytes = TRUE)
  txt <- gsub("\343", "a", txt, useBytes = TRUE)
  txt <- gsub("\344", "a", txt, useBytes = TRUE)
  txt <- gsub("\345", "a", txt, useBytes = TRUE)
  txt <- gsub("\346", "a", txt, useBytes = TRUE)
  txt <- gsub("\347", "c", txt, useBytes = TRUE)
  txt <- gsub("\350", "e", txt, useBytes = TRUE)
  txt <- gsub("\351", "e", txt, useBytes = TRUE)
  txt <- gsub("\352", "e", txt, useBytes = TRUE)
  txt <- gsub("\353", "e", txt, useBytes = TRUE)
  txt <- gsub("\354", "i", txt, useBytes = TRUE)
  txt <- gsub("\355", "i", txt, useBytes = TRUE)
  txt <- gsub("\356", "i", txt, useBytes = TRUE)
  txt <- gsub("\357", "i", txt, useBytes = TRUE)
  txt <- gsub("\360", "e", txt, useBytes = TRUE)
  txt <- gsub("\361", "n", txt, useBytes = TRUE)
  txt <- gsub("\362", "o", txt, useBytes = TRUE)
  txt <- gsub("\363", "o", txt, useBytes = TRUE)
  txt <- gsub("\364", "o", txt, useBytes = TRUE)
  txt <- gsub("\365", "o", txt, useBytes = TRUE)
  txt <- gsub("\366", "o", txt, useBytes = TRUE)
  txt <- gsub("\367", "/", txt, useBytes = TRUE)
  txt <- gsub("\370", "o", txt, useBytes = TRUE)
  txt <- gsub("\371", "u", txt, useBytes = TRUE)
  txt <- gsub("\372", "u", txt, useBytes = TRUE)
  txt <- gsub("\373", "u", txt, useBytes = TRUE)
  txt <- gsub("\374", "u", txt, useBytes = TRUE)
  txt <- gsub("\375", "y", txt, useBytes = TRUE)
  txt <- gsub("\376", "t", txt, useBytes = TRUE)
  txt <- gsub("\377", "y", txt, useBytes = TRUE)
  txt <- gsub("\303\237", "AY", txt, useBytes = TRUE)
  
  # replace any remaining non UTF-8 characters with question marks
  # otherwise they crash the `extract_outcomes` function
  txt <- iconv(txt, from = "", to = "UTF-8", sub = "???")
  
  ##
  txt
}


read.rm5.csv <- function(file, sep = ",", quote = "\"",
                         title, numbers.in.labels = TRUE) {
  ##
  selvar <- function(x, sel, value = NA) {
    res <-
      if (!is.null(x))
        x[sel]
      else value[sel]
    res
  }
  ##
  numchar <- function(x) {
    res <- as.numeric(as.character(x))
    res
  }
  ##
  nachar <- function(x) {
    res <- x
    res[is.na(res)] <- ""
    res
  }
  ##
  if (missing(title)) {
    title <- strsplit(file, "\\.csv$")[[1]]
    tmp <- strsplit(title, "\\/")
    title <- tmp[[1]][length(tmp[[1]])]
  }
  ##
  ## Check number of fields in input data file
  ##
  if (length(unique(count.fields(file, sep = sep, quote = quote,
                                 comment.char = ""))) != 1)
    stop("Input file has different number of elements. ",
         "Probably your export data file contains the field ",
         "'Risk of bias tables'. Please create a new export data file in ",
         "RevMan 5 without this field.")
  ##
  tdata <- read.table(file, header = TRUE,
                      sep = sep, quote = quote,
                      comment.char = "")
  ##
  nam <- names(tdata)
  
  ##
  ## Overcome problems to import files with UTF-8, byte order mark encoding
  ## (see https://en.wikipedia.org/wiki/Byte_order_mark)
  ##
  nam[grep("Comparison.Number$", nam)] <- "Comparison.Number"
  ##
  if (!all(c("Comparison.Number", "Outcome.Number",
             "Subgroup.Number") %in% nam))
    stop("Mandatory fields 'Comparison Number', 'Outcome Number',",
         " and 'Subgroup Number' not included in export file ",
         deparse(substitute(file)),
         " (see help page of function read.rm5).")
  ##
  nam[nam == "Comparison.Number"] <- "comp.no"
  nam[nam == "Outcome.Number"] <- "outcome.no"
  nam[nam == "Subgroup.Number"] <- "group.no"
  nam[nam == "Name"] <- "author"
  ##
  nam[nam == "Data.Type"] <- "type"
  nam[nam == "Statistical.Method"] <- "method"
  nam[nam == "Effect.Measure"] <- "sm"
  nam[nam == "Analysis.Model"] <- "model"
  ##
  ## Variable 'Totals' not included in data frame 
  ## Variable 'Study.Confidence.Interval' not included in data frame 
  ## Variable 'Total.Confidence.Interval' not included in data frame 
  ## Variable 'Test.for.subgroup.differences' not included in data frame 
  ##
  nam[nam == "Swap.event.and.non.event"] <- "swap.events"
  nam[nam == "Entered.data.are.on.log.scale..Generic.Inverse.Variance.only."] <-
    "logscale"
  nam[nam == "Enter.number.of.participants..Generic.Inverse.Variance..for.display.only."] <-
    "enter.n"
  ##
  nam[nam == "Events.1"] <- "event.e"
  nam[nam == "Mean.1"] <- "mean.e"
  nam[nam == "SD.1"] <- "sd.e"
  nam[nam == "Total.1"] <- "n.e"
  ##
  nam[nam == "Events.2"] <- "event.c"
  nam[nam == "Mean.2"] <- "mean.c"
  nam[nam == "SD.2"] <- "sd.c"
  nam[nam == "Total.2"] <- "n.c"
  ##
  nam[nam == "O.E"] <- "O.E"
  nam[nam == "Var"] <- "V"
  ##
  nam[nam == "Effect.Estimate"] <- "TE"
  nam[nam == "SE"] <- "seTE"
  ##
  nam[nam == "CI.Start"] <- "lower"
  nam[nam == "CI.End"]   <- "upper"
  ##
  nam[nam == "Weight"]   <- "weight"
  ##
  ## Variable 'Q' is Q
  nam[nam == "P.Q."]  <- "pval.Q"
  nam[nam == "I..Q."] <- "I2"
  nam[nam == "Tau."]  <- "tau2"
  ##
  ## Variable 'Z' is Z
  nam[nam == "P.Z."]     <- "pval.TE"
  ## Variable 'Qint' is Qint
  nam[nam == "P.Qint."]  <- "pval.Qint"
  nam[nam == "I..Qint."] <- "I2.Qint"
  ## Variable 'df' is df
  ##
  nam[nam == "Group.Label.1"]     <- "label.e"
  nam[nam == "Group.Label.2"]     <- "label.c"
  nam[nam == "Left.Graph.Label"]  <- "label.left"
  nam[nam == "Right.Graph.Label"] <- "label.right"
  ##
  nam[nam == "Year.of.study"]      <- "year"
  nam[nam == "User.defined.order"] <- "order"
  ##
  ## Fancy coding to set variable names of I^2 and tau^2
  ## as "I2" and "tau2" (not necessary for Linux)
  ## I don't know how to do this differently with Windows ...
  ##
  if (!all(c("I2", "tau2") %in% nam)) {
    pos1 <- seq(along=nam)[nam == "pval.Q"]
    pos2 <- seq(along=nam)[nam == "Z"]
    if ((pos2 - pos1) == 3)
      nam[pos1 + 1:2] <- c("I2", "tau2")
  }
  if (!"I2.Qint" %in% nam) {
    pos3 <- seq(along=nam)[nam == "pval.Qint"]
    pos4 <- seq(along=nam)[nam == "df"]
    if ((pos4 - pos3) == 2)
      nam[pos1 + 1] <- c("I2.Qint")
  }
  ##
  names(tdata) <- nam
  ##
  tdata$author <- as.character(tdata$author)
  tdata$type <- as.character(tdata$type)
  tdata$method <- as.character(tdata$method)
  tdata$method[tdata$method == "IV"] <- "Inverse"
  tdata$sm <- as.character(tdata$sm)
  tdata$model <- as.character(tdata$model)
  tdata$common <- tdata$model == "Fixed"
  tdata$random <- tdata$model == "Random"
  ##
  tdata$sm <- em2sm(tdata$sm)
  ##
  sel.oe <- tdata$method == "EXP_O_E_VAR"
  tdata$method[sel.oe] <- "Peto"
  tdata$sm[sel.oe] <- "OR"
  ##
  sel.peto <- tdata$method == "PETO" & tdata$sm == "PETO_OR"
  tdata$method[sel.peto] <- "Peto"
  tdata$sm[sel.peto] <- "OR"
  ##
  tdata$event.e <- as.numeric(gsub(",", "", tdata$event.e))
  tdata$n.e <- as.numeric(gsub(",", "", tdata$n.e))
  tdata$event.c <- as.numeric(gsub(",", "", tdata$event.c))
  tdata$n.c <- as.numeric(gsub(",", "", tdata$n.c))
  tdata$mean.e <- as.numeric(gsub(",", "", tdata$mean.e))
  tdata$sd.e <- as.numeric(gsub(",", "", tdata$sd.e))
  tdata$mean.c <- as.numeric(gsub(",", "", tdata$mean.c))
  tdata$sd.c <- as.numeric(gsub(",", "", tdata$sd.c))
  if (!is.null(tdata$O.E))
    tdata$O.E <- as.numeric(gsub(",", "", tdata$O.E))
  if (!is.null(tdata$V))
    tdata$V <- as.numeric(gsub(",", "", tdata$V))
  tdata$TE <- as.numeric(gsub(",", "", tdata$TE))
  ## No warning concerning infinite values for
  ## seTE, lower and upper CI bound
  tdata$seTE <- suppressWarnings(as.numeric(gsub(",", "", tdata$seTE)))
  tdata$lower <- suppressWarnings(as.numeric(gsub(",", "", tdata$lower)))
  tdata$upper <- suppressWarnings(as.numeric(gsub(",", "", tdata$upper)))
  ##
  tdata$weight <- as.numeric(gsub(",", "", tdata$weight))
  tdata$Q <- as.numeric(gsub(",", "", tdata$Q))
  tdata$tau2 <- as.numeric(gsub(",", "", tdata$tau2))
  tdata$Qint <- as.numeric(gsub(",", "", tdata$Qint))
  ##
  diffcomp <- c(1, diff(tdata$comp.no))
  diffoutc <- c(1, diff(tdata$outcome.no))
  diffgrp  <- c(1, diff(tdata$group.no))
  ##
  sel.comp  <- diffcomp != 0
  sel.outc  <- diffoutc != 0 & !sel.comp
  sel.grp   <- tdata$type != "" & diffgrp  != 0 & !sel.comp & !sel.outc
  sel.study <- !sel.comp & !sel.outc & !sel.grp
  ##
  res <- list(title = title,
              comparison = 
                data.frame(comp.no = selvar(tdata$comp.no, sel.comp),
                           complab = selvar(tdata$author, sel.comp)
                           ),
              outcome = 
                data.frame(comp.no    = selvar(tdata$comp.no, sel.outc),
                           outcome.no = selvar(tdata$outcome.no, sel.outc),
                           ##
                           type        = selvar(tdata$type, sel.outc),
                           method      = selvar(tdata$method, sel.outc),
                           sm          = selvar(tdata$sm, sel.outc),
                           model       = selvar(tdata$model, sel.outc),
                           common = selvar(tdata$common, sel.outc),
                           random = selvar(tdata$random, sel.outc),
                           outclab     = selvar(tdata$author, sel.outc),
                           ##
                           k = selvar(tdata$df, sel.outc) + 1,
                           ##
                           event.e.pooled = selvar(tdata$event.e, sel.outc),
                           n.e.pooled     = selvar(tdata$n.e, sel.outc),
                           event.c.pooled = selvar(tdata$event.c, sel.outc),
                           n.c.pooled     = selvar(tdata$n.c, sel.outc),
                           ##
                           TE.pooled    = selvar(tdata$TE, sel.outc),
                           lower.pooled = selvar(tdata$lower, sel.outc),
                           upper.pooled = selvar(tdata$upper, sel.outc),
                           ##
                           weight.pooled = selvar(tdata$weight, sel.outc),
                           ##
                           Z.pooled       = selvar(tdata$Z, sel.outc),
                           pval.TE.pooled = selvar(tdata$pval.TE, sel.outc),
                           ##
                           Q      = selvar(tdata$Q, sel.outc),
                           pval.Q = selvar(tdata$pval.Q, sel.outc),
                           I2     = selvar(tdata$I2, sel.outc),
                           tau2   = selvar(tdata$tau2, sel.outc),
                           ##
                           Q.w      = selvar(tdata$Qint, sel.outc),
                           pval.Q.w = selvar(tdata$pval.Qint, sel.outc),
                           I2.w     = selvar(tdata$I2.Qint, sel.outc),
                           ##
                           swap.events = selvar(tdata$swap.events, sel.outc),
                           enter.n     = selvar(tdata$enter.n, sel.outc),
                           logscale    = selvar(tdata$logscale, sel.outc),
                           ##
                           label.e     = selvar(tdata$label.e, sel.outc),
                           label.c     = selvar(tdata$label.c, sel.outc),
                           label.left  = selvar(tdata$label.left, sel.outc),
                           label.right = selvar(tdata$label.right, sel.outc)
                           ),
              group =
                if (sum(sel.grp) > 0 ) {
                  data.frame(comp.no    = selvar(tdata$comp.no, sel.grp),
                             outcome.no = selvar(tdata$outcome.no, sel.grp),
                             group.no   = selvar(tdata$group.no, sel.grp),
                             grplab     = selvar(tdata$author, sel.grp)
                             )}
                else {
                  NA
                             },
              study =
                data.frame(comp.no    = selvar(tdata$comp.no, sel.study),
                           outcome.no = selvar(tdata$outcome.no, sel.study),
                           group.no   = selvar(tdata$group.no, sel.study),
                           studlab    = selvar(tdata$author, sel.study),
                           year       = selvar(tdata$year, sel.study),
                           ##
                           event.e = selvar(tdata$event.e, sel.study),
                           n.e     = selvar(tdata$n.e, sel.study),
                           event.c = selvar(tdata$event.c, sel.study),
                           n.c     = selvar(tdata$n.c, sel.study),
                           ##
                           mean.e = selvar(tdata$mean.e, sel.study),
                           sd.e   = selvar(tdata$sd.e, sel.study),
                           mean.c = selvar(tdata$mean.c, sel.study),
                           sd.c   = selvar(tdata$sd.c, sel.study),
                           ##
                           O.E = selvar(tdata$O.E, sel.study),
                           V   = selvar(tdata$V, sel.study),
                           ##
                           TE   = selvar(tdata$TE, sel.study),
                           seTE = selvar(tdata$seTE, sel.study),
                           ##
                           lower.TE = selvar(tdata$lower, sel.study),
                           upper.TE = selvar(tdata$upper, sel.study),
                           weight   = selvar(tdata$weight, sel.study),
                           ##
                           order = selvar(tdata$order, sel.study)
                           )
              )
  ##
  if (is.data.frame(res$group))
    res2 <- merge(res$study, res$group,
                  by = c("comp.no", "outcome.no", "group.no"),
                  all.x = TRUE)
  else
    res2 <- res$study
  res2 <- merge(res2, res$outcome,
                by = c("comp.no", "outcome.no"))
  res2 <- merge(res2, res$comparison,
                by = c("comp.no"))
  ##
  sel.nocont <- (res2$mean.e == 0 & res2$sd.e == 0 &
                 res2$mean.c == 0 & res2$sd.c == 0)
  res2$mean.e[sel.nocont] <- NA
  res2$sd.e[sel.nocont]   <- NA
  res2$mean.c[sel.nocont] <- NA
  res2$sd.c[sel.nocont]   <- NA
  ##
  if (is.factor(res2$O.E))
    res2$O.E <- numchar(res2$O.E)
  if (is.factor(res2$V))
    res2$V <- numchar(res2$V)
  sel.noOE <- res2$O.E == 0 & res2$V == 0
  res2$O.E[sel.noOE] <- NA
  res2$V[sel.noOE]   <- NA
  ##
  res2$complab <- as.character(res2$complab)
  res2$outclab <- nachar(as.character(res2$outclab))
  res2$studlab <- nachar(as.character(res2$studlab))
  if (is.data.frame(res$group)) {
    res2$grplab <- nachar(as.character(res2$grplab))
    ## Replace some XML code:
    res2$grplab <- gsub("\\&lt;", "<", res2$grplab)
    res2$grplab <- gsub("\\&gt;", ">", res2$grplab)
  }
  ##
  res2$sm     <- as.character(res2$sm)
  res2$method <- as.character(res2$method)
  res2$type   <- substring(res2$type, 1, 1)
  res2$model  <- as.character(res2$model)
  ##
  res2$label.e     <- nachar(as.character(res2$label.e))
  res2$label.c     <- nachar(as.character(res2$label.c))
  res2$label.left  <- nachar(as.character(res2$label.left))
  res2$label.right <- nachar(as.character(res2$label.right))
  ##
  if (numbers.in.labels & length(res2$comp.no) > 0) {
    res2$complab <- paste0(res2$comp.no, " ", res2$complab)
    res2$outclab <- paste0(res2$comp.no, ".", res2$outcome.no, " ",
                           res2$outclab)
    if (is.data.frame(res$group))
      res2$grplab <- paste0(res2$comp.no, ".", res2$outcome.no, ".",
                            res2$group.no, " ", res2$grplab)
  }
  ##
  res2 <- res2[order(res2$comp.no,
                     res2$outcome.no,
                     res2$group.no),]
  ##
  attr(res2, "title") <- res$title
  ##
  if (length(res2$comp.no) > 0) {
    for (i in unique(res2$comp.no)) {
      for (j in unique(res2$outcome.no[res2$comp.no == i])) {
        sel2 <- res2$comp.no == i & res2$outcome.no == j
        if (unique(res2$sm[sel2]) == "OTHER")
          warning("Summary measure unclear for outcome '",
                  ##paste0(unique(res2$comp.no[sel2]), ".",
                  ##       unique(res2$outcome.no[sel2])),
                  unique(res2$outclab[sel2]),
                  "'. Please use argument 'sm' in function metacr ",
                  "to choose adequate summary measure.")
      }
    }
  }
  
  
  ##res2$TE[res2$sm == "HR"] <- log(res2$TE[res2$sm == "HR"])
  ##res2$sm[res2$sm == "log HR"] <- "HR"

  iswap.events <- charmatch(tolower(as.character(res2$swap.events)),
                            c("yes", "no"), nomatch = NA)
  res2$swap.events <- ifelse(iswap.events == 1, TRUE,
                      ifelse(iswap.events == 2, FALSE, NA))
  ##
  ienter.n <- charmatch(tolower(as.character(res2$enter.n)),
                        c("yes", "no"), nomatch = NA)
  res2$enter.n <- ifelse(ienter.n == 1, TRUE,
                  ifelse(ienter.n == 2, FALSE, NA))
  ##
  ilogscale <- charmatch(tolower(as.character(res2$logscale)),
                         c("no", "yes"), nomatch = NA)
  res2$logscale <- ifelse(ilogscale == 1, TRUE,
                   ifelse(ilogscale == 2, FALSE, NA))
  
  attr(res2, "version") <- packageDescription("meta")$Version
  
  class(res2) <- c("rm5", "data.frame")
  res2
}


read.rm5.rm5 <- function(file, title, numbers.in.labels = TRUE, debug = 0) {
  
  
  ##
  ## Read data file
  ##
  rdata <- readLines(file, warn = FALSE)
  
  
  ##
  ## Convert extended ASCII codes
  ##
  rdata <- oct2txt(rdata)
  
  
  ## 
  ## Deal with missing line breaks
  ##
  if (length(rdata) == 1 && sum(gregexpr(">", rdata)[[1]] > 0)) {
    rdata <- unlist(strsplit(rdata, "(?<=>)", perl = TRUE))
  }

  
  ##
  ## Extract title
  ##
  if (missing(title)) {
    sel1 <- grep("<TITLE", rdata)
    sel2 <- grep("</TITLE", rdata)
    title.orig <- rdata[unique(c(sel1[1], sel2[1]))]
    title <- xml_text(as_xml_document(paste(title.orig, collapse = " ")))
  }
  if (debug) {
    cat("\n***** TITLE *****\n")
    cat(paste0("Title: ", title, "\n"))
    ##
    if (debug == 3)
      print(title.orig)
  }
  
  
  ##
  ## Determine whether study data are available
  ##
  sel.data1 <- grepl("<ANALYSES_AND_DATA", rdata)
  sel.data2 <- grepl("</ANALYSES_AND_DATA", rdata)
  ##
  if (!any(sel.data2)) {
    warning("Cochrane review does not contain any study data.",
            call. = FALSE)
    return(NULL)
  }
  ##
  sel.data1 <- grep("<ANALYSES_AND_DATA", rdata)
  sel.data2 <- grep("</ANALYSES_AND_DATA", rdata)
  
  
  ##
  ## Extract data (as vector of character strings)
  ##
  txt <- rdata[(sel.data1 + 1):(sel.data2 - 1)]
  ##
  sel.data1 <- grepl("<DICH_DATA", txt)
  sel.data2 <- grepl("<CONT_DATA", txt)
  sel.data3 <- grepl("<IV_DATA", txt)
  sel.data4 <- grepl("<IPD_DATA", txt)
  ##
  if (!any(sel.data1 | sel.data2 | sel.data3 | sel.data4)) {
    warning("Cochrane review does not contain any usable study data.",
            call. = FALSE)
    return(NULL)
  }
  
  
  res <- data.frame(comp.no = NA, outcome.no = NA, group.no = NA,
                    id = "",
                    event.e = NA, n.e = NA, event.c = NA, n.c = NA,
                    mean.e = NA, sd.e = NA, mean.c = NA, sd.c = NA,
                    O.E = NA, V = NA,
                    TE = NA, seTE = NA, lower.TE = NA, upper.TE = NA,
                    weight = NA, order = NA, level = NA,
                    grplab = "",
                    overall = NA,
                    test.subgroup = NA,
                    type = "", method = "", sm = "", model = "",
                    common = NA, random = NA,
                    outclab = "",
                    k = NA,
                    event.e.pooled = NA, n.e.pooled = NA,
                    event.c.pooled = NA, n.c.pooled = NA,
                    TE.pooled = NA, lower.pooled = NA, upper.pooled = NA,
                    level.ma = NA,
                    weight.pooled = NA,
                    Z.pooled = NA, pval.TE.pooled = NA,
                    Q = NA, pval.Q = NA,
                    I2 = NA, tau2 = NA,
                    Q.w = NA, pval.Q.w = NA, I2.w = NA,
                    swap.events = NA, logscale = NA,
                    label.e = "", label.c = "",
                    label.left = "", label.right = "",
                    complab = "")
  
  
  ##
  ##
  ## Comparisons
  ##
  ##
  sel.comp1 <- grep("<COMPARISON", txt)
  sel.comp2 <- grep("</COMPARISON", txt)
  ##
  if (length(sel.comp1) != length(sel.comp2))
    stop("Malformed XML file (tag: COMPARISON)")
  ##
  for (i in seq(along = sel.comp1)) {
    txt.i <- txt[sel.comp1[i]:sel.comp2[i]]
    sel.data.i <- c(grep("</NAME>", txt.i))
    if (length(sel.data.i) > 0)
      sel.i <- c(1:(max(c(2, min(sel.data.i)))), length(txt.i))
    else
      sel.i <- c(1:2, length(txt.i))
    ##
    xml.i <- as_xml_document(paste(txt.i[sel.i], collapse = " "))
    ##
    comp.no <- as.numeric(xml_attr(xml.i, "NO"))
    complab <- xml_text(xml_find_all(xml.i, "//NAME"))
    ##
    if (debug) {
      cat("\n**** COMPARISON ****\n")
      cat(paste0("Comp.no: ", comp.no, "\n"))
      cat(paste0("complab: ", complab, "\n"))
      ##
      if (debug == 3)
        print(txt.i[sel.i])
    }
    ##
    ## Binary outcomes
    ##
    res <- extract_outcomes(txt.i, "DICH", res,
                            comp.no, complab,
                            debug = debug)
    ##
    ## Continuous outcomes
    ##
    res <- extract_outcomes(txt.i, "CONT", res,
                            comp.no, complab,
                            debug = debug)
    ##
    ## IV outcomes
    ##
    res <- extract_outcomes(txt.i, "IV", res,
                            comp.no, complab,
                            debug = debug)
    ##
    ## IPD outcomes
    ##
    res <- extract_outcomes(txt.i, "IPD", res,
                            comp.no, complab,
                            debug = debug)
  }
  ##
  res <- res[-1, ]
  
  
  ##
  ## Determine whether information on included studies are available
  ##
  sel.incl1 <- grepl("<INCLUDED_STUDIES ", rdata) |
    grepl("<INCLUDED_STUDIES>", rdata)
  sel.incl2 <- grepl("</INCLUDED_STUDIES>", rdata)
  ##
  if (any(sel.incl2)) {
    ##
    ## Extract study information (as vector of character strings)
    ##
    sel.incl1 <- c(grep("<INCLUDED_STUDIES ", rdata),
                   grep("<INCLUDED_STUDIES>", rdata))
    sel.incl2 <- grep("</INCLUDED_STUDIES>", rdata)
    ##
    study <- rdata[(sel.incl1 + 1):(sel.incl2 - 1)]
    ##
    sel.study1 <- grep("<STUDY", study)
    sel.study2 <- grep("</STUDY", study)
    ##
    id <- studlab <- year <- rep("", length(sel.study1))
    ##
    if ((length(sel.study1) > 0 | length(sel.study2) > 0) &&
        (length(sel.study1) == length(sel.study2))) {
      if (debug)
        cat(paste0("\n* STUDY *\n"))
      for (s in seq(along = sel.study1)) {
        study.s <- study[sel.study1[s]:sel.study2[s]]
        xml.s <- as_xml_document(paste(study.s, collapse = " "))
        ##
        id[s] <- xml_attr(xml.s, "ID")
        studlab[s] <- xml_attr(xml.s, "NAME")
        year[s] <- xml_attr(xml.s, "YEAR")
        ##
        if (debug >= 2) {
          cat("- id:\n")
          print(id)
          cat("- studlab:\n")
          print(studlab)
          cat("- year:\n")
          print(year)
          ##
          if (debug == 3)
            print(study.s)
        }
      }
      ##
      names.res <- names(res)
      res.study <- data.frame(id, studlab, year)
      res <- merge(res, res.study, by = "id",
                   all.x = TRUE)
      ##
      names.res <- c(names.res[1:4], "studlab", "year",
                     names.res[5:length(names.res)])
      ##
      res <- res[, names.res]
    }
  }
  
  
  res <- res[order(res$comp.no, res$outcome.no, res$group.no), ]
  ##
  if (numbers.in.labels & length(res$comp.no) > 0) {
    res$complab <- paste0(res$comp.no, " ", res$complab)
    res$outclab <- paste0(res$comp.no, ".", res$outcome.no, " ",
                          res$outclab)
    sel.grp <- !is.na(res$group.no)
    if (any(sel.grp))
      res$grplab[sel.grp] <-
        paste0(res$comp.no, ".", res$outcome.no, ".",
               res$group.no, " ", res$grplab)[sel.grp]
  }
  ##
  res$test.subgroup[is.na(res$group.no)] <- NA
  ##
  row.names(res) <- 1:nrow(res)
  ##
  res$id <- as.character(res$id)
  res$studlab <- as.character(res$studlab)
  res$year <- as.character(res$year)
  res$grplab <- as.character(res$grplab)
  res$type <- as.character(res$type)
  res$method <- as.character(res$method)
  res$sm <- as.character(res$sm)
  res$model <- as.character(res$model)
  res$label.e <- as.character(res$label.e)
  res$label.c <- as.character(res$label.c)
  res$label.left <- as.character(res$label.left)
  res$label.right <- as.character(res$label.right)
  ##
  ## sel.rel <- is_relative_effect(res$sm)
  ## res$TE[sel.rel] <- log(res$TE[sel.rel])
  ## res$lower.TE[sel.rel] <- log(res$lower.TE[sel.rel])
  ## res$upper.TE[sel.rel] <- log(res$upper.TE[sel.rel])
  ## res$seTE <- TE.seTE.ci(res$lower.TE, res$upper.TE, res$level)$seTE
  ##
  res$fixed <- res$common
  ##
  attr(res, "title") <- title
  attr(res, "version") <- packageDescription("meta")$Version
  
  class(res) <- c("rm5", "data.frame")
  res
}
