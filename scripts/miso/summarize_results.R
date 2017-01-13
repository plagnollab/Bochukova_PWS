list.samples <- c("Control1", "Control2", "Control3", "Control4", "PWS1", "PWS2", "PWS3", "PWS4")

for (type in c("A3SS", "A5SS", "MXE", "RI", "SE")[1]) {

  message(type)
  list.events <- c()
  for (sample in list.samples) {

    my.tab <- read.table(paste0("MISO/summary/", type, "_", sample, "/summary/", type, ".miso_summary"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    #message(nrow(my.tab), " ", sample)
    list.events <- c(list.events, my.tab$event_name)
  }

  unique.events <- unique(list.events)
  

  res.frame <- data.frame(event = unique.events)
  for (sample in list.samples) {
    my.tab <- read.table(paste0("MISO/summary/", type, "_", sample, "/summary/", type, ".miso_summary"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    res.frame[, sample ] <- my.tab$miso_posterior_mean[ match(res.frame$event, table = my.tab$event) ]
  }

  
  control.count <- as.matrix(res.frame[, grepl(pattern = "Control", names(res.frame))])
  PWS.count <- as.matrix(res.frame[, grepl(pattern = "PWS", names(res.frame))])


  n.controls <-  apply(control.count, MAR = 1, FUN = function(x) {sum(!is.na(x))})
  n.PWS <-  apply(PWS.count, MAR = 1, FUN = function(x) {sum(!is.na(x))})
  good.rows <- n.controls >= 4 & n.PWS >= 4

  ## now we subset for the rows with enough data
  PWS.count <- PWS.count[good.rows,]
  control.count <- control.count[good.rows,]
  n.controls <- n.controls[ good.rows ]
  n.PWS <- n.PWS[ good.rows ]
  event.names <- res.frame$event[ good.rows ]
  
  ### and now we compute
  delta.psi <- rowMeans(PWS.count, na.rm = TRUE) - rowMeans(control.count, na.rm = TRUE)
  sd.controls <- apply( control.count, MAR = 1, FUN = sd, na.rm = TRUE)
  sd.PWS <- apply( PWS.count, MAR = 1, FUN = sd, na.rm = TRUE)

  s.combined <- sqrt( sd.controls^2/n.controls  + sd.PWS^2/n.PWS)
  t <- delta.psi  / s.combined

  sig.results <- which(delta.psi > 0.2 & abs(t) > 3)
  good.results <- data.frame(event = event.names[ sig.results ],
                             delta.psi = delta.psi[ sig.results ],
                             t = t [ sig.results])
  
    
  message("Number of events: ", length(delta.psi), " and median value: ", median(delta.psi))
  
  
}
