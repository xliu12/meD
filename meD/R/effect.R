

effect <- function(Yrt, Yt1only=FALSE, btt=NULL, bm=NULL) {
  effects <- within(Yrt, {
    toD <- (Yt1r1 - Yt0r1) - (Yt1r0 - Yt0r0)
    meD <- (Yt1r1 - Yt0r1 ) - (Yt1r1.Mt1r0 - Yt0r1.Mt0r0)
    reD <- (Yt1r1.Mt1r0 - Yt0r1.Mt0r0) - (Yt1r0 - Yt0r0)
    
    if (length(btt)>0) {
      meD <- btt["R"]*(bm["ttM"]+bm["ttRM"]) + btt["ttR"]*(bm["M"]+bm["RM"]) + btt["ttR"]*(bm["ttM"]+bm["ttRM"])
      # reD <- bm["ttR"] + bm["RM"]*(btt["tt"]+btt["ttR"]*0) + bm["ttRM"]*(btt["intercept"])
    }
  })
  
  
  if (Yt1only) {
    effects <- within(Yrt, {
      toD <- (Yt1r1 - Yt1r0) 
      meD <- (Yt1r1 - Yt1r1.Mt1r0) 
      reD <- (Yt1r1.Mt1r0 - Yt1r0) 
    })
  }
  
  
  effects
}