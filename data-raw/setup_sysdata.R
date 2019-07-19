internal_state <- new.env(parent = emptyenv())

internal_state$NOAAsites <- NULL

save(internal_state, file="../R/sysdata.rda", version=2)
