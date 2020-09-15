model.frame.vglmMRSea <- function (formula, data = NULL, subset = NULL, na.action = na.fail, drop.unused.levels = FALSE, xlev = NULL, splineParams=NULL, ...) 
{

  possible_newdata <- !missing(data) && is.data.frame(data) && 
    identical(substitute(data), quote(newdata)) && (nr <- nrow(data)) > 
    0
  if (!missing(formula) && nargs() == 1 && is.list(formula) && 
      !is.null(m <- formula$model)) 
    return(m)
  if (!missing(formula) && nargs() == 1 && is.list(formula) && 
      all(c("terms", "call") %in% names(formula))) {
    fcall <- formula$call
    print(fcall)
    m <- match(c("formula", "data", "subset", "weights", 
                 "na.action"), names(fcall), 0)
    fcall <- fcall[c(1, m)]
    fcall[[1L]] <- quote(stats::model.frame)
    env <- environment(formula$terms)
    if (is.null(env)) 
      env <- parent.frame()
    return(eval(fcall, env))
  }
  
  if (missing(formula)) {
    if (!missing(data) && inherits(data, "data.frame") && 
        length(attr(data, "terms"))) 
      return(data)
    formula <- as.formula(data)
  }
  else if (missing(data) && inherits(formula, "data.frame")) {
    if (length(attr(formula, "terms"))) 
      return(formula)
    data <- formula
    formula <- as.formula(data)
  }
  if(length(formula)!=length(as.formula(formula))){
    data<-formula$data
    splineParams<-formula$splineParams
  }
  formula <- as.formula(formula)
  
  if(!is.null(splineParams)){
    dists<-splineParams[[1]]$dist
    aR<-splineParams[[1]]$knotPos
    radii<-splineParams[[1]]$radii
    radiusIndices<-splineParams[[1]]$radiusIndices
  }
  
  if (missing(na.action)) {
    if (!is.null(naa <- attr(data, "na.action")) & mode(naa) != 
        "numeric") 
      na.action <- naa
    else if (!is.null(naa <- getOption("na.action"))) 
      na.action <- naa
  }
  
  if (missing(data)) 
    data <- environment(formula)
  else if (!is.data.frame(data) && !is.environment(data) && 
           !is.null(attr(data, "class"))) 
    data <- as.data.frame(data)
  else if (is.array(data)) 
    stop("'data' must be a data.frame, not a matrix or an array")
  if (!inherits(formula, "terms")) 
    formula <- terms(formula, data = data)
  if(!is.null(splineParams)){
    env<-environment()
  }else{
    env <- environment(formula)  
  }
  
  rownames <- .row_names_info(data, 0L)
  vars <- attr(formula, "variables")
  predvars <- attr(formula, "predvars")
  if (is.null(predvars)) 
    predvars <- vars
  varnames <- sapply(vars, function(x) paste(deparse(x, width.cutoff = 500), 
                                             collapse = " "))[-1L]
  variables <- eval(predvars, data, env)
  resp <- attr(formula, "response")
  
  if (is.null(rownames) && resp > 0L) {
    lhs <- variables[[resp]]
    rownames <- if (is.matrix(lhs)) 
      rownames(lhs)
    else names(lhs)
  }
  
  if (possible_newdata && length(variables)) {
    nr2 <- max(sapply(variables, NROW))
    if (nr2 != nr) 
      warning(sprintf(paste0(ngettext(nr, "'newdata' had %d row", 
                                      "'newdata' had %d rows"), " ", ngettext(nr2, 
                                                                              "but variable found had %d row", "but variables found have %d rows")), 
                      nr, nr2), call. = FALSE, domain = NA)
  }
  
  if (is.null(attr(formula, "predvars"))) {
    for (i in seq_along(varnames)) predvars[[i + 1L]] <- makepredictcall(variables[[i]], 
                                                                         vars[[i + 1L]])
    attr(formula, "predvars") <- predvars
  }
  
  extras <- substitute(list(...))
  extranames <- names(extras[-1L])
  extras <- eval(extras, data, env)
  subset <- eval(substitute(subset), data, env)
  data <- .External2(stats:::C_modelframe, formula, rownames, variables, 
                     varnames, extras, extranames, subset, na.action)
  
  if (length(xlev)) {
    for (nm in names(xlev)) if (!is.null(xl <- xlev[[nm]])) {
      xi <- data[[nm]]
      if (is.character(xi)) 
        xi <- as.factor(xi)
      if (!is.factor(xi) || is.null(nxl <- levels(xi))) 
        warning(gettextf("variable '%s' is not a factor", 
                         nm), domain = NA)
      else {
        ctr <- attr(xi, "contrasts")
        xi <- xi[, drop = TRUE]
        nxl <- levels(xi)
        if (any(m <- is.na(match(nxl, xl)))) 
          stop(sprintf(ngettext(length(m), "factor %s has new level %s", 
                                "factor %s has new levels %s"), nm, paste(nxl[m], 
                                                                          collapse = ", ")), domain = NA)
        data[[nm]] <- factor(xi, levels = xl, exclude = NULL)
        if (!identical(attr(data[[nm]], "contrasts"), 
                       ctr)) 
          warning(gettext(sprintf("contrasts dropped from factor %s", 
                                  nm), domain = NA), call. = FALSE)
      }
    }
  }
  else if (drop.unused.levels) {
    for (nm in names(data)) {
      x <- data[[nm]]
      if (is.factor(x) && length(unique(x[!is.na(x)])) < 
          length(levels(x))) {
        ctr <- attr(x, "contrasts")
        data[[nm]] <- x[, drop = TRUE]
        if (!identical(attr(data[[nm]], "contrasts"), 
                       ctr)) 
          warning(gettext(sprintf("contrasts dropped from factor %s due to missing levels", 
                                  nm), domain = NA), call. = FALSE)
      }
    }
  }
  
  attr(formula, "dataClasses") <- vapply(data, .MFclass, "")
  attr(data, "terms") <- formula
  data
}

