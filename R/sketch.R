# Sketches for Linear Regression
#
# Copyright: Leo Geppert, Ludger Sandig
# License: GPL v3

sketch <- function(data, file, epsilon = NULL, obs_sketch = NULL, warn = TRUE,
                   affine = TRUE, method = c("R", "S", "C"), ...) {
    if (missing(data) & missing(file))
        stop ("One of data or file has to be specified.")
    if (!missing(file) & !missing(data))
        warning("Both data and file specified, file ignored.")
    if (!missing(data))
        if (!is.data.frame(data) && !is.matrix(data))
            stop("data must be a data.frame or matrix.")
    if (!missing(file) & missing(data)) {
        if (!file.exists(file))
            stop("File not found.")
        data <- read.table(file, ...)
    }
    if (!is.numeric(epsilon) && !is.numeric(obs_sketch))
        stop ("One of epsilon and obs_sketch must be numeric.")
    if (!is.logical(affine))
        stop ("affine must be logical.")
    method <- match.arg(method)
    if (is.numeric(epsilon) && (epsilon > 0.5 || epsilon <= 0))
        stop("Possible values for epsilon lie in the interval (0, 0.5].")
    if (warn) {
        n <- dim(data)[1]
        d <- dim(data)[2] + affine
        if (is.numeric(epsilon)) {
            s <- if (method == "R" || method == "S") {
                     as.integer(ceiling(d * log(d) / (epsilon * epsilon)))
                 } else if (method == "C") {
                     as.integer(2^ceiling(log(ceiling(d*d/(20*epsilon*epsilon)),
                                              base = 2)))
                 }
            if (s > n)
                warning(paste0("Sketch has more rows than original data set, ",
                               "consider using larger epsilon with method ",
                               method, "."))
        } else {
            if (method == "C") {
                next_power_of_two <- as.integer(2^ceiling(log(obs_sketch,
                                                              base = 2)))
                if (obs_sketch != next_power_of_two) {
                    warning(paste0("Rounding up obs_sketch from ", obs_sketch,
                                   " to ", next_power_of_two,
                                   " for method 'C'."))
                    obs_sketch <- next_power_of_two
                }
            }
            if (obs_sketch > n)
                warning(paste0("Sketch has more rows than original data set, ",
                               "consider using smaller obs_sketch with method ",
                               method, "."))
        }
    }
    if (is.numeric(epsilon) && is.numeric(obs_sketch))
        warning(paste0("Both epsilon and obs_sketch specified.\n",
                       "Used epsilon only, obs_sketch was ignored."))

    # ensure right storage mode and class before calling to C
    data <- as.matrix(data)
    storage.mode(data) <- "double"
    sketch_cols <- ncol(data) + ifelse(affine, 1, 0)
    # add intercept column if requested
    if (affine) {
        data <- cbind(x0 = rep(1, nrow(data)), data)
    }
    # calculate sketches in a .Call() to the shared object
    if (method == "R") {
        if (!is.null(epsilon)) {
            sketch_rows <- ceiling(sketch_cols * log(sketch_cols) / (epsilon^2))
        } else {
            sketch_rows <- obs_sketch
        }
        the_sketch <- .Call("sketch_rad", data, as.integer(sketch_rows))
    } else if (method == "S") {
        if (!is.null(epsilon)) {
            sketch_rows <- ceiling(sketch_cols * log(sketch_cols) / (epsilon^2))
        } else {
            sketch_rows <- obs_sketch
        }
        the_sketch <- .Call("sketch_srht", data, as.integer(sketch_rows))
    } else { # method "C"
        if (!is.null(epsilon)) {
            sketch_rows <- 2^ceiling(log(ceiling(sketch_cols^2/(20*epsilon^2)),
                                         base = 2))
        } else {
            sketch_rows <- obs_sketch
        }
        the_sketch <- .Call("sketch_cw", data, as.integer(sketch_rows))
    }
    colnames(the_sketch) <- colnames(data)
    return(as.data.frame(the_sketch))
}

readinandsketch <- function(file, nrows = 50000L, epsilon = NULL,
                            obs_sketch = NULL, affine = TRUE,
                            method= c("C", "S", "R"), header = FALSE,
                            sep = "", col.names, skip = 0,
                            warn = FALSE, ...) {
    if (!file.exists(file))
        stop("File not found.")
    if (!is.numeric(epsilon) && !is.numeric(obs_sketch))
        stop("One of epsilon and obs_sketch must be numeric.")
    if (!is.logical(affine))
        stop("affine must be logical.")
    method <- match.arg(method)
    if (is.numeric(epsilon) && is.numeric(obs_sketch))
        warning(paste0("Both epsilon and obs_sketch specified.\n",
                       "Used epsilon only, obs_sketch was ignored."))
    if (nrows <= 0)
        stop("nrows has to be a positive integer.")
    epsgiven <- is.numeric(epsilon)
    con <- file(file, open = "r")
    if (header) {
        hdr <- readLines(con = con, n = skip)  # skip junk lines
        hdr <- readLines(con = con, n = 1)     # real header
        hdr <- as.character(read.table(text = hdr, sep = sep,  # parse
                                       stringsAsFactors = FALSE)[1, ])
    }
    # calculate dimensions of sketch
    s_columns <- ifelse(affine, length(hdr) + 1, length(hdr))
    if (!epsgiven) {
        if (method == "C") {
            s_rows <- 2^ceiling(log(obs_sketch, base = 2))
            warning(paste0("Method CW can only generate sketch sizes",
                           "that are powers of two.\n",
                           "Rounding up."))
        } else {
            s_rows <- obs_sketch
        }
    } else {
        if (method == "R" || method == "S")
            s_rows <- ceiling(s_columns * log(s_columns) / (epsilon^2))
        else  # method == "C"
            s_rows <- 2^ceiling(log(ceiling(s_columns^2 / (20 * epsilon^2)),
                                    base = 2))
                  
    }
    # create emtpy data frame, fix names, and accumulate sketch
    skizze <- data.frame(matrix(data = 0, nrow = s_rows, ncol = s_columns))
    names(skizze) <- c(if(affine) "x0", hdr)
    while (length(txt <- readLines(con = con, n = nrows)) > 0) {
        if (method == "S" && length(txt) < s_rows) {
            warning(paste0("Performing SRHT on final block is impossible:",
                           "Fewer rows than in sketch.\n",
                           "Ignoring this block."))
            break
        }
        nextblock <- read.table(text = txt, col.names = hdr,
                                sep = sep, ...)
        skizze <- skizze + 
            if (epsgiven) {
                sketch(data = nextblock, epsilon = epsilon, warn = warn,
                       affine = affine, method = method)
            } else {
                sketch(data = nextblock, obs_sketch = obs_sketch, warn = warn,
                       affine = affine, method = method)
            }
    }
    close(con)
    return(skizze)
}
