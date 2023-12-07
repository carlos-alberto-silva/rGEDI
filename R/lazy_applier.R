lazy_call <- function(func) {
    is_formula <- tryCatch(lazyeval::is_formula(func), error = function(e) FALSE)
    if (!is_formula) func <- lazyeval::f_capture(func)

    func <- lazyeval::f_interp(func)
    call <- lazyeval::as_call(func)

    call
}

lazy_apply_dt_call <- function(dt, call, group.by = "") {
    requireNamespace("data.table")
    inner.env <- new.env()
    func.names.args <- all.names(call)
    fn.name <- func.names.args[1]
    dt.name <- paste0(func.names.args, collapse = "")

    fn.is.function <- tryCatch(is.function(get(fn.name)), error = function(e) FALSE)
    if (! fn.is.function) {
        inner.env[[fn.name]] <- parent.frame(2)[[fn.name]]
        stopifnot("Could not find function" = is.function(inner.env[[fn.name]]))
    }
    inner.env[[dt.name]] <- data.table::as.data.table(dt)
    result <- eval(parse(text = paste0(dt.name, "[, ", deparse(call), ", ", group.by, "]")), inner.env)
    return(result)
}
