
#transformation can be
# NULL
# a single transformation object (that applies to all variables)
# a single character name of a transformation (that applies to all variables)
# a named list of any length or an unnamed list of length n.var whose elements are either
#    NULL, a transformation object, or a character name of a transformation object
# a named character vector of any length or an unnamed character vector of length n.var

#' exporting for now for debugging
#'@export
process.transformations <- function(var.names, transformations, n.var=length(var.names),
                                    var.names.name='var.names')
{
    parsed.transformations = lapply(var.names, function(name){distributions::get.defined.transformation('identity')})
    names(parsed.transformations) = names(transformations)
    transformations = as.list(transformations)

    if (is.null(transformations) || length(transformations)==0)
    {}
    else if (is.null(var.names))
    {
        if (length(transformations)==1)
            parsed.transformations = lapply(1:n.var, function(i){transformations[[1]]})
        else if (length(transformations)==n.var)
            parsed.transformations = transformations
        else
            stop(paste0("If ", var.names.name, " are not specified, the transformations parameter must have either length 1 or the same length as the number of variables"))
    }
    else if (is.null(names(transformations)))
    {
        if (length(transformations)==1)
            parsed.transformations = lapply(var.names, function(name){transformations[[1]]})
        else if (length(transformations) == length(var.names))
        {
            parsed.transformations = transformations
            names(parsed.transformations) = var.names
        }
        else
            stop(paste0("If names are not specified for the 'transformations' parameter, it must have the same length as 'var.names' (",
                        length(var.names), ")"))
    }
    else
    {
        unused.transformation.names = setdiff(names(transformations), var.names)
        if (length(unused.transformation.names) > 0)
            warning(paste0("The following transformation",
                           ifelse(length(unused.transformation.names>1, "s were", 'was')),
                           " given but ",
                           ifelse(length(unused.transformation.names>1, "are not variables", 'is not a variable')),
                           " specified in '", var.names.name, "': ",
                           paste0("'", unused.transformation.names, "'", collapse=", ")
            ))

        parsed.transformations[names(transformations)] = transformations[names(transformations)]
    }

    transformations = lapply(parsed.transformations, function(one.transformation){
        if (is.null(one.transformation))
            distributions::get.defined.transformation('identity')
        else if (is(one.transformation, 'transformation'))
            one.transformation
        else if (is(one.transformation, 'character'))
        {
            if (length(one.transformation)!=1)
                stop("Only single names (ie not vectors) are permitted as elements of 'transformations")

            distributions::get.defined.transformation(one.transformation, throw.error.if.no.match = T)
        }
        else
            stop("transformations must be either NULL, an instance of the 'transformation' class, or a character name of a defined transformation")
    })

    names(transformations) = var.names
    transformations
}

do.transform.parameters <- function(control, parameters)
{
    rv = sapply(1:control@n.var, function(i){
        if (is.null(control@transformations[[i]]))
            parameters[i]
        else
            control@transformations[[i]]@transform(parameters[i])
    })

    if (any(is.na(rv)))
    {
        na.mask = is.na(rv)
        if (is.null(control@var.names))
        {
            stop(paste0("The transformation",
                        ifelse(sum(na.mask)>1, "s", ""),
                        " for the ",
                        get.ordinal.list.string((1:control@n.var)[na.mask]),
                        "parameter",
                        ifelse(sum(na.mask)>1, "s", ""),
                        "produced NA value",
                        ifelse(sum(na.mask)>1, "s", ""),
                        ". The input parameter",
                        ifelse(sum(na.mask)>1, "s were", " was"),
                        ": ",
                        paste0(parameters[na.mask], collapse=", ")))
        }
        else
        {
            stop(paste0("The transformation",
                        ifelse(sum(na.mask)>1, "s", ""),
                        " for the following parameter",
                        ifelse(sum(na.mask)>1, "s", ""),
                        "produced NA value",
                        ifelse(sum(na.mask)>1, "s", ""),
                        ": ",
                        paste0("'", control@var.names[na.mask], "'", collapse=', '),
                        ". The input parameter",
                        ifelse(sum(na.mask)>1, "s were", " was"),
                        ": ",
                        paste0(parameters[na.mask], collapse=", ")))
        }
    }

    if (any(is.infinite(rv)))
    {
        inf.mask = is.infinite(rv)
        if (is.null(control@var.names))
        {
            stop(paste0("The transformation",
                        ifelse(sum(inf.mask)>1, "s", ""),
                        " for the ",
                        get.ordinal.list.string((1:control@n.var)[inf.mask]),
                        "parameter",
                        ifelse(sum(inf.mask)>1, "s", ""),
                        "produced infinite value",
                        ifelse(sum(inf.mask)>1, "s", ""),
                        ". The input parameter",
                        ifelse(sum(inf.mask)>1, "s were", " was"),
                        ": ",
                        paste0(parameters[inf.mask], collapse=", ")))
        }
        else
        {
            stop(paste0("The transformation",
                        ifelse(sum(inf.mask)>1, "s", ""),
                        " for the following parameter",
                        ifelse(sum(inf.mask)>1, "s", ""),
                        "produced infinite value",
                        ifelse(sum(inf.mask)>1, "s", ""),
                        ": ",
                        paste0("'", control@var.names[inf.mask], "'", collapse=', '),
                        ". The input parameter",
                        ifelse(sum(inf.mask)>1, "s were", " was"),
                        ": ",
                        paste0(parameters[inf.mask], collapse=", ")))
        }
    }

    rv
}

do.reverse.transform.parameters <- function(control, transformed.parameters)
{
    rv = sapply(1:control@n.var, function(i){
        if (is.null(control@transformations[[i]]))
            parameters[i]
        else
            control@transformations[[i]]@reverse.transform(transformed.parameters[i])
    })

    if (any(is.na(rv)))
    {
        na.mask = is.na(rv)
        if (is.null(control@var.names))
        {
            stop(paste0("The reverse transformation",
                        ifelse(sum(na.mask)>1, "s", ""),
                        " for the ",
                        get.ordinal.list.string((1:control@n.var)[na.mask]),
                        "parameter",
                        ifelse(sum(na.mask)>1, "s", ""),
                        "produced NA value",
                        ifelse(sum(na.mask)>1, "s", ""),
                        ". The input parameter",
                        ifelse(sum(na.mask)>1, "s were", " was"),
                        ": ",
                        paste0(parameters[na.mask], collapse=", ")))
        }
        else
        {
            stop(paste0("The reverse transformation",
                        ifelse(sum(na.mask)>1, "s", ""),
                        " for the following parameter",
                        ifelse(sum(na.mask)>1, "s", ""),
                        "produced NA value",
                        ifelse(sum(na.mask)>1, "s", ""),
                        ": ",
                        paste0("'", control@var.names[na.mask], "'", collapse=', '),
                        ". The input parameter",
                        ifelse(sum(na.mask)>1, "s were", " was"),
                        ": ",
                        paste0(parameters[na.mask], collapse=", ")))
        }
    }

    rv
}
