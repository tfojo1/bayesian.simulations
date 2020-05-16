
##-- SUBSETTING MCMC SAMPLES --##
##-----------------------------##

do.subset.mcmc <- function(mcmc,
                           chains=1:mcmc@n.chains,
                           var.names=mcmc@var.names,
                           additional.burn=0,
                           additional.thin=1,
                           exact.var.names=F)
{
    #Match var.names
    var.names = match.var.names(mcmc, var.names, exact.var.names)

    #Check arguments
    chains = check.chains(mcmc, chains)

    if (additional.burn < 0 || additional.burn >= mcmc@n.iter)
        stop("'additional.burn' must be >= 0 and < the number of iteration in the MCMC object (",
             mcmc@n.iter, ")")
    if (additional.thin < 1 || additional.thin > (mcmc@n.iter-additional.burn))
        stop("'additional.thin' must be >= 1 and <= the number of iteration in the MCMC object less the number of additional iterations to burn (",
             mcmc@n.iter-additional.burn, ")")


    #Decide what we're going to keep
    keep.mask = ((1:mcmc@n.iter) > additional.burn) &
        ((1:mcmc@n.iter-additional.burn) %% additional.thin == 0)

    if (sum(keep.mask)==0)
        stop(paste0("No simulations are available after applying an additional burn of ",
                    additional.burn, " and an additional thin of ", additional.thin))

    dim.names = list(chain=as.character(chains),
                     iteration=as.character((1:mcmc@n.iter)[keep.mask]),
                     variable=var.names)

    samples = mcmc@samples[chains, keep.mask, var.names]
    dim(samples) = sapply(dim.names, length)
    dimnames(samples) = dim.names

    simulation.indices = mcmc@simulation.indices[chains, keep.mask]
    dim(simulation.indices) = sapply(dim.names[1:2], length)
    dimnames(simulation.indices) = dim.names[1:2]

    list(samples=samples,
         simulation.indices=simulation.indices,
         n.iter=sum(keep.mask))
}

##--------------------##
##-- MATCHING NAMES --##
##--------------------##

match.var.names <- function(mcmc, var.names, exact.var.names)
{
    do.match.names(to.match=var.names,
                   possible.matches=mcmc@var.names,
                   exact=exact.var.names,
                   to.match.name='var.names',
                   category.descriptor='variable')
}

match.block.names <- function(mcmc, block.names, exact.block.names)
{
    do.match.names(to.match=block.names,
                   possible.matches=mcmc@sample.steps,
                   exact=exact.block.names,
                   to.match.name='block.names',
                   category.descriptor='block')
}

#The same as match.block.names but referring to 'steps' instead of 'blocks
match.step.names <- function(mcmc, step.names, exact.step.names)
{
    do.match.names(to.match=step.names,
                   possible.matches=mcmc@sample.steps,
                   exact=exact.step.names,
                   to.match.name='step.names',
                   category.descriptor='step')
}


do.match.names <- function(to.match, possible.matches, exact,
                           to.match.name='var.names',
                           category.descriptor='variable')
{
    if (is(to.match, 'integer') ||
        (is(to.match, 'numeric') && all(round(to.match)==to.match)))
    {
        if (any(is.na(to.match)))
            stop(paste0("'", to.match, "' cannot contain NA values"))
        else if (any(to.match<1) || any(to.match>length(possible.matches)))
            stop("Indices for '", to.match.name, "' must be between 1 and ", length(possible.matches))

        rv = possible.matches[to.match]
    }
    else
    {
        if (!is(to.match, 'character'))
            stop(paste0(to.match.name, " must be either a character/character vector or an integer vector of indices"))
        else if (length(to.match)==0)
            stop(paste0(to.match.name, " cannot have length 0"))
        else if (any(is.na(to.match)))
            stop(paste0("'", to.match, "' cannot contain NA values"))

        if (exact)
            rv = intersect(to.match, possible.matches)
        else
        {
            regexes = gsub("\\*", ".*", to.match)
            regexes = paste0("^", regexes, ".*$")
            if (length(regexes)==1)
                matches = grepl(regexes, possible.matches, ignore.case = T)
            else
                matches = apply(sapply(regexes, function(regex){
                    grepl(regex, possible.matches, ignore.case = T)
                }), 1, any)

            rv = possible.matches[matches]
        }

        if (length(rv)==0)
        {
            if (exact)
                match.str = 'are named'
            else if (length(possible.matches)==1)
                match.str = 'match the pattern'
            else
                match.str = 'match the patterns'

            if (length(to.match)==1)
                stop(paste0("No ", category.descriptor,
                            "s in the given mcmc ",
                            match.str, " '", to.match, "'"))
            else if (length(to.match)==2)
                stop(paste0("No ", category.descriptor,
                            "s in the given mcmc ",
                            match.str, " '", to.match[1],
                            "' or '", to.match[2], "'"))
            else
                stop(paste0("No ", category.descriptor,
                            "s in the given mcmc ",
                            match.str, " ",
                            paste0("'", to.match[-length(to.match)], "'", collapse=", "),
                            ", or '", to.match[length(to.match)], "'"))
        }
    }

    unique(rv)
}

##-----------------------------##
##-- CHECKING VARIABLE NAMES --##
##-----------------------------##

check.variable.names <- function(v, desired.names=NULL, arg.name.for.error=substitute(v),
                                 desired.length=length(desired.names),
                                 allow.expand.scalar=T)
{
    if (all(class(v) != 'numeric') && all(class(v) != 'integer'))
        stop(paste0("'", arg.name.for.error, "' must be a numeric or integer vector"))
    else if (is.null(names(v)) || is.null(desired.names))
    {
        if (allow.expand.scalar && length(v)==1)
            v = rep(v, desired.length)

        if (length(v) == desired.length)
            names(v) = desired.names
        else
        {
            if (is.null(desired.names))
                stop(paste0("'", arg.name.for.error, "' must be a numeric of length ", desired.length))
            else
                stop(paste0("'", arg.name.for.error, "' must either be a numeric of length ", desired.length,
                            " or a named vector that includes the names ",
                            paste0("'", desired.names, "'", collapse=', ')))
        }
    }
    else
    {
        missing.desired.names = sapply(desired.names, function(name){all(names(v)!=name)})
        if (any(missing.desired.names))
            stop(paste0("'", arg.name.for.error, "' is missing ",
                        ifelse(sum(missing.desired.names)==1, "an element", "elements"),
                        " named ",
                        paste0("'", desired.names[missing.desired.names], "'", collapse=", ")))
        else
            v = v[desired.names]
    }

    v
}

check.cov.mat.names <- function(mat, desired.names, arg.name.for.error=substitute(mat))
{
    if (all(class(mat) != 'matrix') && all(class(mat) != 'array'))
        stop(paste0("'", arg.name.for.error, "' must be a matrix or an array with two dimensions"))
    else if (length(dim(mat))!=2)
        stop(paste0("'", arg.name.for.error, "' must be a matrix or an array with two dimensions"))
    else if (is.null(dimnames(mat)))
    {
        if (dim(mat)[1] == length(desired.names) || dim(mat)[2] == length(desired.names))
            dimnames(mat) = list(desired.names, desired.names)
        else
            stop(paste0("'", arg.name.for.error, "' must either be a ", length(desired.names), "x", length(desired.names),
                        "matrix or have dimnames set for both dimensions that include the names ",
                        paste0("'", desired.names, "'", collapse=', ')))
    }
    else
    {
        missing.desired.names = sapply(desired.names, function(name){all(dimnames(mat)[[1]]!=name)})
        if (any(missing.desired.names))
            stop(paste0("The dimnames for the first dimension of '", arg.name.for.error, "' are missing the ",
                        ifelse(sum(missing.desired.names)==1, "name", "names"),
                        paste0("'", desired.names[missing.desired.names], "'", collapse=", ")))
        else
        {
            missing.desired.names = sapply(desired.names, function(name){all(dimnames(mat)[[2]]!=name)})
            if (any(missing.desired.names))
                stop(paste0("The dimnames for the second dimension of '", arg.name.for.error, "' are missing the ",
                            ifelse(sum(missing.desired.names)==1, "name", "names"),
                            paste0("'", desired.names[missing.desired.names], "'", collapse=", ")))

            mat = mat[desired.names,desired.names]
        }
    }

    mat
}

##-----------------------##
##-- UPDATE PARAMETERS --##
##-----------------------##

#'@export
visualize.update.parameters <- function(base.update=1,
                                        update.decay=0.5,
                                        update.prior.iter=100,
                                        iter.to.summarize=c(1,10,50,100,500,1000))
{
    iter = 1:max(iter.to.summarize)
    weight = base.update / (update.prior.iter + iter)^update.decay
    p0 = cumprod(1-weight)
    cum.weight = 1-p0

    raw.summary.mat = rbind(weight[iter.to.summarize],cum.weight[iter.to.summarize])
    summary.mat = paste0(100*round(raw.summary.mat,2), '%')
    summary.dim.names = list(statistic=c('Iteration Weight','Cumulative Weight'), iteration=as.character(iter.to.summarize))
    summary.mat = array(summary.mat, dim=sapply(summary.dim.names, length), dimnames=summary.dim.names)
    print(summary.mat)

    df = data.frame(Iteration=c(iter,iter),
                    Weight=c(weight,cum.weight),
                    Type=rep(c('Iteration Weight','Cumulative Weight'), each=length(iter)))

    print(ggplot(df, aes(Iteration, Weight, color=Type)) + geom_line(size=2))

    invisible(raw.summary.mat)
}

#'@export
estimate.update.parameters <- function(weight1, iter1, weight2, iter2)
{
    if (iter1 >= iter2 || weight1 >= weight2)
        stop("iter2 must be > iter1, and weight2 must be > weight1")

    objective.fn <- function(params){
        C = exp(params[1])
        prior.iter = exp(params[2])
        alpha = 1 / (1+exp(-params[3]))

        gammas = C/(prior.iter+1:iter2)^alpha

        p0.1 = prod(1-gammas[1:iter1])
        p0.2 = prod(1-gammas)

        guess.weight1 = 1-p0.1
        guess.weight2 = 1-p0.2

        (guess.weight1 - weight1)^2 + (guess.weight2 - weight2)^2
    }

    optim.result = optim(par=c(C=log(1), prior.iter=log(100), alpha=0),
                         fn=objective.fn, control=list(maxit=1000))


    list(base.update = as.numeric(exp(optim.result$par[1])),
         update.prior.iter = as.numeric(exp(optim.result$par[2])),
         update.decay = as.numeric(1 / (1+exp(-optim.result$par[3]))))
}


##----------------------------------##
##-- FOR PRETTY PRINTING MESSAGES --##
##----------------------------------##

get.timespan.text <- function(seconds,
                              max.spans.to.list=2,
                              digits.for.last.span=1,
                              allowed.spans = c('week','day','hour','minute','second'),
                              show.zeros=T)
{
    mult.from.previous.span = c(week=7,
                                day=24,
                                hour=60,
                                minute=60,
                                second=1)

    total.mult.for.span = rev(cumprod(rev(mult.from.previous.span)))

    total.mult.for.span = total.mult.for.span[intersect(names(total.mult.for.span), allowed.spans)]

    sapply(seconds, function(one.seconds){

        spans = one.seconds / total.mult.for.span

        first.keep.span = c((1:length(spans))[spans>1], length(spans))[1]
        keep.spans = intersect(1:length(spans), first.keep.span + 1:max.spans.to.list - 1)

        spans[-1] = spans[-1] %% mult.from.previous.span[-length(spans)]

        spans = spans[keep.spans]

        spans[-length(spans)] = floor(spans[-length(spans)])
        round.factor = 10^digits.for.last.span
        spans[length(spans)] = floor(round.factor * spans[length(spans)]) / round.factor

        if (!show.zeros)
            spans = spans[span>0 || (1:length(spans))==length(spans)]

        span.text = sapply(1:length(spans), function(i){
            paste0(format(spans[i], big.mark=','), " ",
                   names(spans)[i],
                   ifelse(spans[i]==1, '', 's')
            )
        })

        paste0(span.text, collapse=", ")
    })

}


get.ordinal <- function(nums)
{
    sapply(nums, function(num){
        last.two = as.numeric(num) %% 100
        last.one = as.numeric(num) %% 10

        if (last.two >= 11 && last.two <= 13)
            paste0(num, 'th')
        else if (last.one==1)
            paste0(num, 'st')
        else if (last.one==2)
            paste0(num, 'nd')
        else if (last.one==3)
            paste0(num, 'rd')
        else
            paste0(num, 'th')
    })
}

get.ordinal.list.string <- function(nums)
{
    if (length(nums)==1)
        get.ordinal(nums)
    else if (length(nums==2))
        paste0(get.ordinal(nums[1]), " and ", get.ordinal(nums[2]))
    else
    {
        len = length(nums)
        paste0(paste0(sapply(nums[-len], get.ordinal), collapse=", "),
               ", and ",
               get.ordinal(nums[len]))
    }
}


##-----------##
##-- OTHER --##
##-----------##



list.equals <- function(x, y)
{
    if (is.null(x) && is.null(y))
        T
    else
        !is.null(x) && !is.null(y) &&
        length(x)==length(y) &&
        sapply(1:length(x), function(i){
            all(x[[i]]==y[[i]])
        })
}

cast.to.integer <- function(value, arg.name.for.error=substitute(value))
{
    if (class(value) != 'integer' && class(value) != 'numeric')
        stop(paste0("'", arg.name.for.error, "' must be an integer"))
    else if (length(value) != 1)
        stop(paste0("'", arg.name.for.error, "' must a scalar, integer value"))
    else if ((value%%1)!=0)
        stop(paste0("'", arg.name.for.error, "' must be an integer"))
    else
        as.integer(value)
}

map.bounds <- function(var.names, bounds, n.var=length(var.names),
                       var.names.names)
{

}
