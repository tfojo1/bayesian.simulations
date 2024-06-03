
create.random.viewer.output.stream <- function(print.chain=T,
                                               use.viewer=T)
{
    html.file = tempfile(pattern='updates', fileext='.html')

    create.output.stream(html.file,
                         use.html = T,
                         use.viewer = use.viewer,
                         print.chain=print.chain)
}

create.output.stream <- function(dst,
                                 use.html,
                                 use.viewer,
                                 print.chain)
{
    if (use.viewer)
        viewer = getOption('viewer')

    if (use.html)
        newline = '<BR>'
    else
        newline = '\n'

    if (use.html)
        tab = '&nbsp;&nbsp;&nbsp;&nbsp;'
    else
        tab = '\t'

    function(..., chain=NA, sep='')
    {
        to.print = paste(..., sep=sep)

        split = unlist(strsplit(to.print,
                        split="\\n|<BR>"))
        ends.with.newline = grepl("\\n|<BR>$", to.print)

        if (print.chain && !is.na(chain))
            split[split != ''] = paste0('CHAIN ', chain, ":", tab, split[split != ''])

        to.print = paste0(split, collapse=newline)

#        if (dst != '')
#            to.print = paste0(newline, to.print)

        if (ends.with.newline)
            to.print = paste0(to.print, newline)

        if (dst=='')
            cat(to.print)
        else
            write(to.print, file=dst, append=T)

        if (use.viewer)
            viewer(dst)
    }
}

wrap.output.stream.for.chain <- function(stream, chain)
{
    function(..., sep='')
    {
        stream(..., chain=chain, sep=sep)
    }
}

get.default.output.stream <- function(user.specified,
                                      n.chains,
                                      initial.message,
                                      print.output.location.if.random.without.viewer=T)
{


    if (is.null(user.specified) || is.na(user.specified) || user.specified=='')
    {
        if (n.chains==1)
        {
            dst=''
            is.random.dst=F
            is.html=F
        }
        else
        {
            dst = tempfile(pattern='updates', fileext='.html')
            is.random.dst=T
            is.html=T
        }
    }
    else
    {
        dst = user.specified
        is.random.dst=F
        is.html = grepl('.html$', user.specified)
    }

    stream=create.output.stream(dst=dst, use.html=is.html, use.viewer=F, print.chain=n.chains>1)
    if (!is.null(initial.message))
        stream(initial.message)

    is.rstudio = Sys.getenv("RSTUDIO")==1 || .Platform$GUI == "RStudio"
    if (is.rstudio && is.html)
    {
        viewer = getOption('viewer')
        viewer(dst)
        launched.viewer=T
    }
    else
        launched.viewer=F

    if (is.random.dst && !launched.viewer)
        cat("Updates will be printed to '", dst, "'", sep='')

    stream
}
