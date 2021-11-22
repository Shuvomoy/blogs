using Dates

# Includes the Package, bringing the lx_functions into scope
using FranklinTheorems

# Includes the custom markdown files, bringing the `\newcommand` and `\newenvironment` definitions into scope.
Franklin.include_external_config(FranklinTheorems.config_path()) 

function hfun_bar(vname)
  val = Meta.parse(vname[1])
  return round(sqrt(val), digits=2)
end

function hfun_m1fill(vname)
  var = vname[1]
  return pagevar("index", var)
end

function lx_baz(com, _)
  # keep this first line
  brace_content = Franklin.content(com.braces[1]) # input string
  # do whatever you want here
  return uppercase(brace_content)
end

function hfun_timestamp_now()
    return string(Dates.now()) * "+00:00"
end


"""
    {{blogposts}}

Plug in the list of blog posts contained in the `/posts/` folder.
"""
function hfun_blogposts()
    today = Dates.today()
    curyear = year(today)
    curmonth = month(today)
    curday = day(today)

    list = readdir("posts")
    filter!(f -> endswith(f, ".md"), list)
    sorter(p) = begin
	    @show "calling sorter"
        ps  = splitext(p)[1]
        url = "/posts/$ps/"
        surl = strip(url, '/')
		println("surl = ", surl)
        pubdate = pagevar(surl, :published)
		@show pubdate
        if isnothing(pubdate)
		    @show "💀 something is not right"
		    @show "nopubdate"
            return Date(Dates.unix2datetime(stat(surl * ".md").ctime))
        end
		@show Date(pubdate, DateFormat("U d, Y"))
        return Date(pubdate, DateFormat("U d, Y"))
    end
    sort!(list, by=sorter, rev=true)

    io = IOBuffer()
    write(io, """<ul class="blog-posts">""")
    for (i, post) in enumerate(list)
        if post == "index.md"
            continue
        end
        ps  = splitext(post)[1]
		@show "ps = $ps"
        write(io, "<li><span><i>")
        url = "/posts/$ps/"
        surl = strip(url, '/')
        title = pagevar(surl, :title)
        pubdate = pagevar(surl, :published)
        if isnothing(pubdate)
            date    = "$curyear-$curmonth-$curday"
        else
            date    = Date(pubdate, DateFormat("U d, Y"))
        end
        write(io, """($date) </i></span><a href="$url">$title</a>""")
    end
    write(io, "</ul>")
    return String(take!(io))
end

"""
    {{custom_taglist}}

Plug in the list of blog posts with the given tag
"""
function hfun_custom_taglist()::String
    tag = locvar(:fd_tag)
    rpaths = globvar("fd_tag_pages")[tag]
    sorter(p) = begin
        pubdate = pagevar(p, :published)
        if isnothing(pubdate)
            return Date(Dates.unix2datetime(stat(p * ".md").ctime))
        end
        return Date(pubdate, DateFormat("U d, Y"))
    end
    sort!(rpaths, by=sorter, rev=true)

    io = IOBuffer()
    write(io, """<ul class="blog-posts">""")
    # go over all paths
    for rpath in rpaths
        write(io, "<li><span><i>")
        url = get_url(rpath)
        title = pagevar(rpath, :title)
        pubdate = pagevar(rpath, :published)
        if isnothing(pubdate)
            date    = "$curyear-$curmonth-$curday"
        else
            date    = Date(pubdate, DateFormat("U d, Y"))
        end
        # write some appropriate HTML
        write(io, """$date</i></span><a href="$url">$title</a>""")
    end
    write(io, "</ul>")
    return String(take!(io))
end


