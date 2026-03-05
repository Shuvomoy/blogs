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

function _escape_html(s::AbstractString)
    s = replace(s, "&" => "&amp;")
    s = replace(s, "<" => "&lt;")
    s = replace(s, ">" => "&gt;")
    s = replace(s, "\"" => "&quot;")
    return replace(s, "'" => "&#39;")
end

function _note_link(filename::AbstractString)
    # Franklin emits HTML files as page routes (.../index.html), so link to
    # the stem path for .html notes and direct file path for other formats.
    stem, ext = splitext(filename)
    if lowercase(ext) == ".html"
        return "/assets/Notes/" * replace(stem, " " => "%20") * "/"
    end
    return "/assets/Notes/" * replace(filename, " " => "%20")
end

function hfun_handwritten_notes()
    notes_dir = joinpath(@__DIR__, "assets", "Notes")
    if !isdir(notes_dir)
        return "<p>No handwritten notes are available yet.</p>"
    end

    pattern = r"^(\d{4})_(\d{2})_(\d{2})_(.+)\.(html|pdf)$"i
    entries = Dict{String, Dict{Symbol, Any}}()

    for filename in readdir(notes_dir)
        fullpath = joinpath(notes_dir, filename)
        isfile(fullpath) || continue

        m = match(pattern, filename)
        m === nothing && continue

        y, mo, d, raw_title, ext = m.captures
        note_date = try
            Date(parse(Int, y), parse(Int, mo), parse(Int, d))
        catch
            continue
        end

        key = "$(y)_$(mo)_$(d)_$(raw_title)"
        if !haskey(entries, key)
            entries[key] = Dict(
                :date => note_date,
                :title => replace(raw_title, "_" => " "),
                :html => nothing,
                :pdf => nothing,
            )
        end

        ext = lowercase(ext)
        if ext == "html"
            entries[key][:html] = filename
        elseif ext == "pdf"
            entries[key][:pdf] = filename
        end
    end

    if isempty(entries)
        return "<p>No handwritten notes are available yet. Add files to <code>assets/Notes/</code> using <code>YYYY_MM_DD_Title.html</code> or <code>YYYY_MM_DD_Title.pdf</code>.</p>"
    end

    sorted_entries = sort(
        collect(values(entries));
        by = entry -> (entry[:date], entry[:title]),
        rev = true,
    )

    io = IOBuffer()
    write(io, "<ul class=\"handwritten-notes-list\">")
    for entry in sorted_entries
        date_str = Dates.format(entry[:date], dateformat"yyyy-mm-dd")
        title_str = _escape_html(string(entry[:title]))

        html_href = isnothing(entry[:html]) ? nothing : _escape_html(_note_link(string(entry[:html])))
        pdf_href = isnothing(entry[:pdf]) ? nothing : _escape_html(_note_link(string(entry[:pdf])))

        primary_href = isnothing(html_href) ? pdf_href : html_href
        primary_href === nothing && continue

        write(io, "<li>($date_str) <a href=\"$primary_href\">$title_str</a>")
        if !isnothing(html_href) && !isnothing(pdf_href)
            write(io, " (<a href=\"$html_href\">HTML</a> | <a href=\"$pdf_href\">PDF</a>)")
        end
        write(io, "</li>")
    end
    write(io, "</ul>")

    return String(take!(io))
end
