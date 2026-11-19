using Documenter
using ResilienceCore

# Source and edit links need a git remote; outside a git checkout, build without them.
ingit = success(pipeline(`git -C $(@__DIR__) rev-parse`; stdout = devnull, stderr = devnull))

makedocs(;
         modules = [ResilienceCore],
         authors = "Simon Frost",
         sitename = "ResilienceCore.jl",
         checkdocs = :exports,
         remotes = ingit ? Dict() : nothing,
         format = Documenter.HTML(;
                                  edit_link = ingit ? "main" : nothing,
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://ecorecipes.github.io/ResilienceCore.jl"),
         pages = ["Home" => "index.md", "Metric catalogue" => "catalogue.md", "Gallery" => "gallery.md", "API Reference" => "api.md"])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/ResilienceCore.jl.git", push_preview = true)
end
