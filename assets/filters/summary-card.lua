-- Engineering summary card for project pages.
--
-- Reads a `summary:` block from page metadata (each project keeps it in its
-- own folder's _metadata.yml so the page source, and therefore its freeze
-- hash, stays untouched) and prepends a compact card to the page body:
--
--   summary:
--     problem: ...
--     data: ...
--     methods: ...
--     result: ...
--     stack: [Python, GeoPandas]
--     code: {text: ..., href: ...}   # optional; defaults to first code-links entry
--     live: "Live · updated nightly"  # optional badge in the card title
--
-- Values are markdown. Any value starting with "TODO" is flagged visually.
--
-- A `reproduce:` block adds a "Reproduce" row of badges (Codespaces, Binder,
-- Environment, Data, and DOI once `doi:` is set) and a "Data availability"
-- appendix section:
--
--   reproduce:
--     env: environment.yml            # relative to the page's folder
--     data:
--       - {name: ..., url: ..., licence: ..., fetch: ..., path: ...}
--     note: ...                       # optional

local stringify = pandoc.utils.stringify

local ROWS = {
  { key = "problem", label = "Problem" },
  { key = "data",    label = "Data" },
  { key = "methods", label = "Methods" },
  { key = "result",  label = "Key result" },
}

local function to_blocks(value)
  local t = pandoc.utils.type(value)
  if t == "Blocks" then
    return value
  elseif t == "Inlines" then
    return pandoc.Blocks({ pandoc.Plain(value) })
  end
  return pandoc.Blocks({ pandoc.Plain({ pandoc.Str(stringify(value)) }) })
end

local function is_todo(value)
  return stringify(value):match("^%s*TODO") ~= nil
end

local function row(label, blocks, extra_class)
  local classes = { "es-row" }
  if extra_class then table.insert(classes, extra_class) end
  return pandoc.Div({
    pandoc.Div({ pandoc.Plain({ pandoc.Str(label) }) }, pandoc.Attr("", { "es-label" })),
    pandoc.Div(blocks, pandoc.Attr("", { "es-value" })),
  }, pandoc.Attr("", classes))
end

local function code_link(meta, summary)
  local link = summary.code
  if link == nil and meta["code-links"] ~= nil and #meta["code-links"] > 0 then
    link = meta["code-links"][1]
  end
  if link == nil or link.href == nil then
    return nil
  end
  local text = link.text and stringify(link.text) or "Source"
  return pandoc.Blocks({ pandoc.Plain({
    pandoc.Link({ pandoc.Str(text) }, stringify(link.href)),
  }) })
end

local REPO = "tilak-pentapati/tilak-pentapati.github.io"
local BRANCH = "main"

-- Folder of the page being rendered, relative to the project root.
local function page_dir()
  local input = quarto.doc.input_file
  local root = quarto.project.directory
  if input == nil or root == nil then return nil end
  local rel = input:sub(#root + 2)
  return rel:match("^(.*)/[^/]+$") or "."
end

local function badge(kind, icon, text, href, title)
  return pandoc.RawInline("html", string.format(
    '<a class="repro-badge repro-%s" href="%s" title="%s"%s>' ..
    '<i class="bi bi-%s" aria-hidden="true"></i><span>%s</span></a>',
    kind, href, title, href:sub(1, 1) == "#" and "" or ' target="_blank" rel="noopener"',
    icon, text))
end

local function repro_badges(meta, repro, dir)
  local env = stringify(repro.env or "environment.yml")
  local inlines = pandoc.Inlines({
    badge("codespaces", "github", "Codespaces",
      "https://codespaces.new/" .. REPO .. "?quickstart=1",
      "Open the repository in GitHub Codespaces, then run scripts/repro/setup-env.sh " .. dir),
    badge("binder", "rocket-takeoff", "Binder",
      "https://mybinder.org/v2/gh/" .. REPO .. "/" .. BRANCH .. "?labpath=" .. dir .. "/index.qmd",
      "Launch this page's source on mybinder.org"),
    badge("env", "box-seam", "Environment",
      "https://github.com/" .. REPO .. "/blob/" .. BRANCH .. "/" .. dir .. "/" .. env,
      "Pinned package versions: " .. env),
    badge("data", "database", "Data", "#data-availability",
      "Where the data comes from and how to fetch it"),
  })
  if meta.doi ~= nil then
    local doi = stringify(meta.doi)
    inlines:insert(badge("doi", "bookmark-check", "DOI " .. doi, "https://doi.org/" .. doi,
      "Archived release on Zenodo"))
  end
  return pandoc.Blocks({ pandoc.Div({ pandoc.Plain(inlines) }, pandoc.Attr("", { "repro-badges" })) })
end

local function field(label, value)
  local content = pandoc.utils.blocks_to_inlines(to_blocks(value))
  local out = pandoc.Inlines({ pandoc.Strong({ pandoc.Str(label .. ":") }), pandoc.Space() })
  out:extend(content)
  return pandoc.Div({ pandoc.Plain(out) })
end

-- "Data availability" appendix: one entry per source, then the rebuild recipe.
local function data_section(repro, dir)
  local blocks = pandoc.Blocks({
    pandoc.Header(2, { pandoc.Str("Data availability") },
      pandoc.Attr("data-availability", { "appendix", "unnumbered" })),
  })
  local items = {}
  for _, src in ipairs(repro.data or {}) do
    local name = pandoc.utils.blocks_to_inlines(to_blocks(src.name))
    local title = src.url and pandoc.Link(name, stringify(src.url)) or pandoc.Span(name)
    local entry = { pandoc.Div({ pandoc.Plain({ pandoc.Strong({ title }) }) }) }
    if src.licence then table.insert(entry, field("Licence", src.licence)) end
    if src.fetch then table.insert(entry, field("How to get it", src.fetch)) end
    if src.path then table.insert(entry, field("Expected path", src.path)) end
    table.insert(items, pandoc.Div(entry, pandoc.Attr("", { "repro-source" })))
  end
  blocks:extend(items)
  if repro.note then blocks:extend(to_blocks(repro.note)) end
  blocks:insert(pandoc.Para({
    pandoc.Str("Raw data is not part of the repository (see"), pandoc.Space(),
    pandoc.Link({ pandoc.Code("data/README.md") },
      "https://github.com/" .. REPO .. "/blob/" .. BRANCH .. "/data/README.md"),
    pandoc.Str(")."), pandoc.Space(), pandoc.Str("To rebuild the environment and re-run this page:"),
    pandoc.Space(), pandoc.Code("scripts/repro/setup-env.sh " .. dir), pandoc.Str("."),
  }))
  return blocks
end

function Pandoc(doc)
  if not quarto.doc.is_format("html") then return nil end
  local summary = doc.meta.summary
  if summary == nil or pandoc.utils.type(summary) ~= "table" then return nil end

  local rows = {}
  for _, r in ipairs(ROWS) do
    local v = summary[r.key]
    if v ~= nil then
      local blocks = to_blocks(v)
      if is_todo(v) then
        blocks = pandoc.Blocks({ pandoc.Div(blocks, pandoc.Attr("", { "es-todo" })) })
      end
      table.insert(rows, row(r.label, blocks, r.key == "result" and "es-result" or nil))
    end
  end

  if summary.stack ~= nil then
    local chips = {}
    for _, item in ipairs(summary.stack) do
      table.insert(chips, pandoc.Span({ pandoc.Str(stringify(item)) }, pandoc.Attr("", { "es-chip" })))
      table.insert(chips, pandoc.Space())
    end
    table.insert(rows, row("Stack", pandoc.Blocks({ pandoc.Plain(chips) })))
  end

  local code = code_link(doc.meta, summary)
  table.insert(rows, row("Code", code or to_blocks(pandoc.Inlines("TODO: add code link"))))

  local repro = doc.meta.reproduce
  local dir = repro and page_dir()
  if repro ~= nil and dir ~= nil then
    table.insert(rows, row("Reproduce", repro_badges(doc.meta, repro, dir)))
    doc.blocks:extend(data_section(repro, dir))
  end

  local title = pandoc.Div({ pandoc.Plain({
    pandoc.RawInline("html", '<i class="bi bi-clipboard-data" aria-hidden="true"></i>'),
    pandoc.Str("At a glance"),
  }) }, pandoc.Attr("", { "es-title" }))
  if summary.live ~= nil then
    title.content[1].content:insert(pandoc.Span({ pandoc.Str(stringify(summary.live)) },
      pandoc.Attr("", { "es-live" })))
  end

  table.insert(rows, 1, title)
  local card = pandoc.Div(rows, pandoc.Attr("project-summary", { "eng-summary" },
    { role = "note", ["aria-label"] = "Project summary" }))

  doc.blocks:insert(1, card)
  return doc
end
