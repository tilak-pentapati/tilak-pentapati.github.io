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
--
-- Values are markdown. Any value starting with "TODO" is flagged visually.

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

  local title = pandoc.Div({ pandoc.Plain({
    pandoc.RawInline("html", '<i class="bi bi-clipboard-data" aria-hidden="true"></i>'),
    pandoc.Str("At a glance"),
  }) }, pandoc.Attr("", { "es-title" }))

  table.insert(rows, 1, title)
  local card = pandoc.Div(rows, pandoc.Attr("project-summary", { "eng-summary" },
    { role = "note", ["aria-label"] = "Project summary" }))

  doc.blocks:insert(1, card)
  return doc
end
