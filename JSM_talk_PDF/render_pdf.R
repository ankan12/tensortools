script_file <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]])))
  }

  source_file <- sys.frames()[[1]]$ofile
  if (!is.null(source_file)) {
    return(normalizePath(source_file))
  }

  normalizePath("render_pdf.R")
}

project_dir <- dirname(script_file())
qmd_file <- file.path(project_dir, "JSM_talk_2026.qmd")
html_file <- file.path(project_dir, "JSM_talk_2026.html")
pdf_file <- file.path(project_dir, "JSM_talk_2026.pdf")

quarto_bin <- Sys.which("quarto")
if (!nzchar(quarto_bin)) {
  stop("Quarto was not found on PATH.")
}

chrome_bin <- Sys.getenv(
  "JSM_CHROME",
  unset = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
)
if (!file.exists(chrome_bin)) {
  stop("Google Chrome was not found at: ", chrome_bin)
}

if (!requireNamespace("processx", quietly = TRUE)) {
  stop("The installed processx package is required.")
}

r_library <- Sys.getenv(
  "JSM_R_LIB",
  unset = paste0(
    "/Users/anupamakannan/Library/CloudStorage/",
    "OneDrive-BaylorUniversity/tensortools/.Rlib"
  )
)

build_dir <- tempfile("jsm-talk-", tmpdir = "/private/tmp")
#dir.create(build_dir, recursive = TRUE)
on.exit(unlink(build_dir, recursive = TRUE, force = TRUE), add = TRUE)

quarto_tmp <- file.path(build_dir, "quarto")
chrome_profile <- file.path(build_dir, "chrome-profile")
#dir.create(quarto_tmp)
#dir.create(chrome_profile)

quarto_env <- Sys.getenv()
quarto_env[["TMPDIR"]] <- quarto_tmp
quarto_env[["R_LIBS_USER"]] <- r_library

quarto_result <- processx::run(
  quarto_bin,
  c("render", qmd_file),
  env = quarto_env,
  echo = TRUE,
  error_on_status = TRUE
)

print_html <- file.path(build_dir, "JSM_talk_2026_print.html")
html_text <- readChar(html_file, nchars = file.info(html_file)$size)
html_text <- sub(
  'id="thank-you-slide" class="slide level2 smaller scrollable"',
  'id="thank-you-slide" class="slide level2"',
  html_text,
  fixed = TRUE
)
html_text <- gsub(
  ' <a href="https://doi\\.org/[^"]+">https://doi\\.org/[^<]+</a>\\.',
  '',
  html_text
)
html_text <- gsub("\\.\\s+\\.", ".", html_text)
html_text <- sub(
  "config: 'TeX-AMS_HTML-full'",
  "config: 'TeX-AMS_SVG-full'",
  html_text,
  fixed = TRUE
)

# Reveal's built-in separate-fragment PDF mode clones slides after MathJax has
# typeset them. Build the static fragment states before Reveal and MathJax start
# instead, so every exported state receives a normal math-typesetting pass.
html_text <- sub(
  "'pdfSeparateFragments': true",
  "'pdfSeparateFragments': false",
  html_text,
  fixed = TRUE
)
html_text <- sub(
  "slideNumber: ('c/t'|true),",
  "slideNumber: false,",
  html_text
)

fragment_export_script <- paste0(
  "<script id=\"jsm-static-pdf-fragments\">\n",
  "(() => {\n",
  "  if (!new URLSearchParams(window.location.search).has('print-pdf')) return;\n",
  "  const slides = document.querySelector('.reveal .slides');\n",
  "  if (!slides) return;\n",
  "  const logicalSlides = Array.from(slides.children).filter((slide) =>\n",
  "    slide.matches('section')\n",
  "  );\n",
  "  const countedSlides = logicalSlides.filter((slide) =>\n",
  "    slide.dataset.visibility !== 'uncounted'\n",
  "  );\n",
  "  const totalSlides = countedSlides.length;\n",
  "  const numberStyle = document.createElement('style');\n",
  "  numberStyle.textContent =\n",
  "    '.jsm-logical-slide-number {' +\n",
  "    'position: absolute; top: 8px; right: 8px; z-index: 1000;' +\n",
  "    'font: 14px/1 Arial, sans-serif; color: #444;' +\n",
  "    '}';\n",
  "  document.head.appendChild(numberStyle);\n",
  "  const addPdfPageNumbers = () => {\n",
  "    document.querySelectorAll('.pdf-page').forEach((page) => {\n",
  "      if (page.querySelector(':scope > .jsm-logical-slide-number')) return;\n",
  "      const slide = page.querySelector(':scope > section');\n",
  "      const label = slide && slide.dataset.jsmLogicalNumber;\n",
  "      if (!label) return;\n",
  "      const number = document.createElement('div');\n",
  "      number.className = 'jsm-logical-slide-number';\n",
  "      number.textContent = label;\n",
  "      page.appendChild(number);\n",
  "    });\n",
  "  };\n",
  "  const pageObserver = new MutationObserver(addPdfPageNumbers);\n",
  "  pageObserver.observe(slides, { childList: true });\n",
  "  let countedIndex = 0;\n",
  "  logicalSlides.forEach((slide) => {\n",
  "    const isUncounted = slide.dataset.visibility === 'uncounted';\n",
  "    if (!isUncounted) countedIndex += 1;\n",
  "    const logicalNumber = isUncounted\n",
  "      ? null\n",
  "      : countedIndex + '/' + totalSlides;\n",
  "    const fragmentCount = slide.querySelectorAll('.fragment').length;\n",
  "    if (fragmentCount === 0) {\n",
  "      if (logicalNumber) slide.dataset.jsmLogicalNumber = logicalNumber;\n",
  "      return;\n",
  "    }\n",
  "    const baseId = slide.id || 'slide';\n",
  "    for (let step = 0; step <= fragmentCount; step += 1) {\n",
  "      const clone = slide.cloneNode(true);\n",
  "      clone.id = \`${baseId}-pdf-step-${step}\`;\n",
  "      Array.from(clone.querySelectorAll('.fragment')).forEach((node, index) => {\n",
  "        if (index < step) {\n",
  "          node.classList.remove('fragment', 'visible', 'current-fragment');\n",
  "          node.removeAttribute('data-fragment-index');\n",
  "        } else {\n",
  "          node.remove();\n",
  "        }\n",
  "      });\n",
  "      if (logicalNumber) clone.dataset.jsmLogicalNumber = logicalNumber;\n",
  "      slides.insertBefore(clone, slide);\n",
  "    }\n",
  "    slide.remove();\n",
  "  });\n",
  "})();\n",
  "</script>\n"
)

reveal_init_marker <- paste0(
  "  <script>\n\n",
  "      // Full list of configuration options available at:"
)
if (!grepl(reveal_init_marker, html_text, fixed = TRUE)) {
  stop("Could not locate Reveal initialization in the rendered HTML.")
}
html_text <- sub(
  reveal_init_marker,
  paste0(fragment_export_script, reveal_init_marker),
  html_text,
  fixed = TRUE
)
writeChar(html_text, print_html, eos = NULL)

temporary_pdf <- file.path(build_dir, "JSM_talk_2026.pdf")
print_url <- utils::URLencode(
  paste0("file://", print_html, "?print-pdf"),
  reserved = FALSE
)

chrome_process <- processx::process$new(
  chrome_bin,
  c(
    "--headless=new",
    "--disable-gpu",
    "--disable-background-networking",
    "--disable-component-update",
    "--disable-sync",
    "--no-first-run",
    "--no-default-browser-check",
    "--allow-file-access-from-files",
    paste0("--user-data-dir=", chrome_profile),
    "--run-all-compositor-stages-before-draw",
    "--virtual-time-budget=25000",
    "--no-pdf-header-footer",
    paste0("--print-to-pdf=", temporary_pdf),
    print_url
  ),
  stdout = "|",
  stderr = "|",
  cleanup_tree = TRUE
)

last_size <- 0
stable_count <- 0
for (attempt in seq_len(60)) {
  Sys.sleep(1)
  if (file.exists(temporary_pdf)) {
    current_size <- file.info(temporary_pdf)$size
    if (!is.na(current_size) && current_size > 0 && current_size == last_size) {
      stable_count <- stable_count + 1
    } else {
      stable_count <- 0
      last_size <- current_size
    }
    if (stable_count >= 2) {
      break
    }
  }
}

if (chrome_process$is_alive()) {
  chrome_process$kill_tree()
}
chrome_process$wait(timeout = 2000)

if (!file.exists(temporary_pdf) || file.info(temporary_pdf)$size == 0) {
  stop("Chrome did not create the PDF within 60 seconds.")
}

if (!file.copy(temporary_pdf, pdf_file, overwrite = TRUE)) {
  stop("Could not copy the completed PDF to: ", pdf_file)
}

message("Created: ", pdf_file)
