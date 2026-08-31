# Format/style code
# Must have air (https://posit-dev.github.io/air);
# if missing, the format step below fails silently
system2("air", c("format", "."))

# Generate documentation
roxygen2md::roxygen2md()
devtools::document()

# The vignette is a Quarto document, so the Quarto CLI must be findable.
# It ships with Positron and RStudio but is not on PATH; point QUARTO_PATH at
# it if quarto::quarto_version() errors.
if (inherits(try(quarto::quarto_version(), silent = TRUE), "try-error")) {
  Sys.setenv(
    QUARTO_PATH =
      "C:/Program Files/Positron/resources/app/quarto/bin/quarto.exe"
  )
}

# The vignette calls the package's own functions, and Quarto renders it in a
# fresh session that loads the INSTALLED package rather than the source tree.
# Install before checking, or the vignette fails on anything added since the
# last install.
devtools::install(upgrade = "never", build_vignettes = FALSE)

# README
devtools::build_readme()

# Build site
pkgdown::build_site()
browseURL("docs/index.html")

# Test & Check
devtools::test()
devtools::check()

rcmdcheck::rcmdcheck(
  args = c("--no-manual", "--as-cran"),
  build_args = "--resave-data=best",
  error_on = "warning"
)

# Code Coverage
covr::report(covr::package_coverage())

# Build Package
devtools::build()
