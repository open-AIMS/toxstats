# The disclaimer is attached to the session rather than to each result. A
# reader of a printed result is nearly always the person who ran it, and a
# message on every print would be ignored within a day. The two artefacts that
# do travel away from the session, the Shiny report and the generated script,
# carry the text in full themselves; see the shinytoxstats package.

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "toxstats was written with generative AI and is for testing and\n",
    "validation only. Do not use it to derive toxicity estimates for\n",
    "regulatory, reporting or any other official purpose. See ?toxstats."
  )
}
