## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE,include=FALSE------------------------------------------

## ---- echo=FALSE,results="asis"-----------------------------------------------
files <- list.files("web", "\\.Rmd$")
for (file in file.path("web", files)) {
  lines <- readLines(file, n = 10)
  title_idx <- grep("^title: ", lines)
  if (length(title_idx) > 0) {
    title <- sub("^title: ", "", lines[title_idx[1]])
    title <- sub('^"', "", title)
    title <- sub('"$', "", title)
    cat("* [",
      title,
      "](https://lindesaysh.github.io/MRSea/articles/",
      sub("\\.Rmd", ".html", file),
      ")\n",
      sep = ""
    )
  }
}

## ---- echo=FALSE,results="asis"-----------------------------------------------
files <- list.files(".", "\\.Rmd$")
for (file in files) {
  lines <- readLines(file, n = 10)
  title_idx <- grep("^title: ", lines)
  if (length(title_idx) > 0) {
    title <- sub("^title: ", "", lines[title_idx[1]])
    title <- sub('^"', "", title)
    title <- sub('"$', "", title)
    cat("* [",
      title,
      "](https://lindesaysh.github.io/MRSea/articles/",
      sub("\\.Rmd", ".html", file),
      ")\n",
      sep = ""
    )
  }
}

