scripts <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (s in scripts) {
  cat("\n#####", s, "\n")
  source(s, echo = FALSE, print.eval = TRUE)
}
