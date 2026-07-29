get_top_modindices <- get_top_mi <- function(sem_fit, min_mi = 10) {
  modindices(sem_fit) %>%
    filter(mi > min_mi, op == "~") %>% # Focus on regressions
    arrange(desc(mi)) %>%
    head(5)
}