# R_code
## purrr::map_dfr
tbl <-
  list.files(path = "./../02_stamp/output", recursive = T, pattern = ".*.xlsx", full.names = TRUE) %>%
  purrr::set_names() %>%
  purrr::map_dfr(readxl::read_excel, .id = "filepath")
