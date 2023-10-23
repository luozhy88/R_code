# R_code
## purrr::map_dfr
tbl <-
  list.files(path = "./../02_stamp/output", recursive = T, pattern = ".*.xlsx", full.names = TRUE) %>%
  purrr::set_names() %>%
  purrr::map_dfr(readxl::read_excel, .id = "filepath")

## 三线表
 library(kableExtra)
data <- data.frame(
  Name = c("Alice", "Bob", "Charlie"),
  Age = c(25, 30, 28),
  Score = c(92, 85, 89)
)

table_q=kbl(data %>% head(), align = "c") %>%
  kable_classic(full_width = F) %>%
  kableExtra::footnote(general = "Here is a general comments of the table. ",
           # number = c("Footnote 1; ", "Footnote 2; "),
           # alphabet = c("Footnote A; ", "Footnote B; "),
           symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2")
           ) #%>% kableExtra::as_image()
# out.name=paste0( dirname(filename), gsub(".csv",".html",basename(filename))  )
# save_kable(table_q, out.name )

table_q
# print(glue::glue("You can get the result of t.test in " ,out.name) )
