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
out.name=paste0( dirname(filename), gsub(".csv",".html",basename(filename))  )
save_kable(table_q, out.name )

table_q
print(glue::glue("You can get the result of t.test in " ,out.name) )

# 获取变量名称并转换为字符串
variable_name <- deparse(substitute(example_immune))

# Deseq2 
// 结果默认返回deseq2中最后一个分组顺序，如Groups是最后一个
ddsFull <- phyloseq::phyloseq_to_deseq2(phy_tmp, design = ~ Sex +Age + BMI + 吸烟 + 荤素饮食习惯 + 食用酸奶制品的频率 +食用腌制类食物的频率 + Groups ) 
ddsReduced <- phyloseq::phyloseq_to_deseq2(phy_tmp, design = ~ Sex +Age  + Groups ) 

ddsFull_dds <- DESeq2::DESeq(ddsFull)
ddsReduced_dds<-DESeq2::DESeq(ddsReduced)
