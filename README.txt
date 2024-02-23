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

# Rmarkdown 

换页：\newpage
引用图片：![](%22../Metaphlan4_TOFU/version3_Metaphlan4_SRP320766_165Subjects_TOFU/output/differential/consensus/differential_consensus_Species.png%22){width="100%"
height="95%"}
字体大小及行距：\fontsize{12}{12}\selectfont

---
Title: Report
author: "YNYK"
date: "`r Sys.Date()`"
bibliography: myref.bib
output:
  pdf_document:
    toc: yes
    toc_depth: '4'
    keep_tex: yes
    latex_engine: xelatex
    number_sections: true
  html_document:
    toc_depth: 4
    toc: yes
    thumbnails: yes
    number_sections: yes
    toc_float: yes
    gallery: yes
header-includes:
- \usepackage{fancyhdr}
- \pagestyle{fancy}
- \fancyhead[CO]{Innoihealth Biopharmaceutical Co.LTD}
editor_options: 
  markdown: 
    wrap: 72
---

# 读取excel文件中的所有sheet名称  
sheets <- readxl::excel_sheets("../../input/20240222 差异基因通路富集cluster.xlsx")  

# 产生一个空表进行保存循环的结果
results <- data.frame(feature=character(), log2fc=numeric(), p_value=numeric(), p_adjust=numeric())  # Create an empty data frame to store the results
results <- rbind(results, data.frame(feature=feature, log2fc=log2fc, p_value=test_result$p.value))  # Store the results in the data frame

# 用pull拉取变量的数据
 group1 <- lcms %>% filter(class == "ineffective") %>% pull(!!sym(feature))  # Extract data for group1

