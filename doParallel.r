


library(dplyr)
library(optparse)

# 描述参数的解析方式 多线程
option_list <- list(
  make_option(c("-f", "--first"), type = "character", default = FALSE,
              action = "store", help = "输入正向model，如/database/Models/Ms2query/KEGG/v20240625/positive_mode/exact_matches_test_sets_splits/test_split_1/models ")
)
# 解析参数
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))
First_arg<-opt$first
SCORE<<-First_arg
print(opt)
print(glue::glue("args f:",SCORE) )
options(stringsAsFactors = F)

#################################################load database#################################################
# SCORE<<-0.95
combine_CID.SID_to_synonyms <- readRDS("output/04_combine_CID.SID_to_synonyms/X04_combine_CID.SID_to_synonyms.rds")
# combine_CID.SID_to_synonyms=combine_CID.SID_to_synonyms[1:1000,1:6]
# head(combine_CID.SID_to_synonyms)[1:2,1:6]

################################################predict data################################################


library(readxl)
library(foreach)
library(doParallel)
library(stringdist)
library(dplyr)
library(purrr)
all_measurement <- read.csv("Test1/Metabolite_data.csv")
all_measurement=all_measurement #%>% head(8)
# all_measurement$Metabolite
# compund.name.mapping=purrr::map_chr( compund.name , mapping.kegg)
# compund.name.mapping=purrr::map_chr( all_measurement$Metabolite , mapping.kegg)


mapping.kegg <- function(compound_name, SCORE, combine_CID.SID_to_synonyms) {
  # compound_name="ferulic acid"
  print(compound_name)
  print(glue::glue("SCORE:{SCORE}"))
  
  # 使用foreach并行处理每一行
  combine_CID.SID_to_synonyms$Score <- foreach(r = 1:nrow(combine_CID.SID_to_synonyms), .combine = c, .packages = c("stringdist", "dplyr")) %dopar% {
    row <- combine_CID.SID_to_synonyms[r, "Syn"]
    row_vector <- unlist(strsplit(row, split = ";"))
    row_vector_max <- max(stringsim(a = compound_name, b = row_vector, method = "cosine"))
    row_vector_max
  }
  # SCORE=0.98
  better_score <- combine_CID.SID_to_synonyms  %>%  dplyr::filter(Score > SCORE)#%>% dplyr::arrange(de"Score")
  # 按Score进行降序排列
  better_score <- better_score[order(-better_score$Score), ] %>% head(1)
  print(glue::glue("nrow:{nrow(better_score)}"))
  better_score_kegg <- paste(better_score$cid, collapse = ";")
  print(better_score_kegg)
  return(better_score_kegg)
}

# qq=mapping.kegg(all_measurement$Metabolite[1], SCORE, combine_CID.SID_to_synonyms)
# 使用foreach并行处理多个输入
# re <- foreach(metabolite = all_measurement$Metabolite, .combine = c, .packages = c("dplyr", "stringdist", "foreach")) %dopar% {
#   mapping.kegg(metabolite, SCORE=SCORE, combine_CID.SID_to_synonyms=combine_CID.SID_to_synonyms)
# }

# 使用tryCatch执行并行任务
re <- tryCatch({
              # 检测核心数并创建集群
              numCores <- detectCores() / 4
              cl <- makeCluster(numCores)
              registerDoParallel(cl)
              print(glue::glue("numCores:{numCores}"))
              
              foreach(metabolite = all_measurement$Metabolite, .combine = c, .packages = c("dplyr", "stringdist", "foreach")) %dopar% {
                    mapping.kegg(metabolite, SCORE=SCORE, combine_CID.SID_to_synonyms=combine_CID.SID_to_synonyms)
                      }
  
          }, error = function(e) {
            message("捕获到错误: ", e$message)
            NULL  # 返回NULL或其他默认值
          }, finally = {
            # 确保集群在任务完成或出现错误时被释放
            stopCluster(cl)
            message("集群已停止")
          })



# print(re)
all_measurement$Score=re

dir.create("output/05_predict", showWarnings = FALSE)
saveRDS(all_measurement, glue::glue("output/05_predict/X05_predict{SCORE}.rds"))
write.csv(all_measurement, glue::glue("output/05_predict/X05_predict{SCORE}.csv"), row.names = FALSE)


###################################################sumarry#####################
anno_count=sum(all_measurement$Score!="")
anno_percentage=anno_count/nrow(all_measurement)

print(glue::glue("anno_count:{anno_count}"))
print(glue::glue("anno_percentage:{anno_percentage}"))


