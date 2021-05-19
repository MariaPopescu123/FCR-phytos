phytos1 <- phytos1 %>%
  mutate(MonthYear = as.numeric(factor(format(Date, "%y-%m"))))
my.combn.may <- combn(c(1,6,11,16),2)
my.combn.june <- combn(c(2,7,12,17),2)
my.combn.july <- combn(c(3,8,13,18),2)
my.combn.aug <- combn(c(4,9,14,19),2)
my.combn.sept <- combn(c(5,10,15,20),2)

m.y <- c(1:20)
m.y.c <- c("May 2016","June 2016","July 2016","August 2016","September 2016",
           "May 2017","June 2017","July 2017","August 2017","September 2017",
           "May 2018","June 2018","July 2018","August 2018","September 2018",
           "May 2019","June 2019","July 2019","August 2019","September 2019")
m.y.key <- data.frame(m.y,m.y.c)
colnames(m.y.key) <- c("MonthYear_A","MonthYear_Char_A")

m.y.table <- left_join(final.m.y,m.y.key,by = "MonthYear_A")

m.y.key <- data.frame(m.y,m.y.c)
colnames(m.y.key) <- c("MonthYear_B","MonthYear_Char_B")

m.y.table2 <- left_join(m.y.table,m.y.key,by = "MonthYear_B") %>%
  mutate(comparison = paste(MonthYear_Char_A,MonthYear_Char_B, sep = "-"),
         R = round(R,2),
         p = round(p,3)) %>%
  select(comparison, R, p)
write.csv(m.y.table2,"./2_Data_analysis/Table_S7.csv",row.names = FALSE)
