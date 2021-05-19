final.2016$Year <- "2016"
final.2017$Year <- "2017"
final.2018$Year <- "2018"
final.2019$Year <- "2019"
final <- bind_rows(final.2016,final.2017,final.2018,final.2019) %>%
  mutate(Month_A = ifelse(Month_A == 5,"May",
                          ifelse(Month_A == 6, "June",
                                 ifelse(Month_A == 7,"July",
                                        ifelse(Month_A == 8,"August",
                                               ifelse(Month_A == 9,"September",""))))),
         Month_B = ifelse(Month_B == 5,"May",
                          ifelse(Month_B == 6, "June",
                                 ifelse(Month_B == 7,"July",
                                        ifelse(Month_B == 8,"August",
                                               ifelse(Month_B == 9,"September",""))))),
         comparison = paste(Month_A,Month_B, sep = "-"),
         R = round(R,2),
         p = round(p,3)) %>%
  select(Year, comparison, R, p)
write.csv(final,"./2_Data_analysis/intra-annual_pairwise_month_ANOSIM.csv",row.names = FALSE)
