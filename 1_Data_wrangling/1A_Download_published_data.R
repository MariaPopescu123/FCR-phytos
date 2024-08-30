#Title: 1A_Download_published_data
#Author: Mary Lofton
#Date: 14JUL21

#Purpose: Download all published data from Environmental Data Initiative (EDI) so
#it can be pulled in for data tidying and analysis

#download chemistry data from EDI 
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.199.6&entityid=2b3dc84ae6b12d10bd5485f1c300af13"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/chemistry.csv", method='libcurl')

#download CTD data from EDI 
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.200.10&entityid=2461524a7da8f1906bfc3806d594f94c"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/CTD.csv", method='libcurl')

#download YSI data from EDI
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.198.7&entityid=25b5e8b7f4291614d5c6d959a08148d8"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/YSI.csv", method='libcurl')

#download SCC thermistor string data from EDI
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.271.4&entityid=151ddb643fbfec06e058f27af837e5ae"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/SCC.csv", method='libcurl')

options(timeout = 300)  # Increase timeout to 300 seconds

#download met station data from EDI
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.389.5&entityid=3d1866fecfb8e17dc902c76436239431"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/met.csv", method='libcurl')

#download Secchi data from EDI 
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.198.8&entityid=336d0a27c4ae396a75f4c07c01652985"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/Secchi.csv", method='libcurl')

#download FluoroProbe data from EDI
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.272.4&entityid=e7e3e6e513985a602d9a5f22687d4efc"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/FP.csv", method='libcurl')

#download inflow stream discharge data from EDI
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.202.7&entityid=f5fa5de4b49bae8373f6e7c1773b026e"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/discharge.csv", method='libcurl')

#download phytoplankton biovolume data from EDI
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.875.1&entityid=00d64d23fc2b75d973f69cc09bb5d083"
destination <- "./0_Data_files"

download.file(data,destfile = "./0_Data_files/phytoplankton.csv", method='libcurl')
