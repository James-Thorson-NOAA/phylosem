
data_dir =  R'(C:\Users\James.Thorson\Desktop\Git\phylosem\data-raw)'

# Load and format
Mlifehist_ver1_0 = read.csv( file.path(data_dir,"Mlifehist_ver1.0.csv") )
Mlifehist_ver1_0 = Mlifehist_ver1_0[,-match(c("GrowthRef","Mref","tmaxRef"),names(Mlifehist_ver1_0))]

# Export
setwd( R'(C:\Users\James.Thorson\Desktop\Git\phylosem)' )
usethis::use_data( Mlifehist_ver1_0, overwrite=TRUE )
