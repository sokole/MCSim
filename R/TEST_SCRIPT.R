# -- old MCSim
devtools::install_github("sokole/MCSim")

fn.make.landscape3()
?fn.make.landscape
require(MCSim)

landscape<-fn.make.landscape3(x=c(1:10),
                   y=c(1:10),
                   management = 1,
                   area_m2 = rep(1,10),
                   Ef = rep(.5,10),
                   JM = 1000)

sim1<-fn.metaSIM3(landscape=landscape )
