#如果是同一个数据框里面的两个变量，可以如下合并，直接table两列
library("MASS")
car.data <- data.frame(Cars93$AirBags, Cars93$Type)
car.data = table(Cars93$AirBags, Cars93$Type) 
print(car.data)
print(chisq.test(car.data))


#如果分别从两个数据框提取数据，可以如下整理数据，分别table两个变量，然后rbind合并
a <- table(group_high$gleason_score)
b <- table(group_low$gleason_score)
c <- rbind(a,b)
c
par(mfrow=c(2,1))
barplot(a)
barplot(b)
chisq.test(c)#if got any warning, better to use fishertest
fisher.test(c)
p <- fisher.test(c)$p.value
format(p,scientific = F)# use non-scientific format, easier to read
