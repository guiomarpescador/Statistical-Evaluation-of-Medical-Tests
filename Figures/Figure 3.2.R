#FIGURE 2
rm(list = ls())
require(ggplot2)
library(gridExtra)

# Densities
seqh1 <- seq(-4,16,by=0.1)
seqd1 <- seq(-4,16,by=0.1)

densd1 <- dnorm(seq(-4,16,by=0.1))
densh1 <- dnorm(seq(-4,16,by=0.1), mean=0.2)

dfden1 <- data.frame(den = c(c(densh1), c(densd1)), 
                     seq    = c(seqh1, seqd1),
                     Status = c(rep("Nondiseased", length(seqh1)),
                                rep("Diseased", length(seqd1))))
densd2 <- dnorm(seq(-4,16,by=0.1), mean=9)
densh2 <- dnorm(seq(-4,16,by=0.1), mean=12)

dfden2 <- data.frame(den = c(c(densh2), c(densd2)), 
                     seq    = c(seqh1, seqd1),
                     Status = c(rep("Nondiseased", length(seqh1)),
                                rep("Diseased", length(seqd1))))


d1 <- ggplot(dfden1, aes(x = seq, y = den, linetype = Status)) +
  labs(title = "", x = "Test outcome", y = "Density") +
  coord_cartesian(xlim = c(range(dfden1$seq)[1], range(dfden1$seq)[2]), ylim =c(0, 0.45)) +
  geom_area(data = subset(dfden1, Status == 'Nondiseased'), fill = "blue", alpha = 0.3, show.legend = FALSE) +
  geom_area(data = subset(dfden1, Status == 'Diseased'), fill = "blue", alpha = 0.3, show.legend = FALSE) +
  geom_area(data = subset(dfden2, Status == 'Nondiseased'), fill = "red", alpha = 0.3, show.legend = FALSE) +
  geom_area(data = subset(dfden2, Status == 'Diseased'), fill = "red", alpha = 0.3, show.legend = FALSE) +
  geom_line() +
  geom_line(aes(x = dfden2$seq, y =dfden2$den, linetype = dfden2$Status)) +
  geom_text(x=0.1, y=0.45, label = "Males", colour = "blue", size = 10) +
  geom_text(x=10.5, y=0.45, label = "Females", colour = "red", size = 10) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  theme_bw() +
  theme_void() +
  theme(legend.position = "left", 
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

# ROC curve
p <- seq(0,1,by=0.01); np <- length(p)
roc1 <- pnorm(0.2 + qnorm(seq(0,1,by=0.01)))
roc2 <-  pnorm(3 + qnorm(seq(0,1,by=0.01)))
dfr1 <- data.frame(p = p, roc = roc1)
dfr2 <- data.frame(p = p, roc = roc2)

r1 <- ggplot(dfr1 ,aes(x = p, y = roc)) +
  labs(title = "", x = "FPF", y = "TPF") +
  geom_line(size = 1, colour = 'blue') +
  geom_abline(intercept = 0, slope =1, linetype=2, colour='grey', size = 1) +
  geom_line(aes(x=dfr2$p, y= dfr2$roc), linetype=1, colour='red', size = 1) +
  geom_text(x= 0.70, y = 0.2, label = "AUC Males: 0.556", colour = 'blue', size = 10) + 
  geom_text(x= 0.70, y = 0.10, label = "AUC Females: 0.982", colour = 'red', size = 10) + 
  theme_bw() +
  theme(aspect.ratio = 1,
       strip.text.x = element_text(size = 20),
       axis.text = element_text(size = 20),
       axis.title = element_text(size = 20),
       legend.position = "",
       plot.margin = unit(c(0, 0, 0, 0), "cm"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       strip.background = element_blank(),
       panel.border = element_rect(colour = "black", fill = NA))

# Arrange and plot
grid.arrange(d1, r1, ncol = 2)

