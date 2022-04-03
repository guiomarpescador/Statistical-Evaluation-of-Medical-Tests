# Figure 1
rm(list = ls())
require(ggplot2)
library(gridExtra)

# Density curve 1
muh1 <- 6
mud1 <- 6

seqh1 <- seq(muh1 - 4, muh1 + 4, len = 200)
seqd1 <- seq(mud1 - 4, mud1 + 4, len = 200)

densd1 <- dnorm(seqd1, mud1, 1)
densh1 <- dnorm(seqh1, muh1, 1)

dfden1 <- data.frame(den = c(c(densh1), c(densd1)), 
                     seq    = c(seqh1, seqd1),
                     Status = c(rep("Nondiseased", length(seqh1)),
                                rep("Diseased", length(seqd1))))

d1 <- ggplot(dfden1, aes(x = seq, y = den, linetype = Status, group = Status)) +
  labs(title = "", x = "Test outcome", y = "Density") +
  coord_cartesian(xlim = c(range(dfden1$seq)[1], range(dfden1$seq)[2]), ylim =c(0, 0.45)) +
  geom_area(stat = 'identity', aes(fill=Status), alpha = 0.4) +
  geom_line() +
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_fill_manual(values=c("#0072B2","#D55E00")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.text.x = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
  
# ROC curve 1
p <- seq(0, 1, len = 100); np <- length(p)
roc1 <- 1 - pnorm(qnorm(1 - p, muh1, 1), mud1, 1)
auc1 <- sum(roc1)/np
dfr1 <- data.frame(p = p, roc = roc1)

r1 <- ggplot(dfr1 ,aes(x = p, y = roc)) +
  labs(title = "", x = "FPF", y = "TPF") +
  geom_line(size = 1) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.text.x = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        legend.position = "",
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))


# Density curve 2
muh2 <- 6
mud2 <- 6.75

seqh2 <- seq(muh2 - 4, muh2 + 4, len = 200)
seqd2 <- seq(mud2 - 4, mud2 + 4, len = 200)

densd2 <- dnorm(seqd2, mud2, 1)
densh2 <- dnorm(seqh2, muh2, 1)

dfden2 <- data.frame(den = c(c(densh2), c(densd2)), 
                     seq    = c(seqh2, seqd2),
                     Status = c(rep("Nondiseased", length(seqh2)),
                                rep("Diseased", length(seqd2))))

d2 <- ggplot(dfden2, aes(x = seq, y = den, linetype = Status)) +
  labs(title = "", x = "Test outcome", y = "Density") +
  coord_cartesian(xlim = c(range(dfden2$seq)[1], range(dfden2$seq)[2]), ylim =c(0, 0.45)) +
  geom_area(data = subset(dfden2, Status == 'Nondiseased'), fill = "#0072B2", alpha = 0.4, show.legend = FALSE) +
  geom_area(data = subset(dfden2, Status == 'Diseased'), fill = "#D55E00", alpha = 0.4, show.legend = FALSE) +
  geom_line() +
  scale_linetype_manual(values=c("dashed", "solid")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.text.x = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        legend.position = "",
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# ROC curve 2
roc2 <- 1 - pnorm(qnorm(1 - p, muh2, 1), mud2, 1)
auc2 <- sum(roc2)/np
dfr2 <- data.frame(p = p, roc = roc2)

r2 <- ggplot(dfr2 ,aes(x = p, y = roc)) +
  labs(title = "", x = "FPF", y = "TPF") +
  geom_line(size = 1) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.text.x = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        legend.position = "",
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# Density curve 3
muh3 <- 6
mud3 <- 11.5

seqh3 <- seq(muh3 - 4, muh3 + 4, len = 200)
seqd3 <- seq(mud3 - 4, mud3 + 4, len = 200)

densd3 <- dnorm(seqd3, mud3, 1)
densh3 <- dnorm(seqh3, muh3, 1)

dfden3 <- data.frame(den = c(c(densh3), c(densd3)), 
                     seq    = c(seqh3, seqd3),
                     Status = c(rep("Nondiseased", length(seqh3)),
                                rep("Diseased", length(seqd3))))
d3 <-ggplot(dfden3, aes(x = seq, y = den, linetype = Status)) +
  labs(title = "", x = "Test outcome", y = "Density") +
  coord_cartesian(xlim = c(range(dfden3$seq)[1], range(dfden3$seq)[2]), ylim =c(0, 0.45)) +
  geom_area(data = subset(dfden3, Status == 'Nondiseased'), fill = "#0072B2", alpha = 0.4, show.legend = FALSE) +
  geom_area(data = subset(dfden3, Status == 'Diseased'), fill = "#D55E00", alpha = 0.4, show.legend = FALSE) +
  geom_line() +
  scale_linetype_manual(values=c("dashed", "solid")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.text.x = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        legend.position = "",
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# ROC curve 3
roc3 <- 1 - pnorm(qnorm(1 - p, muh3, 1), mud3, 1)
auc3 <- sum(roc3)/np
dfr3 <- data.frame(p = p, roc = roc3)

r3 <- ggplot(dfr3 ,aes(x = p, y = roc)) +
  labs(title = "", x = "FPF", y = "TPF") +
  geom_line(size = 1) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.text.x = element_text(size = 15), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        legend.position = "",
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# Arrange and plot
grid.arrange(d1,d2,d3,r1,r2,r3, ncol =3)
