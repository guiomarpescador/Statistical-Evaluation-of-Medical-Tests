Exploratory Data Analysis
================
Guiomar Pescador-Barrios

This file contains the figures of the EDA in Chapter 6.

## Set-up

``` r
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.2, cex.sub=1)
```

``` r
# Test results AD, MCI and CN respectively
yd <- ADNI$tau[ADNI$DX == 3]
ym <- ADNI$tau[ADNI$DX == 2]
yh <- ADNI$tau[ADNI$DX == 1]
```

### Helper function for colours

``` r
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
```

## Test results

### Histograms

``` r
par(mfrow = c(1,3))
hist(yd, freq = T, main = "Histogram AD group", xlab = "Test oucomes")
hist(ym, freq = T, main = "Histogram MCI group",  xlab = "Test oucomes")
hist(yh, freq = T, main = "Histogram CN group", xlab = "Test oucomes")
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Boxplot

``` r
qplot(factor(DX), tau, data = ADNI, geom = c("boxplot"), fill = factor(DX)) +
  theme(legend.position="none") +
  scale_x_discrete(labels = c("CN", "MCI", "AD")) +
  labs(title="",x="Group", y = "Test oucomes")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.position = 'none')
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Age covariate

### Scatter plot

``` r
ggplot(data = ADNI) +
  geom_point(aes(x = age, y = tau, colour = factor(DX))) +
  scale_colour_discrete(name = "Group", labels = c("CN", "MCI", "AD")) +
 labs(title="",x="Age", y = "Test oucomes") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Barplot

``` r
barplot(table(ADNI[,c("DX","age")]),main="",col=ggplotColours(n=3),xlab="age", 
        legend=c("CN","MCI","AD"), args.legend=list(x="topleft",bty="n", title = "Group", cex = 1.5))
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Gender covariate

### Boxplot

``` r
par(mfrow = c(1, 1))
ggplot(data = ADNI) +
  geom_boxplot(aes(x = factor(gender, labels=c("Male", "Female")), y = tau),  
               fill = "light blue") +
  labs(title="",x="Gender", y = "Test oucomes")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.position = 'none')
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Distribution plot

``` r
tableGender <- table(ADNI$gender, ADNI$DX)/nrow(ADNI)
rownames(tableGender) <- c("Male", "Female")
colnames(tableGender) <- c("CN", "MCI", "AD")
plot(tableGender,main="",ylab="Group",xlab="Gender",col=ggplotColours(n=3), cex.axis=1.5, cex.lab = 1.5)
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## APOE4 covariate

### Boxplot

``` r
ggplot(data = ADNI) +
  geom_boxplot(aes(x = factor(APOE4), y = tau),  
               fill = "light blue") +
  labs(title="",x="APOE4", y = "Test oucomes") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.position = 'none')
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Distribution plot

``` r
tableAPOE4 <- table(ADNI$APOE4, ADNI$DX)/nrow(ADNI)
colnames(tableAPOE4) <- c("CN", "MCI", "AD")
plot(tableAPOE4, main="", ylab="Group", xlab="APOE4", col=ggplotColours(n=3), cex.axis=1.5, cex.lab = 1.5)
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## MMSE covariate

### Scatter plot

``` r
ggplot(data = ADNI) +
  geom_point(aes(x = MMSE, y = tau, colour = factor(DX))) +
  scale_colour_discrete(name = "Group", labels = c("CN", "MCI", "AD"))+
  labs(title="",x="MMSE", y = "Test oucomes") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### Barplot

``` r
barplot(table(ADNI[,c("DX","age")]),main="",col=ggplotColours(n=3),xlab="age", 
        legend=c("CN","MCI","AD"), args.legend=list(x="topleft",bty="n", title = "Group", cex = 1.5))
```

![](figures_EDA_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
