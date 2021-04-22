---
layout: post
title:  "Data Visualization in R"
date:   2021-04-06 20:35:58 +0800
usemathjax: true
categories: jekyll update
---

- Packages may useful for visualization, ggplot2, lattice, and others
- The grid system in R - draw any thing you want
- Recreation of potentially useful figures


- The overall apperence of the figure

```R
library(ggplot2)
t <- theme_classic() + # Set to a white theme
  theme(axis.text = element_text( size = 16,family="Arial",color="black"), # x tick label and y tick label
        axis.title = element_text(size =18, face="bold",family="Arial",color="black"), # x label and y label
        legend.text = element_text( size = 14,family="Arial",color="black"),
        legend.title = element_text( size = 16,family="Arial",color="black"),
        axis.line = element_blank(), # remove x axis and y axis
        panel.background = element_rect(fill="#F5F5F5"), # Background color
        panel.border = element_rect(color="black",fill=NA,size=1) # border of the figure
  )  # Control details of the figures appearance
```

- Scatter plot

```R
g <- ggplot(iris, aes(x=Petal.Length, y=Sepal.Length,col=Species)) + 
  geom_point(size=2, shape=22) + # Size and Shape of the dot
  #scale_colour_viridis_d() + # Set color scale
  scale_color_manual(values=c("#104E8B","#BCD2EE","#70B5B3")) + 
  xlab("Petal") + # Change x title
  ylab("Sepal") + # Change y title
  t
```

- Boxplot / Violin plot

```R
g <- ggplot(iris, aes(x=Species, y=Sepal.Length,fill=Species)) + 
  geom_boxplot(width=0.4) + # Size and Shape of the dot
  # geom_violin() + 
  scale_fill_manual(values=c("#104E8B","#BCD2EE","#70B5B3")) + # Set color scale
  xlab("Petal") + # Change x title
  ylab("Sepal") + # Change y title
  t +
  theme(axis.text.x = element_text(angle = 45,hjust = 1) # x tick label and y tick label, right-horizonton alignment 
  ) 
```

- Use `grid` package for more customized plots

- Heatmap

```R
library(pheatmap)
library(pals)
# Generate a demo
exp.mat = matrix(rnorm(200), 20, 10)
exp.mat[1:10, seq(1, 10, 2)] = exp.mat[1:10, seq(1, 10, 2)] + 3
exp.mat[11:20, seq(2, 10, 2)] = exp.mat[11:20, seq(2, 10, 2)] + 2
exp.mat[15:20, seq(2, 10, 2)] = exp.mat[15:20, seq(2, 10, 2)] + 4
colnames(exp.mat) = paste("Sample", 1:10, sep = "")
rownames(exp.mat) = paste("Gene", 1:20, sep = "")
pheatmap(exp.mat,color = pals::coolwarm(100))
```

## Useful links
- <https://www.r-graph-gallery.com/>
- <https://bookdown.org/rdpeng/RProgDA/the-grid-package.html>

