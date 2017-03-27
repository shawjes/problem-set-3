knitr::opts_chunk$set(echo = TRUE, comment = "#>")
library(tidyverse)

# bgh data
bgh_url <- 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE43163&format=file&file=GSE43163_CompleteCountTable_Bgh.txt.gz'

# 2 info and blank lines at top, skip them
raw_data <- read_tsv(bgh_url, skip = 2)

# the header for the first column is "NA", set it manually
names(raw_data)[1] <- 'gene.name'

raw_data <- as_data_frame(raw_data)

raw_data01<-gather(raw_data, plant_fungus_time.point_rep_value, value=value, 2:49)
head(raw_data01)
attach(raw_data01)
cleaned_data<-separate(raw_data01, "plant_fungus_time.point_rep_value",
                       into = c("plant","fungus","time.point","rep"), sep = "_")
# report the cleaned_data by just naming it, uncomment the following line:
cleaned_data

# A tibble: 310,896 × 6
#               gene.name plant fungus time.point   rep value
# *                <chr> <chr>  <chr>      <chr> <chr> <int>
# 1             bgh04079   B12     A6       6hpi     1    13
# 2             bgh01634   B12     A6       6hpi     1    31
# 3  bghG000012000001001   B12     A6       6hpi     1   121
# 4  bghG000012000002001   B12     A6       6hpi     1     3
# 5             bgh00757   B12     A6       6hpi     1   253

attach(cleaned_data)

# Problem 2:
# Which plant has the highest expression of any gene in the 6hpi time point?
filter(cleaned_data,time.point=="6hpi")%>%top_n(1,value) %>% select(plant)

# Which plant / fungus pair has the highest expression in the 18hpi time point?
filter(cleaned_data,time.point=="18hpi") %>% top_n(1,value) %>% select(plant,fungus)


# Problem 3:
library(tidyr)
#install.packages("stringr")
library(stringr)
#install.packages("dplyr")
library(dplyr)

cleaned_data01<-cleaned_data %>% mutate(time.value = str_replace(time.point, 'hpi', ''),
         time.hpi = str_c('hpi.', time.value))
# A tibble: 310,896 Ã— 8
# gene.name plant fungus time.point   rep value time.value time.hpi
# *                <chr> <chr>  <chr>      <chr> <chr> <int>      <chr>    <chr>
#   1             bgh04079   B12     A6       6hpi     1    13          6    hpi.6
# 2             bgh01634   B12     A6       6hpi     1    31          6    hpi.6
# 3  bghG000012000001001   B12     A6       6hpi     1   121          6    hpi.6
# 4  bghG000012000002001   B12     A6       6hpi     1     3          6    hpi.6
# 5             bgh00757   B12     A6       6hpi     1   253          6    hpi.6
# 6             bgh01273   B12     A6       6hpi     1    45          6    hpi.6
# 7             bgh01274   B12     A6       6hpi     1    32          6    hpi.6
# 8             bgh01277   B12     A6       6hpi     1     9          6    hpi.6
# 9             bgh06140   B12     A6       6hpi     1    47          6    hpi.6
# 10            bgh05774   B12     A6       6hpi     1    42          6    hpi.6
# ... with 310,886 more rows

attach(cleaned_data)


# Q3:
# Identify the top 3 most consistently differentially expressed genes between the earliest and latest
# time points for each combination of plant and fungus strains.
# "Differential expression"" is the difference between value (i.e., gene expression level) between time points.
# "Consistency" is the smallest variance in value between replicates.

# Strategy:
# Create a new table from the cleaned data by moving each hpi value to a new column name, with counts for each
# in the column (hint: use a tidyr verb).
# It is helpful to reformat the hpi values by converting from e.g. 6hpi to hpi.6.
# You can use mutate to do this, i.e.:
# Create a new column containing the expression difference between the relevant time points.
# Calculate summary statistics (mean and variance) of the expression differences by grouping (hint) the gene.name, plant, and virus columns.
# Sort by these statistics and use the dplyr verb slice to pull the ones you want (i.e., the top 3). Note you will have to remove gene.name from the grouping so that sorting works.

help(spread)

relative_consistency<-
  cleaned_data01%>%
  mutate(time.value = str_replace(time.point, 'hpi', ''), time.hpi = str_c('hpi.', time.value))%>%
  select(-time.point, -time.value)%>%
  # why does the spread function not work if I don't remove time.point and time.value?
  spread(time.hpi, value)%>%
  mutate(expr.diff = hpi.24 - hpi.6)%>%
  mutate(abs.expr.diff = abs(hpi.24 - hpi.6))%>%
  group_by(gene.name,plant,fungus)%>%
  summarise(mean.abs.diff=mean(abs.expr.diff),var.signed.diff=(0.00001+var(expr.diff,na.rm=TRUE)))%>%
  ungroup()%>%
  mutate(grp=paste(plant,fungus,sep="-")) %>%
  mutate(grp_gene=paste(gene.name,plant,fungus,sep="-")) %>%
  mutate(relative_consistency=var.signed.diff/mean.abs.diff) %>%
  group_by(grp)%>%
  mutate(unique=n_distinct(grp_gene,na.rm=TRUE)) %>%
  # I'm assuming that if variance=0 then we only have one record for that gene within that plant x fungus group. Therefore not a good basis
  # for finding consistently differentially expressed genes.
  #  filter(var.signed.diff>0)%>%
  # arrange(desc(unique))%>% # Each (gene x plant x fungus) has between 6266 and 6301 entries. Therefore, those with
  # variance of zero are actually valid records (i.e., variance is not zero due to only having one sample).
  #  filter(mean.abs.diff,cume_dist(grp)>=0.99)%>%
  #  arrange(var.signed.diff,desc(mean.abs.diff))%>%
  arrange(relative_consistency)%>%
  slice(1:3)


head(relative_consistency)


# Produces same result as above:
#    cleaned_data%>%
#    mutate(time.value = str_replace(time.point, 'hpi', ''), time.hpi = str_c('hpi.', time.value))%>%
#    select(-time.point, -time.value)%>%
#    spread(time.hpi, value)%>%
#    mutate(expr.diff = hpi.24 - hpi.6)%>%
#    mutate(abs.expr.diff = abs(hpi.24 - hpi.6))%>%
#    group_by(gene.name,plant,fungus)%>%
#    summarise(mean.abs.diff=mean(abs.expr.diff),var.signed.diff=(var(expr.diff,na.rm=TRUE)))%>%
#    ungroup()%>%
#    mutate(grp=paste(plant,fungus,sep="-")) %>%
#    mutate(grp_gene=paste(gene.name,plant,fungus,sep="-")) %>%
#    mutate(relative_consistency=var.signed.diff/mean.abs.diff) %>%
#    group_by(grp)%>%
#    mutate(unique=n_distinct(grp_gene,na.rm=TRUE)) %>%
#  I'm assuming that if variance=0 then we only have one record for that gene within that plant x fungus group. Therefore not a good basis
#  for finding consistently differentially expressed genes.
#  arrange(desc(unique))%>% # Each (gene x plant x fungus) has between 6266 and 6301 entries. Therefore, those with
#  variance of zero are actually valid records (i.e., variance is not zero due to only having one sample).
# filter(mean.abs.diff,cume_dist(grp)>=0.99)%>% # This would be useful if, for example, we knew that only differences greater than a certain threshold
# would be clinically significant - e.g., produce a phenotypic change.
#   arrange(var.signed.diff,desc(mean.abs.diff))%>%
#   slice(1:3)

# Q4:
# Now examine the above final data frame above and write a few sentences putting inline code in least 3 places.
# For example, There are r nrow(mtcars) rows of data in mtcars.
The above data frame identifies the top three most consistently differentially expressed genes between the earliest
and latest time points for each of combination of plant and fungus strains. Stratifying simultaneously by both plant and fungus
produces `r group_by(plant,fungus)` groups. Each group contains a different number of genes. Examining all genes found within each grouping of
plant and fungus results in a data frame with `r group_by(gene.name,plant,fungus)` rows. Meanwhile, there are only `r length(unique(gene.name))`
unique genes, which indicates that our plant x fungus groups share overlapping genes.


# Q5:
# Plot the expression (value) by time (hpi) for the above genes. Format the plot as follows:
# Plot each point.
# Connect the points with a smooth line
# Create a facet for each plant / bgh pair
# Assign a different color to each gene
# Assign a different shape to each replicate
# Add a descriptive title to the plot
# Use the "black & white" theme
# Extra credit: add error bars to the plot (use geom_errorbar).
```{r problem_4, warning=FALSE}
library(ggplot2)
help(left_join)

plot_clean<-cleaned_data01 %>%
  mutate(time.value = str_replace(time.point, 'hpi', ''), time.hpi = str_c('hpi.', time.value)) %>%
  select(gene.name, plant, fungus, value, time.value,rep)
plot_clean
plot_rel_cons<-select(relative_consistency,gene.name,plant,fungus)
plot_rel_cons
plot.data<-plot_rel_cons %>%
  ungroup() %>%
  #return all rows from x, and all cols from x and y
  left_join(., plot_clean) %>%
  mutate(time.value = as.integer(time.value))
dim(plot_clean) # Full cleaned data has dimensions 310,896 × 5
dim(plot_rel_cons) # Data frame listing most consistently differentially expressed genes has dimensions 12 x 4
dim(plot.data) # Left-joined data frame has dimensions 144 x 6
head(plot.data)
tail(plot.data,n=100)
relative_consistency

help(ggplot)
library(ggplot2)
# Plot the expression (value) by time (hpi) for the above genes. Format the plot as follows:
# Plot each point.
help(theme_bw)

# Scatterplot of mpg vs. hp for each combination of gears and cylinders
# in each facet, transmittion type is represented by shape and color
help(qplot)
Q4plot<-qplot(time.value, mean.value, data=df, color=gene.name, #shape=rep, 
      facets=plant~fungus, size=I(3),
      xlab="Time", ylab="Mean expression")+
      #main="Top consistently differentialy expressed genes within plant and fungus")+
      labs(title="Top consistently differentialy expressed genes within plant and fungus",
           subtitle="Mean expression values over time")+
      geom_smooth() +
      theme_bw() +
      theme(legend.key = element_blank())
    
help(qplot)
