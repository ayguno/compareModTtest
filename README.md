# compareModTtest
Designed to help us to interpret similarities between various treatments based on the proteins

This function is designed to help us to interpret similarities between various treatments based on the proteins
being co-effected by the individual treatments
Multiple TMT10 experiments can be compared
"ModT test table results" printed by the Shiny server will be provided into the same working directory, along with the class vector files specific for each experiment (both in csv format)
The function will take a customly prepared csv file as input, which matches the data tables with class vector files
The function will generate heatmaps of top 10 up/down regulated proteins for each experiment
The function will also attempt to plot the entire set of experiments into a single scatterplot
ModT_ClassV_match (chr) = csv table that matches ModT table file names (first column), with the respective Class vector files (second column) 

Functional as of 06/21/2016 (last update)
