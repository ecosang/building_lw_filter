
## 

This is accompanies code of *Appendix C* in paper
[“xxx”](https://dx.doi.org/). This is a simple demonstaration of
Liu-West filter for a simple building grey-box model. We generated a
synthetic data and applied Liu-West filter to see if the filter can be
applicable for this problem. All code is written in `R` language. The
main purpose of this document is to provide reproducible example.

The pacakge depedency of this code is managed by `renv` package. You can
look at `renv.lock` file to see the required package. However, for the
simplicity, just run following script on this rproject.

After clone the repository and run `building_lw_filter.Rproj` file (you
must have [Rstudio](https://rstudio.com/products/rstudio/)) with
[`R>3.5.3`](https://www.r-project.org/).

    git clone https://github.com/ecosang/building_lw_filter.git

In R console, run following script.

``` r
# Run this code in Rstudio
install.packages('renv',repos="https://cran.rstudio.com")
renv::equip() #install required software
renv::restore()
```

Technically, this installs all required packages, and you can reproduce
all codes below. However, if this doesn’t work, please report
[issues](https://github.com/ecosang/building_lw_filter/issues).

All functions used in this code is in `code/utility.R`. Also, all
generated data and trained model are stored in `data` folder.


[Complete code notebook](https://ecosang.github.io/building_lw_filter/filter.html)

Render notebook

```{r}
#with .Rproj
rmarkdown::render(input='code/liu_west_filter_for_a_building_grey_box_model.Rmd', 
output_file = '../docs/notebook.html',knit_root_dir = getwd())
```
