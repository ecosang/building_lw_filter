
# run this code to render notebook and publish to github page.
rmarkdown::render(input='code/liu_west_filter_for_a_building_gray_box_model.Rmd')
file.rename(from ="code/liu_west_filter_for_a_building_gray_box_model.html" , to ="code/filter.html")
file.copy(from="code/filter.html",to="docs/filter.html",overwrite =TRUE)
file.remove('code/filter.html')

