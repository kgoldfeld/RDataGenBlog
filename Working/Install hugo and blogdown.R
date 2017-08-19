# From directions:
# https://proquestionasker.github.io/blog/Making_Site/

devtools::install_github("rstudio/blogdown")

# installed homebrew (using terminal):
# /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

blogdown::install_hugo()
install_theme("halogenica/beautifulhugo", theme_example = FALSE, update_config = FALSE)
new_site()

