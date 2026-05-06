install.packages("usethis")
library(usethis)
usethis::use_git_config(user.name = "Jameson Hinkle", user.email = "hinkyisme@gmail.com")

usethis::create_github_token()

install.packages("gitcreds")
library(gitcreds)
gitcreds::gitcreds_set()
# Choose option 2 and paste your token when prompted

usethis::use_git()
# Say yes when asked to commit existing files

usethis::use_github()
