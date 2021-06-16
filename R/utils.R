# Check lgpr version
check_version <- function(pkg_name) {
  cur_ver <- packageVersion("lgpr")
  min_ver <- "1.1"
  if (cur_ver < min_ver) {
    msg <- sprintf(
      "lgpr version at least %s is needed, found %s", min_ver, cur_ver
    )
    stop(msg)
  }
}
