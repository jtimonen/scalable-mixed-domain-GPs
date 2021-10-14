# Create state-region map
create_region_map <- function() {
  state_groups <- list(
    c("ME", "NH", "VT", "MA", "RI", "CT"), c("NY", "NJ", "PA", "MD", "DE"), c("OH", "MI", "IL", "WI", "MN"),
    c("WA", "OR", "CA", "HI"), c("AZ", "CO", "NM", "NV"), c("IA", "NE", "KS", "ND", "SD"), c("KY", "TN", "MO", "WV", "IN"), c("VA", "OK", "FL", "TX", "NC"), c("AL", "MS", "LA", "GA", "SC", "AR"), c("MT", "ID", "WY", "UT", "AK")
  )
  region_names <- c(
    "New England", "Mid-Atlantic", "Midwest", "West Coast", "Southwest", "Plains", "Border South", "Outer South", "Deep South",
    "Mountain West"
  )
  state_region_map <- mapply(
    FUN = function(states, region) {
      data.frame(
        state = states,
        region = rep(region, length(states)),
        stringsAsFactors = F
      )
    }, state_groups, region_names,
    SIMPLIFY = F
  )
  state_region_map <- bind_rows(state_region_map) %>%
    arrange(state) %>%
    mutate(
      region_ind = as.integer(as.factor(region))
    )
  return(state_region_map)
}

# Add regions based on states
create_regions <- function(dat) {
  map <- create_region_map()
  regions <- rep("foo", nrow(dat))
  n_states <- nrow(map)
  for (j in 1:n_states) {
    s <- map$state[j]
    inds <- which(dat$state == s)
    regions[inds] <- map$region[j]
  }
  return(regions)
}

# Load the election data
load_election_data <- function() {
  dat <- readRDS("pres_vote_historical.RDS")
  dat <- dat %>% filter(state != "DC")
  dat$state <- as.factor(dat$state)
  dat$incumbent <- as.factor(dat$incumbent)
  dat$dem <- as.numeric(dat$dem)
  dat$rep <- as.numeric(dat$rep)
  dat$total <- as.numeric(dat$total)
  dat$rep_share <- dat$rep / (dat$dem + dat$total)

  # Add regions
  regions <- create_regions(dat)
  dat$region <- as.factor(regions)
  # dat$rinc <- as.integer(dat$incumbent == "R")
  return(dat)
}

# Create test data
create_test_x <- function(dat, years) {
  states <- levels(dat$state)
  df_state <- rep(states, each = length(years))
  df_year <- rep(years, times = length(states))
  df <- data.frame(df_year, as.factor(df_state))
  colnames(df) <- c("year", "state")
  regions <- create_regions(df)
  df$region <- as.factor(regions)
  return(df)
}
