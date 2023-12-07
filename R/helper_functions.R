## A few helper functions and a ggplot theme

`%notin%`    <- Negate(`%in%`)

## Split up data frame into list with entries defined by 'col'
split_tibble <- function(tibble, col = 'col') {
  temp_list <- tibble %>% split(., .[, col])
  ## allow for multiple grouping columns, but drop those entries with no records
  temp_list[(lapply(temp_list, nrow) %>% unlist() > 0)]
} 

## Ugly brute force function to figure out parameters of positive gamma distribution
gamma_param_solve <- function(a_base, b_base, new_mean, new_var) {
  test_vals <- expand.grid(
    a_diff = seq(0, 12, by = 0.01)
  , b_diff = seq(-4, 4, by = 0.01)
  ) %>% mutate(
    new_a  = a_base + a_diff
  , new_b  = b_base + b_diff
  ) %>% mutate(
    desired_mean = new_mean
  , desired_var  = new_var
  ) %>% mutate(
    calc_mean    = new_a * new_b
  , calc_var     = new_a * new_b^2
  ) %>% mutate(
    diff_mean    = desired_mean - calc_mean
  , diff_var     = desired_var  - calc_var
  ) %>% mutate(
    joint_min    = abs(diff_mean) + abs(diff_var)
  , joint_min_p  = joint_min / min(joint_min)
  )
  
  test_vals %>% arrange(joint_min_p) %>% slice(1) 
  
}

## logistic functions
logit2   <- function(L, b, k, x) {
  L / (
    1 + exp(-k*(x - b))
  )
}
i_logit2 <- function(L, b, k, y) {
  (log((L - y) / y) / -k) + b
}

## ggplot theme
theme_set(theme_bw())
suppressWarnings(
  theme_update(
    axis.text.x = element_text(size = 16)
  , axis.text.y = element_text(size = 16)
  , axis.title.x = element_text(size = 16)
  , axis.title.y = element_text(size = 16)
  , legend.title = element_text(size = 12)
  , panel.grid.major = element_blank()
  , panel.grid.minor = element_blank()
  , strip.background = element_blank()
  , panel.margin = unit(0, "lines")
  , legend.key.size = unit(.55, "cm")
  , legend.key = element_rect(fill = "white")
  , panel.margin.y = unit(0.5, "lines")
  , panel.border = element_rect(colour = "black", fill = NA, size = 1)
  , strip.text.x = element_text(size = 16, colour = "black", face = "bold"))
)
