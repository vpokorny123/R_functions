target_oz <- 32

brew_ratio <- .1 # 1 to 4 ratio (i.e., 1/5)
coffee_portion <- target_oz * brew_ratio

  
costco_3lb_price <- 30
costco_price_per_oz <- costco_3lb_price / (16 * 3) #16 oz in a pound

costco_price <- coffee_portion * costco_price_per_oz 

chameleon_price <- 10.79
money_saved <- chameleon_price - costco_price
money_saved

money_saved_in_month <- money_saved *4
money_saved_in_month
