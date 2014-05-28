riverLabeller <- function(var, value){
  value <- as.character(value)
  if (var=="site") { 
    value[value=="1"] <- "WB"
    value[value=="2"] <- "OL"
    value[value=="3"] <- "OS"
    value[value=="4"] <- "Is"
  }
  return(value)
}

riverLabeller2 <- function(var, value){
  value <- as.character(value)
  if (var=="site") { 
    value[value=="WEST BROOK"] <- "WB"
    value[value=="WB JIMMY"] <- "OL"
    value[value=="WB MITCHELL"] <- "OS"
    value[value=="WB OBEAR"] <- "Is"
  }
  return(value)
}

# Function that returns Root Mean Squared Error
rmse <- function(error) {
  sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error) {
  mean(abs(error))
}
