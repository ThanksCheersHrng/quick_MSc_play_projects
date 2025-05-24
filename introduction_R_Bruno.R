4 + 2

# You should comment your code with the hashtag symbol

### You can use more than one if you prefer, but only one is necessary.

############################################# 
# You can use to mark sessions of your code #
#############################################

a <- 1
a

A

aa <- 2; class(aa)

bb <- 3L; class(bb)

class(aa > bb)

class("Bruno")

my_first_vector <- c(1, 2, 3)

my_second_vector <- c(4, 5, 6)

my_first_vector + my_second_vector

my_first_vector * my_second_vector

object <- c(1, 2, 3, "Any string")

object

numbers_logic <- c(1, 2, 3, TRUE, FALSE)

numbers_logic

my_list <- list(c(1, 2, 3), 
                "my name", 
                c(123, 123, 123), 
                c("A", 2))

my_list

class(my_list[[4]])

my_other_list <- list(l1 = c(1, 2, 3), 
                      l2 = "my name", 
                      l3 = c(123, 123, 123))

my_other_list$l1

my_other_list[[2]]

m1 <- cbind(c(1, 2, 3), 
           c(22, 22, 22))
m1

m2 <- rbind(c(11, 11, 11), 
            c(1, 2, 3))
m2

matrix(1:9, nrow = 3)

matrix(1:9, nrow = 3, byrow = TRUE)

d1 <- data.frame(var1 = c(11, 11, 11), 
                 var2 = c(1, 2, 3),
                 var3 = c("A", "B", "A"))
d1

d1$var2

v1 <- c(123, 456, 789)
v1[2]

c1 <- rbind(1, 2, 3)
c1[2, ]

c2 <- cbind(11, 22, 33)
c2[1, 2]

## beatles <- read.csv("beatles.csv")

## beatles <- readr::read_csv("beatles.csv")

## library(readr)
## beatles <- read_csv("beatles.csv")

## haven::read_sas()
## haven::read_stata()
## haven::read_sav()

hello_world <- function() print("Hello, World")

hello_world <- function() {
  print("Hello, World")
} 

hello_world_v2 <- function() {
  print("Hello, World")
  return(0)
}
a <- hello_world_v2()
a

hello_world_v2 <- function() {
  print("Hello, World")
  0
}

print_message <- function(arg) {
  print(arg)
}
print_message("New message here")

f1 <- function(arg1, arg2) {
  print(paste0("1st argument: ", arg1, "; 2nd argument ", arg2, "."))
}
f1(1, 2)

f1(arg2 = 1, arg1 = 2)

f1(arg1 = 1)

f2 <- function(arg1 = 1, arg2 = 2) {
  print(paste0("1st argument: ", arg1, "; 2nd argument ", arg2, "."))
}
f2()

f2(arg2 = 3)

f3 <- function(vector_number, ...){
  matrix(vector_number, ...)
}

f3(1:6, ncol = 2)

f3(1:6, nrow = 2)

f3(c(1:6), nrow = 2, byrow = TRUE)

args(summary.lm)

test_function <- function(){
  obj_inside <- 2
  print(paste("Value object = ", obj_inside))
}

test_function()

obj_inside

## plot(x = 1, y = 1)

plot(x = 1, y = 1)

plot(x = c(1, 2), y = c(1, 2))

x2 <- c(1, 2, 3)
y2 <- c(2, 4, 6)
plot(x = x2, y = y2)

plot(x = x2, y = y2, 
     main = "My first plot", xlab = "X label", 
     ylab = "Y label", col = 2, type = 'l', xlim = c(1.5, 3))

plot(x = x2, y = y2, 
     main = "Plot with horizontal lines", 
     col = "blue", type = 'h')

## ?dpois

(exp(-3) * 3 ^2) / 2

dpois(x = 2, lambda = 3)

x <- 0:10
dpois(x = x, lambda = 3)

prob_values <- dpois(x = x, lambda = 3)

## plot(x, prob_values, type = 'h',
##      main = "Probability function \n Random variable with Poisson distribution",
##      ylab = expression("P(X = x)"))

plot(x, prob_values, type = 'h', 
     main = "Probability function \n Random variable with Poisson distribution", 
     ylab = expression("P(X = x)"))

plot(0:20, dpois(0:20, 10), type = 'h', 
     main = "Probability function \n Random variable with Poisson distribution", 
     ylab = expression("P(X = x)"))

n <- 23
1 - prod(365:(365-n+1))/(365^n)

birthday_problem <- function(initial_value, final_value){
  seq_values <- initial_value:final_value
  prob_function <- function(n) 1 - prod(365:(365-n+1))/(365^n)
  values_prob <- sapply(seq_values, prob_function)
  plot(seq_values, values_prob, type = "l", 
       ylab = "Probabilities", xlab = "# People")  
}

## birthday_problem(10, 60)

birthday_problem(10, 60)

birthday_problem(60, 10)

birthday_problem <- function(initial_value, final_value, 
                             makes_plot = TRUE){
  seq_values <- initial_value:final_value
  prob_function <- function(n) 1 - prod(365:(365-n+1))/(365^n)
  values_prob <- sapply(seq_values, prob_function)
  if (makes_plot){
    plot(seq_values, values_prob, type = "l", 
       ylab = "Probabilities", xlab = "# People")  
  }
  values_prob
}

probs <- birthday_problem(10, 30, makes_plot = FALSE)

head(probs)

probs <- birthday_problem(10, 30)
head(probs)

birthday_problem <- function(chosen_n, initial_value, final_value, 
                             makes_plot = TRUE){
  seq_values <- initial_value:final_value
  prob_function <- function(n) 1 - prod(365:(365-n+1))/(365^n)
  values_prob <- sapply(seq_values, prob_function) 
  if (makes_plot){
    plot(seq_values, values_prob, type = "l", 
       ylab = "Probabilities", xlab = "# People")  
  }
  chosen_prob <- values_prob[chosen_n == seq_values]
  print(paste0("The probability for n =  ", chosen_n, 
               " is equal to ", round(100 * chosen_prob, 3), "%"))
  values_prob
}

new_probs <- birthday_problem(23, 10, 30, makes_plot = FALSE)

head(new_probs)

birthday_problem <- function(chosen_n, initial_value, final_value, 
                             makes_plot = TRUE){
  seq_values <- initial_value:final_value
  prob_function <- function(n) 1 - prod(365:(365-n+1))/(365^n)
  values_prob <- sapply(seq_values, prob_function) 
  if (makes_plot){
    plot(seq_values, values_prob, type = "l", ylab = "Probabilities", xlab = "# People")  
    plot_made <- recordPlot()   
  }
  else plot_made <- NULL   
  chosen_prob <- values_prob[chosen_n == seq_values]
  print(paste0("The probability for n =  ", chosen_n, 
               " is equal to ", round(100 * chosen_prob, 3), "%"))
  
  list(seq_values = seq_values, 
       probabilies_plot = plot_made,
       chosen_n = chosen_n,
       chosen_probability = chosen_prob, 
       values_probability = values_prob)
}

list_prob <- birthday_problem(42, 10, 60, makes_plot = FALSE)

list_prob$chosen_n

list_prob$chosen_probability

head(list_prob$values_prob)

list_prob <- birthday_problem(42, 10, 60, makes_plot = TRUE)

x <- list_prob$chosen_n
y <- list_prob$chosen_probability
list_prob$probabilies_plot
abline(v = x, lty = 2)
abline(h = y, lty = 2)
text(x, y, adj = c(0, 1), 
     labels =  paste0("n = ", x, "\n ", 
             "Prob = ", round(y, 3)))



beatles <- readr::read_csv("beatles.csv")

library(dplyr)
beatles1 <- dplyr::select(beatles, 
                          track_name, album_release_year, 
                          danceability, duration_ms)

beatles2 <- mutate(beatles1, 
                   dance_binary = ifelse(danceability > 0.5, 1, 0))

beatles3 <- filter(beatles2, 
                   album_release_year < 1971)

beatles4 <- arrange(beatles3, 
                    duration_ms)

beatles_final <- arrange(filter(mutate(dplyr::select(beatles, track_name, 
                                                     album_release_year,
                                                     danceability, 
                                                     duration_ms), 
                                       dance_binary = ifelse(
                                         danceability > 0.5, 1, 0)), 
                                album_release_year < 1971), duration_ms)

beatles_final <- 
  arrange(
    filter(
      mutate(dplyr::select(beatles, 
                           track_name, album_release_year, 
                           danceability, duration_ms),
             dance_binary = ifelse(danceability > 0.5, 1, 0)),
      album_release_year < 1971), 
    duration_ms)

beatles_final <- beatles %>%
  dplyr::select(track_name, album_release_year, 
                danceability, duration_ms) %>% 
  mutate(dance_binary = ifelse(danceability > 0.5, 1, 0)) %>% 
  filter(album_release_year < 1971) %>% 
  arrange(duration_ms)

head(beatles_final)

1:10 %>% "^"(2, .) 

list_info <- beatles %>% 
  dplyr::select(track_name, duration_ms) %>% 
  mutate(duration_min = (duration_ms/1000)/60) %>% 
  filter(duration_min > 4) %>% 
  list(all_data = beatles, 
       subset_data = .)

dim(list_info$all_data)

dim(list_info$subset_data)
