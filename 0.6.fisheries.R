library(tidyverse)
library(parallel)
library(doParallel)

## pull in the data
fish_catch <- read.csv("~/GitHubRepos/SC.Iceland/fisheries/Fisheries.csv")
cod <- read.csv("fisheries/coddetprop.csv")
age <- read.csv("metadata/AgeOut.csv")


cod$age <-  1950-age$mean[match(cod$X,gsub("-","_",age$ID))]


###parallel function to calcualte bootstraps

bootstrap_yearly_catch_parallel <- function(df, n_simulations = 10000) {
  
  # Define confidence_levels *inside* the function
  confidence_levels <- c("red" = 0.20, "amber" = 0.10, "green" = 0.05)
  
  # Identify columns for each fishery
  fishery_cols <- c("NorFish17_Iceland", 
                    "NorFish18_Dutch", 
                    "NorFish19_French", 
                    "NorFish20_English")
  
  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  results_matrix <- foreach(i = seq_len(n_simulations), .combine = rbind) %dopar% {
    yearly_totals <- numeric(nrow(df))
    
    for (j in seq_len(nrow(df))) {
      total_catch <- 0
      for (fishery in fishery_cols) {
        conf_suffix <- gsub("NorFish[0-9]+_", "", fishery)
        conf_col <- paste0("Conf_", conf_suffix)
        
        catch_val <- df[[fishery]][j]
        conf_val <- df[[conf_col]][j]
        sd_percent <- confidence_levels[[conf_val]]
        
        sampled_catch <- rnorm(1, mean = catch_val, sd = catch_val * sd_percent)
        sampled_catch <- max(0, sampled_catch)
        
        total_catch <- total_catch + sampled_catch
      }
      yearly_totals[j] <- total_catch
    }
    yearly_totals
  }
  
  stopCluster(cl)
  colnames(results_matrix) <- df$Year
  results_matrix
}


# --- Run the parallelized bootstrap ---
set.seed(123)
sim_results <- bootstrap_yearly_catch_parallel(fish_catch, n_simulations = 10000)

# Convert to a tidy data frame for plotting
df_plot <- as_tibble(sim_results) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Year",
    values_to = "SimulatedCatch"
  ) %>%
  group_by(Year) %>%
  summarise(
    mean_catch = mean(SimulatedCatch),
    lower_95   = quantile(SimulatedCatch, 0.025),
    upper_95   = quantile(SimulatedCatch, 0.975)
  ) %>%
  mutate(Year = as.numeric(Year))


plot(df_plot$Year, df_plot$mean_catch, type = "l", col = "blue", lwd = 2,
     xlab = "Year", ylab = "Catch", main = "Bootstrap Mean Catch with 95% CI",
     ylim = range(df_plot$lower_95, df_plot$upper_95))

# Add shaded confidence interval
polygon(c(df_plot$Year, rev(df_plot$Year)), 
        c(df_plot$upper_95, rev(df_plot$lower_95)), 
        col = rgb(0, 0, 1, 0.3), border = NA)

points(fish_catch$Year,fish_catch$Sum,type="l",lwd=2,col="darkred")

####great now lets blur the signal to give a number that is a bit realistic

# Gaussian weighting function with edge correction
get_gaussian_weights <- function(window_size = 10, actual_size = NULL) {
  x <- seq(-window_size, window_size, length.out = 2 * window_size-1)  # Properly symmetric window
  weights <- dnorm(x, mean = 0, sd = window_size / 2)  # Normal distribution weights
  
  # If the window is smaller than expected (at dataset edges)
  if (!is.null(actual_size)) {
    center_idx <- ceiling(length(weights) / 2)  # Center index
    half_size <- floor(actual_size / 2)  # How much to take from each side
    
    start_idx <- max(1, center_idx - half_size)  # Avoid going below index 1
    end_idx <- min(length(weights), center_idx + half_size)  # Avoid exceeding length
    
    weights <- weights[start_idx:end_idx]  # Take a **centered** subset
  }
  
  weights <- weights / sum(weights)  # Normalize so weights sum to 1
  return(weights)
}

# Parallel bootstrap with moving average smoothing (corrected edges)
bootstrap_yearly_catch_parallel_blurred <- function(df, n_simulations = 10000) {
  
  confidence_levels <- c("red" = 0.50, "amber" = 0.10, "green" = 0.05)
  fishery_cols <- c("NorFish17_Iceland", 
                    "NorFish18_Dutch", 
                    "NorFish19_French", 
                    "NorFish20_English")
  
  num_cores <- min(detectCores() - 1, 8)  # Limit cores for stability
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  results_matrix <- tryCatch({
    
    foreach(i = seq_len(n_simulations), .combine = rbind,
            .export = c("get_gaussian_weights")) %dopar% {
              
              yearly_totals <- numeric(nrow(df))
              
              for (j in seq_len(nrow(df))) {
                total_catch <- 0
                for (fishery in fishery_cols) {
                  conf_suffix <- gsub("NorFish[0-9]+_", "", fishery)
                  conf_col <- paste0("Conf_", conf_suffix)
                  
                  catch_val <- df[[fishery]][j]
                  conf_val <- df[[conf_col]][j]
                  sd_percent <- confidence_levels[[conf_val]]
                  
                  sampled_catch <- rnorm(1, mean = catch_val, sd = catch_val * sd_percent)
                  sampled_catch <- max(0, sampled_catch)
                  
                  total_catch <- total_catch + sampled_catch
                }
                yearly_totals[j] <- total_catch
              }
              
              # Apply Gaussian-weighted moving average smoothing with edge handling
              smoothed_totals <- yearly_totals  
              
              for (j in seq_len(nrow(df))) {
                # Determine valid range for this window
                min_idx <- max(1, j - 5)
                max_idx <- min(nrow(df), j + 5)
                
                actual_size <- max_idx - min_idx + 1  # Actual number of valid years in window
                weights <- get_gaussian_weights(10, actual_size)  # Adjust weights safely
                
                smoothed_totals[j] <- sum(yearly_totals[min_idx:max_idx] * weights)
              }
              
              smoothed_totals
            }
  }, finally = {
    stopCluster(cl)  # Ensure cluster shuts down
  })
  
  colnames(results_matrix) <- df$Year
  results_matrix
}

# Run the modified parallel bootstrap with edge handling
set.seed(123)
sim_results_blurred <- bootstrap_yearly_catch_parallel_blurred(fish_catch, n_simulations = 10000)
# Convert to a tidy data frame for plotting
df_plot <- as_tibble(sim_results_blurred ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Year",
    values_to = "SimulatedCatch"
  ) %>%
  group_by(Year) %>%
  summarise(
    mean_catch = mean(SimulatedCatch),
    lower_95   = quantile(SimulatedCatch, 0.025),
    upper_95   = quantile(SimulatedCatch, 0.975)
  ) %>%
  mutate(Year = as.numeric(Year))


plot(df_plot$Year, df_plot$mean_catch, type = "l", col = "blue", lwd = 2,
     xlab = "Year", ylab = "Catch (Metric Tonnes)", main = "",
     ylim = range(12000, df_plot$upper_95),xlim = c(1500,1875))

# Add shaded confidence interval
polygon(c(df_plot$Year, rev(df_plot$Year)), 
        c(df_plot$upper_95, rev(df_plot$lower_95)), 
        col = rgb(0, 0, 1, 0.3), border = NA)

points(fish_catch$Year,fish_catch$Sum,type="l",lwd=2,col="darkred")

###Let's plot cod and fisheries adjacent 

dev.off()
# ---- Define y-axis scaling for cod detections ----
cod_y_base <- 7000  # Set base of cod vertical lines
cod_max_height <- 13500  # Maximum height for cod lines

# Scale cod detection proportions (`cod$x`) into the cod height range
cod_heights <- scales::rescale(cod$x, to = c(cod_y_base, cod_max_height))


pdf("fisheries/fishery.pdf",height = 6,width = 9)
# ---- Main plot: Fisheries & Cod Detection ----
plot(df_plot$Year, df_plot$mean_catch, type = "l", col = "blue", lwd = 2,
     xlab = "Year", ylab = "Catch (Mt)", main = "Bootstrap Mean Catch with 95% CI",
     ylim = c(8000, 47000),  # Extend y-axis for cod detections
     xlim = c(min(df_plot$Year), max(df_plot$Year)))  # Keep x-axis consistent

# Add shaded confidence interval
polygon(c(df_plot$Year, rev(df_plot$Year)), 
        c(df_plot$upper_95, rev(df_plot$lower_95)), 
        col = rgb(0, 0, 1, 0.3), border = NA)


# Add fishery sum data (red line for total catch)
#lines(fish_catch$Year, fish_catch$Sum, lwd = 2, col = "darkred")

# ---- Overlay cod detections as vertical line segments ----
segments(x0 = cod$age,  # X position (historical years)
         y0 = rep(cod_y_base, length(cod$age)),  # Base of the vertical line
         x1 = cod$age,  # Same X position for verticality
         y1 = cod_heights,  # Scaled height
         col = "grey45", lwd = 6)

# Label for cod detection lines
text(1825, 10000, "Cod Detections", col = "grey45")

# Add a reference line for 13300
abline(h = 13300, col = "black", lty = 2)

#optional horrible cods
#par(new = TRUE)  # Allows overlaying a second plot
#plot(cod$age, cod$x, col = "purple", pch = 16, axes = FALSE, xlab = "", ylab = "",
#     xlim = c(min(df_plot$Year), max(df_plot$Year)))  # Keep native scale
#axis(side = 4)  # Add secondary y-axis
#mtext("Cod Data", side = 4, line = 3, col = "purple")  # Label for second y-axis

dev.off()

# ---- Cod Catch vs Detection Proportion (Scatter Plot) ----
# Align cod$catch with df_plot$mean_catch based on matching years
cod$catch <- df_plot$mean_catch[match(cod$age, df_plot$Year)]

pdf("fisheries/cod.fishery.pdf",height = 4,width = 4)
# Scatter plot of cod detection proportion vs Icelandic Catch
plot(cod$x, cod$catch, pch = 16, col = "darkgreen",
     xlab = "Gadus Detection Proportion", 
     ylab = "Yearly Icelandic Catch (Mt)",
     main = "Cod Detection Proportion vs Icelandic Catch")

dev.off()

# Linear regression model
model <- lm(cod$catch~cod$x)
summary(model)



##### Now lets look at the data from Campana et al.!

campanaCod <- read.csv("fisheries/cod_estimates_by_century_Campana.csv")


# Define a function to classify centuries
classify_century <- function(year) {
  if (year > 0) {
    century <- (year %/% 100) + 1
    return(century)
  } else {
    century <- abs(year) %/% 100 + 1
    return(-century)
  }
}

# Apply function to the age column
cod$century <- sapply(cod$age, classify_century)

mean_values <- tapply(cod$x, cod$century, mean, na.rm = TRUE)
sd_values <- tapply(cod$x, cod$century, sd, na.rm = TRUE)


# Combine results into a data frame
result_table <- data.frame(
  Century = names(mean_values),
  Mean = mean_values,
  SD = sd_values
)


# Convert Century columns to character for consistent matching
result_table$Century <- as.character(result_table$Century)

# Filter result_table to only keep centuries present in campanaCod
filtered_result_table <- result_table[result_table$Century %in% as.character(campanaCod$Century), ]

campanaCod$meanMTBcod <- result_table$Mean[match(as.character(campanaCod$Century),result_table$Century,)]
campanaCod$sdMTBcod <- result_table$SD[match(as.character(campanaCod$Century),result_table$Century,)]


campanaCod$Century
cod$century

cod$match(cod$century,campanaCod$Century)

cod$x[match(cod$century,campanaCod$Century)]

pdf("fisheries/CampanaVsMTB.pdf",height = 3,width=6)
par(mfrow=c(1,2),mar=c(5.1,4.1,1.1,1.1))
plot(jitter(cod$x,0.4)~campanaCod$Adult_abundance[match(cod$century,campanaCod$Century)],pch=16,cex=0.8,ylab="Cod detection proportion",xlab="Adult (6 yr+) Abundance Index")
plot(jitter(cod$x,0.4)~campanaCod$Z[match(cod$century,campanaCod$Century)],pch=16,cex=0.8,xlab="Instantaneous mortality rate (Z)",ylab="")
dev.off()


campanaCod$Z[match(cod$century,campanaCod$Century)]

summary(lm(campanaCod$TotalCatch_000t~campanaCod$meanMTBcod))
summary(lm(campanaCod$Adult_abundance~campanaCod$meanMTBcod))
summary(lm(campanaCod$N~campanaCod$meanMTBcod))
summary(lm(campanaCod$G_index~campanaCod$meanMTBcod))
summary(lm(campanaCod$Z~campanaCod$meanMTBcod))
summary(lm(campanaCod$CC_Z~campanaCod$meanMTBcod))


plot(campanaCod$meanMTBcod,campanaCod$Z,pch=16)
plot(campanaCod$meanMTBcod,campanaCod$Z,pch=16)


campanaCod$meanMTBcod
plot(campanaCod$Century,campanaCod$meanMTBcod,pch="-",cex=3,ylim=c(0,8))
points(jitter(as.numeric(cod$century)),cod$x,pch=16)
points(campanaCod$Century,campanaCod$meanMTBcod,pch="-",cex=3,col="darkred")






# Lets run an ANOVA
shapiro.test(resid(aov(cod$x[cod$century>9]~as.factor(cod$century[cod$century>9]))))

# Homogeneity of variance
leveneTest(cod$x[cod$century>9]~as.factor(cod$century[cod$century>9]))


kruskal.test(cod$x[cod$century > 9] ~ as.factor(cod$century[cod$century > 9]))

pairwise.wilcox.test(cod$x[cod$century > 9], as.factor(cod$century[cod$century > 9]), p.adjust.method = "none")
pairwise.wilcox.test(cod$x, as.factor(cod$century))

pairwise.wilcox.test(cod$x, as.factor(cod$century), p.adjust.method = "none")
points(campanaCod$Century,campanaCod$TotalCatch_000t/10)

plot(campanaCod$meanMTBcod,campanaCod$TotalCatch_000t,ylim=c(0,50))

cod$centFact <- as.ordered(cod$century)
smallCod <- cod[cod$century > 9,]


library(brms)
bayes_model <- brm(x ~ centFact, data = smallCod, family = cumulative())
summary(bayes_model)

mean_values <- tapply(df_plot$mean_catch, sapply(df_plot$Year, classify_century), mean, na.rm = TRUE)


