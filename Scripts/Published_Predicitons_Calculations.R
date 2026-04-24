# SPAN-100 --------------------------------------

# Formula
# SPAN-100 = Age + NIHSS  
# SPAN-100 positive = SPAN-100 ≥ 100

# Calculation
# SPAN-100 calculation
span100 <- age + nihss
span_positive <- span100 >= 100

# Output
cat("SPAN-100 Score:", span100, "\nSPAN-100 Positive:", span_positive)

# iScore (Canadian Score) --------------------------------------

# Formula / Calculation online
# https://www.sorcan.ca/iscore


# ASTRAL Score (logistic regression model) --------------------------------------

# Calculation
# Input values
age <- 78
nihss <- 13
time_admission <- 300  # in minutes
glucose <- 7.5
loc <- TRUE
visual_loss <- TRUE
prestroke_mrs_gt1 <- TRUE

# Calculate linear predictor (logit)
logit <- -5.5 +
    0.4 * (age / 10) +
    0.12 * nihss +
    0.003 * time_admission +
    ifelse(loc, 0.75, 0) +
    ifelse(visual_loss, 0.62, 0) +
    0.1 * glucose +
    ifelse(prestroke_mrs_gt1, 0.8, 0)

# Convert logit to probability
p_poor_outcome <- 1 / (1 + exp(-logit))

# Output
cat("ASTRAL predicted probability of poor outcome:", round(p_poor_outcome, 3))

# PREP2 Algorithm (decision tree) --------------------------------------

# Calculation
predict_PREP2 <- function(safe_score, mep_present, age) {
    if (safe_score >= 5) {
        if (age < 80) {
            return("Excellent")
        } else {
            return("Good")
        }
    } else {
        if (mep_present) {
            return("Good")
        } else {
            if (age < 80) {
                return("Limited")
            } else {
                return("Poor")
            }
        }
    }
}

# Example usage
safe_score <- 3
mep_present <- FALSE
age <- 82

result <- predict_PREP2(safe_score, mep_present, age)
cat("PREP2 Prediction:", result)
