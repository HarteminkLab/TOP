test_that("predict_TOP_mean_coef()", {
  data <- data.frame(pwm_score = c(5,10),
                     bin1 = c(1,5),
                     bin2 = c(20,30),
                     bin3 = c(1,0),
                     bin4 = c(10,20),
                     bin5 = c(2,3))

  mean_coef <- c(intercept = 5, pwm_score = 3, bin = c(1, 3, -1, 2, 1))

  res <- predict_TOP_mean_coef(data, mean_coef, transform = 'asinh')
  TOP_predictions <- res$predicted

  data_full <- as.matrix(data.frame(intercept = 1, data, check.names = F))
  coefficients <- as.matrix(mean_coef, ncol = 1)
  correct_predictions <- sinh(as.numeric(data_full %*% coefficients))

  expect_equal(TOP_predictions, correct_predictions)
})
