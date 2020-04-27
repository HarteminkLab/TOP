test_that("predict_coef_BH_mean()", {
  data <- data.frame(pwm_score = c(5,10),
                     DNase_1 = c(1,5),
                     DNase_2 = c(20,30),
                     DNase_3 = c(1,0),
                     DNase_4 = c(10,20),
                     DNase_5 = c(2,3))

  coef_mean <- c(intercept = 5, pwm_score = 3, c(1, 3, -1, 2, 1))
  predicted <- predict_coef_BH_mean(data, coef_mean, transform = 'asinh')

  data_full <- as.matrix(data.frame(intercept = 1, data, check.names = F))
  coefficients <- as.matrix(coef_mean, ncol = 1)
  predictions <- as.numeric(data_full %*% coefficients)
  predictions <- sinh(predictions)

  expect_equal(predicted, predictions)
})
