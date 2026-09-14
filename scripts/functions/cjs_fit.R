##########################################
# cjs_fit for c_hat model fitting
##########################################
#purpose: Calculate goodness-of-fit for CJS models using Test 2 and Test 3
#returns: c-hat (vif), chi-square statistic, and degrees of freedom
#inputs:
#ng=1 number of "groups"
#ic=capture history (ch)
#ig=dimensions of groups
#min_expected= minimum expected count for chi-square


cjs_fit <- function(ch, beta, cap_X, surv_X, ints, min_expected = 1) {
  
  #calculate fitted probabilities from your model
  prob_matrices <- fill_prob_matrices(ch, beta, cap_X, surv_X, ints)
  p_hat <- prob_matrices$p_hat
  s_hat <- prob_matrices$s_hat
  
  nan <- nrow(ch)
  ns <- ncol(ch)
  
  # ------------------------------------------------------------------
  # TEST 2: Compare observed vs expected m-array
  # ------------------------------------------------------------------
  
  # Create OBSERVED m-array
  m_obs <- matrix(0, nrow = ns, ncol = ns)
  releases <- numeric(ns)
  
  for (i in 1:nan) {
    for (j in 1:(ns-1)) {
      if (ch[i, j] >= 1) {
        releases[j] <- releases[j] + 1
        found_next <- FALSE
        for (k in (j+1):ns) {
          if (ch[i, k] >= 1) {
            m_obs[j, k] <- m_obs[j, k] + 1
            found_next <- TRUE
            break
          }
        }
      }
    }
  }
  
  #create EXPECTED m-array from model predictions
  m_exp <- matrix(0, nrow = ns, ncol = ns)
  
  for (j in 1:(ns-1)) {
    if (releases[j] > 0) {
      for (i in 1:nan) {
        if (ch[i, j] >= 1) {
          # Probability of surviving from j to k and being captured at k
          # but not captured between j and k
          for (k in (j+1):ns) {
            surv_prob <- 1
            for (t in j:(k-1)) {
              surv_prob <- surv_prob * s_hat[i, t]
            }
            cap_prob <- p_hat[i, k]
            not_cap_prob <- 1
            for (t in (j+1):(k-1)) {
              not_cap_prob <- not_cap_prob * (1 - p_hat[i, t])
            }
            m_exp[j, k] <- m_exp[j, k] + surv_prob * cap_prob * not_cap_prob
          }
        }
      }
    }
  }
  
  #calculate Test 2 chi-square
  chi_sq_2 <- 0
  df_2 <- 0
  
  if (ns >= 4) {
    for (j in 2:(ns-2)) {
      # Compare newly marked vs previously marked
      prev_marked_obs <- sum(m_obs[1:(j-1), (j+1):ns])
      new_marked_obs <- sum(m_obs[j, (j+1):ns])
      prev_marked_exp <- sum(m_exp[1:(j-1), (j+1):ns])
      new_marked_exp <- sum(m_exp[j, (j+1):ns])
      
      # Skip if expected counts too small
      if (prev_marked_exp < min_expected || new_marked_exp < min_expected) {
        next
      }
      
      # Chi-square contribution (simplified - actual would be more detailed)
      chi_sq_2 <- chi_sq_2 + (prev_marked_obs - prev_marked_exp)^2 / prev_marked_exp
      chi_sq_2 <- chi_sq_2 + (new_marked_obs - new_marked_exp)^2 / new_marked_exp
      df_2 <- df_2 + 1
    }
  }
  
  # ------------------------------------------------------------------
  # TEST 3: Compare observed vs expected capture histories
  # ------------------------------------------------------------------
  
  chi_sq_3 <- 0
  df_3 <- 0
  
  for (j in 2:(ns-1)) {
    # Create 2x2 table: captured at j vs captured at j+1
    obs_table <- matrix(0, nrow = 2, ncol = 2)
    exp_table <- matrix(0, nrow = 2, ncol = 2)
    
    for (i in 1:nan) {
      # Only animals alive at j
      alive_at_j <- FALSE
      for (t in j:ns) {
        if (ch[i, t] >= 1) {
          alive_at_j <- TRUE
          break
        }
      }
      
      if (alive_at_j) {
        captured_j <- as.numeric(ch[i, j] >= 1)
        captured_j1 <- as.numeric(ch[i, j+1] >= 1)
        
        obs_table[captured_j + 1, captured_j1 + 1] <- 
          obs_table[captured_j + 1, captured_j1 + 1] + 1
        
        # Expected probability
        if (captured_j == 1) {
          # Animal captured at j
          exp_prob <- s_hat[i, j] * p_hat[i, j+1]
        } else {
          # Animal not captured at j but could be captured at j+1
          exp_prob <- s_hat[i, j] * p_hat[i, j+1] * (1 - p_hat[i, j])
        }
        
        # Simplified expected counts
        exp_table[captured_j + 1, captured_j1 + 1] <- 
          exp_table[captured_j + 1, captured_j1 + 1] + exp_prob
      }
    }
    
    # Calculate chi-square if expected counts sufficient
    if (sum(exp_table) > min_expected * 4) {
      for (r in 1:2) {
        for (c in 1:2) {
          if (exp_table[r, c] > min_expected) {
            chi_sq_3 <- chi_sq_3 + (obs_table[r, c] - exp_table[r, c])^2 / exp_table[r, c]
          }
        }
      }
      df_3 <- df_3 + 1
    }
  }
  
  # ------------------------------------------------------------------
  # Calculate c-hat
  # ------------------------------------------------------------------
  
  total_chi_sq <- chi_sq_2 + chi_sq_3
  total_df <- df_2 + df_3
  
  if (total_df > 0) {
    c_hat <- total_chi_sq / total_df
    c_hat <- max(1.0, c_hat)  # c-hat should be at least 1
  } else {
    c_hat <- 1.0
  }
  
  return(list(
    vif = c_hat,
    chigt = total_chi_sq,
    idfgt = total_df,
    diagnostic = list(
      test2_chi_sq = chi_sq_2,
      test2_df = df_2,
      test3_chi_sq = chi_sq_3,
      test3_df = df_3
    )
  ))
}

# ------------------------------------------------------------------
# Helper function for TEST 2
# ------------------------------------------------------------------
calculate_test2 <- function(ns, m, releases, min_expected, use_chat_rot, chat_rot) {
  
  chi_sq <- 0
  df <- 0
  
  # For each release occasion j from 2 to ns-2
  for (j in 2:(ns-2)) {
    
    # Create contingency table for release occasion j
    # Rows: Previously marked (releases before j) vs Newly marked (releases at j)
    # Columns: Recaptured at occasions j+1 to ns
    
    # Previously marked: sum of releases before j that were recaptured
    prev_marked <- numeric(ns - j)
    new_marked <- numeric(ns - j)
    
    for (k in (j+1):ns) {
      # Previously marked: animals released before j and recaptured at k
      prev_marked[k-j] <- sum(m[1:(j-1), k])
      # Newly marked: animals released at j and recaptured at k
      new_marked[k-j] <- m[j, k]
    }
    
    # Row totals
    prev_total <- sum(prev_marked)
    new_total <- sum(new_marked)
    col_totals <- prev_marked + new_marked
    
    # Column totals (total recaptures at each occasion)
    total_animals <- prev_total + new_total
    
    # Skip if insufficient data
    if (total_animals < min_expected) {
      next
    }
    
    # Calculate expected values and chi-square
    for (k in 1:length(prev_marked)) {
      # Expected counts under independence
      exp_prev <- prev_total * col_totals[k] / total_animals
      exp_new <- new_total * col_totals[k] / total_animals
      
      # Skip if expected count too small
      if (exp_prev < min_expected || exp_new < min_expected) {
        next
      }
      
      # Chi-square contribution
      chi_sq <- chi_sq + ((prev_marked[k] - exp_prev)^2 / exp_prev)
      chi_sq <- chi_sq + ((new_marked[k] - exp_new)^2 / exp_new)
      df <- df + 1
    }
    
    # Adjust degrees of freedom (for each column we lose 1 df)
    df <- df - 1  # Because column totals sum to row totals
  }
  
  return(list(chi_sq = chi_sq, df = max(0, df)))
}

# ------------------------------------------------------------------
# Helper function for TEST 3
# ------------------------------------------------------------------
calculate_test3 <- function(nan, ns, ic, min_expected, use_chat_rot, chat_rot) {
  
  chi_sq <- 0
  df <- 0
  
  # For each occasion j from 2 to ns-1
  for (j in 2:(ns-1)) {
    
    # Create 2x2 contingency table for occasion j
    # Rows: Captured at j (yes/no)
    # Columns: Captured at j+1 (yes/no)
    # This tests if capture probability depends on previous capture
    
    contingency <- matrix(0, nrow = 2, ncol = 2)
    
    for (i in 1:nan) {
      # Only consider animals alive at j (captured at or after j)
      if (any(ic[i, j:ns] >= 1)) {
        captured_j <- as.numeric(ic[i, j] >= 1)
        captured_j1 <- as.numeric(ic[i, j+1] >= 1)
        
        contingency[captured_j + 1, captured_j1 + 1] <- 
          contingency[captured_j + 1, captured_j1 + 1] + 1
      }
    }
    
    # Skip if any cell has insufficient counts
    if (any(contingency < min_expected)) {
      next
    }
    
    # Calculate expected values
    row_totals <- rowSums(contingency)
    col_totals <- colSums(contingency)
    total <- sum(contingency)
    
    # Skip if totals are too small
    if (total < min_expected * 4) {  # All 4 cells need minimum
      next
    }
    
    # Calculate chi-square
    for (r in 1:2) {
      for (c in 1:2) {
        expected <- row_totals[r] * col_totals[c] / total
        if (expected >= min_expected) {
          chi_sq <- chi_sq + ((contingency[r, c] - expected)^2 / expected)
        }
      }
    }
    
    df <- df + 1  # Each 2x2 table contributes 1 degree of freedom
  }
  
  return(list(chi_sq = chi_sq, df = max(0, df)))
}

# ------------------------------------------------------------------
# Simplified wrapper function for model fitting
# ------------------------------------------------------------------
cjs_fit_simple <- function(ch, beta, cap_X, surv_X, ints) {
  # Simple wrapper that matches how you're calling it in CJS_run()
  result <- calculate_simple_c_hat(
    ch = as.matrix(ch),
    beta = beta,
    cap_X = as.matrix(cap_X),
    surv_X = as.matrix(surv_X),
    ints = as.matrix(ints)
  )
  return(result)
}

calculate_simple_c_hat <- function(ch, beta, cap_X, surv_X, ints) {
  nan <- nrow(ch)
  ns <- ncol(ch)
  
  # Get predicted probabilities
  probs <- fill_prob_matrices(ch, beta, cap_X, surv_X, ints)
  p_hat <- probs$p_hat
  
  # Calculate Pearson residuals
  residuals <- matrix(0, nan, ns)
  for (i in 1:nan) {
    for (j in 1:ns) {
      if (!is.na(p_hat[i, j]) && p_hat[i, j] > 0 && p_hat[i, j] < 1) {
        observed <- as.numeric(ch[i, j] >= 1)
        expected <- p_hat[i, j]
        residuals[i, j] <- (observed - expected) / sqrt(expected * (1 - expected))
      }
    }
  }
  
  # Simple c-hat estimate
  c_hat <- sum(residuals^2, na.rm = TRUE) / (nan * ns - length(beta))
  return(max(1.0, c_hat))
}